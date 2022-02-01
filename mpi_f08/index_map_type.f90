!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2022  Neil N. Carlson
!!
!! Permission is hereby granted, free of charge, to any person obtaining a
!! copy of this software and associated documentation files (the "Software"),
!! to deal in the Software without restriction, including without limitation
!! the rights to use, copy, modify, merge, publish, distribute, sublicense,
!! and/or sell copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following conditions:
!!
!! The above copyright notice and this permission notice shall be included
!! in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!! DEALINGS IN THE SOFTWARE.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module index_map_type

  use mpi_f08
  implicit none
  private

  type, public :: index_map
    type(mpi_comm) :: comm
    integer :: onp_size = 0    ! number of indices assigned to this process (on-process)
    integer :: offp_size = 0   ! number of off-process indices referenced from this process
    integer :: local_size = 0  ! number of local indices (on and off-process)
    integer :: global_size = 0 ! size of the global index set
    integer :: first        ! first global index of the range assigned to this process
    integer :: last         ! last global index of the range assigned to this process
    integer, allocatable :: offp_index(:), send_index(:)
    integer, allocatable :: send_counts(:), send_displs(:)
    integer, allocatable :: recv_counts(:), recv_displs(:)
  contains
    procedure :: init
    procedure :: gather
    procedure :: global_index
  end type

contains

  subroutine init(this, comm, bsize, offp_index)

    class(index_map), intent(out) :: this
    type(mpi_comm), intent(in) :: comm
    integer, intent(in) :: bsize, offp_index(:)

    integer :: nproc, ierr

    this%comm = comm

    call MPI_Comm_size(comm, nproc)

    this%onp_size = bsize
    this%offp_size = 0
    this%local_size = this%onp_size + this%offp_size
    call MPI_Scan(bsize, this%last, 1, MPI_INTEGER, MPI_SUM, this%comm)
    this%first = this%last - this%onp_size + 1
    this%global_size = this%last
    call MPI_Bcast(this%global_size, 1, MPI_INTEGER, nproc-1, this%comm)

    call add_offp_index(this, offp_index)

  end subroutine init


  subroutine add_offp_index(this, offp_index)

    class(index_map), intent(inout) :: this
    integer, intent(in) :: offp_index(:)

    integer :: ierr, nproc, my_rank, rank, i, j, j1, n
    integer, allocatable :: last(:), offp_rank(:), send_index(:)
    integer, allocatable :: recv_counts(:), recv_ranks(:)
    integer, allocatable :: send_counts(:), send_ranks(:)
    type(mpi_comm) :: new_comm

    !TODO: ensure offp_index is strictly increasing

    this%offp_index = offp_index
    this%offp_size  = size(this%offp_index)
    this%local_size = this%onp_size + this%offp_size

    call MPI_Comm_size(this%comm, nproc)
    call MPI_Comm_rank(this%comm, my_rank)

    !! Determine which rank owns each off-process index (OFFP_RANK).
    !! OFFP_RANK will be ordered if OFFP_INDEX is ordered (we need this later).
    allocate(last(0:nproc-1), offp_rank(this%offp_size))
    call MPI_Allgather(this%last, 1, MPI_INTEGER, last, 1, MPI_INTEGER, this%comm)
    rank = 0
    do j = 1, size(this%offp_index)
      do while (this%offp_index(j) > last(rank))
        rank = rank + 1
        !ASSERT(rank < nproc)
      end do
      if (rank == my_rank) stop 1 !TODO: replace by assertion
      offp_rank(j) = rank
    end do
    deallocate(last)

    !! We know which ranks we will need to receive indexed data from, but not
    !! the corresponding info of which ranks we will need to send indexed data.
    !! The next step is to communicate that info.

    !! Count the number of off-process indices for each rank (RECV_COUNTS).
    !! This relies on the OFFP_RANK array being ordered.
    allocate(recv_counts(0:nproc-1))
    j1 = 1
    do rank = 0, nproc-1
      do j = j1, size(offp_rank)
        if (offp_rank(j) > rank) exit
      end do
      recv_counts(rank) = j - j1
      j1 = j
    end do
    deallocate(offp_rank)

    !! Distribute the number of indices received from each rank to that rank.
    !! The result is the number of indices sent to each rank (SEND_COUNTS).
    allocate(send_counts(0:nproc-1))
    call MPI_Alltoall(recv_counts, 1, MPI_INTEGER, send_counts, 1, MPI_INTEGER, this%comm)

    !! We now know which ranks we will be receiving indexed data from, and which
    !! ranks we will be sending indexed data to.  This all that is needed to
    !! define a virtual process topology.

    !! Compress the RECV_COUNTS array and generate the list of ranks that
    !! we receive data from: RECV_RANKS is the list of incoming neighbors.
    n = count(recv_counts > 0)
    allocate(recv_ranks(n), this%recv_counts(n))
    n = 0
    do rank = 0, nproc-1
      if (recv_counts(rank) > 0) then
        n = n + 1
        recv_ranks(n) = rank
        this%recv_counts(n) = recv_counts(rank)
      end if
    end do
    deallocate(recv_counts)

    !! Compress the SEND_COUNTS array and generate the list of ranks that
    !! we send data to: SEND_RANKS is the list of outgoing neighbors.
    n = count(send_counts > 0)
    allocate(send_ranks(n), this%send_counts(n))
    n = 0
    do rank = 0, nproc-1
      if (send_counts(rank) > 0) then
        n = n + 1
        send_ranks(n) = rank
        this%send_counts(n) = send_counts(rank)
      end if
    end do
    deallocate(send_counts)

    !! Create the virtual topology. Here we use the RECV_COUNTS and SEND_COUNTS
    !! as the edge weights (used if rank reordering is enabled). An alternative
    !! is to replace them with MPI_UNWEIGHTED.
    call MPI_Dist_graph_create_adjacent(this%comm, &
        size(recv_ranks), recv_ranks, this%recv_counts, &
        size(send_ranks), send_ranks, this%send_counts, &
        MPI_INFO_NULL, .false., new_comm)
    this%comm = new_comm

    !! The components %RECV_COUNTS, %RECV_DISPLS, %SEND_COUNTS and %SEND_DIPSLS
    !! initialized here are meant for use with MPI_Neighbor_alltoallv. The
    !! component %SEND_INDEX will be used to fill the send buffer; the receive
    !! buffer will be the off-process data array itself.

    !! Generate displacements into the send buffer for the start of the data
    !! for each neighbor rank. This is the exclusive sum scan of %SEND_COUNT.
    allocate(this%send_displs, mold=this%send_counts)
    if (size(this%send_displs) > 0) then
      this%send_displs(1) = 0
      do i = 2, size(this%send_displs)
        this%send_displs(i) = this%send_displs(i-1) + this%send_counts(i-1)
      end do
    end if

    !! Generate displacements into the receive buffer for the start of the data
    !! for each neighbor rank. This is the exclusive sum scan of %RECV_COUNT.
    allocate(this%recv_displs, mold=this%recv_counts)
    if (size(this%recv_displs) > 0) then
      this%recv_displs(1) = 0
      do i = 2, size(this%recv_displs)
        this%recv_displs(i) = this%recv_displs(i-1) + this%recv_counts(i-1)
      end do
    end if

    !! Communicate the global off-process indices to their owning ranks.
    !! NB: We send the "receive index" data and receive the "send index" data.
    allocate(this%send_index(sum(this%send_counts)))
    call MPI_Neighbor_alltoallv(this%offp_index, this%recv_counts, this%recv_displs, MPI_INTEGER, &
        this%send_index, this%send_counts, this%send_displs, MPI_INTEGER, this%comm, ierr)
    this%send_index = this%send_index - this%first + 1  ! map to local indices
    !ASSERT(all(this%send_index >= 1 .and. this%send_index <= this%onp_size))

  end subroutine add_offp_index

  elemental function global_index(this, n) result(gid)
    class(index_map), intent(in) :: this
    integer, intent(in) :: n
    integer :: gid
    gid = -1
    if (n < 1) return
    if (n <= this%onp_size) then
      gid = this%first + n - 1
    else if (n <= this%local_size) then
      gid = this%offp_index(n-this%onp_size)
    end if
  end function

  subroutine gather(this, local_data)
    class(index_map), intent(in) :: this
    integer, intent(inout) :: local_data(:)
    call gather_aux(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  subroutine gather_aux(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer, intent(in)  :: onp_data(:)
    integer, intent(out) :: offp_data(:)
    integer :: ierr
    integer, allocatable :: send_buf(:)
    send_buf = onp_data(this%send_index)
    call MPI_Neighbor_alltoallv(send_buf, this%send_counts, this%send_displs, MPI_INTEGER, &
        offp_data, this%recv_counts, this%recv_displs, MPI_INTEGER, this%comm)
  end subroutine

end module index_map_type
