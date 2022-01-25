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
    integer :: onP_size = 0    ! number of indices assigned to this process (on-process)
    integer :: offP_size = 0   ! number of off-process indices referenced from this process
    integer :: local_size = 0  ! number of local indices (on and off-process)
    integer :: global_size = 0 ! size of the global index set
    integer :: first        ! first global index of the range assigned to this process
    integer :: last         ! last global index of the range assigned to this process
    integer, allocatable :: offP_index(:)
    integer, allocatable :: send_counts(:), sdispls(:), send_index(:)
    integer, allocatable :: recv_counts(:), rdispls(:)
  contains
    procedure :: init
    procedure :: gather
    procedure :: global_index
  end type

contains

  subroutine init(this, comm, bsize, offP_index)

    type(mpi_comm), intent(in) :: comm
    class(index_map), intent(out) :: this
    integer, intent(in) :: bsize, offP_index(:)

    integer :: nproc, my_rank, pe, ierr
    integer, allocatable :: bsizes(:)

    this%comm = comm

    call MPI_Comm_size(comm, nproc)
    call MPI_Comm_rank(comm, my_rank)

    allocate(bsizes(0:nproc-1))
    call MPI_Allgather(bsize, 1, MPI_INTEGER, bsizes, 1, MPI_INTEGER, this%comm)

    this%onP_size = bsizes(my_rank)
    this%offP_size = 0
    this%last = sum(bsizes(0:my_rank))
    this%first = this%last - this%onP_size + 1
    this%local_size = this%onP_size
    this%global_size = sum(bsizes)

    call add_offP_index(this, offP_index)

  end subroutine init


  subroutine add_offP_index(this, offP_index)

    class(index_map), intent(inout) :: this
    integer, intent(in) :: offP_index(:)

    integer :: ierr, nproc, my_rank, rank, i, j, j1, n
    integer, allocatable :: last(:), offP_rank(:), send_index(:)
    integer, allocatable :: recv_counts(:), recv_displs(:), recv_ranks(:)
    integer, allocatable :: send_counts(:), send_displs(:), send_ranks(:)
    type(mpi_comm) :: new_comm

    this%offP_index = offP_index
    !TODO: ensure offP_index is strictly increasing

    this%offP_size  = size(this%offP_index)
    this%local_size = this%onP_size + this%offP_size

    call MPI_Comm_size(this%comm, nproc)
    call MPI_Comm_rank(this%comm, my_rank)
    
    !! Determine which rank owns each off-process index (OFFP_RANK).
    !! OFFP_RANK will be ordered if OFFP_INDEX is ordered (we need this later).
    allocate(last(0:nproc-1), offP_rank(this%offP_size))
    call MPI_Allgather(this%last, 1, MPI_INTEGER, last, 1, MPI_INTEGER, this%comm)
    rank = 0
    do j = 1, size(this%offP_index)
      do while (this%offP_index(j) > last(rank))
        rank = rank + 1
        !ASSERT(rank < nproc)
      end do
      if (rank == my_rank) stop 1 !TODO: replace by assertion
      offP_rank(j) = rank
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
      do j = j1, size(offP_rank)
        if (offP_rank(j) > rank) exit
      end do
      recv_counts(rank) = j - j1
      j1 = j
    end do
    deallocate(offP_rank)

    !! Distribute the number of indices received from each rank to that rank.
    !! The result is the number of indices sent to each rank (SEND_COUNTS).
    allocate(send_counts(0:nproc-1))
    call MPI_Alltoall(recv_counts, 1, MPI_INTEGER, send_counts, 1, MPI_INTEGER, this%comm)

    !! Now that we know how many indices are being sent to/received from each
    !! rank, we can communicate the indices themselves.
    !! Generate the displacements into the index data arrays for the start of
    !! the data for each rank. These are exclusive sum scans of the count arrays
    allocate(recv_displs(0:nproc-1), send_displs(0:nproc-1))
    recv_displs(0) = 0
    send_displs(0) = 0
    do rank = 1, nproc-1
      recv_displs(rank) = recv_displs(rank-1) + recv_counts(rank-1)
      send_displs(rank) = send_displs(rank-1) + send_counts(rank-1)
    end do
    !! NB: We send the "receive index" data and receive the "send index" data.
    allocate(send_index(sum(send_counts)))
    call MPI_Alltoallv(this%offp_index, recv_counts, recv_displs, MPI_INTEGER, &
                            send_index, send_counts, send_displs, MPI_INTEGER, &
                       this%comm)
    deallocate(recv_displs, send_displs)

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

    !! The components %RECV_COUNTS, %RDISPLS, %SEND_COUNTS, and %SDIPSLS
    !! initialized here are meant for use with MPI_Neighbor_alltoallv.
    !! The component %SEND_INDEX will be used to fill the send buffer;
    !! the receive buffer will be the off-process array itself.

    !! Convert the global SEND_INDEX values to process-local indices.
    !ASSERT(all(send_index >= this%first))
    !ASSERT(all(send_index <= this%last))
    this%send_index = send_index - this%first + 1

    !! Generate the displacements into the send buffer for the start of the
    !! data for each neighbor. This is the exclusive sum scan of %SEND_COUNT.
    allocate(this%sdispls, mold=this%send_counts)
    if (size(this%sdispls) > 0) then
      this%sdispls(1) = 0
      do i = 2, size(this%sdispls)
        this%sdispls(i) = this%sdispls(i-1) + this%send_counts(i-1)
      end do
    end if

    !! Generate the displacements into the receive buffer for the start of the
    !! data for each neighbor. This is the exclusive sum scan of %RECV_COUNT.
    allocate(this%rdispls, mold=this%recv_counts)
    if (size(this%rdispls) > 0) then
      this%rdispls(1) = 0
      do i = 2, size(this%rdispls)
        this%rdispls(i) = this%rdispls(i-1) + this%recv_counts(i-1)
      end do
    end if

  end subroutine add_offP_index

  elemental function global_index(this, n) result(gid)
    class(index_map), intent(in) :: this
    integer, intent(in) :: n
    integer :: gid
    gid = -1
    if (n < 1) return
    if (n <= this%onP_size) then
      gid = this%first + n - 1
    else if (n <= this%local_size) then
      gid = this%offP_index(n-this%onP_size)
    end if
  end function

  subroutine gather(this, local_data)
    class(index_map), intent(in) :: this
    integer, intent(inout) :: local_data(:)
    call gather_aux(this, local_data(:this%onP_size), local_data(this%onP_size+1:))
  end subroutine

  subroutine gather_aux(this, onP_data, offP_data)
    class(index_map), intent(in) :: this
    integer, intent(in)  :: onP_data(:)
    integer, intent(out) :: offP_data(:)
    integer :: ierr
    integer, allocatable :: send_buf(:)
    send_buf = onP_data(this%send_index)
    call MPI_Neighbor_alltoallv(send_buf, this%send_counts, this%sdispls, MPI_INTEGER, &
        offp_data, this%recv_counts, this%rdispls, MPI_INTEGER, this%comm)
  end subroutine

end module index_map_type
