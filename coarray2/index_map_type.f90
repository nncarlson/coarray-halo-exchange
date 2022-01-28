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

  implicit none
  private

  type, public :: index_map
    integer :: onP_size = 0    ! number of indices assigned to this process (on-process)
    integer :: offP_size = 0   ! number of off-process indices referenced from this process
    integer :: local_size = 0  ! number of local indices (on and off-process)
    integer :: global_size = 0 ! size of the global index set
    integer :: first        ! first global index of the range assigned to this process
    integer :: last         ! last global index of the range assigned to this process
    integer, allocatable :: offP_index(:)
    integer, allocatable :: recv_displs(:), recv_counts(:)
    integer, allocatable :: send_displs(:), send_image(:), send_index(:)
    integer, allocatable :: src_pe(:), src_id(:)
    integer, allocatable :: offset(:) !
  contains
    procedure :: init
    procedure :: gather
    procedure :: global_index
  end type

contains

  subroutine init(this, bsize, offP_index)

    class(index_map), intent(out) :: this
    integer, intent(in) :: bsize, offP_index(:)

    integer :: nPE, this_PE
    integer, allocatable :: bsizes(:)[:]

    nPE = num_images()
    this_PE = this_image()

    !! MPI_Allgather of bsize
    allocate(bsizes(nPE)[*])
    bsizes(this_PE)[1] = bsize
    sync all
    call co_broadcast(bsizes, 1)
    sync all

    this%onP_size = bsizes(this_PE)
    this%last = sum(bsizes(1:this_PE))
    this%first = this%last - this%onP_size + 1
    this%local_size = this%onP_size
    this%global_size = sum(bsizes)

    call add_offP_index(this, offP_index)

    sync all

  end subroutine init

  subroutine add_offP_index(this, offP_index)

    class(index_map), intent(inout), target :: this
    integer, intent(in) :: offP_index(:)

    integer :: npe, i, j, j1, n, k
    integer, allocatable :: offP_image(:), send_counts(:), recv_displs(:)
    integer, allocatable :: last(:)[:], recv_counts(:)[:], send_displs(:)[:]

    type box
      integer, pointer :: data(:)
    end type
    type(box), allocatable :: buffer[:]
    

    !TODO: ensure offP_index is strictly increasing

    this%offP_index = offP_index
    this%offP_size  = size(this%offP_index)
    this%local_size = this%onP_size + this%offP_size

    npe = num_images()

    !! Determine the image that owns each off-process index (OFFP_IMAGE).
    !! OFFP_IMAGE will be ordered if OFFP_INDEX is ordered (exploited later).
    allocate(last(npe)[*], offP_image(this%offP_size))
    last(this_image())[1] = this%last
    sync all
    call co_broadcast(last, 1)
    sync all
    i = 1
    do j = 1, size(this%offP_index)
      do while (this%offP_index(j) > last(i))
        i = i + 1
        !ASSERT(i <= num_images())
      end do
      if (i == this_image()) stop 1 !TODO: replace by assertion
      offp_image(j) = i
    end do
    
    !! We know which images we will need to receive indexed data from, but not
    !! the corresponding info of which images we will need to send indexed data.
    !! The next step is to communicate that info.

    !! Count the number of off-process indices received from each image
    !! (RECV_COUNTS). This relies on the OFFP_IMAGE array being ordered.
    call move_alloc(last, recv_counts)  ! reuse coarray
    j1 = 1
    do i = 1, npe
      do j = j1, size(offp_image)
        if (offp_image(j) > i) exit
      end do
      recv_counts(i) = j - j1
      j1 = j
    end do
    sync all
    
    !! Distribute the number of indices received from each image to that image.
    !! The result is the number of indices sent to each image (SEND_COUNTS).
    allocate(send_counts(npe))
    do i = 1, npe
      send_counts(i) = recv_counts(this_image())[i]
    end do
    sync all
    n = sum(send_counts)
    allocate(this%send_index(n))

    !! Now that we know how many indices are being sent to/received from each
    !! image, we can communicate the indices themselves.

    !! Generate the displacements into the index data arrays for the start of
    !! the data for each rank. These are exclusive sum scans of the count arrays
    allocate(recv_displs(npe), send_displs(npe)[*])
    recv_displs(1) = 0
    send_displs(1) = 0
    do i = 2, npe
      recv_displs(i) = recv_displs(i-1) + recv_counts(i-1)
      send_displs(i) = send_displs(i-1) + send_counts(i-1)
    end do
    deallocate(send_counts)
    sync all

    !! Compress the RECV_COUNTS array and generate the list of ranks that
    !! we receive data from: RECV_RANKS is the list of incoming neighbors.
    n = count(recv_counts > 0)
    allocate(this%send_image(n), this%recv_counts(n), this%recv_displs(n), this%send_displs(n))
    n = 0
    do i = 1, npe
      if (recv_counts(i) > 0) then
        n = n + 1
        this%send_image(n) = i
        this%recv_counts(n) = recv_counts(i)
        this%recv_displs(n) = recv_displs(i)
        this%send_displs(n) = send_displs(this_image())[i]
      end if
    end do
    deallocate(send_displs, recv_displs, recv_counts)
    
    allocate(buffer[*])
    buffer%data => this%send_index
    sync all
    do j = 1, size(this%recv_counts)
      associate (i => this%send_image(j), n => this%recv_counts(j), &
                 s0 => this%send_displs(j), r0 => this%recv_displs(j))
        buffer[i]%data(s0+1:s0+n) = this%offP_index(r0+1:r0+n)
      end associate
    end do
    sync all
    this%send_index = this%send_index - this%first + 1
      
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
    integer, intent(in), target :: onP_data(:)
    integer, intent(out) :: offP_data(:)

    integer :: j

    type box
      !integer, allocatable :: data(:)
      integer, pointer :: data(:)
    end type
    type(box), allocatable :: send_buf[:]
    
    allocate(send_buf[*])
    allocate(send_buf%data(size(this%send_index)))
    send_buf%data = onP_data(this%send_index)
    sync all
    
    do j = 1, size(this%recv_counts)
      associate (i  => this%send_image(j), n => this%recv_counts(j), &
                 r0 => this%recv_displs(j), s0 => this%send_displs(j))
        offP_data(r0+1:r0+n) = send_buf[i]%data(s0+1:s0+n)
      end associate
    end do
    sync all

  end subroutine gather_aux

end module index_map_type
