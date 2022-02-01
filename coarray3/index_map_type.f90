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

  use coarray_collectives
  implicit none
  private

  type, public :: index_map
    integer :: onp_size = 0    ! number of indices assigned to this process (on-process)
    integer :: offp_size = 0   ! number of off-process indices referenced from this process
    integer :: local_size = 0  ! number of local indices (on and off-process)
    integer :: global_size = 0 ! size of the global index set
    integer :: first        ! first global index of the range assigned to this process
    integer :: last         ! last global index of the range assigned to this process
    integer, allocatable :: offp_index(:), send_index(:)
    integer, allocatable :: dest_image(:), dest_index(:)
  contains
    procedure :: init
    procedure :: gather
    procedure :: global_index
  end type

contains

  subroutine init(this, bsize, offp_index)

    class(index_map), intent(out) :: this
    integer, intent(in) :: bsize, offp_index(:)

    integer :: nproc

    nproc = num_images()

    this%onp_size = bsize
    this%offp_size = 0
    this%local_size = this%onp_size + this%offp_size
    call co_sum_scan(bsize, this%last)
    this%first = this%last - this%onp_size + 1
    this%global_size = this%last
    call co_broadcast(this%global_size, nproc)

    call add_offp_index(this, offp_index)

  end subroutine init

  subroutine add_offp_index(this, offp_index)

    class(index_map), intent(inout), target :: this
    integer, intent(in) :: offp_index(:)

    integer :: nproc, i, j, j1, n, offset
    integer, allocatable :: offp_image(:), send_counts(:)
    integer, allocatable :: last(:)[:], recv_counts(:)[:], recv_displs(:)[:]

    type box
      integer, pointer :: data(:)
    end type
    type(box), allocatable :: buffer[:]

    !TODO: ensure offp_index is strictly increasing

    this%offp_index = offp_index
    this%offp_size  = size(this%offp_index)
    this%local_size = this%onp_size + this%offp_size

    nproc = num_images()

    !! Determine the image that owns each off-process index (OFFP_IMAGE).
    !! OFFP_IMAGE will be ordered if OFFP_INDEX is ordered (exploited later).
    allocate(last(nproc)[*], offp_image(this%offp_size))
    last(this_image())[1] = this%last
    sync all
    call co_broadcast(last, 1)
    i = 1
    do j = 1, size(this%offp_index)
      do while (this%offp_index(j) > last(i))
        i = i + 1
        !ASSERT(i <= nproc)
      end do
      if (i == this_image()) stop 1 !TODO: replace by assertion
      offp_image(j) = i
    end do

    !! Count the number of off-process indices received from each image
    !! (RECV_COUNTS). This relies on the OFFP_IMAGE array being ordered.
    call move_alloc(last, recv_counts)  ! reuse coarray
    j1 = 1
    do i = 1, nproc
      do j = j1, size(offp_image)
        if (offp_image(j) > i) exit
      end do
      recv_counts(i) = j - j1
      j1 = j
    end do
    deallocate(offp_image)

    !! Generate the displacements into the off-process data arrays for the
    !! start of the data coming from each image.
    allocate(recv_displs(nproc)[*])
    recv_displs(1) = 0
    do i = 2, nproc
      recv_displs(i) = recv_displs(i-1) + recv_counts(i-1)
    end do

    !! Distribute the number of indices received from each image to that image.
    !! The result is the number of indices sent to each image (SEND_COUNTS).
    allocate(send_counts(nproc))
    sync all
    do i = 1, nproc
      send_counts(i) = recv_counts(this_image())[i]
    end do

    !! The components %DEST_IMAGE and %DEST_INDEX initialized here establish
    !! the coarray communication pattern used by GATHER_AUX. The component
    !! %SEND_INDEX generated next will be used to access elements of the
    !! on-process data array.

    n = sum(send_counts)
    allocate(this%send_index(n), this%dest_image(n), this%dest_index(n))
    deallocate(recv_counts)
    sync all

    n = 0
    do i = 1, nproc
      offset = recv_displs(this_image())[i]
      do j = 1, send_counts(i)
        n = n + 1
        this%dest_image(n) = i
        this%dest_index(n) = offset + j
      end do
    end do
    deallocate(recv_displs, send_counts)

    !! Communicate the global off-process indices to their owning images.
    allocate(buffer[*])
    buffer%data => this%offp_index
    sync all
    do j = 1, size(this%send_index)
      this%send_index(j) = buffer[this%dest_image(j)]%data(this%dest_index(j))
    end do
    sync all
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
    integer, intent(in) :: onp_data(:)
    integer, intent(out), target :: offp_data(:)

    integer :: j

    type box
      integer, pointer :: data(:)
    end type
    type(box), allocatable :: dest[:]
    allocate(dest[*])

    dest%data => offp_data
    sync all
    do j = 1, size(this%send_index)
      dest[this%dest_image(j)]%data(this%dest_index(j)) = onp_data(this%send_index(j))
    end do
    sync all

  end subroutine gather_aux

end module index_map_type
