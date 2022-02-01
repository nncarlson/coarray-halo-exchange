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
    integer, allocatable :: offp_index(:)
    integer, allocatable :: src_image(:), src_index(:)
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

    class(index_map), intent(inout) :: this
    integer, intent(in) :: offp_index(:)

    integer :: nproc, i, j
    integer, allocatable :: last(:)[:]

    !TODO: ensure offp_index is strictly increasing

    this%offp_index = offp_index
    this%offp_size  = size(this%offp_index)
    this%local_size = this%onp_size + this%offp_size

    nproc = num_images()

    !! Gather the last global index owned by each image
    allocate(last(0:nproc)[*])
    if (this_image() == 1) last(0) = 0
    last(this_image())[1] = this%last
    sync all
    call co_broadcast(last, 1)

    !! Determine the image that owns each off-process index (SRC_IMAGE)
    !! and the corresponding local index in that image (SRC_INDEX).
    allocate(this%src_image(this%offp_size), this%src_index(this%offp_size))
    i = 1
    do j = 1, size(this%offp_index)
      do while (this%offp_index(j) > last(i))
        i = i + 1
        !ASSERT(i <= nproc)
      end do
      if (i == this_image()) stop 1 !TODO: replace by assertion
      this%src_image(j) = i
      this%src_index(j) = this%offp_index(j) - last(i-1)
    end do

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
    integer, intent(in), target :: onp_data(:)
    integer, intent(out) :: offp_data(:)

    integer :: j

    type box
      integer, pointer :: data(:)
    end type
    type(box), allocatable :: src[:]
    allocate(src[*])

    src%data => onp_data
    sync all
    do j = 1, this%offp_size
      offp_data(j) = src[this%src_image(j)]%data(this%src_index(j))
    end do

  end subroutine gather_aux

end module index_map_type
