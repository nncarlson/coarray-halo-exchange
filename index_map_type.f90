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

    integer :: nPE, this_PE, pe
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

    !! offset is exclusive prefix sum of bsizes
    allocate(this%offset(nPE+1))
    this%offset(1) = 0
    do pe = 1, nPE
      this%offset(pe+1) = this%offset(pe) + bsizes(pe)
    end do

    call add_offP_index(this, offP_index)

    sync all

  end subroutine init

  subroutine add_offP_index(this, offP_index)

    class(index_map), intent(inout) :: this
    integer, intent(in) :: offP_index(:)

    integer :: j, p

    this%offP_index = offP_index
    !TODO: ensure offP_index is strictly increasing

    this%offP_size  = size(this%offP_index)
    this%local_size = this%onP_size + this%offP_size

    allocate(this%src_pe(this%offP_size), this%src_id(this%offP_size))

    p = 1
    do j = 1, this%offP_size
      do while (this%offP_index(j) > this%offset(p+1))
        p = p + 1
      end do
      if (p == this_image()) stop 1 !TODO: replace by assertion
      this%src_pe(j) = p
      this%src_id(j) = this%offP_index(j) - this%offset(p)
    end do

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
      integer, pointer :: data(:)
    end type
    type(box), allocatable :: src[:]
    allocate(src[*])

    src%data => onP_data
    sync all
    do j = 1, this%offP_size
      offP_data(j) = src[this%src_pe(j)]%data(this%src_id(j))
    end do
    sync all

  end subroutine gather_aux

end module index_map_type
