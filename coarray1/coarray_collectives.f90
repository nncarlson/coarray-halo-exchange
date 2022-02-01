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

module coarray_collectives

  implicit none
  private

  public :: co_sum_scan

  interface co_sum_scan
    procedure :: co_sum_scan1, co_sum_scan2
  end interface

contains

  subroutine co_sum_scan1(x)
    integer, intent(inout) :: x
    integer, allocatable :: c[:]
    allocate(c[*])
    c = x
    call sum_scan_co(c)
    x = c
  end subroutine

  subroutine co_sum_scan2(x, y)
    integer, intent(in)  :: x
    integer, intent(out) :: y
    integer, allocatable :: c[:]
    allocate(c[*])
    c = x
    call sum_scan_co(c)
    y = c
  end subroutine

  !! This version of a prefix sum implements a "work-efficient" algorithm [1,2]
  !! that requires 2 log2 N - 2 (sequential) steps, but only twice the work of
  !! the sequential algorithm. Here "work" is the summing of two numbers, one
  !! of which is read from another image. Thus this parallel algorithm is very
  !! efficient compared to other parallel algorithms with respect to the amount
  !! of remote communication. This algorithm also presents a very interesting
  !! exercise/challenge in the use of SYNC IMAGES.
  !!
  !! [1] Sendgupta, Lefohn, and Owens, "A work-efficient step-efficient
  !!     prefix-sum algorithm," 2006.
  !! [2] https://en.wikipedia.org/wiki/Prefix_sum

  subroutine sum_scan_co(x)

    integer, intent(inout) :: x[*]

    integer :: me, m

    me = this_image()
    m = 1
    do while (2*m <= num_images())
      if (modulo(me,2*m) == 0) then ! replaceable by some faster bit-level operation?
        sync images (me-m)
        x = x + x[me-m]
      else
        if (me+m <= num_images()) sync images (me+m)
        exit
      end if
      m = 2*m
    end do

    if (me > m) then
      sync images (me-m)
      x = x + x[me-m]
    end if

    do while (m > 1)
      m = m/2
      if (me+m <= num_images()) sync images (me+m)
    end do

  end subroutine sum_scan_co

end module coarray_collectives
