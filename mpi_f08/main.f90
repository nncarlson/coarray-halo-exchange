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

program main

  use,intrinsic :: iso_fortran_env, only: int64, error_unit
  use index_map_type
  use mpi_f08
  implicit none

  logical :: is_IOP
  integer :: num_arg, repeat, this_PE, nPE, lun, bsize, noffP, j, n
  integer, allocatable :: offP_index(:), array(:)
  type(index_map) :: imap
  character(63) :: prog, datadir, filename, arg
  integer(int64) :: t1, t2, rate

  call MPI_Init()
  call MPI_Comm_size(MPI_COMM_WORLD, nPE)
  call MPI_Comm_rank(MPI_COMM_WORLD, this_PE)
  is_IOP = (this_PE == 0)

  call get_command_argument(0, prog)
  num_arg = command_argument_count()

  if (num_arg < 1) then
    write(error_unit,'(a)') 'Usage: ' // trim(prog) // ' DATADIR [NUM-REPEATS]'
    error stop
  end if

  call get_command_argument(1, datadir)

  if (num_arg > 1) then
    call get_command_argument(2, arg)
    read(arg,*) repeat
  else
    repeat = 1
  end if

  write(filename,'(a,i3.3)') trim(datadir) // '/data', this_PE+1
  open(newunit=lun,file=filename,access='stream',form='unformatted',action='read')

  read(lun) bsize, noffP
  allocate(offP_index(noffP))
  read(lun) offP_index
  close(lun)

  call imap%init(MPI_COMM_WORLD, bsize, offP_index)

  call MPI_Allreduce(imap%offP_size, n, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD)
  if (is_IOP) then
    write(*,'(a,i0,a)') 'Timing gather of ', n, ' off-process data elements'
    write(*,'(i0,a,i0,a)') imap%global_size, ' elements distributed across ', nPE, ' processes'
  end if

  !! Fill local on-process array elements with their GID; off-process elements with -1
  allocate(array(imap%local_size))
  array = -1
  do j = 1, imap%onP_size
    array(j) = imap%global_index(j)
  end do

  !! Gather the off-process values
  if (is_IOP) call system_clock(t1)
  do j = 1, repeat
    call imap%gather(array)
  end do
  if (is_IOP) then
    call system_clock(t2, rate)
    write(*,'(a,g0,a)') 'Wall time: ', real(t2-t1)/real(repeat*rate), ' sec'
  end if

  !! Verify the gathered off-process values
  if (any(array(1+imap%onP_size:) /= imap%offP_index)) error stop

  call MPI_Finalize()

end program
