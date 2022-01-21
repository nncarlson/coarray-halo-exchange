program main

  use,intrinsic :: iso_fortran_env, only: int64, error_unit
  use index_map_type
  implicit none

  integer :: num_arg, repeat, this_PE, nPE, lun, bsize, noffP, j
  integer, allocatable :: offP_index(:), array(:)
  type(index_map) :: imap
  character(63) :: prog, datadir, filename, arg
  integer(int64) :: t1, t2, rate
  
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
  
  this_PE = this_image()
  nPE = num_images()

  write(filename,'(a,i3.3)') trim(datadir) // '/data', this_PE
  open(newunit=lun,file=filename,access='stream',form='unformatted',action='read')

  read(lun) bsize, noffP
  allocate(offP_index(noffP))
  read(lun) offP_index
  close(lun)

  call imap%init(bsize, offP_index)
  
  print *, this_PE, ':', imap%onP_size, imap%offP_size
  
  if (this_PE == 1) then
    write(*,'(a)') 'Timing gather of off-process data elements'
    write(*,'(i0,a,i0,a)') imap%global_size, ' elements distributed across ', nPE, ' processes'
  end if
  
  !! Fill local on-process array elements with their GID; off-process elements with -1
  allocate(array(imap%local_size))
  array = -1
  do j = 1, imap%onP_size
    array(j) = imap%global_index(j)
  end do

  !! Gather the off-process values
  sync all
  if (this_PE == 1) call system_clock(t1)
  do j = 1, repeat
    call imap%gather(array)
  end do
  if (this_PE == 1) then
    call system_clock(t2, rate)
    print *, 'Wall time:', real(t2-t1)/real(repeat*rate)
  end if

  !! Verify the gathered off-process values
  if (any(array(1+imap%onP_size:) /= imap%offP_index)) error stop

end program
