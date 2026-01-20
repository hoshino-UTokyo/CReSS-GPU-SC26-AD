program test_params
  implicit none
  character(len=108) :: exbvar
  integer :: ape, unit_num
  character(len=256) :: line
  
  unit_num = 11
  open(unit_num, file='data/params.txt', status='old')
  read(unit_num, '(A)') line
  exbvar = adjustl(line(index(line,'=')+1:))
  close(unit_num)
  
  ape = 9
  print *, 'exbvar = "', trim(exbvar), '"'
  print *, 'len = ', len_trim(exbvar)
  print *, 'ape = ', ape
  print *, 'exbvar(ape:ape) = "', exbvar(ape:ape), '"'
  print *, 'check x: ', exbvar(ape:ape) == 'x'
end program
