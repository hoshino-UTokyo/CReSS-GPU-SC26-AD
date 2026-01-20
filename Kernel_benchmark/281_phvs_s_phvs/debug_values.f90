program debug_values
  implicit none
  integer :: ni, nj, nk
  integer :: exbopt, wbc, ebc, sbc, nbc, advopt, mpopt, mfcopt, ape
  integer :: ebe, ebn, ebs, ebw, isub, jsub, nisub, njsub
  real :: gdxdt, gdxdtn, gdydt, gdydtn
  character(len=108) :: exbvar
  real, allocatable :: scpx_ref(:,:,:), scpy_ref(:,:,:)
  integer :: unit_num
  character(len=256) :: line

  unit_num = 11
  open(unit_num, file='data/params.txt', status='old')
  read(unit_num, '(A)') line; exbvar = adjustl(line(index(line,'=')+1:))
  read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) exbopt
  read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) wbc
  read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) ebc
  read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) sbc
  read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) nbc
  read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) advopt
  read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) mpopt
  read(unit_num, '(A)') line; read(line(index(line,'=')+1:), *) mfcopt
  do unit_num = 1, 3
    read(11, '(A)') line
  end do
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) ape
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) ni
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) nj
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) nk
  do unit_num = 1, 4
    read(11, '(A)') line
  end do
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) ebe
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) ebn
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) ebs
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) ebw
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) gdxdt
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) gdxdtn
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) gdydt
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) gdydtn
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) isub
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) jsub
  do unit_num = 1, 4
    read(11, '(A)') line
  end do
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) nisub
  do unit_num = 1, 4
    read(11, '(A)') line
  end do
  read(11, '(A)') line; read(line(index(line,'=')+1:), *) njsub
  close(11)

  print *, 'ni, nj, nk = ', ni, nj, nk
  print *, 'gdxdt = ', gdxdt
  print *, 'gdxdtn = ', gdxdtn
  print *, 'gdydt = ', gdydt
  print *, 'gdydtn = ', gdydtn
  print *, 'wbc,ebc,sbc,nbc = ', wbc, ebc, sbc, nbc
  print *, 'ebw,ebe,ebs,ebn = ', ebw, ebe, ebs, ebn
  print *, 'isub,jsub = ', isub, jsub
  print *, 'nisub,njsub = ', nisub, njsub
  print *, 'ape = ', ape
  print *, 'exbvar(ape:ape) = "', exbvar(ape:ape), '"'

  allocate(scpx_ref(1:nj, 1:nk, 1:2))
  allocate(scpy_ref(1:ni, 1:nk, 1:2))

  open(12, file='data/scpx_ref.bin', status='old', access='stream', form='unformatted')
  read(12) scpx_ref
  close(12)
  open(12, file='data/scpy_ref.bin', status='old', access='stream', form='unformatted')
  read(12) scpy_ref
  close(12)

  print *, 'scpx_ref(1,2,1) = ', scpx_ref(1,2,1)
  print *, 'scpx_ref(1,2,2) = ', scpx_ref(1,2,2)
  print *, 'scpy_ref(1,2,1) = ', scpy_ref(1,2,1)
  print *, 'scpy_ref(1,2,2) = ', scpy_ref(1,2,2)
end program
