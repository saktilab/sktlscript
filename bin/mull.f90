program mul2charge
implicit none
character(len=100) :: mulliken, label
integer :: iAtom, nAtom, nOrbs, iStep,iOrbs, nStep
double precision, allocatable :: c(:), mul(:)
integer :: npt, nbr, ns, no, nn, nc, nh

read(5,*) mulliken 
read(5,*) nStep
read(5,*) nOrbs
read(5,*) npt
read(5,*) nbr
read(5,*) ns
read(5,*) no
read(5,*) nn
read(5,*) nc
read(5,*) nh

open (11,file=mulliken, action='read', status='old')
allocate(c(nOrbs))
allocate(mul(nOrbs))

do iStep = 1, nStep
  read(11,*)
  read(11,*)
  do iOrbs=1,nOrbs
     read(11,*) label, mul(iOrbs)
  end do
  write(6,*)"STEP ", iStep

!!! Calculate the charge for each Pt atom
do iOrbs=1, (3*npt+1)-3, 3
   c(iOrbs)=2-(mul(iOrbs)+mul(iOrbs+1)+mul(iOrbs+2))
  ! write(6,*) "Pt", c(iOrbs) 
end do

!!!Calculate the charge for each Br atom
do iOrbs=(3*npt+1), (3*nbr+3*npt+1)-3, 3
   c(iOrbs)=7-(mul(iOrbs)+mul(iOrbs+1)+mul(iOrbs+2))
  ! write(6,*) "Br", c(iOrbs)
end do

!!!Calculate the charge for each S atom
do iOrbs=(3*nbr+3*npt+1), (3*ns+3*nbr+3*npt+1)-3, 3
   c(iOrbs)=6-(mul(iOrbs)+mul(iOrbs+1)+mul(iOrbs+2))
  ! write(6,*) "S", c(iOrbs)
end do

!!!Calculate the charge for each O atom
do iOrbs=(3*ns+3*nbr+3*npt+1), (2*no+3*ns+3*nbr+3*npt+1)-2, 2
   c(iOrbs)=6-(mul(iOrbs)+mul(iOrbs+1))
   write(6,*) "O", c(iOrbs)
end do

!!!Calculate the charge for each N atom
do iOrbs=(2*no+3*ns+3*nbr+3*npt+1),(2*nn+2*no+3*ns+3*nbr+3*npt+1)-2, 2
   c(iOrbs)=5-(mul(iOrbs)+mul(iOrbs+1))
  ! write(6,*) "N", c(iOrbs)
end do

!!!Calculate the charge for each C atom
do iOrbs=(2*nn+2*no+3*ns+3*nbr+3*npt+1),(2*nc+2*nn+2*no+3*ns+3*nbr+3*npt+1)-2,2
   c(iOrbs)=4-(mul(iOrbs)+mul(iOrbs+1))
  ! write(6,*) "C", c(iOrbs)
end do

!!!Calculate the charge for each H atom
do iOrbs=(2*nc+2*nn+2*no+3*ns+3*nbr+3*npt+1), (nh+2*nc+2*nn+2*no+3*ns+3*nbr+3*npt+1)-1,1
  c(iOrbs)=1-mul(iOrbs)
 ! write(6,*) "H", c(iOrbs)
end do
end do

deallocate(c)
deallocate(mul)



end program mul2charge
