program centerofmass
implicit none

integer :: iStep, NStep, NAtom, iAtom 
character(len=10) :: traject, label
double precision, allocatable :: rx(:), ry(:), rz(:)
double precision, parameter :: mH=1.00794
double precision, parameter :: mC=12.0107
double precision, parameter :: mO=15.9990
double precision, parameter :: mLi=15.9990


read(5,*) traject
read(5,*) NStep
read(5,*) NAtom

allocate(rx(NAtom))
allocate(ry(NAtom))
allocate(rz(NAtom))


open(10, file=traject)

do iStep=1,NStep
read(10,*)
read(10,*)
do iAtom=1, NAtom 
   read(10,*) label, rx(iAtom), ry(iAtom), rz(iAtom) 

end do

end do

deallocate(rx)
deallocate(ry)
deallocate(rz)

end program centerofmass