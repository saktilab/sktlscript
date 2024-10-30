program print_charge
implicit none
    integer :: NAtom,iStep, iAtom, NStep, serial
    character(len=100), allocatable :: dummy(:)
    double precision, allocatable :: c(:)
    double precision :: dt
    character(len=100) :: nac
integer :: ppos
character(len=4)  :: new_ext="dat"

ppos = scan(trim(nac),".", BACK= .true.)
if ( ppos > 0 ) nac = nac(1:ppos)//new_ext


read(5,*) nac
read(5,*) NAtom
read(5,*) NStep
read(5,*) dt
read(5,*) serial

open(10, file=nac, action="read", status="old")
open(20, file="time_charge.dat", action="write")

allocate(c(NAtom))
allocate(dummy(NAtom))

do iStep = 1, NStep
    read(10,*) 
    read(10,*)
   do iAtom = 1, NAtom
      read(10,*) dummy(iAtom), c(iAtom)
if (iAtom == serial) write(20,*) iStep*dt/100, c(iAtom), iAtom
   end do
end do



deallocate(c)
deallocate(dummy)


end program
