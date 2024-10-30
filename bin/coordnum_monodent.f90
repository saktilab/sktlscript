program coordnum
implicit none
character(len=20) :: coordfile
double precision, parameter :: ONE = 1.0d+00
integer :: iStep, nStep, nAtom, iAtom, ind1, ind2, ind3, ind4
character(len=2), allocatable :: label(:)
double precision, allocatable :: rx(:), ry(:), rz(:)
double precision, allocatable :: rx12(:), rx14(:), rx23(:), rx24(:), ry12(:), ry14(:), ry23(:), ry24(:)
double precision, allocatable :: rz12(:), rz14(:), rz23(:), rz24(:)
double precision, allocatable :: r12(:), r14(:), r23(:), r24(:), rat12(:), rat14(:), rat23(:), rat24(:), n(:)
double precision :: Lx, Ly, Lz, L(3), Linv(3)
double precision, parameter :: r0 = 2.50
integer :: ppos
character(len=4)  :: new_ext="traj"



read(5,*) coordfile
read(5,*) ind1, ind2
read(5,*) nAtom
read(5,*) nStep
read(5,*) Lx, Ly, Lz

ppos = scan(trim(coordfile),".", BACK= .true.)
if ( ppos > 0 ) coordfile = coordfile(1:ppos)//new_ext

L(1) = Lx
L(2) = Ly
L(3) = Lz
Linv(1:3) = ONE/L(1:3)

open(10, file=coordfile, action="read", status="old")
open(20, file="time_n.dat", action="write")
allocate(label(nAtom))
allocate(rx(nAtom))
allocate(ry(nAtom))
allocate(rz(nAtom))
allocate(rx12(nStep))
allocate(rx14(nStep))
allocate(rx23(nStep))
allocate(rx24(nStep))
allocate(ry12(nStep))
allocate(ry14(nStep))
allocate(ry23(nStep))
allocate(ry24(nStep))
allocate(rz12(nStep))
allocate(rz14(nStep))
allocate(rz23(nStep))
allocate(rz24(nStep))
allocate(r12(nStep))
allocate(r14(nStep))
allocate(r23(nStep))
allocate(r24(nStep))
allocate(rat12(nStep))
allocate(rat14(nStep))
allocate(rat23(nStep))
allocate(rat24(nStep))
allocate(n(nStep))

do iStep = 1, nStep
    read(10,*) 
    read(10,*)
    do iAtom = 1, nAtom
        read(10,*) label(iAtom), rx(iAtom), ry(iAtom), rz(iAtom)
        rx12(iStep) = rx(ind1)-rx(ind2)
        ry12(iStep) = ry(ind1)-ry(ind2)
        rz12(iStep) = rz(ind1)-rz(ind2)
        
        rx12(iStep) = rx12(iStep) - dnint(rx12(iStep)*Linv(1))*L(1)
        
        ry12(iStep) = ry12(iStep) - dnint(ry12(iStep)*Linv(2))*L(2)
        
        rz12(iStep) = rz12(iStep) - dnint(rz12(iStep)*Linv(3))*L(3)
        
        r12(iStep) = (rx12(iStep)**2 + ry12(iStep)**2 + rz12(iStep)**2)**0.5
        
        rat12(iStep) = (1-(r12(iStep)/r0)**6)/(1-(r12(iStep)/r0)**12)
        
        n(iStep) = rat12(iStep)
    end do

write(20, *) iStep, n(iStep)

end do



deallocate(label)
deallocate(rx)
deallocate(ry)
deallocate(rz)
deallocate(rx12)
deallocate(rx14)
deallocate(rx23)
deallocate(rx24)
deallocate(ry12)
deallocate(ry14)
deallocate(ry23)
deallocate(ry24)
deallocate(rz12)
deallocate(rz14)
deallocate(rz23)
deallocate(rz24)
deallocate(r12)
deallocate(r14)
deallocate(r23)
deallocate(r24)
deallocate(rat12)
deallocate(rat14)
deallocate(rat23)
deallocate(rat24)
end program coordnum