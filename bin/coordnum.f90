program coordnum
implicit none
character(len=20) :: coordfile
double precision, parameter :: ONE = 1.0d+00
integer :: iStep, nStep, nAtom, iAtom, ind1, ind2, ind3, ind4
character(len=2), allocatable :: label(:)
double precision, allocatable :: rx(:), ry(:), rz(:)
double precision, allocatable :: rx13(:), rx14(:), rx23(:), rx24(:), ry13(:), ry14(:), ry23(:), ry24(:)
double precision, allocatable :: rz13(:), rz14(:), rz23(:), rz24(:)
double precision, allocatable :: r13(:), r14(:), r23(:), r24(:), rat13(:), rat14(:), rat23(:), rat24(:), n(:)
double precision :: Lx, Ly, Lz, L(3), Linv(3)
double precision, parameter :: r0 = 2.50
integer :: ppos
character(len=4)  :: new_ext="traj"



read(5,*) coordfile
read(5,*) ind1, ind2, ind3, ind4
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
allocate(rx13(nStep))
allocate(rx14(nStep))
allocate(rx23(nStep))
allocate(rx24(nStep))
allocate(ry13(nStep))
allocate(ry14(nStep))
allocate(ry23(nStep))
allocate(ry24(nStep))
allocate(rz13(nStep))
allocate(rz14(nStep))
allocate(rz23(nStep))
allocate(rz24(nStep))
allocate(r13(nStep))
allocate(r14(nStep))
allocate(r23(nStep))
allocate(r24(nStep))
allocate(rat13(nStep))
allocate(rat14(nStep))
allocate(rat23(nStep))
allocate(rat24(nStep))
allocate(n(nStep))

do iStep = 1, nStep
    read(10,*) 
    read(10,*)
    do iAtom = 1, nAtom
        read(10,*) label(iAtom), rx(iAtom), ry(iAtom), rz(iAtom)
        rx13(iStep) = rx(ind1)-rx(ind3)
        rx14(iStep) = rx(ind1)-rx(ind4)
        rx23(iStep) = rx(ind2)-rx(ind3)
        rx24(iStep) = rx(ind2)-rx(ind4)
        ry13(iStep) = ry(ind1)-ry(ind3)
        ry14(iStep) = ry(ind1)-ry(ind4)
        ry23(iStep) = ry(ind2)-ry(ind3)
        ry24(iStep) = ry(ind2)-ry(ind4)
        rz13(iStep) = rz(ind1)-rz(ind3)
        rz14(iStep) = rz(ind1)-rz(ind4)
        rz23(iStep) = rz(ind2)-rz(ind3)
        rz24(iStep) = rz(ind2)-rz(ind4)

        rx13(iStep) = rx13(iStep) - dnint(rx13(iStep)*Linv(1))*L(1)
        rx14(iStep) = rx14(iStep) - dnint(rx14(iStep)*Linv(1))*L(1)
        rx23(iStep) = rx23(iStep) - dnint(rx23(iStep)*Linv(1))*L(1)
        rx24(iStep) = rx24(iStep) - dnint(rx24(iStep)*Linv(1))*L(1)
        ry13(iStep) = ry13(iStep) - dnint(ry13(iStep)*Linv(2))*L(2)
        ry14(iStep) = ry14(iStep) - dnint(ry14(iStep)*Linv(2))*L(2)
        ry23(iStep) = ry23(iStep) - dnint(ry23(iStep)*Linv(2))*L(2)
        ry24(iStep) = ry24(iStep) - dnint(ry24(iStep)*Linv(2))*L(2)
        rz13(iStep) = rz13(iStep) - dnint(rz13(iStep)*Linv(3))*L(3)
        rz14(iStep) = rz14(iStep) - dnint(rz14(iStep)*Linv(3))*L(3)
        rz23(iStep) = rz23(iStep) - dnint(rz23(iStep)*Linv(3))*L(3)
        rz24(iStep) = rz24(iStep) - dnint(rz24(iStep)*Linv(3))*L(3)
        
        r13(iStep) = (rx13(iStep)**2 + ry13(iStep)**2 + rz13(iStep)**2)**0.5
        r14(iStep) = (rx14(iStep)**2 + ry14(iStep)**2 + rz14(iStep)**2)**0.5
        r23(iStep) = (rx23(iStep)**2 + ry23(iStep)**2 + rz23(iStep)**2)**0.5
        r24(iStep) = (rx24(iStep)**2 + ry24(iStep)**2 + rz24(iStep)**2)**0.5
        
        rat13(iStep) = (1-(r13(iStep)/r0)**6)/(1-(r13(iStep)/r0)**12)
        rat14(iStep) = (1-(r14(iStep)/r0)**6)/(1-(r14(iStep)/r0)**12)
        rat23(iStep) = (1-(r23(iStep)/r0)**6)/(1-(r23(iStep)/r0)**12)
        rat24(iStep) = (1-(r24(iStep)/r0)**6)/(1-(r24(iStep)/r0)**12)
        n(iStep) = 0.5*(rat13(iStep) + rat14(iStep) + rat23(iStep) + rat24(iStep))
    end do

write(20, *) iStep, n(iStep)

end do



deallocate(label)
deallocate(rx)
deallocate(ry)
deallocate(rz)
deallocate(rx13)
deallocate(rx14)
deallocate(rx23)
deallocate(rx24)
deallocate(ry13)
deallocate(ry14)
deallocate(ry23)
deallocate(ry24)
deallocate(rz13)
deallocate(rz14)
deallocate(rz23)
deallocate(rz24)
deallocate(r13)
deallocate(r14)
deallocate(r23)
deallocate(r24)
deallocate(rat13)
deallocate(rat14)
deallocate(rat23)
deallocate(rat24)
end program coordnum