program coordnum
implicit none
character(len=20) :: coordfile
double precision, parameter :: ONE = 1.0d+00
integer :: iStep, nStep, nAtom, iAtom, ind1, ind2, ind3, ind4, ind5, ind6
character(len=2), allocatable :: label(:)
double precision, allocatable :: rx(:), ry(:), rz(:)
double precision, allocatable :: rx14(:), rx15(:), rx16(:), ry14(:), ry15(:), ry16(:)
double precision, allocatable :: rz14(:), rz15(:), rz16(:)
double precision, allocatable :: rx24(:), rx25(:), rx26(:), ry24(:), ry25(:), ry26(:)
double precision, allocatable :: rz24(:), rz25(:), rz26(:)
double precision, allocatable :: rx34(:), rx35(:), rx36(:), ry34(:), ry35(:), ry36(:)
double precision, allocatable :: rz34(:), rz35(:), rz36(:)
double precision, allocatable :: r14(:), r15(:), r16(:), rat14(:), rat15(:), rat16(:), n(:)
double precision, allocatable :: r24(:), r25(:), r26(:), rat24(:), rat25(:), rat26(:)
double precision, allocatable :: r34(:), r35(:), r36(:), rat34(:), rat35(:), rat36(:)
double precision :: Lx, Ly, Lz, L(3), Linv(3)
double precision, parameter :: r0 = 2.50
integer :: ppos
character(len=4)  :: new_ext="traj"



read(5,*) coordfile
read(5,*) ind1, ind2, ind3, ind4, ind5, ind6
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
allocate(rx14(nStep))
allocate(rx15(nStep))
allocate(rx16(nStep))
allocate(rx24(nStep))
allocate(rx25(nStep))
allocate(rx26(nStep))
allocate(rx34(nStep))
allocate(rx35(nStep))
allocate(rx36(nStep))
allocate(ry14(nStep))
allocate(ry15(nStep))
allocate(ry16(nStep))
allocate(ry24(nStep))
allocate(ry25(nStep))
allocate(ry26(nStep))
allocate(ry34(nStep))
allocate(ry35(nStep))
allocate(ry36(nStep))
allocate(rz14(nStep))
allocate(rz15(nStep))
allocate(rz16(nStep))
allocate(rz24(nStep))
allocate(rz25(nStep))
allocate(rz26(nStep))
allocate(rz34(nStep))
allocate(rz35(nStep))
allocate(rz36(nStep))
allocate(r14(nStep))
allocate(r15(nStep))
allocate(r16(nStep))
allocate(r24(nStep))
allocate(r25(nStep))
allocate(r26(nStep))
allocate(r34(nStep))
allocate(r35(nStep))
allocate(r36(nStep))
allocate(rat14(nStep))
allocate(rat15(nStep))
allocate(rat16(nStep))
allocate(rat24(nStep))
allocate(rat25(nStep))
allocate(rat26(nStep))
allocate(rat34(nStep))
allocate(rat35(nStep))
allocate(rat36(nStep))

allocate(n(nStep))

do iStep = 1, nStep
    read(10,*) 
    read(10,*)
    do iAtom = 1, nAtom
        read(10,*) label(iAtom), rx(iAtom), ry(iAtom), rz(iAtom)
        rx14(iStep) = rx(ind1)-rx(ind4)
        rx15(iStep) = rx(ind1)-rx(ind5)
        rx16(iStep) = rx(ind1)-rx(ind6)
        
        rx24(iStep) = rx(ind2)-rx(ind4)
        rx25(iStep) = rx(ind2)-rx(ind5)
        rx26(iStep) = rx(ind2)-rx(ind6)
        
        rx34(iStep) = rx(ind3)-rx(ind4)
        rx35(iStep) = rx(ind3)-rx(ind5)
        rx36(iStep) = rx(ind3)-rx(ind6)
        
        ry14(iStep) = ry(ind1)-ry(ind4)
        ry15(iStep) = ry(ind1)-ry(ind5)
        ry16(iStep) = ry(ind1)-ry(ind6)
        
        ry24(iStep) = ry(ind2)-ry(ind4)
        ry25(iStep) = ry(ind2)-ry(ind5)
        ry26(iStep) = ry(ind2)-ry(ind6)
        
        ry34(iStep) = ry(ind3)-ry(ind4)
        ry35(iStep) = ry(ind3)-ry(ind5)
        ry36(iStep) = ry(ind3)-ry(ind6)
        
        rz14(iStep) = rz(ind1)-rz(ind4)
        rz15(iStep) = rz(ind1)-rz(ind5)
        rz16(iStep) = rz(ind1)-rz(ind6)
        
        rz24(iStep) = rz(ind2)-rz(ind4)
        rz25(iStep) = rz(ind2)-rz(ind5)
        rz26(iStep) = rz(ind2)-rz(ind6)
        
        rz34(iStep) = rz(ind3)-rz(ind4)
        rz35(iStep) = rz(ind3)-rz(ind5)
        rz36(iStep) = rz(ind3)-rz(ind6)
        

        rx14(iStep) = rx14(iStep) - dnint(rx14(iStep)*Linv(1))*L(1)
        rx15(iStep) = rx15(iStep) - dnint(rx15(iStep)*Linv(1))*L(1)
        rx16(iStep) = rx16(iStep) - dnint(rx16(iStep)*Linv(1))*L(1)
        rx24(iStep) = rx24(iStep) - dnint(rx24(iStep)*Linv(1))*L(1)
        rx25(iStep) = rx25(iStep) - dnint(rx25(iStep)*Linv(1))*L(1)
        rx26(iStep) = rx26(iStep) - dnint(rx26(iStep)*Linv(1))*L(1)
        rx34(iStep) = rx34(iStep) - dnint(rx34(iStep)*Linv(1))*L(1)
        rx35(iStep) = rx35(iStep) - dnint(rx35(iStep)*Linv(1))*L(1)
        rx36(iStep) = rx36(iStep) - dnint(rx36(iStep)*Linv(1))*L(1)
        
        ry14(iStep) = ry14(iStep) - dnint(ry14(iStep)*Linv(2))*L(2)
        ry15(iStep) = ry15(iStep) - dnint(ry15(iStep)*Linv(2))*L(2)
        ry16(iStep) = ry16(iStep) - dnint(ry16(iStep)*Linv(2))*L(2)
        ry24(iStep) = ry24(iStep) - dnint(ry24(iStep)*Linv(2))*L(2)
        ry25(iStep) = ry25(iStep) - dnint(ry25(iStep)*Linv(2))*L(2)
        ry26(iStep) = ry26(iStep) - dnint(ry26(iStep)*Linv(2))*L(2)
        ry34(iStep) = ry34(iStep) - dnint(ry34(iStep)*Linv(2))*L(2)
        ry35(iStep) = ry35(iStep) - dnint(ry35(iStep)*Linv(2))*L(2)
        ry36(iStep) = ry36(iStep) - dnint(ry36(iStep)*Linv(2))*L(2)
        
        rz14(iStep) = rz14(iStep) - dnint(rz14(iStep)*Linv(3))*L(3)
        rz15(iStep) = rz15(iStep) - dnint(rz15(iStep)*Linv(3))*L(3)
        rz16(iStep) = rz16(iStep) - dnint(rz16(iStep)*Linv(3))*L(3)
        rz24(iStep) = rz24(iStep) - dnint(rz24(iStep)*Linv(3))*L(3)
        rz25(iStep) = rz25(iStep) - dnint(rz25(iStep)*Linv(3))*L(3)
        rz26(iStep) = rz26(iStep) - dnint(rz26(iStep)*Linv(3))*L(3)
        rz34(iStep) = rz34(iStep) - dnint(rz34(iStep)*Linv(3))*L(3)
        rz35(iStep) = rz35(iStep) - dnint(rz35(iStep)*Linv(3))*L(3)
        rz36(iStep) = rz36(iStep) - dnint(rz36(iStep)*Linv(3))*L(3)
        
        r14(iStep) = (rx14(iStep)**2 + ry14(iStep)**2 + rz14(iStep)**2)**0.5
        r15(iStep) = (rx15(iStep)**2 + ry15(iStep)**2 + rz15(iStep)**2)**0.5
        r16(iStep) = (rx16(iStep)**2 + ry16(iStep)**2 + rz16(iStep)**2)**0.5
        
        r24(iStep) = (rx24(iStep)**2 + ry24(iStep)**2 + rz24(iStep)**2)**0.5
        r25(iStep) = (rx25(iStep)**2 + ry25(iStep)**2 + rz25(iStep)**2)**0.5
        r26(iStep) = (rx26(iStep)**2 + ry26(iStep)**2 + rz26(iStep)**2)**0.5
        
        r34(iStep) = (rx34(iStep)**2 + ry34(iStep)**2 + rz34(iStep)**2)**0.5
        r35(iStep) = (rx35(iStep)**2 + ry35(iStep)**2 + rz35(iStep)**2)**0.5
        r36(iStep) = (rx36(iStep)**2 + ry36(iStep)**2 + rz36(iStep)**2)**0.5
        

        rat14(iStep) = (1-(r14(iStep)/r0)**6)/(1-(r14(iStep)/r0)**12)
        rat15(iStep) = (1-(r15(iStep)/r0)**6)/(1-(r15(iStep)/r0)**12)
        rat16(iStep) = (1-(r16(iStep)/r0)**6)/(1-(r16(iStep)/r0)**12)
        
        rat24(iStep) = (1-(r24(iStep)/r0)**6)/(1-(r24(iStep)/r0)**12)
        rat25(iStep) = (1-(r25(iStep)/r0)**6)/(1-(r25(iStep)/r0)**12)
        rat26(iStep) = (1-(r26(iStep)/r0)**6)/(1-(r26(iStep)/r0)**12)
        
        rat34(iStep) = (1-(r34(iStep)/r0)**6)/(1-(r34(iStep)/r0)**12)
        rat35(iStep) = (1-(r35(iStep)/r0)**6)/(1-(r35(iStep)/r0)**12)
        rat36(iStep) = (1-(r36(iStep)/r0)**6)/(1-(r36(iStep)/r0)**12)
        
        n(iStep) = ((rat14(iStep) + rat15(iStep) + rat16(iStep) + rat24(iStep) + rat25(iStep) &
        + rat26(iStep) + rat34(iStep) + rat35(iStep) + rat36(iStep)))/3.0
    end do

write(20, *) iStep, n(iStep)

end do



deallocate(label)
deallocate(rx)
deallocate(ry)
deallocate(rz)
deallocate(rx14)
deallocate(rx15)
deallocate(rx16)
deallocate(rx24)
deallocate(rx25)
deallocate(rx26)
deallocate(rx34)
deallocate(rx35)
deallocate(rx36)
deallocate(ry14)
deallocate(ry15)
deallocate(ry16)
deallocate(ry24)
deallocate(ry25)
deallocate(ry26)
deallocate(ry34)
deallocate(ry35)
deallocate(ry36)
deallocate(rz14)
deallocate(rz15)
deallocate(rz16)
deallocate(rz24)
deallocate(rz25)
deallocate(rz26)
deallocate(rz34)
deallocate(rz35)
deallocate(rz36)
deallocate(r14)
deallocate(r15)
deallocate(r16)
deallocate(r24)
deallocate(r25)
deallocate(r26)
deallocate(r34)
deallocate(r35)
deallocate(r36)
deallocate(rat14)
deallocate(rat15)
deallocate(rat16)
deallocate(rat24)
deallocate(rat25)
deallocate(rat26)
deallocate(rat34)
deallocate(rat35)
deallocate(rat36)

end program coordnum