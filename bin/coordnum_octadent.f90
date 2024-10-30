program coordnum
implicit none
character(len=20) :: coordfile
double precision, parameter :: ONE = 1.0d+00
integer :: iStep, nStep, nAtom, iAtom, ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8
character(len=2), allocatable :: label(:)
double precision, allocatable :: rx(:), ry(:), rz(:)
double precision, allocatable :: rx15(:), rx16(:), rx17(:), rx18(:), ry15(:), ry16(:), ry17(:), ry18(:)
double precision, allocatable :: rz15(:), rz16(:), rz17(:), rz18(:)
double precision, allocatable :: rx25(:), rx26(:), rx27(:), rx28(:), ry25(:), ry26(:), ry27(:), ry28(:)
double precision, allocatable :: rz25(:), rz26(:), rz27(:), rz28(:)
double precision, allocatable :: rx35(:), rx36(:), rx37(:), rx38(:), ry35(:), ry36(:), ry37(:), ry38(:)
double precision, allocatable :: rz35(:), rz36(:), rz37(:), rz38(:)
double precision, allocatable :: rx45(:), rx46(:), rx47(:), rx48(:), ry45(:), ry46(:), ry47(:), ry48(:)
double precision, allocatable :: rz45(:), rz46(:), rz47(:), rz48(:)
double precision, allocatable :: r15(:), r16(:), r17(:), r18(:), rat15(:), rat16(:), rat17(:), rat18(:), n(:)
double precision, allocatable :: r25(:), r26(:), r27(:), r28(:), rat25(:), rat26(:), rat27(:), rat28(:)
double precision, allocatable :: r35(:), r36(:), r37(:), r38(:), rat35(:), rat36(:), rat37(:), rat38(:)
double precision, allocatable :: r45(:), r46(:), r47(:), r48(:), rat45(:), rat46(:), rat47(:), rat48(:)
double precision :: Lx, Ly, Lz, L(3), Linv(3)
double precision, parameter :: r0 = 2.50
integer :: ppos
character(len=4)  :: new_ext="traj"



read(5,*) coordfile
read(5,*) ind1, ind2, ind3, ind4, ind5, ind6, ind7, ind8
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
allocate(rx15(nStep))
allocate(rx16(nStep))
allocate(rx17(nStep))
allocate(rx18(nStep))
allocate(rx25(nStep))
allocate(rx26(nStep))
allocate(rx27(nStep))
allocate(rx28(nStep))
allocate(rx35(nStep))
allocate(rx36(nStep))
allocate(rx37(nStep))
allocate(rx38(nStep))
allocate(rx45(nStep))
allocate(rx46(nStep))
allocate(rx47(nStep))
allocate(rx48(nStep))
allocate(ry15(nStep))
allocate(ry16(nStep))
allocate(ry17(nStep))
allocate(ry18(nStep))
allocate(ry25(nStep))
allocate(ry26(nStep))
allocate(ry27(nStep))
allocate(ry28(nStep))
allocate(ry35(nStep))
allocate(ry36(nStep))
allocate(ry37(nStep))
allocate(ry38(nStep))
allocate(ry45(nStep))
allocate(ry46(nStep))
allocate(ry47(nStep))
allocate(ry48(nStep))
allocate(rz15(nStep))
allocate(rz16(nStep))
allocate(rz17(nStep))
allocate(rz18(nStep))
allocate(rz25(nStep))
allocate(rz26(nStep))
allocate(rz27(nStep))
allocate(rz28(nStep))
allocate(rz35(nStep))
allocate(rz36(nStep))
allocate(rz37(nStep))
allocate(rz38(nStep))
allocate(rz45(nStep))
allocate(rz46(nStep))
allocate(rz47(nStep))
allocate(rz48(nStep))
allocate(r15(nStep))
allocate(r16(nStep))
allocate(r17(nStep))
allocate(r18(nStep))
allocate(r25(nStep))
allocate(r26(nStep))
allocate(r27(nStep))
allocate(r28(nStep))
allocate(r35(nStep))
allocate(r36(nStep))
allocate(r37(nStep))
allocate(r38(nStep))
allocate(r45(nStep))
allocate(r46(nStep))
allocate(r47(nStep))
allocate(r48(nStep))
allocate(rat15(nStep))
allocate(rat16(nStep))
allocate(rat17(nStep))
allocate(rat18(nStep))
allocate(rat25(nStep))
allocate(rat26(nStep))
allocate(rat27(nStep))
allocate(rat28(nStep))
allocate(rat35(nStep))
allocate(rat36(nStep))
allocate(rat37(nStep))
allocate(rat38(nStep))
allocate(rat45(nStep))
allocate(rat46(nStep))
allocate(rat47(nStep))
allocate(rat48(nStep))
allocate(n(nStep))

do iStep = 1, nStep
    read(10,*) 
    read(10,*)
    do iAtom = 1, nAtom
        read(10,*) label(iAtom), rx(iAtom), ry(iAtom), rz(iAtom)
        rx15(iStep) = rx(ind1)-rx(ind5)
        rx16(iStep) = rx(ind1)-rx(ind6)
        rx17(iStep) = rx(ind1)-rx(ind7)
        rx18(iStep) = rx(ind1)-rx(ind8)
 
        rx25(iStep) = rx(ind2)-rx(ind5)
        rx26(iStep) = rx(ind2)-rx(ind6)
        rx27(iStep) = rx(ind2)-rx(ind7)
        rx28(iStep) = rx(ind2)-rx(ind8)
 
        rx35(iStep) = rx(ind3)-rx(ind5)
        rx36(iStep) = rx(ind3)-rx(ind6)
        rx37(iStep) = rx(ind3)-rx(ind7)
        rx38(iStep) = rx(ind3)-rx(ind8)
        
        rx45(iStep) = rx(ind4)-rx(ind5)
        rx46(iStep) = rx(ind4)-rx(ind6)
        rx47(iStep) = rx(ind4)-rx(ind7)
        rx48(iStep) = rx(ind4)-rx(ind8)
 
        ry15(iStep) = ry(ind1)-ry(ind5)
        ry16(iStep) = ry(ind1)-ry(ind6)
        ry17(iStep) = ry(ind1)-ry(ind7)
        ry18(iStep) = ry(ind1)-ry(ind8)
 
        ry25(iStep) = ry(ind2)-ry(ind5)
        ry26(iStep) = ry(ind2)-ry(ind6)
        ry27(iStep) = ry(ind2)-ry(ind7)
        ry28(iStep) = ry(ind2)-ry(ind8)
 
        ry35(iStep) = ry(ind3)-ry(ind5)
        ry36(iStep) = ry(ind3)-ry(ind6)
        ry37(iStep) = ry(ind3)-ry(ind7)
        ry38(iStep) = ry(ind3)-ry(ind8)
        
        ry45(iStep) = ry(ind4)-ry(ind5)
        ry46(iStep) = ry(ind4)-ry(ind6)
        ry47(iStep) = ry(ind4)-ry(ind7)
        ry48(iStep) = ry(ind4)-ry(ind8)

        rz15(iStep) = rz(ind1)-rz(ind5)
        rz16(iStep) = rz(ind1)-rz(ind6)
        rz17(iStep) = rz(ind1)-rz(ind7)
        rz18(iStep) = rz(ind1)-rz(ind8)
 
        rz25(iStep) = rz(ind2)-rz(ind5)
        rz26(iStep) = rz(ind2)-rz(ind6)
        rz27(iStep) = rz(ind2)-rz(ind7)
        rz28(iStep) = rz(ind2)-rz(ind8)
 
        rz35(iStep) = rz(ind3)-rz(ind5)
        rz36(iStep) = rz(ind3)-rz(ind6)
        rz37(iStep) = rz(ind3)-rz(ind7)
        rz38(iStep) = rz(ind3)-rz(ind8)
        
        rz45(iStep) = rz(ind4)-rz(ind5)
        rz46(iStep) = rz(ind4)-rz(ind6)
        rz47(iStep) = rz(ind4)-rz(ind7)
        rz48(iStep) = rz(ind4)-rz(ind8)

        rx15(iStep) = rx15(iStep) - dnint(rx15(iStep)*Linv(1))*L(1)
        rx16(iStep) = rx16(iStep) - dnint(rx16(iStep)*Linv(1))*L(1)
        rx17(iStep) = rx17(iStep) - dnint(rx17(iStep)*Linv(1))*L(1)
        rx18(iStep) = rx18(iStep) - dnint(rx18(iStep)*Linv(1))*L(1)
        
        rx25(iStep) = rx25(iStep) - dnint(rx25(iStep)*Linv(1))*L(1)
        rx26(iStep) = rx26(iStep) - dnint(rx26(iStep)*Linv(1))*L(1)
        rx27(iStep) = rx27(iStep) - dnint(rx27(iStep)*Linv(1))*L(1)
        rx28(iStep) = rx28(iStep) - dnint(rx28(iStep)*Linv(1))*L(1)
        
        rx35(iStep) = rx35(iStep) - dnint(rx35(iStep)*Linv(1))*L(1)
        rx36(iStep) = rx36(iStep) - dnint(rx36(iStep)*Linv(1))*L(1)
        rx37(iStep) = rx37(iStep) - dnint(rx37(iStep)*Linv(1))*L(1)
        rx38(iStep) = rx38(iStep) - dnint(rx38(iStep)*Linv(1))*L(1)
        
        rx45(iStep) = rx45(iStep) - dnint(rx45(iStep)*Linv(1))*L(1)
        rx46(iStep) = rx46(iStep) - dnint(rx46(iStep)*Linv(1))*L(1)
        rx47(iStep) = rx47(iStep) - dnint(rx47(iStep)*Linv(1))*L(1)
        rx48(iStep) = rx48(iStep) - dnint(rx48(iStep)*Linv(1))*L(1)
        

        ry15(iStep) = ry15(iStep) - dnint(ry15(iStep)*Linv(2))*L(2)
        ry16(iStep) = ry16(iStep) - dnint(ry16(iStep)*Linv(2))*L(2)
        ry17(iStep) = ry17(iStep) - dnint(ry17(iStep)*Linv(2))*L(2)
        ry18(iStep) = ry18(iStep) - dnint(ry18(iStep)*Linv(2))*L(2)
        
        ry25(iStep) = ry25(iStep) - dnint(ry25(iStep)*Linv(2))*L(2)
        ry26(iStep) = ry26(iStep) - dnint(ry26(iStep)*Linv(2))*L(2)
        ry27(iStep) = ry27(iStep) - dnint(ry27(iStep)*Linv(2))*L(2)
        ry28(iStep) = ry28(iStep) - dnint(ry28(iStep)*Linv(2))*L(2)
        
        ry35(iStep) = ry35(iStep) - dnint(ry35(iStep)*Linv(2))*L(2)
        ry36(iStep) = ry36(iStep) - dnint(ry36(iStep)*Linv(2))*L(2)
        ry37(iStep) = ry37(iStep) - dnint(ry37(iStep)*Linv(2))*L(2)
        ry38(iStep) = ry38(iStep) - dnint(ry38(iStep)*Linv(2))*L(2)
        
        ry45(iStep) = ry45(iStep) - dnint(ry45(iStep)*Linv(2))*L(2)
        ry46(iStep) = ry46(iStep) - dnint(ry46(iStep)*Linv(2))*L(2)
        ry47(iStep) = ry47(iStep) - dnint(ry47(iStep)*Linv(2))*L(2)
        ry48(iStep) = ry48(iStep) - dnint(ry48(iStep)*Linv(2))*L(2)
    
        rz15(iStep) = rz15(iStep) - dnint(rz15(iStep)*Linv(3))*L(3)
        rz16(iStep) = rz16(iStep) - dnint(rz16(iStep)*Linv(3))*L(3)
        rz17(iStep) = rz17(iStep) - dnint(rz17(iStep)*Linv(3))*L(3)
        rz18(iStep) = rz18(iStep) - dnint(rz18(iStep)*Linv(3))*L(3)
        
        rz25(iStep) = rz25(iStep) - dnint(rz25(iStep)*Linv(3))*L(3)
        rz26(iStep) = rz26(iStep) - dnint(rz26(iStep)*Linv(3))*L(3)
        rz27(iStep) = rz27(iStep) - dnint(rz27(iStep)*Linv(3))*L(3)
        rz28(iStep) = rz28(iStep) - dnint(rz28(iStep)*Linv(3))*L(3)
        
        rz35(iStep) = rz35(iStep) - dnint(rz35(iStep)*Linv(3))*L(3)
        rz36(iStep) = rz36(iStep) - dnint(rz36(iStep)*Linv(3))*L(3)
        rz37(iStep) = rz37(iStep) - dnint(rz37(iStep)*Linv(3))*L(3)
        rz38(iStep) = rz38(iStep) - dnint(rz38(iStep)*Linv(3))*L(3)
        
        rz45(iStep) = rz45(iStep) - dnint(rz45(iStep)*Linv(3))*L(3)
        rz46(iStep) = rz46(iStep) - dnint(rz46(iStep)*Linv(3))*L(3)
        rz47(iStep) = rz47(iStep) - dnint(rz47(iStep)*Linv(3))*L(3)
        rz48(iStep) = rz48(iStep) - dnint(rz48(iStep)*Linv(3))*L(3)


        r15(iStep) = (rx15(iStep)**2 + ry15(iStep)**2 + rz15(iStep)**2)**0.5
        r16(iStep) = (rx16(iStep)**2 + ry16(iStep)**2 + rz16(iStep)**2)**0.5
        r17(iStep) = (rx17(iStep)**2 + ry17(iStep)**2 + rz17(iStep)**2)**0.5
        r18(iStep) = (rx18(iStep)**2 + ry18(iStep)**2 + rz18(iStep)**2)**0.5
        
        r25(iStep) = (rx25(iStep)**2 + ry25(iStep)**2 + rz25(iStep)**2)**0.5
        r26(iStep) = (rx26(iStep)**2 + ry26(iStep)**2 + rz26(iStep)**2)**0.5
        r27(iStep) = (rx27(iStep)**2 + ry27(iStep)**2 + rz27(iStep)**2)**0.5
        r28(iStep) = (rx28(iStep)**2 + ry28(iStep)**2 + rz28(iStep)**2)**0.5
        
        r35(iStep) = (rx35(iStep)**2 + ry35(iStep)**2 + rz35(iStep)**2)**0.5
        r36(iStep) = (rx36(iStep)**2 + ry36(iStep)**2 + rz36(iStep)**2)**0.5
        r37(iStep) = (rx37(iStep)**2 + ry37(iStep)**2 + rz37(iStep)**2)**0.5
        r38(iStep) = (rx38(iStep)**2 + ry38(iStep)**2 + rz38(iStep)**2)**0.5
        
        r45(iStep) = (rx45(iStep)**2 + ry45(iStep)**2 + rz45(iStep)**2)**0.5
        r46(iStep) = (rx46(iStep)**2 + ry46(iStep)**2 + rz46(iStep)**2)**0.5
        r47(iStep) = (rx47(iStep)**2 + ry47(iStep)**2 + rz47(iStep)**2)**0.5
        r48(iStep) = (rx48(iStep)**2 + ry48(iStep)**2 + rz48(iStep)**2)**0.5
        
        rat15(iStep) = (1-(r15(iStep)/r0)**6)/(1-(r15(iStep)/r0)**12)
        rat16(iStep) = (1-(r16(iStep)/r0)**6)/(1-(r16(iStep)/r0)**12)
        rat17(iStep) = (1-(r17(iStep)/r0)**6)/(1-(r17(iStep)/r0)**12)
        rat18(iStep) = (1-(r18(iStep)/r0)**6)/(1-(r18(iStep)/r0)**12)

        rat25(iStep) = (1-(r25(iStep)/r0)**6)/(1-(r25(iStep)/r0)**12)
        rat26(iStep) = (1-(r26(iStep)/r0)**6)/(1-(r26(iStep)/r0)**12)
        rat27(iStep) = (1-(r27(iStep)/r0)**6)/(1-(r27(iStep)/r0)**12)
        rat28(iStep) = (1-(r28(iStep)/r0)**6)/(1-(r28(iStep)/r0)**12)

        rat35(iStep) = (1-(r35(iStep)/r0)**6)/(1-(r35(iStep)/r0)**12)
        rat36(iStep) = (1-(r36(iStep)/r0)**6)/(1-(r36(iStep)/r0)**12)
        rat37(iStep) = (1-(r37(iStep)/r0)**6)/(1-(r37(iStep)/r0)**12)
        rat38(iStep) = (1-(r38(iStep)/r0)**6)/(1-(r38(iStep)/r0)**12)

        rat45(iStep) = (1-(r45(iStep)/r0)**6)/(1-(r45(iStep)/r0)**12)
        rat46(iStep) = (1-(r46(iStep)/r0)**6)/(1-(r46(iStep)/r0)**12)
        rat47(iStep) = (1-(r47(iStep)/r0)**6)/(1-(r47(iStep)/r0)**12)
        rat48(iStep) = (1-(r48(iStep)/r0)**6)/(1-(r48(iStep)/r0)**12)

        n(iStep) = 0.25*(rat15(iStep) + rat16(iStep) + rat17(iStep) + rat18(iStep) + rat25(iStep) &
        + rat26(iStep) + rat27(iStep) + rat28(iStep) + rat35(iStep) + rat36(iStep) + rat37(iStep) &
        + rat38(iStep)+ rat45(iStep) + rat46(iStep) + rat47(iStep) + rat48(iStep))
    end do

write(20, *) iStep, n(iStep)

end do



deallocate(label)
deallocate(rx)
deallocate(ry)
deallocate(rz)
deallocate(rx15)
deallocate(rx16)
deallocate(rx17)
deallocate(rx18)
deallocate(rx25)
deallocate(rx26)
deallocate(rx27)
deallocate(rx28)
deallocate(rx35)
deallocate(rx36)
deallocate(rx37)
deallocate(rx38)
deallocate(rx45)
deallocate(rx46)
deallocate(rx47)
deallocate(rx48)
deallocate(ry15)
deallocate(ry16)
deallocate(ry17)
deallocate(ry18)
deallocate(ry25)
deallocate(ry26)
deallocate(ry27)
deallocate(ry28)
deallocate(ry35)
deallocate(ry36)
deallocate(ry37)
deallocate(ry38)
deallocate(ry45)
deallocate(ry46)
deallocate(ry47)
deallocate(ry48)
deallocate(rz15)
deallocate(rz16)
deallocate(rz17)
deallocate(rz18)
deallocate(rz25)
deallocate(rz26)
deallocate(rz27)
deallocate(rz28)
deallocate(rz35)
deallocate(rz36)
deallocate(rz37)
deallocate(rz38)
deallocate(rz45)
deallocate(rz46)
deallocate(rz47)
deallocate(rz48)
deallocate(r15)
deallocate(r16)
deallocate(r17)
deallocate(r18)
deallocate(r25)
deallocate(r26)
deallocate(r27)
deallocate(r28)
deallocate(r35)
deallocate(r36)
deallocate(r37)
deallocate(r38)
deallocate(r45)
deallocate(r46)
deallocate(r47)
deallocate(r48)
deallocate(rat15)
deallocate(rat16)
deallocate(rat17)
deallocate(rat18)
deallocate(rat25)
deallocate(rat26)
deallocate(rat27)
deallocate(rat28)
deallocate(rat35)
deallocate(rat36)
deallocate(rat37)
deallocate(rat38)
deallocate(rat45)
deallocate(rat46)
deallocate(rat47)
deallocate(rat48)
end program coordnum