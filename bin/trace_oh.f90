program trace_oh
    implicit none
    double precision, parameter :: ONE = 1.0d+00
    double precision :: Lx, Ly, Lz, Linv(3), L(3)
    character(len=20) :: coordfile
    integer :: NAtom, NStep, iStep, iAtom, jAtom, coord
    character(len=20), allocatable :: label(:)
    double precision, allocatable :: x(:), y(:), z(:), rx(:), ry(:), rz(:), r(:)
    integer, allocatable :: cn(:)
    integer :: n
read(5,*) coordfile
read(5,*) NAtom
read(5,*) NStep
read(5,*) Lx, Ly, Lz


L(1) = Lx
L(2) = Ly
L(3) = Lz
Linv(1:3) = ONE/L(1:3)

open(10, file=coordfile, action="read", status="old")
open(20, file="id_oh.dat", action="write")

allocate(label(NAtom))
allocate(x(NAtom))
allocate(y(NAtom))
allocate(z(NAtom))
allocate(rx(NAtom))
allocate(ry(NAtom))
allocate(rz(NAtom))
allocate(r(NAtom))
allocate(cn(NAtom))
do iStep = 1, NStep
    read(10,*)
    read(10,*)
    
    do iAtom = 1, NAtom
        read(10,*) label(iAtom), x(iAtom), y(iAtom), z(iAtom)
    end do
    
    do iAtom = 1, NAtom
        n = 0
        do jAtom = 1, NAtom
            if (label(iAtom) == "O" .and. label(jAtom) == "H") then
                rx(iAtom) = x(iAtom) - x(jAtom)
                ry(iAtom) = y(iAtom) - y(jAtom)
                rz(iAtom) = z(iAtom) - z(jAtom)
                
                rx(iAtom) = rx(iAtom) - dnint(rx(iAtom)*Linv(1))*L(1)
                ry(iAtom) = ry(iAtom) - dnint(ry(iAtom)*Linv(2))*L(2)
                rz(iAtom) = rz(iAtom) - dnint(rz(iAtom)*Linv(3))*L(3)
                
                r(iAtom) = (rx(iAtom)**2 + ry(iAtom)**2 + rz(iAtom)**2)**0.5
                ! cn(iAtom) = count(r .lt. 1.00 .and. r.gt.0.00)                        
                if (r(iAtom) .lt. 1.5 .and. r(iAtom) .gt. 0) then
                    n = n + 1
                    cn(iAtom) = n
                      if (maxval(cn)==1) write(20,*) iStep, iAtom, jAtom, r(iAtom), cn(iAtom)              
                end if
                
                end if
                    
            end do
            
        
        ! write(20,*) iStep, iAtom, cn(iAtom)
                ! if (cn(iAtom) == 1) then
                !     write(20,*) iStep, iAtom, cn(iAtom)
                ! end if
    end do
    
end do

deallocate(label)
deallocate(x)
deallocate(y)
deallocate(z)
deallocate(rx)
deallocate(ry)
deallocate(rz)
deallocate(r)
deallocate(cn)
end program trace_oh