program edl
!****************************************************************
!  Electric Double-Layer
!      Calculate Electrostatic Potential from .xyz file
!****************************************************************
implicit none
  character(len=40) filename, filename_xyz, filename_mul, filename_out, step_cha, atom
  integer i, j, k, l, Natom, max, len, leng, step, skip
  integer Xdelta, Ydelta, Zdelta
  real min, delta, length, width, distance, potential, pot
  real Xmax, Xmin, Ymax, Ymin, Zmax, Zmin, rX, rY, rZ
  real :: element = 1.60217657
  real :: coulomb = 8.9876
  real,allocatable,dimension(:,:) :: r
  real,allocatable,dimension(:) :: mul
  real,allocatable,dimension(:) :: Xpotential, Ypotential, Zpotential
!===============================================================
! get argument
!===============================================================
  
  ! call getarg(1, filename_xyz)
  ! call getarg(2, filename_mul)
  ! call getarg(3, step_cha)
  ! write(filename_xyz,*) trim(filename)//'.xyz'
  ! write(filename_mul,*) trim(filename)//'.mul'
  ! write(filename_out,*) trim(filename)//'.csv'
  ! read(step_cha,*) NAtom
  
  ! read(10,*) Natom
  read(5,*) filename_xyz
  read(5,*) filename_mul
  read(5,*) NAtom
  allocate( r(Natom,3) )
  allocate( mul(Natom) )
  ! rewind(10)
open(10,file=filename_xyz, action="read")
open(20,file=filename_mul, action="read")
open(30,file="edl.dat", action="write")
write(*,*) 'A'
!   skip = ( Natom + 2 ) * ( step )
write(*,*) Natom
! write(*,*) step
! write(*,*) skip
!===============================================================
! read xyz file
!===============================================================
  ! if (skip > 0) then
  !    do i = 1, skip
  !       read(10,*)
  !    end do
  ! end if
write(*,*) 'A'
  read(10,*)
  read(10,*)
  do k = 1, Natom
     read(10,*) atom, r(k,1), r(k,2), r(k,3)
  end do
!===============================================================
! read mul file
!===============================================================
  ! if (skip > 0) then
  !    do i = 1, skip
  !       read(20,*)
  !    end do
  ! end if
  read(20,*)
  read(20,*)
  do k = 1, Natom
     read(20,*) atom, mul(k)
  end do
!===============================================================
! calculate distance
!===============================================================
  Xmax = maxval(r(:,1))+1
  Xmin = minval(r(:,1))-1
  Ymax = maxval(r(:,2))+1
  Ymin = minval(r(:,2))-1
  Zmax = maxval(r(:,3))+1
  Zmin = minval(r(:,3))-1
  Xdelta = ( Xmax - Xmin ) * 10
  Ydelta = ( Ymax - Ymin ) * 10
  Zdelta = ( Zmax - Zmin ) * 10
  allocate( Xpotential(Xdelta+1) )
  allocate( Ypotential(Ydelta+1) )
  allocate( Zpotential(Zdelta+1) )
write(*,*) 'A'
!===============================================================
! calculate potential
!===============================================================
  do i = 1, Xdelta+1
     rX = 0
     rX = Xmin + ( i - 1 ) * 0.1
     Ypotential = 0
     do j= 1, Ydelta+1
        rY = 0
        rY = Ymin + ( j - 1 ) * 0.1
        Zpotential = 0
        do k = 1, Zdelta+1
           rZ = 0
           rZ = Zmin + ( k - 1 ) * 0.1
           do l = 1, Natom
              distance = 0
              pot = 0
              distance = ( ( rX - r(l,1) )**2 + ( rY - r(l,2) )**2 + ( rZ - r(l,3) )**2 )**0.5
              if ( distance /= 0 ) then
                 pot = mul(l) / distance
                 potential = potential + pot
              end if
              if(l == Natom) then
                 Zpotential(k) = potential
                 potential = 0
              end if
           end do
           if ( k == Zdelta + 1 ) then
              Ypotential(j) = sum( Zpotential )
           end if
        end do
        if ( j == Ydelta + 1 ) then
           Xpotential(i) = sum( Ypotential ) / ( ( Ydelta + 1 ) * ( Zdelta + 1 ) )
        end if
     end do
     write(30,*) i, Xpotential(i)
  end do

  deallocate(r)
  deallocate(mul)
  deallocate(Xpotential)
  deallocate(Ypotential)
  deallocate(Zpotential)
end program edl
