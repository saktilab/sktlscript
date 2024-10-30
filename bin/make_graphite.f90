program make_graphite
!****************************************************************
!  Input Maker
!      graphite structutre for .xyz file
!****************************************************************
implicit none
  character(len=8) filename
  character(len=3) temp_vert, temp_hori, temp_layer
  integer i, j, k, l, Natom, vert, hori, layer, num, mai
  real,allocatable,dimension(:,:) :: r
  character,allocatable,dimension(:) :: atom
!===============================================================
! get argument
!===============================================================
  call getarg(1, temp_vert)
  call getarg(2, temp_hori)
  call getarg(3, temp_layer)
  read(temp_vert,*) vert
  read(temp_hori,*) hori
  read(temp_layer,*) layer
  open(10,file="gra.xyz", action="write")
  if (mod(vert,2) == 0) then
     Natom = vert * (2 * hori + 1) * layer + (vert * 2 + hori) * layer
  else
     Natom = vert * (2 * hori + 1) * layer + (vert * 2 + hori + 1) * layer
  end if
  allocate( r(Natom,3) )
  allocate( atom(Natom) )
!===============================================================
! calculate
!===============================================================
  r(:,:) = 0.0
  atom(:) = 'C'
  r(1,:) = (/0, 0, 0/)
  do i = 2, vert
     if (mod(i,2) == 0) then
        r(i,:) = (/r(i-1,1)+1.3915, r(i-1,2), r(i-1,3)/)
     else
        r(i,:) = (/r(i-1,1)+2.783, r(i-1,2), r(i-1,3)/)
     end if
  end do
  r(vert+1,:) = (/r(1,1)-1.3915*0.5, r(1,2)+1.3915*3**0.5*0.5, r(1,3)/)
  do i = 1, vert-1
     if (mod(i,2) == 0) then
        r(vert+i+1,:) = (/r(vert+i,1)+1.3915, r(vert+i,2), r(vert+i,3)/)
     else
        r(vert+i+1,:) = (/r(vert+i,1)+2.783, r(vert+i,2), r(vert+i,3)/)
     end if
  end do
!
  if (hori > 1) then
     k = 0
     do i = 1, hori - 1
        do j = 1, 2 * vert
           k = k + 1
           r(2*vert + k,:) = (/r(k,1), r(k,2)+1.3915*3**0.5, r(k,3)/)
        end do
     end do
  end if
!
  k = 0
  do i = 1, vert
     k = k + 1
     r(2*hori*vert + k,:) = (/r(2*vert*(hori-1)+k,1), r(2*vert*(hori-1)+k,2)+1.3915*3**0.5, r(2*vert*(hori-1)+k,3)/)
  end do
!
  if (mod(vert,2) == 0) then
     atom(vert*(2*hori+1)+1:vert*(2*hori+3)+hori) = 'H'
  else
     atom(vert*(2*hori+1)+1:vert*(2*hori+3)+hori+1) = 'H'
  end if
!
  r(vert*(2*hori+1)+1,:) = (/r(1,1)-1.09*0.5, r(1,2)-1.09*3**0.5*0.5, r(1,3)/)
  do i = 1, vert-1
     if (mod(i,2) == 0) then
        r(vert*(2*hori+1)+1+i,:) = (/r(vert*(2*hori+1)+i,1)+1.693, r(vert*(2*hori+1)+i,2), r(vert*(2*hori+1)+i,3)/)
     else
        r(vert*(2*hori+1)+1+i,:) = (/r(vert*(2*hori+1)+i,1)+2.4815, r(vert*(2*hori+1)+i,2), r(vert*(2*hori+1)+i,3)/)
     end if
  end do
  r(vert*(2*hori+2)+1,:) = (/r(vert*2*hori+1,1)-1.09*0.5, r(vert*2*hori+1,2)+1.09*3**0.5*0.5, r(vert*2*hori+1,3)/)
  do i = 1, vert-1
     if (mod(i,2) == 0) then
        r(vert*(2*hori+2)+1+i,:) = (/r(vert*(2*hori+2)+i,1)+1.693, r(vert*(2*hori+2)+i,2), r(vert*(2*hori+2)+i,3)/)
     else
        r(vert*(2*hori+2)+1+i,:) = (/r(vert*(2*hori+2)+i,1)+2.4815, r(vert*(2*hori+2)+i,2), r(vert*(2*hori+2)+i,3)/)
     end if
  end do
!
  if (mod(vert,2) == 0) then
     do i= 1, hori
        r(vert*(2*hori+3)+i,:) = (/r(2*vert*i,1)+1.09, r(2*vert*i,2), r(2*vert*i,3)/)
     end do
  else
     do i = 1, hori+1
        r(vert*(2*hori+3)+i,:) = (/r(2*vert*i-vert,1)+1.09, r(2*vert*i-vert,2), r(2*vert*i-vert,3)/)
     end do
  end if
!===============================================================
!ここまでで1枚分
!===============================================================
  mai = Natom / layer
  k = 0
  do i = 1, layer - 1
     do j = 1, mai
        k = k + 1
        if (mod(i,2) == 0) then
           r(mai+k,:) = (/r(k,1)+1.3915, r(k,2), r(k,3)+3.3/)
           atom(mai+k) = atom(k)
        else
           r(mai+k,:) = (/r(k,1)-1.3915, r(k,2), r(k,3)+3.3/)
           atom(mai+k) = atom(k)
        end if
     end do
  end do
!===============================================================
!write
!===============================================================
  write(10,*) NAtom
  write(10,'(A)')
  do i = 1, Natom
     write(10,'(A,3f15.7)') atom(i), r(i,1), r(i,2), r(i,3)
  end do
end program make_graphite
