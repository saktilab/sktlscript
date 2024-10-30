program analysis

implicit none
character*63 :: char63
character*2, allocatable :: element(:)

real*8, allocatable :: weight(:)
real*8, allocatable :: rx(:,:), ry(:,:), rz(:,:)
real*8, allocatable :: rcomx(:,:), rcomy(:,:), rcomz(:,:)
integer, allocatable :: nx(:,:), ny(:,:), nz(:,:)
integer, allocatable :: ncomx(:,:), ncomy(:,:), ncomz(:,:)
real*8, allocatable :: rsqs(:), rsqa(:), rsqna(:)
integer :: nstep, ndme, nat_dme, nsalt, nat_salt, i, j, k,  natom
real*8 :: alattice, wc, wo, wh, wn, wf, wna, wli, wtemps, wtempa, ws
real*8 :: rsqstemp, rsqatemp, rsqnatemp, rmax
real*8 :: rytemp1, rytemp2, rxtemp1, rxtemp2, rztemp1, rztemp2
integer :: ntemp1, ntemp2, ntemp
real*8 :: rxtemp, rytemp, rztemp, aversqs, aversqa, aversqna
integer, parameter :: print_step=1000
double precision, parameter :: time_step=1
read (*, '(a63)') char63
read (char63, *) nstep
print *, "nsteps: ", nstep
read (*, '(a63)') char63
read (char63, *) ndme, nat_dme
print *, "# of solvents: ", ndme, " , # of atom in solvent: ", nat_dme
read (*, '(a63)') char63
read (char63, *) nsalt, nat_salt
print *, "# of salt: ", nsalt, " , # of atom in salt: ", nat_salt
read (*, '(a63)') char63
read (char63, *) alattice

print *, "lattice: ", alattice
!1000 format(x,a2,3(f20.10))

!
! Set parameters
!

!nstep = 40001
!ndme = 62
!nsalt = 0
!alattice = 22.036835398173
!nat_dme=16
!nat_salt=10

!
! File open
!
!open(10,file="unfolded_traject.xyz")
!open(20,file="msd_unfolded_okoshi.txt")
open(10,file="dftb.traj")
open(20,file="msd.txt")

!
! Depending parameters
!
rmax = alattice / 2.0

natom = ndme*nat_dme+nsalt*nat_salt
wc = 12.01
wo = 16.01
wh = 1.008
wn = 14.007
wf = 19.00
ws = 32.0
!wna = 22.99
wna = 6.941
!
! Allocate memory
!
allocate(element(natom),weight(natom))
allocate(rx(nstep,natom),ry(nstep,natom),rz(nstep,natom))
allocate(rcomx(nstep,ndme+2*nsalt),rcomy(nstep,ndme+2*nsalt),rcomz(nstep,ndme+2*nsalt))
allocate(nx(nstep,natom),ny(nstep,natom),nz(nstep,natom))
allocate(ncomx(nstep,ndme+2*nsalt),ncomy(nstep,ndme+2*nsalt),ncomz(nstep,ndme+2*nsalt))
allocate(rsqs(nstep),rsqa(nstep),rsqna(nstep))

!
! Zero clear
!
element = '0'
weight = 0.0
rx = 0.0
ry = 0.0
rz = 0.0
rcomx = 0.0
rcomy = 0.0
rcomz = 0.0
nx = 0
ny = 0
nz = 0
ncomx = 0
ncomy = 0
ncomz = 0
rsqs = 0.0
rsqa = 0.0
rsqna = 0.0

!
! Read geometries
!
do i=1,nstep
write(*,*)i,"out of",nstep
   do j=1,2
      read(10,'(a63)')char63
   enddo
   do j=1,natom
      read(10,*)element(j),rx(i,j),ry(i,j),rz(i,j)
   enddo
enddo

!
! Weight assignment
!
do i=1,natom
   if(element(i).eq."C")weight(i)=wc
   if(element(i).eq."O")weight(i)=wo
   if(element(i).eq."H")weight(i)=wh
   if(element(i).eq."N")weight(i)=wn
   if(element(i).eq."F")weight(i)=wf
   if(element(i).eq."S")weight(i)=ws
   if(element(i).eq."Na")weight(i)=wna
enddo

!
! Minimal extraction under PBC
!
do i=1,nstep
   do j=1,ndme
      do k=2,nat_dme
         ntemp1 = (j-1)*nat_dme+1
         ntemp2 = (j-1)*nat_dme+k
         rxtemp1 = rx(i,ntemp1)
         rxtemp2 = rx(i,ntemp2)
         if(rxtemp2-rxtemp1.gt.rmax)nx(i,ntemp2)=nx(i,ntemp2)+1
         if(rxtemp1-rxtemp2.gt.rmax)nx(i,ntemp2)=nx(i,ntemp2)-1
         rytemp1 = ry(i,ntemp1)
         rytemp2 = ry(i,ntemp2)
         if(rytemp2-rytemp1.gt.rmax)ny(i,ntemp2)=ny(i,ntemp2)+1
         if(rytemp1-rytemp2.gt.rmax)ny(i,ntemp2)=ny(i,ntemp2)-1
         rztemp1 = rz(i,ntemp1)
         rztemp2 = rz(i,ntemp2)
         if(rztemp2-rztemp1.gt.rmax)nz(i,ntemp2)=nz(i,ntemp2)+1
         if(rztemp1-rztemp2.gt.rmax)nz(i,ntemp2)=nz(i,ntemp2)-1
      enddo
   enddo

   do j=1,nsalt
      do k=2,nat_salt-1
         ntemp1 = nat_dme*ndme+(j-1)*(nat_salt-1)+1
         ntemp2 = nat_dme*ndme+(j-1)*(nat_salt-1)+k
         rxtemp1 = rx(i,ntemp1)
         rxtemp2 = rx(i,ntemp2)
         if(rxtemp2-rxtemp1.gt.rmax)nx(i,ntemp2)=nx(i,ntemp2)+1
         if(rxtemp1-rxtemp2.gt.rmax)nx(i,ntemp2)=nx(i,ntemp2)-1
         rytemp1 = ry(i,ntemp1)
         rytemp2 = ry(i,ntemp2)
         if(rytemp2-rytemp1.gt.rmax)ny(i,ntemp2)=ny(i,ntemp2)+1
         if(rytemp1-rytemp2.gt.rmax)ny(i,ntemp2)=ny(i,ntemp2)-1
         rztemp1 = rz(i,ntemp1)
         rztemp2 = rz(i,ntemp2)
         if(rztemp2-rztemp1.gt.rmax)nz(i,ntemp2)=nz(i,ntemp2)+1
         if(rztemp1-rztemp2.gt.rmax)nz(i,ntemp2)=nz(i,ntemp2)-1
      enddo   
   enddo
enddo

!
! Center of mass
!
do i=1,nstep
   do j=1,ndme
      do k=1,nat_dme
         ntemp = (j-1)*nat_dme+k
         rxtemp = rx(i,ntemp)-alattice*nx(i,ntemp)
         rytemp = ry(i,ntemp)-alattice*ny(i,ntemp)
         rztemp = rz(i,ntemp)-alattice*nz(i,ntemp)
         rcomx(i,j) = rcomx(i,j)+weight(ntemp)*rxtemp
         rcomy(i,j) = rcomy(i,j)+weight(ntemp)*rytemp
         rcomz(i,j) = rcomz(i,j)+weight(ntemp)*rztemp
      enddo
   enddo

   do j=1,nsalt
      do k=1, (nat_salt-1)
         ntemp = ndme*nat_dme+(j-1)*(nat_salt-1)+k
         rxtemp = rx(i,ntemp)-alattice*nx(i,ntemp)
         rytemp = ry(i,ntemp)-alattice*ny(i,ntemp)
         rztemp = rz(i,ntemp)-alattice*nz(i,ntemp)
         rcomx(i,j+ndme) = rcomx(i,j+ndme)+weight(ntemp)*rxtemp
         rcomy(i,j+ndme) = rcomy(i,j+ndme)+weight(ntemp)*rytemp
         rcomz(i,j+ndme) = rcomz(i,j+ndme)+weight(ntemp)*rztemp
      enddo
   enddo

   do j=1,nsalt
      ntemp = ndme*nat_dme+nsalt*(nat_salt-1)+j
      rcomx(i,j+ndme+nsalt) = rx(i,ntemp)-alattice*nx(i,ntemp)
      rcomy(i,j+ndme+nsalt) = ry(i,ntemp)-alattice*ny(i,ntemp)
      rcomz(i,j+ndme+nsalt) = rz(i,ntemp)-alattice*nz(i,ntemp)
   enddo

enddo   

wtemps = 0.0
wtempa = 0.0

do i=1,nat_dme
   wtemps = wtemps+weight(i)
enddo
do i=ndme*nat_dme+1,ndme*nat_dme+(nat_salt-1)
   wtempa = wtempa+weight(i)
enddo

do i=1,ndme
   rcomx(:,i)=rcomx(:,i)/wtemps
   rcomy(:,i)=rcomy(:,i)/wtemps
   rcomz(:,i)=rcomz(:,i)/wtemps
enddo
do i=ndme+1,ndme+nsalt
   rcomx(:,i)=rcomx(:,i)/wtempa
   rcomy(:,i)=rcomy(:,i)/wtempa
   rcomz(:,i)=rcomz(:,i)/wtempa
enddo

!
! Unwrap CoMs
!
do i=2,nstep
   do j=1,ndme+2*nsalt
      rxtemp1=rcomx(i,j)-rcomx(i-1,j)
      rxtemp2=rcomx(i-1,j)-rcomx(i,j)
      rytemp1=rcomy(i,j)-rcomy(i-1,j)
      rytemp2=rcomy(i-1,j)-rcomy(i,j)
      rztemp1=rcomz(i,j)-rcomz(i-1,j)
      rztemp2=rcomz(i-1,j)-rcomz(i,j)
      ncomx(i,j)=ncomx(i-1,j)
      ncomy(i,j)=ncomy(i-1,j)
      ncomz(i,j)=ncomz(i-1,j)
      if(rxtemp1.gt.rmax)ncomx(i,j)=ncomx(i,j)+1
      if(rxtemp2.gt.rmax)ncomx(i,j)=ncomx(i,j)-1
      if(rytemp1.gt.rmax)ncomy(i,j)=ncomy(i,j)+1
      if(rytemp2.gt.rmax)ncomy(i,j)=ncomy(i,j)-1
      if(rztemp1.gt.rmax)ncomz(i,j)=ncomz(i,j)+1
      if(rztemp2.gt.rmax)ncomz(i,j)=ncomz(i,j)-1
   enddo
enddo   

!
! Calculate MSDs
!
do i=2,nstep

   do j=1,ndme
      rxtemp=rcomx(i,j)-alattice*ncomx(i,j)-rcomx(1,j)
      rytemp=rcomy(i,j)-alattice*ncomy(i,j)-rcomy(1,j)
      rztemp=rcomz(i,j)-alattice*ncomz(i,j)-rcomz(1,j)
      rsqs(i)=rsqs(i)+rxtemp**2+rytemp**2+rztemp**2
   enddo
   rsqs(i)=rsqs(i)/ndme

   do j=ndme+1,ndme+nsalt
      rxtemp=rcomx(i,j)-alattice*ncomx(i,j)-rcomx(1,j)
      rytemp=rcomy(i,j)-alattice*ncomy(i,j)-rcomy(1,j)
      rztemp=rcomz(i,j)-alattice*ncomz(i,j)-rcomz(1,j)
      rsqa(i)=rsqa(i)+rxtemp**2+rytemp**2+rztemp**2
   enddo
   if (nsalt == 0) then
       rsqa(i) = 0.0
   else
       rsqa(i)=rsqa(i)/nsalt
   end if

   do j=ndme+nsalt+1,ndme+nsalt+nsalt
      rxtemp=rcomx(i,j)-alattice*ncomx(i,j)-rcomx(1,j)
      rytemp=rcomy(i,j)-alattice*ncomy(i,j)-rcomy(1,j)
      rztemp=rcomz(i,j)-alattice*ncomz(i,j)-rcomz(1,j)
      rsqna(i)=rsqna(i)+rxtemp**2+rytemp**2+rztemp**2
   enddo 
   if (nsalt == 0) then
        rsqna(i) = 0
   else
       rsqna(i)=rsqna(i)/nsalt
   end if

enddo

!
! Output MSDs
!
write(20,'(a)')"#MSDs of solvents, anions, and sodium ions (ang-sq)"
write(20,'(a)')"#Step number is summed over 20 steps, i.e., 10 fs with 0.5 fs of time step"
k=1
do i=1,nstep/k
   rsqstemp=0.0
   rsqatemp=0.0
   rsqnatemp=0.0
   do j=(i-1)*k+1,(i-1)*k+k
      rsqstemp=rsqstemp+rsqs(j)
      rsqatemp=rsqatemp+rsqa(j)
      rsqnatemp=rsqnatemp+rsqna(j)
   enddo   
   aversqs=rsqstemp/k
   aversqa=rsqatemp/k
   aversqna=rsqnatemp/k

   write(20,'(f10.5,3x,3f10.5)') dble(i)*dble(k)*time_step/1000*print_step,aversqs,aversqa,aversqna
enddo   

end program
