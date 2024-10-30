program test
 implicit none
 integer :: iStep, NGaus, NStep, iGaus, icv, freq, init
 double precision :: height, width
!  integer, parameter :: NStep = 70000
!  integer, parameter :: iGaus = 1
!  double precision, parameter :: height = 5e-4
!  integer, parameter :: iCV = 1
!  double precision, parameter :: width = 0.100000000000000
 double precision, allocatable :: cv(:)
 character (len=10) :: restart 
!  open(5, filename = "biaspot", action='read')
 open(6, file = "colvar.dat", action = 'read')
 open(10, file = "bias.dat", action = 'write')
 read(5,*) restart
 read(5,*) init
 read(5,*) NStep
 read(5,*) iGaus
 read(5,*) height
 read(5,*) iCV
 read(5,*) width
 read(5,*) freq

 
 allocate(cv(500))
 
 do iStep = 1, (NStep-init)/freq
    read(6,*) cv(iStep)
 end do

do iStep = init+freq, NStep, freq
 write(10, '(1x, "GAUSSIAN BIAS POTENTIAL:", 1x, i8)') (iStep-init)/freq
 if (restart == "restart") then 
    write(10, '(1x, "*** RECOVERED FROM RESTART FILE")')
 else
 write(10, '(1x, "*** AT T=",f14.2," FSEC, THIS RUN''S STEP NO.=",i8)') &
 &      iStep/2.0 , iStep
 end if
 
 write(10, '(5x, "Gaussian height     =", f20.10, " a.u.")') height
 write(10, '(5x, "Collective variable =", i20)') iCV
 write(10, '(7x, "Coordinate        =", f20.10)') cv((iStep-init)/freq) 
 write(10, '(7x, "Gaussian width    =", f20.10)') width
end do
 deallocate(cv)
end program test
