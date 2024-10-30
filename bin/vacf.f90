program autocorrelation
 implicit none
 !!! This program contains the following variables
 !!! z: autocorrelation function
 !!! i, j, k, l: counter
 !!! m: number of atoms
 !!! t: total MD steps
 !!! ERRV: status variable
 !!! v: Velocity field of trajectory
 !!! AUX: Auxillary variable

 double precision, allocatable :: z(:), v(:,:,:)
 integer :: n, m, t, i, j, k, l, ERRV, mdprint
 character (len=10) :: dummy
 double precision :: AUX, dif, dt
 print *, "Enter the number of MD steps"
 read *, t
 print *, "Enter the number of atoms"
 read *, m
 print *, "Enter your MD time step in femtoseconds"
 read *, dt
 print *, "Enter your deposited frequency of MD step"
 read *, mdprint
 allocate(v(t,m,3), z(t), STAT=ERRV)
 if (ERRV/=0) STOP "Error in allocation"

 open (11, file="oveloc.dat", status='old', action='read')
 open (6, file="vacf.dat", action='write')
 !!! Now you will read your velocities 
 do i=1, t
 read(11,*)
 read(11,*)
     do k=1, m
       read(11,*) dummy, v(i,k,1:3)
     end do
 end do

 !!! sum along the velocities
    do i=1,t
       !!! Sum of the autocorrelation function
     do j=1,t-i+1
    !!! Sum over atoms
      print *, (j+i)
      do k=1,m
       !!! Dot products of velocities
         do l=1,3
            z(i)=z(i)+v(j+i,k,l)*v(j,k,l)
         end do
       end do
     end do
    end do
  !!! Scaling
     AUX=z(1)
     do i=1,t
        z(i)=z(i)/(AUX*(t-i+1))
     end do
    
 !!! Calculate diffusion coefficient
     dif=sum(z)/3
 !!! Write results along the time in picoseconds
    do i=1, t
       write(6,*) i*dt*mdprint/1000, z(i), dif
    end do
end program autocorrelation
