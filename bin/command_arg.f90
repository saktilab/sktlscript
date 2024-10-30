program TestArg
    integer::narg,cptArg !#of arg & counter of arg
    character(len=20)::name !Arg name
    
   !Check if any arguments are found
    narg=command_argument_count()
   !Loop over the arguments
    if(narg>0)then
   !loop across options
    do cptArg=1,narg
     call get_command_argument(cptArg,name)
      select case(adjustl(name))
       case("--help","-h")
        write(*,*)"This is program TestArg : Version 0.1"
       case default
        write(*,*)"Option ",adjustl(name),"unknown"
      end select
    end do
    end if
   end program TestArg