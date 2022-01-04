	module UserCtrl
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	! This module defines the derived data type structure with all filenames
	implicit none
	! -------------------------------------------------------------
	! FOWARD MODE:		INPUT:	.CNTL	.MESH	(.msh.vtk,time slice,all the same)
	!					OUTPUT:	.FIELD	.field2 (.msh.vtk,time slice,all the same)
	character*1, parameter	:: FORWARD = 'f'
	character(1)	:: job
	character(80)	:: file_root
	logical			:: output_level
	! file used in forward process
	character(100) :: inm_file,ctl_file,opt_file,opt2_file,vtk_file,msh_file
	integer,parameter :: inm_unit = 101, ctl_unit = 102, opt_unit = 103, opt2_unit = 111, vtk_unit = 104, msh_unit = 105
	! files used in inversion process
	character(100) :: inv_mesh_file,inv_msh_file,inv_model_file,inv_data_file,inv_file
	integer,parameter :: inv_mesh_unit = 201, inv_msh_unit = 205, inv_model_unit = 202, inv_data_unit=203, inv_unit = 211
	character(100) :: inv_del_file,inv_del2_file,inv_frechet_file,inv_fr_file,inv_ftr_file,inv_wm_file
	integer,parameter :: inv_del_unit = 206, inv_del2_unit = 207, inv_frechet_unit = 231, inv_fr_unit = 232, inv_ftr_unit = 233, inv_wm_unit = 221
    character(100) :: inv_freecenter_file
    integer,parameter :: inv_freecenter_unit = 231
    contains
    
	subroutine parseArgs
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	USE IFPOSIX, only : PXFGETCWD
	implicit none
	integer :: iargc,narg,k,ilen,ierror
	character(80)  :: arg, verbose,buf
	narg = iargc()
	verbose = 'regular'
	if(narg .gt. 1) then
		k=1
		search_arg: do
			k=k+1
			if (k .eq. narg) exit  search_arg
			call getarg ( k, arg )
			if(arg(1:1).eq.'-') then
				if (arg(2:2).eq.'v') then
					call getarg ( k+1, verbose )
				end if
				narg=k-1
				exit search_arg
			end if
		end do search_arg
	end if
	! set the level of output based on the user input
	select case (verbose)
	case ('full')
		print *,'Output all information to screen and results to files.'
		output_level = .true.
	case ('none')
		!print*,'Output nothing to screen but only results to files.'
		output_level = .false.
		case default
		print *,'Output all information to screen and results to files.'
		output_level = .true.
	end select
	if(narg==2) then
		call getarg(1,arg)
		if(arg(1:1).eq.'-') job = arg(2:2)
    else
        write(*,*) 'TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang'
        write(*,*) 'This program comes with ABSOLUTELY NO WARRANTY.'
        write(*,*) 'This is free software, and you are welcome to redistribute it under certain conditions.'
		write(*,*) 'Usage: ./FETD -[job] FILE_ROOT'
		write(*,*) '[FORWARD]		-f	Calculates the predicted data and saves the EM solution'
		write(*,*) 'Optional		-v	[full|none]'
		stop
	endif
	! extract all following command line arguments
	call getarg(2,file_root)
	call PXFGETCWD (buf,ilen,ierror)
	write(*,*), 'TEMF3DT: Run job -',job,' with root name: ',TRIM(file_root)
	!write(*,*), 'FETD: Run job -',job,' in folder ',buf(1:ilen),' with root name: ',TRIM(file_root)
	! set file names
	inm_file = trim(file_root)//".mesh"
	ctl_file = trim(file_root)//".cntl"
	opt_file = trim(file_root)//".field"
	opt2_file = trim(file_root)//".field2"
	!character(20) :: vtk_file = trim(file_root)//".vtk"
	msh_file = trim(file_root)//".msh"
	!character(20) :: inv_mesh_file = "info.mesh"
	inv_msh_file = trim(file_root)//"inv"//".msh"
	inv_model_file = trim(file_root)//".model"
	inv_data_file = trim(file_root)//".data"
	inv_del_file = trim(file_root)//".del"
	inv_del2_file = trim(file_root)//".del2"
	inv_frechet_file = trim(file_root)//".frechet"
	inv_fr_file = trim(file_root)//".fr"
	inv_ftr_file = trim(file_root)//".ftr"
	inv_wm_file = trim(file_root)//".wm"
    inv_freecenter_file = trim(file_root)//".freecenter"
	end subroutine parseArgs
	end module UserCtrl

	module mpi_parameters
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	use mpi
	implicit none
	integer :: rank,size,mpi_stat(MPI_STATUS_SIZE),ierr
	end module mpi_parameters

	program main
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	USE IFPORT
	use fwd_modeling
	!use inv_modeling
	use mkl_service
	use omp_lib
	use system_lib
	use UserCtrl
	use mpi
	use mpi_parameters
	implicit none
	!include 'mpif.h'
	type(model) :: fwd_mod
	integer :: irc,cpuInfo(3)
	! Initialize task data on master MPI process
	call MPI_INIT (ierr)
	call MPI_COMM_SIZE (MPI_COMM_WORLD, size, ierr)
	call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)
	call kmp_set_warnings_off()
	! get file_root
	if(rank==0) call parseArgs
	! reading mesh
	if(rank==0) elapsed_time = TIMEF()
	if(output_level.and.rank==0) print*, "============================ building mesh ..."
	if(rank==0) call model_mesh (fwd_mod)
	if(rank==0) elapsed_time = TIMEF()
	if(output_level.and.rank==0) write(*,"(a,f12.2,a)")," Time passes: ",elapsed_time," seconds."
	! reading control information
	if(output_level.and.rank==0) print*, "============================ building control info ..."
	if(rank==0) call model_ctl (fwd_mod)
	if(rank==0) elapsed_time = TIMEF()
	if(output_level.and.rank==0) write(*,"(a,f12.2,a)")," Time passes: ",elapsed_time," seconds."
	! inversion control information
	call MPI_BCAST(job, 1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
	! mesh calculation
	if(output_level.and.rank==0) print*, "============================ mesh processing ..."
	if(rank==0) call mesh_process (fwd_mod)
	if(rank==0) elapsed_time = TIMEF()
	if(output_level.and.rank==0) write(*,"(a,f12.2,a)")," Time passes: ",elapsed_time," seconds."
	! time loop
	if(output_level.and.rank==0) print*, "============================ time loop ..."
	call time_solver (fwd_mod)
	if(rank==0) elapsed_time = TIMEF()
	if(output_level.and.rank==0) write(*,"(a,f12.2,a)")," Time passes: ",elapsed_time," seconds."
	! clear allocated arrays
	if(rank==0) call clear_model (fwd_mod)
	if(rank==0) CALL mkl_free_buffers()
	!if(rank==0) call clear_inv
	if(rank==0) elapsed_time = TIMEF()
	if(rank==0) write(*,"(a,f12.2,a)") " TEMF3DT: End of TEMF3DT simulation with time period ",elapsed_time," seconds."
	!call MPI_FINALIZE (ierr)
	call MPI_Abort(MPI_COMM_WORLD,0,ierr)
	end program main