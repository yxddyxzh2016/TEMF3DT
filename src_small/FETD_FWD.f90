    module fwd_modeling
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    use UserCtrl
    USE IFCORE, only : commitqq
    use MKL_SPBLAS
    !USE MKL_SOLVERS_EE,only : mkl_sparse_ee_init,mkl_sparse_d_gv
    implicit none
    private
    ! time variables
    REAL(8),public :: elapsed_time
    ! --------------------- mathematical and physical parameters
    real(8),parameter :: pi  = 3.141592653589793238462643383279502884197_8
    real(8),parameter :: mu = pi * 4.d-7
    real(8),parameter :: epsilon = 8.8541878176d-12
    real(8),parameter :: c0 = 2.99792458d+8
    real(8),parameter :: r2d = 180.0_8/pi
    real(8),parameter :: d2r = pi/180.0_8
    real(8),parameter :: beta = 0.25
    ! --------------------- SYSTEM-LOCAL INDEX variables
    integer,parameter :: ind_t2d(3,2)=reshape([1,2,3,2,3,1],[3,2])
    integer,parameter :: ind_tet(6,2)=reshape([1,1,1,2,4,3,2,3,4,3,2,4],[6,2])
    integer,parameter :: ind_tri(4,3)=reshape([1,4,1,1,3,2,4,2,2,3,3,4],[4,3])
    integer,parameter :: ind_edg(4,3)=reshape([4,4,6,5,1,6,2,3,2,5,3,1],[4,3])
    integer,parameter :: ind_dir(4,3)=reshape([-1,1,-1,-1,-1,1,-1,-1,1,1,1,1],[4,3])
    type(MATRIX_DESCR) :: DESCR
    ! time variables
    real(8) :: time_assemble,time_frac,time_stepping
    real(8) :: time_1,time_2 ! used in outer statistics
    real(8) :: time_3,time_4 ! used in inner statistics
    ! --------------------- MODEL TYPE node, edge and tet information
    type model
        ! node
        integer :: nod_num
        real(8),allocatable :: nod_cor(:,:)
        ! source
        integer :: seg_num
        integer,allocatable :: src_vex(:,:)
        integer,allocatable :: src_phy(:)
        integer,allocatable :: src_edg(:)
        integer,allocatable :: src_drc(:)
        ! source info from the control file
        integer :: src_num
        integer :: src_sourcetype,src_pulsetype
        real(8) :: src_pulsewidth,src_current
        integer,allocatable :: src_tp(:),src_count_edges(:)
        real(8),allocatable :: src_end_1(:,:),src_end_2(:,:)
        ! tet
        integer :: tet_num
        integer,allocatable :: tet_vex(:,:)
        integer,allocatable :: tet_phy(:)
        integer,allocatable :: tet_drc(:,:)
        real(8),allocatable :: tet_vol(:)
        integer,allocatable :: tet_edg(:,:)
        real(8),allocatable :: tet_len(:,:)
        real(8),allocatable :: a(:,:),b(:,:),c(:,:),d(:,:)
        ! surface edges
        integer :: sur_edg_num
        logical,allocatable :: sur_edg_list(:)
        integer,allocatable :: sur_edg_index(:)
        ! edge
        integer :: edg_num
        integer,allocatable :: edg_vex(:,:)
        real(8),allocatable :: edg_len(:)
        ! receiver
        integer :: rec_num
        real(8) :: rec_charac
        real(8),allocatable :: rec_cor(:,:)
        integer,allocatable :: rec_tet(:)
        real(8),allocatable :: n_arr(:,:,:),d_arr(:,:,:)
        ! domain character
        logical :: iso
        integer :: dom_num
        integer,allocatable :: dom_ind(:)
        real(8),allocatable :: sig_3d(:),per_3d(:),mag_3d(:)
        real(8),allocatable :: sig_a3d(:,:,:),per_a3d(:,:,:),mag_a3d(:,:,:)
        real(8) :: sig_parameter
        ! system matrix
        type(SPARSE_MATRIX_T) :: a_mkl,b_mkl,c_mkl
        ! surf_kind: 0 - PEC boundary ; 1 - 1st order ABC
        integer :: bound_cond
        ! time variables
        integer :: step_mthd,step_dvsi,step_dble,step_dblesize
        real(8) :: step_init,time_maxm
        ! field variables
        real(8),allocatable :: a_f(:,:,:),a_f_channel(:,:,:)
        integer :: ind_f,init_size_F
        ! output variables: ex ey ez dx dy dz bx by bz
        logical :: o_f(9)
        integer,allocatable :: o_v(:)
        ! vtk output variables
        integer :: vtk_num
        real(8),allocatable :: vtk_times(:)
        logical,allocatable :: vtk_done(:)
        real(8),allocatable :: n_arr_global(:,:,:),d_arr_global(:,:,:)
        real(8),allocatable :: e_field_global(:,:),d_field_global(:,:),b_field_global(:,:)
        ! inversion variables
        logical :: deriv
        real(8),allocatable :: a_prime(:,:) ! the matrix that change the current model edge response e field into receiver point field
        integer :: total_step_record
        real(8),allocatable :: step_record(:),time_record(:),e_record(:,:)
        ! time channels
        integer :: time_channel_num
        real(8),allocatable :: time_channels(:)
        logical,allocatable :: time_channel_ok(:)
    END TYPE model
    public :: model
    public :: model_mesh,model_ctl,mesh_process,time_solver,clear_model
    public :: mkl_add_3,f_kind_iso,mkl_to_csr,msh_file
    contains

    subroutine model_mesh(cmod)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    use utility
    implicit none
    ! in/out variable
    type(model),intent(out) :: cmod
    ! local variables
    character(len=200) :: buffer
    character(len=5) :: word
    integer :: ierr,i,j,k,tri_num,vex1,vex2,tmp
    integer,allocatable :: edge_1(:),edge_2(:),tri_tmp(:,:),src_vex_tmp(:,:),src_phy_tmp(:)
    logical,allocatable :: tri_on_bound(:)
    real(8) :: x(4),y(4),z(4),plane(4),vec1(3),vec2(3)
    real(8),allocatable :: plane_judge(:)
    ! open file
    if(output_level) print*, "--- reading mesh file ..."
    open(inm_unit,file=trim(inm_file))
    ! nod_num nod_cor seg_num src_vex src_phy tet_num tet_vex tet_phy
    do
        read (inm_unit,"(a)",iostat=ierr) buffer
        if (ierr /= 0) exit
        read (buffer,*) word
        if (word == "Verti") then
            read (inm_unit,*) cmod%nod_num
            allocate(cmod%nod_cor(cmod%nod_num,3))
            do i=1,cmod%nod_num
                read (inm_unit,*) cmod%nod_cor(i,1),cmod%nod_cor(i,2),cmod%nod_cor(i,3)
            end do
        else if (word == "Edges") then
            read (inm_unit,*) cmod%seg_num
            allocate(src_vex_tmp(cmod%seg_num,2),src_phy_tmp(cmod%seg_num))
            do i=1,cmod%seg_num
                read (inm_unit,*) src_vex_tmp(i,1),src_vex_tmp(i,2),src_phy_tmp(i)
            end do
            cmod%seg_num = count(src_phy_tmp/=0)
            allocate(cmod%src_vex(cmod%seg_num,2),cmod%src_phy(cmod%seg_num))
            cmod%src_vex(:,1) = pack(src_vex_tmp(:,1),src_phy_tmp/=0)
            cmod%src_vex(:,2) = pack(src_vex_tmp(:,2),src_phy_tmp/=0)
            cmod%src_phy = pack(src_phy_tmp,src_phy_tmp/=0)
        else if (word == "Trian") then
            read (inm_unit,*) tri_num
            allocate(tri_tmp(tri_num,3))
            do i=1,tri_num
                read (inm_unit,*) tri_tmp(i,1),tri_tmp(i,2),tri_tmp(i,3)
            end do
        else if (word == "Tetra") then
            read (inm_unit,*) cmod%tet_num
            allocate(cmod%tet_vex(cmod%tet_num,4),cmod%tet_phy(cmod%tet_num))
            do i=1,cmod%tet_num
                read (inm_unit,*) cmod%tet_vex(i,1),cmod%tet_vex(i,2),cmod%tet_vex(i,3),cmod%tet_vex(i,4),cmod%tet_phy(i)
            end do
        end if
    end do
    close(inm_unit)
    ! edg_num edg_vex
    if(output_level) print*, "--- rebuild edge info ..."
    allocate(edge_1(cmod%tet_num*6),edge_2(cmod%tet_num*6))
    k=0
    do i=1,cmod%tet_num
        do j=1,6
            k=k+1
            edge_1(k) = cmod%tet_vex(i,ind_tet(j,1))
            edge_2(k) = cmod%tet_vex(i,ind_tet(j,2))
            call swap(edge_1(k),edge_2(k),edge_1(k)>edge_2(k))
        end do
    end do
    call unique_ii(edge_1,edge_2)
    cmod%edg_num = size(edge_1)
    allocate(cmod%edg_vex(cmod%edg_num,2))
    cmod%edg_vex(:,1) = edge_1
    cmod%edg_vex(:,2) = edge_2
    ! edg_len
    allocate(cmod%edg_len(cmod%edg_num))
    do i=1,cmod%edg_num
        cmod%edg_len(i) = dist ((/cmod%nod_cor(cmod%edg_vex(i,1),1),cmod%nod_cor(cmod%edg_vex(i,1),2),cmod%nod_cor(cmod%edg_vex(i,1),3)/),&
            &(/cmod%nod_cor(cmod%edg_vex(i,2),1),cmod%nod_cor(cmod%edg_vex(i,2),2),cmod%nod_cor(cmod%edg_vex(i,2),3)/))
    end do
    if(output_level) print*, "--- rebuild tetrahedron info ..."
    ALLOCATE(cmod%tet_edg(cmod%tet_NUM,6),cmod%tet_drc(cmod%tet_NUM,6))
    !$OMP PARALLEL DO private(i,j)
    DO I=1,cmod%tet_NUM
        do j=1,6
            cmod%tet_edg(I,j) = bSearch_Ii ( cmod%EDG_VEX(:,1), cmod%EDG_VEX(:,2), &
                min(cmod%tet_VEX(I,ind_tet(j,1)),cmod%tet_VEX(I,ind_tet(j,2))), max(cmod%tet_VEX(I,ind_tet(j,1)),cmod%tet_VEX(I,ind_tet(j,2))) )
            if ( cmod%tet_VEX(I,ind_tet(j,1)) < cmod%tet_VEX(I,ind_tet(j,2)) ) then
                cmod%tet_drc(i,j) = 1
            else
                cmod%tet_drc(i,j) =-1
            end if
        end do
    END DO
    !$OMP end PARALLEL DO
    allocate(cmod%a(cmod%tet_num,4),cmod%b(cmod%tet_num,4),cmod%c(cmod%tet_num,4),cmod%d(cmod%tet_num,4))
    allocate(cmod%tet_vol(cmod%tet_num),cmod%tet_len(cmod%tet_NUm,6))
    !$OMP PARALLEL DO private(i,j,x,y,z)
    do i=1,cmod%tet_num
        do j=1,4
            x(j) = cmod%nod_cor(cmod%tet_VEX(i,j),1)
            y(j) = cmod%nod_cor(cmod%tet_VEX(i,j),2)
            z(j) = cmod%nod_cor(cmod%tet_VEX(i,j),3)
        end do
        ! a
        cmod%a(i,1) = det_mkl_3(RESHAPE([x(2),x(3),x(4),y(2),y(3),y(4),z(2),z(3),z(4)],[3,3])) * ( 1.d0)
        cmod%a(i,2) = det_mkl_3(RESHAPE([x(1),x(3),x(4),y(1),y(3),y(4),z(1),z(3),z(4)],[3,3])) * (-1.d0)
        cmod%a(i,3) = det_mkl_3(RESHAPE([x(1),x(2),x(4),y(1),y(2),y(4),z(1),z(2),z(4)],[3,3])) * ( 1.d0)
        cmod%a(i,4) = det_mkl_3(RESHAPE([x(1),x(2),x(3),y(1),y(2),y(3),z(1),z(2),z(3)],[3,3])) * (-1.d0)
        ! b
        cmod%b(i,1) = det_mkl_3(RESHAPE([1.d0,1.d0,1.d0,y(2),y(3),y(4),z(2),z(3),z(4)],[3,3])) * (-1.d0)
        cmod%b(i,2) = det_mkl_3(RESHAPE([1.d0,1.d0,1.d0,y(1),y(3),y(4),z(1),z(3),z(4)],[3,3])) * ( 1.d0)
        cmod%b(i,3) = det_mkl_3(RESHAPE([1.d0,1.d0,1.d0,y(1),y(2),y(4),z(1),z(2),z(4)],[3,3])) * (-1.d0)
        cmod%b(i,4) = det_mkl_3(RESHAPE([1.d0,1.d0,1.d0,y(1),y(2),y(3),z(1),z(2),z(3)],[3,3])) * ( 1.d0)
        ! c
        cmod%c(i,1) = det_mkl_3(RESHAPE([1.d0,1.d0,1.d0,x(2),x(3),x(4),z(2),z(3),z(4)],[3,3])) * ( 1.d0)
        cmod%c(i,2) = det_mkl_3(RESHAPE([1.d0,1.d0,1.d0,x(1),x(3),x(4),z(1),z(3),z(4)],[3,3])) * (-1.d0)
        cmod%c(i,3) = det_mkl_3(RESHAPE([1.d0,1.d0,1.d0,x(1),x(2),x(4),z(1),z(2),z(4)],[3,3])) * ( 1.d0)
        cmod%c(i,4) = det_mkl_3(RESHAPE([1.d0,1.d0,1.d0,x(1),x(2),x(3),z(1),z(2),z(3)],[3,3])) * (-1.d0)
        ! d
        cmod%d(i,1) = det_mkl_3(RESHAPE([1.d0,1.d0,1.d0,x(2),x(3),x(4),y(2),y(3),y(4)],[3,3])) * (-1.d0)
        cmod%d(i,2) = det_mkl_3(RESHAPE([1.d0,1.d0,1.d0,x(1),x(3),x(4),y(1),y(3),y(4)],[3,3])) * ( 1.d0)
        cmod%d(i,3) = det_mkl_3(RESHAPE([1.d0,1.d0,1.d0,x(1),x(2),x(4),y(1),y(2),y(4)],[3,3])) * (-1.d0)
        cmod%d(i,4) = det_mkl_3(RESHAPE([1.d0,1.d0,1.d0,x(1),x(2),x(3),y(1),y(2),y(3)],[3,3])) * ( 1.d0)
        ! tet volume
        cmod%tet_vol(i) = tet_volume (X,Y,Z)
        ! tet edge length
        do j=1,6
            cmod%tet_len(i,j) = dist ( [x(ind_tet(j,1)),y(ind_tet(j,1)),z(ind_tet(j,1))],[x(ind_tet(j,2)),y(ind_tet(j,2)),z(ind_tet(j,2))] )
        end do
    end do
    !$OMP end PARALLEL DO
    end subroutine model_mesh

    subroutine model_ctl (cmod)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    use utility
    implicit none
    type(model),intent(inout) :: cmod
    real(8),allocatable :: sig_iso(:),per_iso(:),mag_iso(:)
    real(8),allocatable :: sig_ani(:,:),sig_ang(:,:),per_ani(:,:),per_ang(:,:),mag_ani(:,:),mag_ang(:,:)
    character(180) ::sLine, sCode, sValue
    logical ::bComment
    integer :: ierr,i,j,tmp(9)
    real(8) :: dx(3,3),dy(3,3),dz(3,3),dm(3,3),pm(3,3),direc(3)
    logical :: find_phy
    ! init variables
    cmod%rec_num = 0
    cmod%dom_num = 0
    cmod%bound_cond = 0
    cmod%src_num = 0
    cmod%step_mthd=0
    cmod%step_dvsi = 0
    cmod%time_maxm =0.d0
    cmod%step_dble = 0
    cmod%o_f = .false.
    cmod%vtk_num = 0
    cmod%rec_charac = 10.d0
    ! read in control file
    if(output_level) print*, "--- reading control file ..."
    open(ctl_unit,file=trim(ctl_file))
    do while (.true.)
        read( ctl_unit, '(A180)', iostat = ierr ) sLine
        if (ierr /= 0) exit
        call parseCode( len(sLine), sLine, sCode, sValue, bComment )
        if( bComment ) cycle
        select case (trim(sCode))
        case('time_points')
            read(sValue,*) cmod%time_channel_num
            allocate(cmod%time_channels(cmod%time_channel_num),cmod%time_channel_ok(cmod%time_channel_num))
            do i=1,cmod%time_channel_num
                read (ctl_unit,*) cmod%time_channels(i)
            end do
            cmod%time_channel_ok = .false.
        case ('rec_charac')
            read(sValue,*) cmod%rec_charac
        case ('receivers')
            read(sValue,*) cmod%rec_num
            allocate(cmod%rec_cor(cmod%rec_num,3))
            do i=1,cmod%rec_num
                read (ctl_unit,*) cmod%rec_cor(i,1),cmod%rec_cor(i,2),cmod%rec_cor(i,3)
            end do
            !! filter receivers
            !allocate(rec_filt(cmod%rec_num))
            !rec_filt=.false.
            !do i=1,cmod%rec_num
            !    if (cmod%rec_cor(i,2)<=-0.4) rec_filt(i)=.true.
            !end do
            !cmod%rec_num = count(rec_filt)
            !allocate (rec_tmp(cmod%rec_num,3))
            !rec_tmp(:,1) = pack(cmod%rec_cor(:,1),rec_filt)
            !rec_tmp(:,2) = pack(cmod%rec_cor(:,2),rec_filt)
            !rec_tmp(:,3) = pack(cmod%rec_cor(:,3),rec_filt)
            !deallocate(cmod%rec_cor)
            !allocate(cmod%rec_cor,source=rec_tmp)
        case ('isotropic')
            cmod%iso = .true.
            read(sValue,*) cmod%dom_num
            allocate(cmod%dom_ind(cmod%dom_num))
            allocate(sig_iso(cmod%dom_num),per_iso(cmod%dom_num),mag_iso(cmod%dom_num))
            cmod%sig_parameter = 1.d0
            do i=1,cmod%dom_num
                read (ctl_unit,*) cmod%dom_ind(i),sig_iso(i),per_iso(i),mag_iso(i)
                if (i>1) then
                    cmod%sig_parameter = sig_iso(i)
                elseif (i>2) then
                    if (sig_iso(i)<cmod%sig_parameter) cmod%sig_parameter = sig_iso(i)
                end if
            end do
            cmod%sig_parameter = 1.d0/cmod%sig_parameter
            if (cmod%dom_num>1) then
                do i=2,cmod%dom_num
                    sig_iso(i) = sig_iso(i) * cmod%sig_parameter
                end do
            end if
            allocate(cmod%sig_3d(cmod%tet_num),cmod%per_3d(cmod%tet_num),cmod%mag_3d(cmod%tet_num))
            do i=1,cmod%tet_num
                find_phy = .false.
                do j=1,cmod%dom_num
                    if (cmod%tet_phy(i)==cmod%dom_ind(j)) then
                        find_phy = .true.
                        cmod%sig_3d(i) = sig_iso(j)
                        cmod%per_3d(i) = per_iso(j) * epsilon
                        cmod%mag_3d(i) = mag_iso(j) * mu
                    end if
                end do
                if (.not.find_phy) stop("error - physical domain not found!")
            end do
        case ('anisotropic')
            cmod%iso = .false.
            read(sValue,*) cmod%dom_num
            allocate(cmod%dom_ind(cmod%dom_num))
            allocate(sig_ani(cmod%dom_num,3),sig_ang(cmod%dom_num,3),per_ani(cmod%dom_num,3),per_ang(cmod%dom_num,3),mag_ani(cmod%dom_num,3),mag_ang(cmod%dom_num,3))
            do i=1,cmod%dom_num
                read (ctl_unit,*) cmod%dom_ind(i),sig_ani(i,1),sig_ani(i,2),sig_ani(i,3),sig_ang(i,1),sig_ang(i,2),sig_ang(i,3)
                read (ctl_unit,*) cmod%dom_ind(i),per_ani(i,1),per_ani(i,2),per_ani(i,3),per_ang(i,1),per_ang(i,2),per_ang(i,3)
                read (ctl_unit,*) cmod%dom_ind(i),mag_ani(i,1),mag_ani(i,2),mag_ani(i,3),mag_ang(i,1),mag_ang(i,2),mag_ang(i,3)
            end do
            allocate(cmod%sig_a3d(cmod%tet_num,3,3),cmod%per_a3d(cmod%tet_num,3,3),cmod%mag_a3d(cmod%tet_num,3,3))
            ! convert all deg to reg
            sig_ang = sig_ang * d2r; per_ang = per_ang * d2r; mag_ang = mag_ang * d2r
            do i=1,cmod%tet_num
                find_phy = .false.
                do j=1,cmod%dom_num
                    if (cmod%tet_phy(i)==cmod%dom_ind(j)) then
                        find_phy = .true.
                        ! sig
                        dx = transpose(reshape((/1.d0,0.d0,0.d0,0.d0,cos(sig_ang(j,1)),-sin(sig_ang(j,1)),0.d0,sin(sig_ang(j,1)),cos(sig_ang(j,1))/),(/3,3/)))
                        dy = transpose(reshape((/cos(sig_ang(j,2)),0.d0,sin(sig_ang(j,2)),0.d0,1.d0,0.d0,-sin(sig_ang(j,2)),0.d0,cos(sig_ang(j,2))/),(/3,3/)))
                        dz = transpose(reshape((/cos(sig_ang(j,3)),-sin(sig_ang(j,3)),0.d0,sin(sig_ang(j,3)),cos(sig_ang(j,2)),0.d0,0.d0,0.d0,1.d0/),(/3,3/)))
                        dm = matmul(matmul(dx,dy),dz)
                        pm = transpose(reshape((/sig_ani(j,1),0.d0,0.d0,0.d0,sig_ani(j,2),0.d0,0.d0,0.d0,sig_ani(j,3)/),(/3,3/)))
                        cmod%sig_a3d(i,:,:) = matmul(matmul(dm,pm),transpose(dm))
                        ! mag
                        dx = transpose(reshape((/1.d0,0.d0,0.d0,0.d0,cos(mag_ang(j,1)),-sin(mag_ang(j,1)),0.d0,sin(mag_ang(j,1)),cos(mag_ang(j,1))/),(/3,3/)))
                        dy = transpose(reshape((/cos(mag_ang(j,2)),0.d0,sin(mag_ang(j,2)),0.d0,1.d0,0.d0,-sin(mag_ang(j,2)),0.d0,cos(mag_ang(j,2))/),(/3,3/)))
                        dz = transpose(reshape((/cos(mag_ang(j,3)),-sin(mag_ang(j,3)),0.d0,sin(mag_ang(j,3)),cos(mag_ang(j,2)),0.d0,0.d0,0.d0,1.d0/),(/3,3/)))
                        dm = matmul(matmul(dx,dy),dz)
                        pm = transpose(reshape((/mag_ani(j,1),0.d0,0.d0,0.d0,mag_ani(j,2),0.d0,0.d0,0.d0,mag_ani(j,3)/),(/3,3/)))
                        cmod%mag_a3d(i,:,:) = matmul(matmul(dm,pm),transpose(dm)) * mu
                        ! per
                        dx = transpose(reshape((/1.d0,0.d0,0.d0,0.d0,cos(per_ang(j,1)),-sin(per_ang(j,1)),0.d0,sin(per_ang(j,1)),cos(per_ang(j,1))/),(/3,3/)))
                        dy = transpose(reshape((/cos(per_ang(j,2)),0.d0,sin(per_ang(j,2)),0.d0,1.d0,0.d0,-sin(per_ang(j,2)),0.d0,cos(per_ang(j,2))/),(/3,3/)))
                        dz = transpose(reshape((/cos(per_ang(j,3)),-sin(per_ang(j,3)),0.d0,sin(per_ang(j,3)),cos(per_ang(j,2)),0.d0,0.d0,0.d0,1.d0/),(/3,3/)))
                        dm = matmul(matmul(dx,dy),dz)
                        pm = transpose(reshape((/per_ani(j,1),0.d0,0.d0,0.d0,per_ani(j,2),0.d0,0.d0,0.d0,per_ani(j,3)/),(/3,3/)))
                        cmod%per_a3d(i,:,:) = matmul(matmul(dm,pm),transpose(dm)) * epsilon
                    end if
                end do
                if (.not.find_phy) stop("error - physical domain not found!")
            end do
        case ('boundary condition')
            read(sValue,*) cmod%bound_cond
        case ('source')
            !source: ind source_type(1 electric, 2 magnetic) pulse_type (1 Turn-on step 2 Gaussian pulse 3 delta pulse) t0 current
            !   1    1    1.0e-7    1.0
            !   51
            !   x_start y_start z_start
            !   x_end y_end z_end
            !   ... ...
            read(sValue,*) cmod%src_num
            read(ctl_unit,*) cmod%src_sourcetype,cmod%src_pulsetype,cmod%src_pulsewidth,cmod%src_current
            !!allocate(cmod%src_tp(cmod%src_num,3),cmod%src_t0(cmod%src_num),cmod%src_t1(cmod%src_num),cmod%src_cu(cmod%src_num))
            !allocate(cmod%src_tp(cmod%src_num),cmod%src_end_1(cmod%src_num,3),cmod%src_end_2(cmod%src_num,3))
            !do i=1,cmod%src_num
            !	read (ctl_unit,*)	cmod%src_tp(i)
            !	read (ctl_unit,*)	cmod%src_end_1(i,1),cmod%src_end_1(i,2),cmod%src_end_1(i,3)
            !	read (ctl_unit,*)	cmod%src_end_2(i,1),cmod%src_end_2(i,2),cmod%src_end_2(i,3)
            !end do
            !! src_count_edges
            !allocate(cmod%src_count_edges(cmod%src_num))
            !do i=1,cmod%src_num
            !	cmod%src_count_edges(i) = count(cmod%src_phy==cmod%src_tp(i))
            !end do
            ! src_edg
            if(output_level) print*, "--- rebuild source info ..."
            allocate(cmod%src_edg(cmod%seg_num))
            do i=1,cmod%seg_num
                do j=1,cmod%edg_num
                    if (min(cmod%edg_vex(j,1),cmod%edg_vex(j,2))==min(cmod%src_vex(i,1),cmod%src_vex(i,2)).and.&
                        &max(cmod%edg_vex(j,1),cmod%edg_vex(j,2))==max(cmod%src_vex(i,1),cmod%src_vex(i,2))) then
                        cmod%src_edg(i) = j
                    end if
                end do
            end do
            ! src_drc
            allocate(cmod%src_drc(cmod%seg_num))
            cmod%src_drc = 1
            cmod%src_drc(cmod%seg_num) = -1
            !!!!!!!!cmod%src_drc(1:cmod%seg_num/2-1) = -1
            !!!!!!!!cmod%src_drc(cmod%seg_num/2) = 1
            !!!!!!!!cmod%src_drc(cmod%seg_num/2+1:cmod%seg_num-1) = 1
            !!!!!!!!cmod%src_drc(cmod%seg_num) = -1
            !cmod%src_drc = 1
            !cmod%src_drc(19) = -1
            !cmod%src_drc(77) = -1
            !cmod%src_drc(1:553) = 250
            !cmod%src_drc(554) = -250
            !cmod%src_drc(555:609) = 560
            !cmod%src_drc(610) = -560
        case ('sp_mthd')
            read(sValue,*) cmod%step_mthd
        case ('sp_division')
            read(sValue,*) cmod%step_dvsi
        case ('time_maximum')
            read(sValue,*) cmod%time_maxm
        case ('sp_double')
            read(sValue,*) cmod%step_dble
        case ('sp_dblesize')
            read(sValue,*) cmod%step_dblesize
        case ('output')
            read (ctl_unit,*) tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6),tmp(7),tmp(8),tmp(9)
            do i=1,9
                if (tmp(i)==1) then
                    cmod%o_f(i)=.true.
                else
                    cmod%o_f(i)=.false.
                end if
            end do
            allocate(cmod%o_v(0))
            do i=1,9
                if(cmod%o_f(i)) cmod%o_v=[cmod%o_v,i]
            end do
        case ('vtkout')
            read(sValue,*) cmod%vtk_num
            if (cmod%vtk_num>0) then
                allocate(cmod%vtk_times(cmod%vtk_num),cmod%vtk_done(cmod%vtk_num))
                do i=1,cmod%vtk_num
                    read (ctl_unit,*) cmod%vtk_times(i)
                end do
                cmod%vtk_done = .false.
            end if
            case default
        end select
    end do
    close(ctl_unit)
    ! check if all readed
    if (cmod%rec_num==0) stop "cmod%rec_num error"
    if (cmod%dom_num==0) stop "cmod%dom_num error"
    if (cmod%src_num==0) stop "cmod%src_num error"
    if (cmod%step_mthd==0) stop "cmod%step_mthd error"
    if (cmod%step_dvsi==0) stop "cmod%step_dvsi error"
    if (abs(cmod%time_maxm)<1.d-20) stop "cmod%time_maxm error"
    if (cmod%step_dble==0) stop "cmod%step_dble error"
    if (all(.not.cmod%o_f)) stop "cmod%o_f error"
    ! rec_tet
    if(output_level) print*, "--- rebuild receiver info ..."
    allocate(cmod%rec_tet(cmod%rec_num))
    cmod%rec_tet = 0
    !$OMP PARALLEL DO private(i)
    do i = 1,cmod%rec_num
        call p_in_tet (cmod,[cmod%rec_cor(i,1),cmod%rec_cor(i,2),cmod%rec_cor(i,3)],cmod%rec_tet(i))
        if (cmod%rec_tet(i)==0) then
            print*, "no rec tet No.",i
            stop
        end if
        !if (cmod%rec_tet(i)==0) then
        !	cmod%rec_cor(i,1)=cmod%rec_cor(i,1)+1
        !	cmod%rec_cor(i,2)=cmod%rec_cor(i,2)+1
        !	cmod%rec_cor(i,3)=cmod%rec_cor(i,3)+1
        !	call p_in_tet (cmod,[cmod%rec_cor(i,1),cmod%rec_cor(i,2),cmod%rec_cor(i,3)],cmod%rec_tet(i))
        !	if (cmod%rec_tet(i)==0) then
        !		print*, "no rec tet No.",i
        !		stop
        !	end if
        !end if
    end do
    !$OMP END PARALLEL DO
    ! print out meshing info
    if(output_level) print*, "--- print out meshing info ..."
    if(output_level) print*,"number of nodes          :", cmod%nod_num
    if(output_level) print*,"number of edges (unknows):", cmod%edg_num
    if(output_level) print*,"number of tets           :", cmod%tet_num
    if(output_level) write(*,"(a,f12.3,8x,a,f12.3)")," upper bound x:",maxval(cmod%nod_cor(:,1))," lower bound x:",minval(cmod%nod_cor(:,1))
    if(output_level) write(*,"(a,f12.3,8x,a,f12.3)")," upper bound y:",maxval(cmod%nod_cor(:,2))," lower bound y:",minval(cmod%nod_cor(:,2))
    if(output_level) write(*,"(a,f12.3,8x,a,f12.3)")," upper bound z:",maxval(cmod%nod_cor(:,3))," lower bound z:",minval(cmod%nod_cor(:,3))
    if(output_level) print*, "--- print out source info ..."
    if(output_level) print*,"number of source edges   :", cmod%seg_num
    do i=1,cmod%src_num
        !source: ind source_type(1 electric, 2 magnetic) pulse_type (1 Turn-on step 2 Gaussian pulse 3 delta pulse) t0 t1 current
        !if(output_level) write(*,"(a,i3,a,i4,a,i4)"), " Source ",i,", Electric Turn-on step source,   Edg Num:",cmod%src_count_edges(i)
        !select case (cmod%src_tp(i,2))
        !case(1)
        !    select case(cmod%src_tp(i,3))
        !    case(1)
        !        write(*,"(a,i3,a,i4,a,i4)"), " Source ",i,": Phy num",cmod%src_tp(i,1),", Electric Turn-on step source,   Edg Num:",cmod%src_count_edges(i)
        !    case(2)
        !        write(*,"(a,i3,a,i4,a,i4)"), " Source ",i,": Phy num",cmod%src_tp(i,1),", Electric Gaussian pulse source, Edg Num:",cmod%src_count_edges(i)
        !    case(3)
        !        write(*,"(a,i3,a,i4,a,i4)"), " Source ",i,": Phy num",cmod%src_tp(i,1),", Electric delta pulse source,    Edg Num:",cmod%src_count_edges(i)
        !    end select
        !case(2)
        !    select case(cmod%src_tp(i,3))
        !    case(1)
        !        write(*,"(a,i3,a,i4,a,i4)"), " Source ",i,": Phy num",cmod%src_tp(i,1),", Magnetic Turn-on step source,   Edg Num:",cmod%src_count_edges(i)
        !    case(2)
        !        write(*,"(a,i3,a,i4,a,i4)"), " Source ",i,": Phy num",cmod%src_tp(i,1),", Magnetic Gaussian pulse source, Edg Num:",cmod%src_count_edges(i)
        !    case(3)
        !        write(*,"(a,i3,a,i4,a,i4)"), " Source ",i,": Phy num",cmod%src_tp(i,1),", Magnetic delta pulse source,    Edg Num:",cmod%src_count_edges(i)
        !    end select
        !end select
    end do
    if(output_level) print*, "--- print out receiver info ..."
    if(output_level) print*,"number of receivers      :", cmod%rec_num
    if(cmod%rec_num>5) then
        do i=1,5
            if(output_level) write(*,"(a,3f16.2,a,i10)") " Coor:",cmod%rec_cor(i,1),cmod%rec_cor(i,2),cmod%rec_cor(i,3), " in tet num:",cmod%rec_tet(i)
        end do
    else
        do i=1,cmod%rec_num
            if(output_level) write(*,"(a,3f16.2,a,i10)") " Coor:",cmod%rec_cor(i,1),cmod%rec_cor(i,2),cmod%rec_cor(i,3), " in tet num:",cmod%rec_tet(i)
        end do
    end if
    if(output_level) write(*,*) "Following field components will be output:",cmod%o_f
    if(output_level) print*, "--- print out domain info ..."
    if (cmod%bound_cond==0) then
        if(output_level) print*, "PEC boundary condition applied"
    else
        if(output_level) print*, "1st absorbing boundary condition applied"
    end if
    do i = 1,cmod%dom_num
        if (cmod%iso) then
            if(output_level) write(*,"(a,i4,a,i8,a)"), " Isotropic domain",cmod%dom_ind(i)," is discretized into", count(cmod%tet_phy==cmod%dom_ind(i)), " tets."
            do j=1,cmod%tet_num
                if (cmod%dom_ind(i)==cmod%tet_phy(j)) then
                    if(output_level) write(*,"(a,e10.3,3x,a,e10.3,3x,a,e10.3)"), " Conductivity:",cmod%sig_3d(j),"Permittivity:",cmod%per_3d(j), "Magnetic permeability:",cmod%mag_3d(j)
                    exit
                end if
            end do
        else
            if(output_level) write(*,"(a,i4,a,i8,a)"), " Anisotropic domain",cmod%dom_ind(i)," is discretized into", count(cmod%tet_phy==cmod%dom_ind(i)), " tets."
        end if
    end do
    end subroutine model_ctl

    subroutine mesh_process (cmod)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    use utility
    USE IFPORT
    implicit none
    type(model),intent(inout) :: cmod
    ! coo format of system matrix
    INTEGER :: a_SIZE,b_SIZE,c_SIZE
    integer,allocatable :: a_ROW(:),a_COL(:),b_ROW(:),B_COL(:),c_ROW(:),c_COL(:)
    REAL(8),ALLOCATABLE :: a_VAL(:),B_VAL(:),c_VAL(:)
    ! other variables
    integer :: i,j,k,info,itet
    real(8) :: a_mat(6,6),b_mat(6,6),c_mat(6,6),w_mat(6,6)
    type(SPARSE_MATRIX_T) :: at_mkl,bt_mkl,ct_mkl
    real(8) :: tmp_charac(3,3)
    ! system matrix size initialization
    ALLOCATE(a_row(36*cmod%tet_NUM),a_col(36*cmod%tet_NUM),a_val(36*cmod%tet_NUM))
    a_SIZE=0; a_row=0; a_col=0; a_val=0.d0
    ALLOCATE(b_row(36*cmod%tet_NUM),B_col(36*cmod%tet_NUM),B_val(36*cmod%tet_NUM))
    b_SIZE=0; b_row=0; b_col=0; b_val=0.d0
    ALLOCATE(c_row(36*cmod%tet_NUM),c_col(36*cmod%tet_NUM),c_val(36*cmod%tet_NUM))
    c_SIZE=0; c_row=0; c_col=0; c_val=0.d0
    ! loop over meshes
    if(output_level) print*, "--- calculate tet matrix ..."
    time_1 = TIMEF()
    do i=1,cmod%tet_NUM
        ! element matrix a_mat, b_mat and c_mat
        if (cmod%iso) then
            call f_kind_iso (a_mat,CMOD,cmod%per_3d(i),i)
            call f_kind_iso (b_mat,CMOD,cmod%sig_3d(i),i)
            call e_kind_iso (c_mat,CMOD,1.d0/cmod%mag_3d(i),i)
        else
            call f_kind_ani (a_mat,CMOD,cmod%per_a3d(i,:,:),1.d0,i)
            call f_kind_ani (b_mat,CMOD,cmod%sig_a3d(i,:,:),1.d0,i)
            do j=1,3
                do k=1,3
                    if (cmod%mag_a3d(i,j,k)==0.d0) then
                        tmp_charac(j,k) = 0.d0
                    else
                        tmp_charac(j,k) = 1.d0/cmod%mag_a3d(i,j,k)
                    end if
                end do
            end do
            call e_kind_ani (c_mat,CMOD,tmp_charac,1.d0,i)
        end if
        if (cmod%bound_cond==0) then
            w_mat = 0.d0
        else
            call surface_integral (w_mat,CMOD,i)
        end if
        b_mat = b_mat + w_mat
        ! assign to global matrix
        DO J=1,6
            DO K=1,6
                ! a_SYS(INDEX_M(I,J),INDEX_M(I,K))=a_SYS(INDEX_M(I,J),INDEX_M(I,K))+a_ELE(j,k)
                if (cmod%tet_edg(I,J)<=cmod%tet_edg(I,K)) then
                    a_SIZE = a_SIZE+1
                    a_ROW(a_SIZE) = cmod%tet_edg(I,J)
                    a_COL(a_SIZE) = cmod%tet_edg(I,K)
                    a_VAL(a_SIZE) = a_mat(j,k) * cmod%tet_drc(i,j) * cmod%tet_drc(i,k)
                    ! B_SYS(INDEX_M(I,J),INDEX_M(I,K))=B_SYS(INDEX_M(I,J),INDEX_M(I,K))+B_ELE(j,k)
                    b_SIZE = b_SIZE+1
                    B_ROW(b_SIZE) = cmod%tet_edg(I,J)
                    B_COL(b_SIZE) = cmod%tet_edg(I,K)
                    B_VAL(b_SIZE) = B_mat(j,k) * cmod%tet_drc(i,j) * cmod%tet_drc(i,k)
                    ! c_SYS(INDEX_M(I,J),INDEX_M(I,K))=c_SYS(INDEX_M(I,J),INDEX_M(I,K))+c_ELE(j,k)
                    c_SIZE = c_SIZE+1
                    c_ROW(c_SIZE) = cmod%tet_edg(I,J)
                    c_COL(c_SIZE) = cmod%tet_edg(I,K)
                    c_VAL(c_SIZE) = c_mat(j,k) * cmod%tet_drc(i,j) * cmod%tet_drc(i,k)
                end if
            END DO
        END DO
    end do
    !! boundary condition
    !if(output_level) print*, "--- boundary condition of matrix A ..."
    !!$omp parallel do private(i,j)
    !do j=1,cmod%sur_edg_num
    !	do i=1,a_SIZE
    !		if ((a_ROW(i)==cmod%sur_edg_index(j).or.a_col(i)==cmod%sur_edg_index(j)).and.(a_ROW(i)/=a_col(i))) a_val(i) = 0.d0
    !	end do
    !end do
    !!$omp end parallel do
    !if(output_level) print*, "--- boundary condition of matrix B ..."
    !!$omp parallel do private(i,j)
    !do j=1,cmod%sur_edg_num
    !	do i=1,b_SIZE
    !		if ((b_ROW(i)==cmod%sur_edg_index(j).or.b_col(i)==cmod%sur_edg_index(j)).and.(b_ROW(i)/=b_col(i))) b_val(i) = 0.d0
    !	end do
    !end do
    !!$omp end parallel do
    !if(output_level) print*, "--- boundary condition of matrix C ..."
    !!$omp parallel do private(i,j)
    !do j=1,cmod%sur_edg_num
    !	do i=1,c_SIZE
    !		if ((c_ROW(i)==cmod%sur_edg_index(j).or.c_col(i)==cmod%sur_edg_index(j)).and.(c_ROW(i)/=c_col(i))) c_val(i) = 0.d0
    !	end do
    !end do
    !!$omp end parallel do
    if(output_level) print*, "--- handle global matrix A ..."
    !call coo_handle(a_ROW,a_COL,a_VAL)
    !call coo_to_sym(.true.,a_ROW,a_cOL,a_VAL)
    !if (.not. symmed) stop'not symmed'
    !na=edge_num
    !nnza=size(a_SYS_VAL)
    !call coocsr ( na, nnza, a_SYS_VAL,a_SYS_ROW,a_SYS_COL, va, ja, ia )
    info = mkl_sparse_d_create_coo (at_mkl, SPARSE_INDEX_BASE_one, cmod%edg_num, cmod%edg_num, a_SIZE, a_ROW(1:a_SIZE), a_COL(1:a_SIZE), a_VAL(1:a_SIZE))
    info = mkl_sparse_convert_csr (at_mkl, SPARSE_OPERATION_NON_TRANSPOSE, cmod%a_mkl)
    info = MKL_SPARSE_OPTIMIZE(cmod%a_mkl)
    info = mkl_sparse_destroy(at_mkl)
    if(output_level) print*, "--- handle global matrix B ..."
    info = mkl_sparse_d_create_coo (bt_mkl, SPARSE_INDEX_BASE_one, cmod%edg_num, cmod%edg_num, b_SIZE, b_ROW(1:b_SIZE), b_COL(1:b_SIZE), b_VAL(1:b_SIZE))
    info = mkl_sparse_convert_csr (bt_mkl, SPARSE_OPERATION_NON_TRANSPOSE, cmod%b_mkl)
    info = MKL_SPARSE_OPTIMIZE(cmod%b_mkl)
    info = mkl_sparse_destroy(bt_mkl)
    if(output_level) print*, "--- handle global matrix C ..."
    info = mkl_sparse_d_create_coo (ct_mkl, SPARSE_INDEX_BASE_one, cmod%edg_num, cmod%edg_num, c_SIZE, c_ROW(1:c_SIZE), c_COL(1:c_SIZE), c_VAL(1:c_SIZE))
    info = mkl_sparse_convert_csr (ct_mkl, SPARSE_OPERATION_NON_TRANSPOSE, cmod%c_mkl)
    info = MKL_SPARSE_OPTIMIZE(cmod%c_mkl)
    info = mkl_sparse_destroy(ct_mkl)
    time_2 = TIMEF()
    time_assemble = time_2-time_1
    ! calculate N and D vector
    allocate(cmod%n_arr(cmod%rec_num,6,3),cmod%d_arr(cmod%rec_num,6,3))
    do i=1,cmod%rec_num
        itet = cmod%rec_tet(i)
        do j=1,6
            ! N vector
            cmod%n_arr(i,j,1) = (cmod%a(itet,ind_tet(j,1))+cmod%b(itet,ind_tet(j,1))*cmod%rec_cor(i,1)+cmod%c(itet,ind_tet(j,1))*cmod%rec_cor(i,2)+cmod%d(itet,ind_tet(j,1))*cmod%rec_cor(i,3))*cmod%b(itet,ind_tet(j,2)) &
                & - (cmod%a(itet,ind_tet(j,2))+cmod%b(itet,ind_tet(j,2))*cmod%rec_cor(i,1)+cmod%c(itet,ind_tet(j,2))*cmod%rec_cor(i,2)+cmod%d(itet,ind_tet(j,2))*cmod%rec_cor(i,3))*cmod%b(itet,ind_tet(j,1))
            cmod%n_arr(i,j,2) = (cmod%a(itet,ind_tet(j,1))+cmod%b(itet,ind_tet(j,1))*cmod%rec_cor(i,1)+cmod%c(itet,ind_tet(j,1))*cmod%rec_cor(i,2)+cmod%d(itet,ind_tet(j,1))*cmod%rec_cor(i,3))*cmod%c(itet,ind_tet(j,2)) &
                & - (cmod%a(itet,ind_tet(j,2))+cmod%b(itet,ind_tet(j,2))*cmod%rec_cor(i,1)+cmod%c(itet,ind_tet(j,2))*cmod%rec_cor(i,2)+cmod%d(itet,ind_tet(j,2))*cmod%rec_cor(i,3)) * cmod%c(itet,ind_tet(j,1))
            cmod%n_arr(i,j,3) = (cmod%a(itet,ind_tet(j,1))+cmod%b(itet,ind_tet(j,1))*cmod%rec_cor(i,1)+cmod%c(itet,ind_tet(j,1))*cmod%rec_cor(i,2)+cmod%d(itet,ind_tet(j,1))*cmod%rec_cor(i,3))*cmod%d(itet,ind_tet(j,2)) &
                & - (cmod%a(itet,ind_tet(j,2))+cmod%b(itet,ind_tet(j,2))*cmod%rec_cor(i,1)+cmod%c(itet,ind_tet(j,2))*cmod%rec_cor(i,2)+cmod%d(itet,ind_tet(j,2))*cmod%rec_cor(i,3))*cmod%d(itet,ind_tet(j,1))
            cmod%n_arr(i,j,:) = cmod%n_arr(i,j,:) * cmod%tet_len(itet,j) / (6.d0*cmod%tet_vol(itet))**2.d0
            ! D vector
            cmod%d_arr(i,j,1) = cmod%c(itet,ind_tet(j,1))*cmod%d(itet,ind_tet(j,2)) - cmod%c(itet,ind_tet(j,2))*cmod%d(itet,ind_tet(j,1))
            cmod%d_arr(i,j,2) = cmod%b(itet,ind_tet(j,2))*cmod%d(itet,ind_tet(j,1)) - cmod%b(itet,ind_tet(j,1))*cmod%d(itet,ind_tet(j,2))
            cmod%d_arr(i,j,3) = cmod%b(itet,ind_tet(j,1))*cmod%c(itet,ind_tet(j,2)) - cmod%b(itet,ind_tet(j,2))*cmod%c(itet,ind_tet(j,1))
            cmod%d_arr(i,j,:) = (-2.d0) * cmod%d_arr(i,j,:) * cmod%tet_len(itet,j) / (6.d0*cmod%tet_vol(itet))**2.d0
        end do
    end do
    end subroutine mesh_process
    
    subroutine d_kind (D_VEC,CMOD,scalar,V_VEC,ind)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    ! calculate D-kind integral Dij = scalar * int (Ni*V) dV of element matrix in a single element
    implicit none
    real(8),intent(out) :: D_VEC(6)
    type(model),intent(in) :: cmod
    real(8),intent(in) :: scalar,V_VEC(3)
    integer,intent(in) :: ind
    INTEGER :: I
    real(8) :: d_x,d_y,d_z
    do i=1,6
        d_x = CMOD%tet_len(ind,i) / (36.d0*CMOD%tet_vol(ind)) * (CMOD%a(ind,ind_tet(i,1))*CMOD%b(ind,ind_tet(i,2)) - CMOD%a(ind,ind_tet(i,2))*CMOD%b(ind,ind_tet(i,1)))
        d_y = CMOD%tet_len(ind,i) / (36.d0*CMOD%tet_vol(ind)) * (CMOD%a(ind,ind_tet(i,1))*CMOD%c(ind,ind_tet(i,2)) - CMOD%a(ind,ind_tet(i,2))*CMOD%c(ind,ind_tet(i,1)))
        d_z = CMOD%tet_len(ind,i) / (36.d0*CMOD%tet_vol(ind)) * (CMOD%a(ind,ind_tet(i,1))*CMOD%d(ind,ind_tet(i,2)) - CMOD%a(ind,ind_tet(i,2))*CMOD%d(ind,ind_tet(i,1)))
        D_VEC(i) = d_x*v_VEC(1) + d_y*v_VEC(2) + d_z*v_VEC(3)
    end do
    D_VEC = scalar * D_VEC
    end subroutine d_kind
    
    subroutine intetvec_n (D_VEC,CMOD,scalar,V_VEC,ind)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    ! calculate D-kind integral Dij = scalar * int (Ni*V) dV of element matrix in a single element
    implicit none
    real(8),intent(out) :: D_VEC(6)
    type(model),intent(in) :: cmod
    real(8),intent(in) :: scalar,V_VEC(3)
    integer,intent(in) :: ind
    INTEGER :: I
    real(8) :: d_x,d_y,d_z
    do i=1,6
        d_x = CMOD%tet_len(ind,i) / (36.d0*CMOD%tet_vol(ind)*CMOD%tet_vol(ind)) * (CMOD%a(ind,ind_tet(i,1))*CMOD%b(ind,ind_tet(i,2)) - CMOD%a(ind,ind_tet(i,2))*CMOD%b(ind,ind_tet(i,1)))
        d_y = CMOD%tet_len(ind,i) / (36.d0*CMOD%tet_vol(ind)*CMOD%tet_vol(ind)) * (CMOD%a(ind,ind_tet(i,1))*CMOD%c(ind,ind_tet(i,2)) - CMOD%a(ind,ind_tet(i,2))*CMOD%c(ind,ind_tet(i,1)))
        d_z = CMOD%tet_len(ind,i) / (36.d0*CMOD%tet_vol(ind)*CMOD%tet_vol(ind)) * (CMOD%a(ind,ind_tet(i,1))*CMOD%d(ind,ind_tet(i,2)) - CMOD%a(ind,ind_tet(i,2))*CMOD%d(ind,ind_tet(i,1)))
        D_VEC(i) = d_x*v_VEC(1) + d_y*v_VEC(2) + d_z*v_VEC(3)
    end do
    D_VEC = scalar * D_VEC
    end subroutine intetvec_n
    
    subroutine intetvec_curln (D_VEC,CMOD,scalar,V_VEC,ind)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    ! calculate D-kind integral Dij = scalar * int (Ni*V) dV of element matrix in a single element
    implicit none
    real(8),intent(out) :: D_VEC(6)
    type(model),intent(in) :: cmod
    real(8),intent(in) :: scalar,V_VEC(3)
    integer,intent(in) :: ind
    INTEGER :: I
    real(8) :: d_x,d_y,d_z
    do i=1,6
        d_x = 2.d0 * CMOD%tet_len(ind,i) / (36.d0*CMOD%tet_vol(ind)*CMOD%tet_vol(ind)) * (CMOD%c(ind,ind_tet(i,1))*CMOD%d(ind,ind_tet(i,2)) - CMOD%c(ind,ind_tet(i,2))*CMOD%d(ind,ind_tet(i,1)))
        d_y = 2.d0 * CMOD%tet_len(ind,i) / (36.d0*CMOD%tet_vol(ind)*CMOD%tet_vol(ind)) * (CMOD%b(ind,ind_tet(i,1))*CMOD%d(ind,ind_tet(i,2)) - CMOD%b(ind,ind_tet(i,2))*CMOD%d(ind,ind_tet(i,1)))
        d_z = 2.d0 * CMOD%tet_len(ind,i) / (36.d0*CMOD%tet_vol(ind)*CMOD%tet_vol(ind)) * (CMOD%b(ind,ind_tet(i,1))*CMOD%c(ind,ind_tet(i,2)) - CMOD%b(ind,ind_tet(i,2))*CMOD%c(ind,ind_tet(i,1)))
        D_VEC(i) = d_x*v_VEC(1) + d_y*v_VEC(2) + d_z*v_VEC(3)
    end do
    D_VEC = scalar * D_VEC
    end subroutine intetvec_curln
    
    subroutine e_kind_iso (e_mat,cmod,scalar,ind)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    ! calculate isotropic E-kind integral E ij = scalar * ∫(▽×Ni)·(▽×Nj)dV of element matrix in a single element
    implicit none
    real(8),intent(out) :: e_mat(6,6)
    type(model),intent(in) :: cmod
    real(8),intent(in) :: scalar
    integer,intent(in) :: ind
    real(8) :: tmp_e1(3),tmp_e2(3)
    INTEGER :: I,J
    do i=1,6
        do j=1,6
            ! first part of Eij
            tmp_e1(1) = CMOD%c(ind,ind_tet(i,1))*CMOD%d(ind,ind_tet(i,2)) - CMOD%c(ind,ind_tet(i,2))*CMOD%d(ind,ind_tet(i,1))
            tmp_e1(2) = CMOD%b(ind,ind_tet(i,2))*CMOD%d(ind,ind_tet(i,1)) - CMOD%b(ind,ind_tet(i,1))*CMOD%d(ind,ind_tet(i,2))
            tmp_e1(3) = CMOD%b(ind,ind_tet(i,1))*CMOD%c(ind,ind_tet(i,2)) - CMOD%b(ind,ind_tet(i,2))*CMOD%c(ind,ind_tet(i,1))
            ! second part of Eij
            tmp_e2(1) = CMOD%c(ind,ind_tet(j,1))*CMOD%d(ind,ind_tet(j,2)) - CMOD%c(ind,ind_tet(j,2))*CMOD%d(ind,ind_tet(j,1))
            tmp_e2(2) = CMOD%b(ind,ind_tet(j,2))*CMOD%d(ind,ind_tet(j,1)) - CMOD%b(ind,ind_tet(j,1))*CMOD%d(ind,ind_tet(j,2))
            tmp_e2(3) = CMOD%b(ind,ind_tet(j,1))*CMOD%c(ind,ind_tet(j,2)) - CMOD%b(ind,ind_tet(j,2))*CMOD%c(ind,ind_tet(j,1))
            ! the integration
            e_mat(i,j) = 4.d0*CMOD%tet_len(ind,i)*CMOD%tet_len(ind,j)*CMOD%tet_vol(ind)/(6.d0*CMOD%tet_vol(ind))**4.d0 * dot_product(tmp_e1,tmp_e2)
        end do
    end do
    e_mat = e_mat * scalar
    end subroutine e_kind_iso
    
    subroutine f_kind_iso (f_mat,cmod,scalar,ind)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    ! calculate isotropic f-kind integral Fij = scalar * ∫ Ni·Nj dV of element matrix in a single element
    use utility
    implicit none
    real(8),intent(out) :: f_mat(6,6)
    type(model),intent(in) :: cmod
    real(8),intent(in) :: scalar
    integer,intent(in) :: ind
    INTEGER :: I,J,m,n
    real(8) :: f(2,2)
    do i=1,6
        do j=1,6
            do m=1,2
                do n=1,2
                    f(m,n) = CMOD%b(ind,ind_tet(i,m))*CMOD%b(ind,ind_tet(j,n))+CMOD%c(ind,ind_tet(i,m))*CMOD%c(ind,ind_tet(j,n))+CMOD%d(ind,ind_tet(i,m))*CMOD%d(ind,ind_tet(j,n))
                end do
            end do
            f_mat(i,j) = f(2,2)*(1+dirac(ind_tet(i,1),ind_tet(j,1)))-f(1,2)*(1+dirac(ind_tet(i,2),ind_tet(j,1)))-f(2,1)*(1+dirac(ind_tet(i,1),ind_tet(j,2)))+f(1,1)*(1+dirac(ind_tet(i,2),ind_tet(j,2)))
            f_mat(i,j) = CMOD%tet_len(ind,i)*CMOD%tet_len(ind,j)/(720.d0*CMOD%tet_vol(ind)) * f_mat(i,j)
        end do
    end do
    f_mat = f_mat * scalar
    end subroutine f_kind_iso
    
    subroutine e_kind_ani (e_mat,CMOD,ani_m,scalar,ind)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    ! calculate anisotropic E-kind integral E ij = scalar * ∫(▽×Ni)·ani_m·(▽×Nj)dV of element matrix in a single element
    use utility
    implicit none
    real(8),intent(out) :: e_mat(6,6)
    type(model),intent(in) :: cmod
    real(8),intent(in) :: ani_m(3,3),scalar
    integer,intent(in) :: ind
    real(8) :: tmp_e1(3),tmp_e2(3)
    INTEGER :: I,J
    do i=1,6
        do j=1,6
            ! first part of Eij
            tmp_e1(1) = (CMOD%c(ind,ind_tet(i,1))*CMOD%d(ind,ind_tet(i,2))-CMOD%c(ind,ind_tet(i,2))*CMOD%d(ind,ind_tet(i,1))) * 2.d0*CMOD%tet_len(ind,i)/(6.d0*CMOD%tet_vol(ind))**2.d0
            tmp_e1(2) = (CMOD%b(ind,ind_tet(i,2))*CMOD%d(ind,ind_tet(i,1))-CMOD%b(ind,ind_tet(i,1))*CMOD%d(ind,ind_tet(i,2))) * 2.d0*CMOD%tet_len(ind,i)/(6.d0*CMOD%tet_vol(ind))**2.d0
            tmp_e1(3) = (CMOD%b(ind,ind_tet(i,1))*CMOD%c(ind,ind_tet(i,2))-CMOD%b(ind,ind_tet(i,2))*CMOD%c(ind,ind_tet(i,1))) * 2.d0*CMOD%tet_len(ind,i)/(6.d0*CMOD%tet_vol(ind))**2.d0
            ! second part of Eij
            tmp_e2(1) = (CMOD%c(ind,ind_tet(j,1))*CMOD%d(ind,ind_tet(j,2))-CMOD%c(ind,ind_tet(j,2))*CMOD%d(ind,ind_tet(j,1))) * 2.d0*CMOD%tet_len(ind,j)/(6.d0*CMOD%tet_vol(ind))**2.d0
            tmp_e2(2) = (CMOD%b(ind,ind_tet(j,2))*CMOD%d(ind,ind_tet(j,1))-CMOD%b(ind,ind_tet(j,1))*CMOD%d(ind,ind_tet(j,2))) * 2.d0*CMOD%tet_len(ind,j)/(6.d0*CMOD%tet_vol(ind))**2.d0
            tmp_e2(3) = (CMOD%b(ind,ind_tet(j,1))*CMOD%c(ind,ind_tet(j,2))-CMOD%b(ind,ind_tet(j,2))*CMOD%c(ind,ind_tet(j,1))) * 2.d0*CMOD%tet_len(ind,j)/(6.d0*CMOD%tet_vol(ind))**2.d0
            tmp_e2 = matmul(ani_m,tmp_e2)
            ! the integration
            e_mat(i,j) = CMOD%tet_vol(ind) * dot_product(tmp_e1,tmp_e2)
        end do
    end do
    e_mat = e_mat * scalar
    end subroutine e_kind_ani
    
    subroutine f_kind_ani (f_mat,CMOD,ani_m,scalar,ind)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    ! calculate anisotropic f-kind integral F ij = scalar * (∫ Ni·ani_m·Nj dV) of element matrix in a single element
    use utility
    implicit none
    real(8),intent(out) :: f_mat(6,6)
    type(model),intent(in) :: cmod
    real(8),intent(in) :: ani_m(3,3),scalar
    integer,intent(in) :: ind
    INTEGER :: I,J,m,n
    real(8) :: f(2,2),tmp_v1(3),tmp_v2(3)
    do i=1,6
        do j=1,6
            do m=1,2
                do n=1,2
                    tmp_v1(1) = CMOD%b(ind,ind_tet(i,m)); tmp_v1(2) = CMOD%c(ind,ind_tet(i,m)); tmp_v1(3) = CMOD%d(ind,ind_tet(i,m))
                    tmp_v2(1) = CMOD%b(ind,ind_tet(j,n)); tmp_v2(2) = CMOD%c(ind,ind_tet(j,n)); tmp_v2(3) = CMOD%d(ind,ind_tet(j,n))
                    tmp_v2 = matmul(ani_m,tmp_v2)
                    f(m,n) = dot_product(tmp_v1,tmp_v2)
                end do
            end do
            f_mat(i,j) = f(2,2)*(1+dirac(ind_tet(i,1),ind_tet(j,1)))-f(1,2)*(1+dirac(ind_tet(i,2),ind_tet(j,1)))-f(2,1)*(1+dirac(ind_tet(i,1),ind_tet(j,2)))+f(1,1)*(1+dirac(ind_tet(i,2),ind_tet(j,2)))
            f_mat(i,j) = CMOD%tet_len(ind,i)*CMOD%tet_len(ind,j)/(CMOD%tet_vol(ind)*720.d0) * f_mat(i,j)
        end do
    end do
    f_mat = f_mat * scalar
    end subroutine f_kind_ani
    
    subroutine surface_integral (w_mat,cmod,ind)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    implicit none
    ! in/out variables
    real(8),intent(out) :: w_mat(6,6)
    type(model),intent(in) :: cmod
    integer,intent(in) :: ind
    ! local viariables
    real(8) :: x(3),y(3),z(3),p_vec(4),f_tri(3,3),test_value
    logical :: on_edge,test_p,test_m
    integer :: i,j,k
    w_mat=0.d0
    ! loop over 4 tris
    do i=1,4
        ! assign coors - loop over 3 points
        do j=1,3
            x(j) = cmod%nod_cor(cmod%tet_vex(ind,ind_tri(i,j)),1)
            y(j) = cmod%nod_cor(cmod%tet_vex(ind,ind_tri(i,j)),2)
            z(j) = cmod%nod_cor(cmod%tet_vex(ind,ind_tri(i,j)),3)
        end do
        ! test if any triangle is on the edge of the model
        call plane(x,y,z,p_vec)
        on_edge= .true.
        test_p = .false.
        test_m = .false.
        do j = 1,cmod%nod_num
            test_value = cmod%nod_cor(j,1)*p_vec(1) + cmod%nod_cor(j,2)*p_vec(2) + cmod%nod_cor(j,3)*p_vec(3) + p_vec(4)
            if (test_value>0.d0) test_p = .true.
            if (test_value<0.d0) test_m = .true.
            if (test_p.and.test_m) then
                on_edge = .false.
                exit
            end if
        end do
        if (.not.on_edge) cycle
        call triangle_int (x,y,z,f_tri)
        do j=1,3
            do k=1,3
                w_mat(ind_edg(i,j),ind_edg(i,k)) = w_mat(ind_edg(i,j),ind_edg(i,k)) + dble(ind_dir(i,j)*ind_dir(i,k)) * f_tri(j,k)
            end do
        end do
    end do
    w_mat = w_mat * sqrt(1.257e-6 / 8.85e-12)
    end subroutine surface_integral
    
    subroutine plane(x,y,z,p_vec)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    ! determine a plane according to given 3 points, p_vec stores 4 coefs of (a,b,c,d) of a plane which has the form of ax+by+cx+d=0
    implicit none
    real(8),intent(in) :: x(3),y(3),z(3)
    real(8),intent(out) :: p_vec(4)
    real(8) :: t_vec1(3),t_vec2(3)
    t_vec1(1) = x(2)-x(1); t_vec1(2) = y(2)-y(1); t_vec1(3) = z(2)-z(1)
    t_vec2(1) = x(3)-x(1); t_vec2(2) = y(3)-y(1); t_vec2(3) = z(3)-z(1)
    p_vec (1) = t_vec1(2)*t_vec2(3) - t_vec1(3)*t_vec2(2)
    p_vec (2) = - ( t_vec1(1)*t_vec2(3) - t_vec1(3)*t_vec2(1) )
    p_vec (3) = t_vec1(1)*t_vec2(2) - t_vec1(2)*t_vec2(1)
    p_vec (4) = - ( p_vec(1)*x(1)+p_vec(2)*y(1)+p_vec(3)*z(1) )
    end subroutine plane
    
    subroutine triangle_int (x,y,z,f_mat)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    implicit none
    real(8),intent(in) :: x(3),y(3),z(3)
    real(8),intent(out) :: f_mat(3,3)
    real(8) :: u(3),v(3)
    u(1) = 0.0
    v(1) = 0.0
    u(2) = sqrt ( (x(2)-x(1))**2.0 + (y(2)-y(1))**2.0 + (z(2)-z(1))**2.0 )
    v(2) = 0.0
    u(3) = ( (x(3)-x(1))*(x(2)-x(1)) + (y(3)-y(1))*(y(2)-y(1)) + (z(3)-z(1))*(z(2)-z(1)) ) / u(2)
    v(3) = sqrt ( 1.0 - ( u(3) / sqrt ( (x(3)-x(1))**2.0 + (y(3)-y(1))**2.0 + (z(3)-z(1))**2.0 ) ) ** 2.0 ) * sqrt ( (x(3)-x(1))**2.0 + (y(3)-y(1))**2.0 + (z(3)-z(1))**2.0 )
    call f_int_2d (u,v,f_mat)
    end subroutine triangle_int
    
    subroutine f_int_2d (x,y,f_mat)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    implicit none
    real(8),intent(in) :: x(3),y(3)
    real(8),intent(out) :: f_mat(3,3)
    integer :: i,j
    real(8) :: area,length(3)
    length(1) = sqrt((x(1)-x(2))**2.d0+(y(1)-y(2))**2.d0)
    length(2) = sqrt((x(2)-x(3))**2.d0+(y(2)-y(3))**2.d0)
    length(3) = sqrt((x(3)-x(1))**2.d0+(y(3)-y(1))**2.d0)
    area = 0.5*(x(1)*(y(2)-y(3))-y(1)*(x(2)-x(3))+1.0*(x(2)*y(3)-x(3)*y(2)))
    do i=1,3
        do j=1,3
            f_mat(i,j) = 6.d0*(y(i)-y(3))*(y(j)-y(3)) - 2.d0*(y(i)+y(j)-2.d0*y(3))*(y(1)+y(2)-2.d0*y(3)) + (y(1)-y(3))**2.d0 + (y(2)-y(3))**2.d0 + (y(1)-y(3))*(y(2)-y(3)) + &
                & 6.d0*(x(i)-x(3))*(x(j)-x(3)) - 2.d0*(x(i)+x(j)-2.d0*x(3))*(x(1)+x(2)-2.d0*x(3)) + (x(1)-x(3))**2.d0 + (x(2)-x(3))**2.d0 + (x(1)-x(3))*(x(2)-x(3))
            f_mat(i,j) = f_mat(i,j) * length(i) * length(j) / 24.d0 / area
        end do
    end do
    end subroutine f_int_2d
    
    subroutine p_in_tet (cmod,p,in_which)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    ! Let the tetrahedron have vertices V1=(x1,y1,z1) V2=(x2,y2,z2) V3=(x3,y3,z3) V4=(x4,y4,z4) and your test point be P = (x, y, z).
    ! Then the point P is in the tetrahedron if following five determinants all have the same sign.
    !   |x1 y1 z1 1|    |x  y  z  1|    |x1 y1 z1 1|    |x1 y1 z1 1|    |x1 y1 z1 1|
    !D0=|x2 y2 z2 1| D1=|x2 y2 z2 1| D2=|x  y  z  1| D3=|x2 y2 z2 1| D4=|x2 y2 z2 1|
    !   |x3 y3 z3 1|    |x3 y3 z3 1|    |x3 y3 z3 1|    |x  y  z  1|    |x3 y3 z3 1|
    !   |x4 y4 z4 1|    |x4 y4 z4 1|    |x4 y4 z4 1|    |x4 y4 z4 1|    |x  y  z  1|
    use utility
    implicit none
    type(model),intent(in) :: cmod
    real(8),intent(in) :: p(:)
    integer,intent(out) :: in_which
    real(8) :: x(4),y(4),z(4),d(5),mid(3)
    integer :: i,j
    ! process
    in_which = 0
    do i = 1,cmod%tet_num
        do j=1,4
            x(j) = cmod%nod_cor(cmod%tet_vex(i,j),1)
            y(j) = cmod%nod_cor(cmod%tet_vex(i,j),2)
            z(j) = cmod%nod_cor(cmod%tet_vex(i,j),3)
        end do
        mid(1)=(x(1)+x(2)+x(3)+x(4))/4.d0
        mid(2)=(y(1)+y(2)+y(3)+y(4))/4.d0
        mid(3)=(z(1)+z(2)+z(3)+z(4))/4.d0
        if (sqrt((mid(1)-p(1))**2.d0+(mid(2)-p(2))**2.d0+(mid(3)-p(3))**2.d0)>cmod%rec_charac) cycle
        d(1) = det_mkl_4(reshape([x(1),x(2),x(3),x(4),y(1),y(2),y(3),y(4),z(1),z(2),z(3),z(4),1.d0,1.d0,1.d0,1.d0],[4,4]))
        d(2) = det_mkl_4(reshape([p(1),x(2),x(3),x(4),p(2),y(2),y(3),y(4),p(3),z(2),z(3),z(4),1.d0,1.d0,1.d0,1.d0],[4,4]))
        d(3) = det_mkl_4(reshape([x(1),p(1),x(3),x(4),y(1),p(2),y(3),y(4),z(1),p(3),z(3),z(4),1.d0,1.d0,1.d0,1.d0],[4,4]))
        d(4) = det_mkl_4(reshape([x(1),x(2),p(1),x(4),y(1),y(2),p(2),y(4),z(1),z(2),p(3),z(4),1.d0,1.d0,1.d0,1.d0],[4,4]))
        d(5) = det_mkl_4(reshape([x(1),x(2),x(3),p(1),y(1),y(2),y(3),p(2),z(1),z(2),z(3),p(3),1.d0,1.d0,1.d0,1.d0],[4,4]))
        if ( (d(1)>=0.d0.and.d(2)>=0.d0.and.d(3)>=0.d0.and.d(4)>=0.d0.and.d(5)>=0.d0).or.(d(1)<=0.d0.and.d(2)<=0.d0.and.d(3)<=0.d0.and.d(4)<=0.d0.and.d(5)<=0.d0) ) then
            in_which = i
            exit
        end if
    end do
    !if (in_which==0) stop 'no rec tet'
    end subroutine p_in_tet
    !
    !subroutine init_calculation (cmod,init_e)
    !use pardiso_direct_real
    !implicit none
    !type(model),intent(in) :: cmod
    !real(8),ALLOCATABLE,intent(out) :: init_e(:)
    !real(8),ALLOCATABLE :: init_FAI(:)
    !integer :: i,j,k,INFO,itet
    !real(8) :: loc_m(4,4),loc_vec(3)
    !integer :: m_size/ cmod%tet_vol(i) * (cmod%b(i,j)*loc_vec(1)+cmod%c(i,j)*loc_vec(2)+cmod%d(i,j)*loc_vec(3))
    !			end if
    !		end do
    !	end do
    !	! assign loc_m to global matrix
    !	do j=1,4
    !		do k=1,4
    !			if (cmod%tet_vex(I,J)<=cmod%tet_vex(I,K)) then
    !				m_SIZE = m_SIZE+1
    !				m_ROW(m_SIZE) = cmod%tet_vex(I,J)
    !				m_COL(m_SIZE) = cmod%tet_vex(I,K)
    !				m_VAL(m_SIZE) = loc_m(j,k)
    !			end if
    !		end do
    !	end dotet_vex(i,ind_tet(j,1)),cmod%tet_vex(i,ind_tet(j,2))]),2) - cmod%nod_cor(minval([cmod%tet_vex(i,ind_tet(j,1)),cmod%tet_vex(i,ind_tet(j,2))]),2)
    !!        loc_vec(3) = cmod%nod_cor(maxval([cmod%tet_vex(i,ind_tet(j,1)),cmod%tet_vex(i,ind_tet(j,2))]),3) - cmod%nod_cor(minval([cmod%tet_vex(i,ind_tet(j,1)),cmod%tet_vex(i,ind_tet(j,2))]),3)
    !!        phi_edge(i,j) = d
    !info = mkl_sparse_destroy(mt_mkl)
    !info = mkl_sparse_destroy(m_mkl)
    !!if(output_level) print*, "    finished generating DC solution at step 0"
    !end subroutine init_calculation

    subroutine time_solver (cmod)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    USE IFPORT
    use mpi
    use mpi_parameters
    implicit none
    type(model),intent(inout) :: cmod
    real(8),allocatable :: res_e(:)
    integer :: i
    if(rank==0)	then
        ! output
        if(output_level.and.rank==0) print*, "--- init time variables ..."
        ! initialize step_init
        cmod%step_init = cmod%src_pulsewidth / cmod%step_dvsi
        ! how many time step should be performed
        DESCR%TYPE = SPARSE_MATRIX_TYPE_SYMMETRIC
        DESCR%mode = SPARSE_FILL_MODE_UPPER
        DESCR%diag = SPARSE_DIAG_NON_UNIT
        time_frac = 0.d0
        time_stepping = 0.d0
        time_1 = timef()
    end if
    ! send cmod%step_mthd to all MPI processes
    !call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !if(rank==0)	then
    !	do i = 1, size - 1
    !		call MPI_SEND (cmod%step_mthd, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, ierr)
    !	enddo
    !else
    !	call MPI_RECV (cmod%step_mthd, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpi_stat,ierr)
    !end if
    call MPI_BCAST(cmod%step_mthd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(cmod%src_pulsetype, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    ! time loop subroutine calling
    if (cmod%step_mthd==1) then
        !call unif_newmark (cmod)
    else if (cmod%step_mthd==2 .and. cmod%src_pulsetype==3) then
        ! .true. for output res_e, .false. for input
        call adap_euler (cmod,res_e,.true.)
        call adap_euler (cmod,res_e,.false.)
    else if (cmod%step_mthd==2) then
        call adap_euler (cmod)
    else
        stop "error!"
    end if
    if(rank==0)	then
        time_2 = timef()
        if(output_level) write(*,"(a,f12.2,a)")," time of assembling: ",time_assemble," seconds."
        if(output_level) write(*,"(a,f12.2,a)")," time of fractorization: ",time_frac," seconds."
        if(output_level) write(*,"(a,f12.2,a)")," time of stepping: ",time_stepping," seconds."
    end if
    end subroutine time_solver

    !subroutine unif_newmark (cmod)
    !USE IFPORT
    !use pardiso_direct_real
    !use utility
    !implicit none
    !type(model),intent(inout) :: cmod
    !! system matrix

    subroutine adap_euler (cmod,res_e,inout_judge)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    USE IFPORT
    use pardiso_multi_real_mpi
    use utility
    use mpi
    use mpi_parameters
    implicit none
    type(model),intent(inout) :: cmod
    real(8),allocatable,intent(inout),optional :: res_e(:)
    logical,optional :: inout_judge
    ! system matrix
    integer :: nsys,nnzsys,tmp_nsys,tmp_nnzsys
    integer,allocatable :: isys(:),jsys(:),tmp_isys(:),tmp_jsys(:)
    real(8),allocatable :: vsys(:),tmp_vsys(:)
    ! mkl matrix
    type(SPARSE_MATRIX_T) :: p_mkl,tmp_p_mkl
    ! loop variables
    integer :: info,ind_pardiso,tmp_pardiso
    ! e temp vectors
    real(8),allocatable :: e1(:),e2(:),e3(:),e4(:),init_e(:),r_h(:),r_h_1(:)
    real(8),allocatable :: e_veri_1(:),e_veri_2(:), tmp_e1(:),tmp_e2(:),tmp_e3(:),e_record_1(:)
    ! time variables
    real(8) :: time_cur,step_cur,tmp_time_cur,tmp_step_cur,error,time_limit
    integer :: step,not_changed,i,all_e_index
    real(8),allocatable :: all_e(:,:)
    integer :: pardiso_deallo_index
    ! variable initiation
    if(rank==0) then
        step = 0
        time_cur = 0.d0
        not_changed= 0
        all_e_index = 0
        step_cur = cmod%step_init
        open (opt_unit,name=trim(opt_file))
        if (allocated(cmod%time_channels)) open (opt2_unit,name=trim(opt2_file))
        if(output_level) print*, "--- start time loop ..."
        !if (present(res_e).and.inout_judge) then
        !    time_limit = 500.d0
        !else
        !    time_limit = cmod%time_maxm
        !end if
        time_limit = cmod%time_maxm
    end if
    ! time loop start
    call MPI_BCAST(step, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(time_cur, 1, MPI_real8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(time_limit, 1, MPI_real8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(cmod%step_dble, 1, MPI_real8, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(cmod%step_dblesize, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
    do while (time_cur<time_limit)
        step = step + 1
        time_cur = time_cur + step_cur
        if (step==1) then
            if(rank==0) then
                allocate(all_e(cmod%step_dble+10,cmod%edg_num))
                ! SYSTEM MATRIX
                call mkl_add_3 (cmod%edg_num, 0.d0,cmod%a_mkl, 3.d0, cmod%b_mkl, 2.d0*step_cur, cmod%c_mkl, p_mkl)
                call mkl_to_csr (p_mkl,nsys,nnzsys,vsys,isys,jsys)
                info = mkl_sparse_destroy (p_mkl)
                ! e and r_h variables
                allocate(e1(nsys),e2(nsys),e3(nsys),e4(nsys),r_h(nsys),r_h_1(nsys))
                e1=0.d0; e2=0.d0; e3=0.d0; e4=0.d0; r_h=0.d0; r_h_1=0.d0
                ! output at 0 step
                if (cmod%src_pulsetype==3 .and. (.not.inout_judge)) then
                    e2 = res_e
                    e3 = res_e
                    e4 = res_e
                end if
                if (cmod%src_pulsetype/=3) call solution_output (cmod,0,time_cur-step_cur,step_cur,e2)
                if (cmod%src_pulsetype==3 .and. (.not.inout_judge)) call solution_output (cmod,0,time_cur-step_cur,step_cur,e2)
            end if
            ! PARDISO AND system matrix
            call pardiso_alloc(2)
            ind_pardiso = 1
            if(rank==0) time_3=timef()
            call pardiso_direct_1 (ind_pardiso,.true.,nsys,nnzsys,isys,jsys,vsys)
            if(rank==0) time_4=timef()
            if(rank==0) time_frac = time_frac+time_4-time_3
            if(output_level.and.rank==0) write(*,"(a,i6,a,d10.3,a,d10.3)") "--- fractoratization finished at step: ",step,", current time: ",time_cur,", current step: ",step_cur
        end if
        if ( mod(step,cmod%step_dble)==0 ) then
            if(rank==0)	then
                if (allocated(e_veri_1)) deallocate(e_veri_1,e_veri_2,e_record_1)
                allocate(e_veri_1(nsys),e_veri_2(nsys),e_record_1(nsys))
            end if
            ! calculate e_veri_1 (1deltat)
            if(rank==0)	then
                if (allocated(tmp_e1)) deallocate(tmp_e1,tmp_e2,tmp_e3)
                allocate(tmp_e1(nsys),tmp_e2(nsys),tmp_e3(nsys))
                tmp_e2 = e2
                tmp_e3 = e3
                tmp_time_cur = time_cur
                tmp_step_cur = step_cur
            end if
            do i=1,cmod%step_dblesize
                if(rank==0)	then
                    call sys_right (step,tmp_time_cur,cmod,r_h)
                    if (cmod%src_pulsetype==3 .and. inout_judge) r_h = -r_h
                    info = mkl_sparse_d_mv (SPARSE_OPERATION_NON_TRANSPOSE, 1.d0, cmod%b_mkl, DESCR, 4.d0*tmp_e2-tmp_e3 , 0.d0, r_h_1)
                    r_H = r_h_1 - 2.d0 * tmp_step_cur * r_h
                    time_3=timef()
                end if
                call pardiso_direct_2 (ind_pardiso,r_h,tmp_e1,isys,jsys,vsys)
                if(rank==0) then
                    time_4=timef()
                    time_stepping=time_stepping+time_4-time_3
                end if
                if (i==cmod%step_dblesize) exit
                if(rank==0)	then
                    if (i==1) e_record_1 = tmp_e1
                    tmp_e3 = tmp_e2
                    tmp_e2 = tmp_e1
                    tmp_time_cur = tmp_time_cur + tmp_step_cur
                end if
            end do
            if(rank==0)	then
                e_veri_1 = tmp_e1
            end if
            ! calculate e_veri_2 (2deltat)
            if(rank==0)	then
                tmp_step_cur = step_cur * cmod%step_dblesize
                tmp_time_cur = time_cur - STEP_CUR + tmp_step_cur
                call mkl_add_3 (cmod%edg_num, 0.d0,cmod%a_mkl, 3.d0, cmod%b_mkl, 2.d0*tmp_step_cur, cmod%c_mkl, tmp_p_mkl)
                call mkl_to_csr (tmp_p_mkl,tmp_nsys,tmp_nnzsys,tmp_vsys,tmp_isys,tmp_jsys)
                info = mkl_sparse_destroy (tmp_p_mkl)
                time_3=timef()
            end if
            if (ind_pardiso==1) tmp_pardiso=2
            if (ind_pardiso==2) tmp_pardiso=1
            call pardiso_direct_1 (tmp_pardiso,.true.,tmp_nsys,tmp_nnzsys,tmp_isys,tmp_jsys,tmp_vsys)
            if(rank==0)	then
                time_4=timef()
                time_frac = time_frac+time_4-time_3
                call sys_right (step,tmp_time_cur,cmod,r_h)
                if (cmod%src_pulsetype==3 .and. inout_judge) r_h = -r_h
                info = mkl_sparse_d_mv (SPARSE_OPERATION_NON_TRANSPOSE, 1.d0, cmod%b_mkl, DESCR, 4.d0*all_e(all_e_index,:)-all_e(all_e_index-cmod%step_dblesize,:) , 0.d0, r_h_1)
                r_H = r_h_1 - 2.d0 * tmp_step_cur * r_h
                time_3=timef()
            end if
            call pardiso_direct_2 (tmp_pardiso,r_h,e_veri_2,tmp_isys,tmp_jsys,tmp_vsys)
            ! error
            if(rank==0)	then
                time_4=timef()
                time_stepping=time_stepping+time_4-time_3
                error = two_error (e_veri_2,e_veri_1)
                if(output_level) write(*,"(a,i6,a,d10.3,a,d10.3)") "--- fractoratization finished at step: ",step,", current time: ",time_cur,", current step: ",step_cur
            end if
            ! judge and update
            if(rank==0)	then
                !if (error>1.0e-4.and.not_changed<10.and.job=="f") then
                if (error>1.0e-4.and.not_changed<2.and.job=="f") then
                    if(output_level) print*,"error:",error,"keep step"
                    deallocate(tmp_isys,tmp_jsys,tmp_vsys)
                    pardiso_deallo_index = tmp_pardiso
                    !call pardiso_direct_3(tmp_pardiso)
                    e1 = e_record_1
                    if (cmod%src_pulsetype/=3) call solution_output (cmod,step,time_cur,step_cur,e1)
                    if (cmod%src_pulsetype==3 .and. (.not.inout_judge)) call solution_output (cmod,step,time_cur,step_cur,e1)
                    e4 = e3
                    e3 = e2
                    e2 = e1
                    not_changed = not_changed + 1
                    all_e_index = 0
                    !cycle
                else
                    if(output_level) print*,"error:",error,"change step"
                    step_cur = tmp_step_cur
                    time_cur = tmp_time_cur
                    pardiso_deallo_index = ind_pardiso
                    !call pardiso_direct_3(ind_pardiso)
                    ind_pardiso = tmp_pardiso
                    deallocate(isys,jsys,vsys)
                    allocate(isys,source=tmp_isys)
                    allocate(jsys,source=tmp_jsys)
                    allocate(vsys,source=tmp_vsys)
                    e4 = all_e(all_e_index-cmod%step_dblesize,:)
                    e3 = all_e(all_e_index,:)
                    e2 = e_veri_1
                    if (cmod%src_pulsetype/=3) call solution_output (cmod,step,time_cur,step_cur,e2)
                    if (cmod%src_pulsetype==3 .and. (.not.inout_judge)) call solution_output (cmod,step,time_cur,step_cur,e2)
                    not_changed = 0
                    all_e_index = 0
                    !cycle
                end if
            end if
            call MPI_BCAST(ind_pardiso, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(pardiso_deallo_index, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)
            call pardiso_direct_3(pardiso_deallo_index)
            cycle
        end if
        ! sys_right and solve
        if (rank==0) then
            call sys_right (step,time_cur,cmod,r_h)
            if (cmod%src_pulsetype==3 .and. inout_judge) r_h = -r_h
            info = mkl_sparse_d_mv (SPARSE_OPERATION_NON_TRANSPOSE, 1.d0, cmod%b_mkl, DESCR, 4.d0*e2-e3 , 0.d0, r_h_1)
            r_H = r_h_1 - 2.d0 * step_cur * r_h
            time_3=timef()
        end if
        call pardiso_direct_2 (ind_pardiso,r_h,e1,isys,jsys,vsys)
        ! output
        if (rank==0) then
            time_4=timef()
            time_stepping=time_stepping+time_4-time_3
            if (cmod%src_pulsetype/=3) call solution_output (cmod,step,time_cur,step_cur,e1)
            if (cmod%src_pulsetype==3 .and. (.not.inout_judge)) call solution_output (cmod,step,time_cur,step_cur,e1)
            all_e_index = all_e_index+1
            all_e(all_e_index,:) = e1
            e4 = e3
            e3 = e2
            e2 = e1
        end if
    end do
    ! deallocate Pardiso variables
    call pardiso_direct_3(ind_pardiso)
    call pardiso_dealloc
    ! post-iteration processes
    if (rank==0) then
        close(opt_unit)
        if (allocated(cmod%time_channels)) close(opt2_unit)
        if(output_level) print*,"total steps:",step
        if (cmod%src_pulsetype==3 .and. inout_judge) then
            !error = two_error (e1,e2)
            !print*,step,time_cur,error
            !if (error<1.d-7.and.time_cur>cmod%src_pulsewidth) then
            if(output_level) write(*,"(a,e10.2,a,i6)"),"finished DC solution calculation at time",time_cur," and step",step
            allocate(res_e(cmod%edg_num))
            res_e = e1
            !return
            !end if
        end if
    end if
    end subroutine adap_euler

    subroutine sys_right (step,t,cmod,r_h)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    use source_mod
    implicit none
    ! in/out variables
    integer,intent(in) :: step
    real(8),intent(in) :: t
    type(model),intent(in) :: cmod
    real(8),intent(out),allocatable :: r_h(:)
    ! local variables
    real(8) :: didt
    integer :: i,j
    logical :: find_phy
    ! process
    allocate(r_h(cmod%edg_num))
    r_h = 0.d0
    if ( step >= cmod%step_dvsi ) return
    !do i=1,cmod%seg_num
    !    find_phy = .false.
    !    do j=1,cmod%src_num
    !        if (cmod%src_phy(i)==cmod%src_tp(j)) then
    !            find_phy = .true.
    !            call time_const_p (cmod%src_sourcetype,cmod%src_pulsetype,t,cmod%src_pulsewidth,cmod%src_current,didt)
    !            ed_num = cmod%src_count_edges(j)
    !        end if
    !    end do
    !    print*,didt * cmod%src_drc(i) / ed_num
    !    r_h(cmod%src_edg(i)) = didt * cmod%src_drc(i) / ed_num !* cmod%edg_len(cmod%src_edg(i))
    !    if (.not.find_phy) stop("error - physical domain not found!")
    !end do
    do i=1,cmod%seg_num
        !find_phy = .false.
        !do j=1,cmod%src_num
        !	if (cmod%src_phy(i)==cmod%src_tp(j)) then
        !		find_phy = .true.
        call time_const_p (cmod%src_sourcetype,cmod%src_pulsetype,t,cmod%src_pulsewidth,cmod%src_current,didt)
        !	end if
        !end do
        !if (.not.find_phy) stop("error - physical domain not found!")
        r_h(cmod%src_edg(i)) = didt * cmod%src_drc(i) * cmod%edg_len(cmod%src_edg(i)) !/ cmod%seg_num !* cmod%edg_len(cmod%src_edg(i))
    end do
    !r_h = r_h * 8000.d0
    !if(cmod%src_pulsetype==2) r_h = r_h / cmod%src_pulsewidth / cmod%src_num * 3.d0
    end subroutine sys_right

    subroutine solution_output (cmod,step,t,step_cur,res_e)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    use system_lib
    use utility
    implicit none
    ! in/out variables
    type(model),intent(inout) :: cmod
    integer,intent(in) :: step
    real(8),intent(in) :: t,step_cur
    real(8),allocatable,intent(in) :: res_e(:)
    ! local viariables
    integer :: i,j,k,itet
    real(8) :: e_p(cmod%rec_num,3),d_p(cmod%rec_num,3)
    real(8),allocatable :: a_ft(:,:,:),step_record_t(:),e_record_t(:,:),time_record_t(:)
    ! loop over receivers
    e_p = 0.d0
    d_p = 0.d0
    do i=1,cmod%rec_num
        itet = cmod%rec_tet(i)
        do j=1,6
            do k=1,3
                e_p(i,k) = e_p(i,k) + cmod%n_arr(i,j,k) * real(cmod%tet_drc(itet,j)) * res_e(cmod%tet_edg(itet,j))
                d_p(i,k) = d_p(i,k) + cmod%d_arr(i,j,k) * cmod%tet_drc(itet,j) * res_e(cmod%tet_edg(itet,j))
            end do
        end do
    end do
    if (step==0) d_p=0.d0
    ! save field to global field variables - allocate stage
    if (.not.allocated(cmod%a_f)) then
        cmod%init_size_f = 100
        cmod%ind_f = 0
        allocate (cmod%a_f_channel(cmod%time_channel_num,cmod%rec_num,9))
        allocate (cmod%a_f(cmod%init_size_f,cmod%rec_num,9))
        allocate (cmod%time_record(cmod%init_size_f))
        allocate (cmod%step_record(cmod%init_size_f),cmod%e_record(cmod%init_size_f,cmod%edg_num))
        ! a_prime
        allocate(cmod%a_prime(6*cmod%rec_num,cmod%edg_num))
        cmod%a_prime = 0.d0
        do i=1,cmod%rec_num
            itet = cmod%rec_tet(i)
            do j=1,6
                cmod%a_prime((i-1)*6+1,cmod%tet_edg(itet,j)) = cmod%n_arr(i,j,1) * cmod%tet_drc(itet,j)
                cmod%a_prime((i-1)*6+2,cmod%tet_edg(itet,j)) = cmod%n_arr(i,j,2) * cmod%tet_drc(itet,j)
                cmod%a_prime((i-1)*6+3,cmod%tet_edg(itet,j)) = cmod%n_arr(i,j,3) * cmod%tet_drc(itet,j)
                cmod%a_prime((i-1)*6+4,cmod%tet_edg(itet,j)) = cmod%d_arr(i,j,1) * cmod%tet_drc(itet,j)
                cmod%a_prime((i-1)*6+5,cmod%tet_edg(itet,j)) = cmod%d_arr(i,j,2) * cmod%tet_drc(itet,j)
                cmod%a_prime((i-1)*6+6,cmod%tet_edg(itet,j)) = cmod%d_arr(i,j,3) * cmod%tet_drc(itet,j)
            end do
        end do
    end if
    if (cmod%ind_f==cmod%init_size_f) then
        allocate (a_ft(cmod%init_size_f,cmod%rec_num,9))
        a_ft = cmod%a_f
        deallocate(cmod%a_f)
        allocate (step_record_t(cmod%init_size_f))
        step_record_t = cmod%step_record
        deallocate(cmod%step_record)
        allocate (time_record_t(cmod%init_size_f))
        time_record_t = cmod%time_record
        deallocate(cmod%time_record)
        allocate (e_record_t(cmod%init_size_f,cmod%edg_num))
        e_record_t = cmod%e_record
        deallocate(cmod%e_record)
        cmod%init_size_f = cmod%init_size_f*2
        allocate (cmod%a_f(cmod%init_size_f,cmod%rec_num,9))
        cmod%a_f(1:cmod%ind_f,:,:) = a_ft
        allocate(cmod%step_record(cmod%init_size_f))
        cmod%step_record(1:cmod%ind_f) = step_record_t
        allocate(cmod%time_record(cmod%init_size_f))
        cmod%time_record(1:cmod%ind_f) = time_record_t
        allocate(cmod%e_record(cmod%init_size_f,cmod%edg_num))
        cmod%e_record(1:cmod%ind_f,:) = e_record_t
    end if
    ! save field to global field variables - save stage
    cmod%ind_f = cmod%ind_f + 1
    ! step record and edge electric field record
    cmod%time_record(cmod%ind_f) = t
    cmod%step_record(cmod%ind_f) = step_cur
    cmod%e_record(cmod%ind_f,:) = res_e
    !if (job=='m'.or.job=='j'.or.job=='t') then
    !    open(unit=1024,file="tmp_res_e_"//itoa(step))
    !    do i=1,cmod%edg_num
    !	    write(1024,*) res_e(i)
    !    end do
    !    close(1024)
    !end if
    ! field record
    do i=1,cmod%rec_num
        ! e field
        cmod%a_f(cmod%ind_f,i,1) = e_p(i,1)
        cmod%a_f(cmod%ind_f,i,2) = e_p(i,2)
        cmod%a_f(cmod%ind_f,i,3) = e_p(i,3)
        ! d field
        cmod%a_f(cmod%ind_f,i,4) = d_p(i,1)
        cmod%a_f(cmod%ind_f,i,5) = d_p(i,2)
        cmod%a_f(cmod%ind_f,i,6) = d_p(i,3)
        ! b field
        if (cmod%ind_f==1) then
            cmod%a_f(cmod%ind_f,i,7) = d_p(i,1)*step_cur
            cmod%a_f(cmod%ind_f,i,8) = d_p(i,2)*step_cur
            cmod%a_f(cmod%ind_f,i,9) = d_p(i,3)*step_cur
        else
            cmod%a_f(cmod%ind_f,i,7) = cmod%a_f(cmod%ind_f-1,i,7) + d_p(i,1)*step_cur
            cmod%a_f(cmod%ind_f,i,8) = cmod%a_f(cmod%ind_f-1,i,8) + d_p(i,2)*step_cur
            cmod%a_f(cmod%ind_f,i,9) = cmod%a_f(cmod%ind_f-1,i,9) + d_p(i,3)*step_cur
        end if
    end do
    ! output to opt_unit
    write(opt_unit,"(i6,e15.5,2x,"//itoa(cmod%rec_num)//"("//itoa(count(cmod%o_f))//"e15.5,2x))") step,t,(cmod%a_f(cmod%ind_f,i,cmod%o_v)/cmod%sig_parameter**2.d0*pi**log10(cmod%sig_parameter),i=1,cmod%rec_num)
    ! extract cmod%a_f to cmod%a_f_channel, and output to opt2_unit
    if (allocated(cmod%time_channels)) then
        do i=1,cmod%time_channel_num
            if (t>=cmod%time_channels(i).and.cmod%time_channel_ok(i)==.false.) then
                if (cmod%ind_f==1) then
                    cmod%a_f_channel(i,:,:) = cmod%a_f(cmod%ind_f,:,:)
                    exit
                end if
                do j=1,cmod%rec_num
                    do k=1,9
                        cmod%a_f_channel(i,j,k) = interp_lin(cmod%time_record(cmod%ind_f-1),cmod%a_f(cmod%ind_f-1,j,k),cmod%time_record(cmod%ind_f),cmod%a_f(cmod%ind_f,j,k),cmod%time_channels(i))
                    end do
                end do
                write(opt2_unit,"(i6,e15.5,2x,"//itoa(cmod%rec_num)//"("//itoa(count(cmod%o_f))//"e15.5,2x))") i,cmod%time_channels(i),(cmod%a_f_channel(i,j,cmod%o_v)/cmod%sig_parameter**2.d0*pi**log10(cmod%sig_parameter),j=1,cmod%rec_num)
                cmod%time_channel_ok(i) = .true.
            end if
        end do
    end if
    ! vtk_output controls
    if (cmod%vtk_num>0) then
        if (.not.all(cmod%vtk_done)) then
            call vtk_field_save (cmod,step,step_cur,res_e)
            do i=1,cmod%vtk_num
                if (cmod%vtk_done(i)/=.true. .and. t>=cmod%vtk_times(i)) then
                    if(output_level) print*, "generate vtk output file at time: ",t
                    call msh_output (cmod,i)
                    !call vtk_output (cmod,step)
                    if(output_level) print*, "generate vtk output file done at time: ",t
                    cmod%vtk_done(i) = .true.
                end if
            end do
        end if
    end if
    end subroutine solution_output

    subroutine vtk_field_save (cmod,step,step_cur,res_e)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    implicit none
    type(model) :: cmod
    integer,intent(in) :: step
    real(8),intent(in) :: step_cur
    real(8),allocatable,intent(in) :: res_e(:)
    integer :: i,j,k
    real(8) :: cp(3)
    if (step==0) then
        allocate(cmod%n_arr_global(cmod%tet_num,6,3),cmod%d_arr_global(cmod%tet_num,6,3))
        allocate(cmod%e_field_global(cmod%tet_num,3),cmod%d_field_global(cmod%tet_num,3),cmod%b_field_global(cmod%tet_num,3))
        cmod%b_field_global = 0.d0
        do i=1,cmod%tet_NUM
            ! assign coordinates
            cp(1) = (cmod%nod_cor(cmod%tet_VEX(i,1),1)+cmod%nod_cor(cmod%tet_VEX(i,2),1)+cmod%nod_cor(cmod%tet_VEX(i,3),1)+cmod%nod_cor(cmod%tet_VEX(i,4),1))/4.0
            cp(2) = (cmod%nod_cor(cmod%tet_VEX(i,1),2)+cmod%nod_cor(cmod%tet_VEX(i,2),2)+cmod%nod_cor(cmod%tet_VEX(i,3),2)+cmod%nod_cor(cmod%tet_VEX(i,4),2))/4.0
            cp(3) = (cmod%nod_cor(cmod%tet_VEX(i,1),3)+cmod%nod_cor(cmod%tet_VEX(i,2),3)+cmod%nod_cor(cmod%tet_VEX(i,3),3)+cmod%nod_cor(cmod%tet_VEX(i,4),3))/4.0
            ! calculate n and d array of this tet
            do j=1,6
                cmod%n_arr_global(i,j,1) = (cmod%a(i,ind_tet(j,1))+cmod%b(i,ind_tet(j,1))*cp(1)+cmod%c(i,ind_tet(j,1))*cp(2)+cmod%d(i,ind_tet(j,1))*cp(3))*cmod%b(i,ind_tet(j,2)) &
                    & - (cmod%a(i,ind_tet(j,2))+cmod%b(i,ind_tet(j,2))*cp(1)+cmod%c(i,ind_tet(j,2))*cp(2)+cmod%d(i,ind_tet(j,2))*cp(3))*cmod%b(i,ind_tet(j,1))
                cmod%n_arr_global(i,j,2) = (cmod%a(i,ind_tet(j,1))+cmod%b(i,ind_tet(j,1))*cp(1)+cmod%c(i,ind_tet(j,1))*cp(2)+cmod%d(i,ind_tet(j,1))*cp(3))*cmod%c(i,ind_tet(j,2)) &
                    & - (cmod%a(i,ind_tet(j,2))+cmod%b(i,ind_tet(j,2))*cp(1)+cmod%c(i,ind_tet(j,2))*cp(2)+cmod%d(i,ind_tet(j,2))*cp(3)) * cmod%c(i,ind_tet(j,1))
                cmod%n_arr_global(i,j,3) = (cmod%a(i,ind_tet(j,1))+cmod%b(i,ind_tet(j,1))*cp(1)+cmod%c(i,ind_tet(j,1))*cp(2)+cmod%d(i,ind_tet(j,1))*cp(3))*cmod%d(i,ind_tet(j,2)) &
                    & - (cmod%a(i,ind_tet(j,2))+cmod%b(i,ind_tet(j,2))*cp(1)+cmod%c(i,ind_tet(j,2))*cp(2)+cmod%d(i,ind_tet(j,2))*cp(3))*cmod%d(i,ind_tet(j,1))
                cmod%n_arr_global(i,j,:) = cmod%n_arr_global(i,j,:) * cmod%tet_len(i,j) / (6.d0*cmod%tet_vol(i))**2.d0
                cmod%d_arr_global(i,j,1) = cmod%c(i,ind_tet(j,1))*cmod%d(i,ind_tet(j,2)) - cmod%c(i,ind_tet(j,2))*cmod%d(i,ind_tet(j,1))
                cmod%d_arr_global(i,j,2) = cmod%b(i,ind_tet(j,2))*cmod%d(i,ind_tet(j,1)) - cmod%b(i,ind_tet(j,1))*cmod%d(i,ind_tet(j,2))
                cmod%d_arr_global(i,j,3) = cmod%b(i,ind_tet(j,1))*cmod%c(i,ind_tet(j,2)) - cmod%b(i,ind_tet(j,2))*cmod%c(i,ind_tet(j,1))
                cmod%d_arr_global(i,j,:) = 2.d0 * cmod%d_arr_global(i,j,:) * cmod%tet_len(i,j) / (6.d0*cmod%tet_vol(i))**2.d0
            end do
        end do
    end if
    ! calculate e d and b output
    cmod%e_field_global = 0.d0
    cmod%d_field_global = 0.d0
    do i=1,cmod%tet_NUM
        do j=1,6
            do k=1,3
                cmod%e_field_global(i,k) = cmod%e_field_global(i,k) + cmod%n_arr_global(i,j,k) * cmod%tet_drc(i,j) * res_e(cmod%tet_edg(i,j))
                cmod%d_field_global(i,k) = cmod%d_field_global(i,k) + cmod%d_arr_global(i,j,k) * cmod%tet_drc(i,j) * res_e(cmod%tet_edg(i,j))
            end do
        end do
    end do
    do i=1,cmod%tet_NUM
        do j=1,3
            cmod%b_field_global(i,j) = cmod%b_field_global(i,j) + cmod%d_field_global(i,j)*step_cur
        end do
    end do
    end subroutine vtk_field_save

    subroutine msh_output (cmod,step)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    USE ifport
    implicit none
    ! in/out variables
    type(model) :: cmod
    integer,intent(in) :: step
    ! local variables
    ! file info
    character(100) :: filename,tmpfile,fnam
    integer :: len1,len2,l
    logical :: logical_result
    integer :: tmp1,tmp2,tmp3,tmp4,tmp5
    ! cell info
    character(500) :: CHAR_list
    integer :: cell_num
    integer,allocatable :: cell_list(:)
    ! temp variables
    integer :: i,j,k,ierr
    ! n array and field variables
    !real(8) :: n_arr(6,3), d_arr(6,3)
    !real(8),allocatable :: e_field(:,:), d_field(:,:)
    ! process
    ! copies a file mshfile_name to filename and open it
    len1 = len_trim(msh_file)
    write (tmpfile, *) step
    filename = trim(file_root)//'_'//trim(adjustl(tmpfile))//".msh"
    len2 = len_trim(filename)
    fnam = 'copy ' //msh_file(1:len1) //' '//filename(1:len2)//' > nul'
    !fnam = 'cp ' //msh_file(1:len1) //' '//filename(1:len2)//' > nul'
    l = len_trim(fnam)
    logical_result = system(fnam(1:l))
    open(msh_unit,file=trim(filename),status="old")
    ! read in cell num and type
    DO
        READ(msh_unit,FMT="(A)",IOSTAT=IERR) CHAR_list
        if (CHAR_list(1:9)=="$Elements") exit
        if (IERR.NE.0) stop("cell number not found")
    end do
    READ(msh_unit,*) tmp1
    cell_num = 0
    do i=1,tmp1
        READ(msh_unit,*) tmp2,tmp3,tmp4,tmp5
        if (tmp4==4) then
            cell_num = cell_num+tmp5
            do j=1,tmp5
                READ(msh_unit,*) tmp2
            end do
        else
            do j=1,tmp5
                READ(msh_unit,*) tmp2
            end do
        end if
    end do
    allocate (cell_list(cell_num))
    rewind(msh_unit)
    DO
        READ(msh_unit,FMT="(A)",IOSTAT=IERR) CHAR_list
        if (CHAR_list(1:9)=="$Elements") exit
        if (IERR.NE.0) stop("cell number not found")
    end do
    k=0
    READ(msh_unit,*) tmp1
    do i=1,tmp1
        READ(msh_unit,*) tmp2,tmp3,tmp4,tmp5
        if (tmp4==4) then
            do j=1,tmp5
                k=k+1
                read(msh_unit,*) tmp2
                cell_list(k) = tmp2
            end do
        else
            do j=1,tmp5
                READ(msh_unit,*) tmp2
            end do
        end if
    end do
    close(msh_unit)
    open(msh_unit,file=trim(filename),status="old",position="append",action="write")
    do i=1,9
        if (cmod%o_f(i)) then
            write(msh_unit,"(a)") "$ElementData"
            write(msh_unit,"(a)") "1"
            SELECT CASE (i)
            CASE (1)
                write(msh_unit,"(a)") '"Ex data"'
            CASE (2)
                write(msh_unit,"(a)") '"Ey data"'
            CASE (3)
                write(msh_unit,"(a)") '"Ez data"'
            CASE (4)
                write(msh_unit,"(a)") '"Dx data"'
            CASE (5)
                write(msh_unit,"(a)") '"Dy data"'
            CASE (6)
                write(msh_unit,"(a)") '"Dz data"'
            CASE (7)
                write(msh_unit,"(a)") '"Bx data"'
            CASE (8)
                write(msh_unit,"(a)") '"By data"'
            CASE (9)
                write(msh_unit,"(a)") '"Bz data"'
            END SELECT
            write(msh_unit,"(a)") "1"
            write(msh_unit,"(a)") "0.0"
            write(msh_unit,"(a)") "3"
            write(msh_unit,"(a)") "0"
            write(msh_unit,"(a)") "1"
            write(msh_unit,*) cell_num
            SELECT CASE (i)
            CASE (1)
                do k=1,cell_num
                    write(msh_unit,"(i10,e16.8)") cell_list(k),abs(cmod%e_field_global(K,1))
                end do
            CASE (2)
                do k=1,cell_num
                    write(msh_unit,"(i10,e16.8)") cell_list(k),abs(cmod%e_field_global(K,2))
                end do
            CASE (3)
                do k=1,cell_num
                    write(msh_unit,"(i10,e16.8)") cell_list(k),abs(cmod%e_field_global(K,3))
                end do
            CASE (4)
                do k=1,cell_num
                    write(msh_unit,"(i10,e16.8)") cell_list(k),abs(cmod%d_field_global(K,1))
                end do
            CASE (5)
                do k=1,cell_num
                    write(msh_unit,"(i10,e16.8)") cell_list(k),abs(cmod%d_field_global(K,2))
                end do
            CASE (6)
                do k=1,cell_num
                    write(msh_unit,"(i10,e16.8)") cell_list(k),abs(cmod%d_field_global(K,3))
                end do
            CASE (7)
                do k=1,cell_num
                    write(msh_unit,"(i10,e16.8)") cell_list(k),abs(cmod%b_field_global(K,1))
                end do
            CASE (8)
                do k=1,cell_num
                    write(msh_unit,"(i10,e16.8)") cell_list(k),abs(cmod%b_field_global(K,2))
                end do
            CASE (9)
                do k=1,cell_num
                    write(msh_unit,"(i10,e16.8)") cell_list(k),abs(cmod%b_field_global(K,3))
                end do
            END SELECT
            write(msh_unit,"(a)") "$EndElementData"
        end if
    end do
    IF (.NOT. COMMITQQ(msh_unit)) stop 'msh_unit writing failed, stop!'
    close(msh_unit)
    end subroutine msh_output

    subroutine vtk_output (cmod,step)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    USE ifport
    implicit none
    ! in/out variables
    type(model),intent(in) :: cmod
    integer,intent(in) :: step
    ! local variables
    ! file info
    character(100) :: filename,tmpfile,fnam
    integer :: len1,len2,l
    logical :: logical_result
    ! cell info
    character(5000) :: CHAR_list,char_temp
    integer :: IERR
    integer :: cell_num
    integer,allocatable :: cell_type(:)
    ! temp variables
    integer :: i,j,k
    ! process
    ! copies a file vtkfile_name to filename and open it
    len1 = len_trim(vtk_file)
    write (tmpfile,*) step
    filename = trim(vtk_file)//'_'//trim(adjustl(tmpfile))//".vtk"
    len2 = len_trim(filename)
    fnam = 'copy ' //vtk_file(1:len1) //' '//filename(1:len2)//' > nul'
    !fnam = 'cp ' //vtk_file(1:len1) //' '//filename(1:len2)//' > nul'
    l = len_trim(fnam)
    logical_result = system(fnam(1:l))
    open(vtk_unit,file=trim(filename))
    ! read in cell num and type
    DO
        READ(vtk_unit,FMT="(A)",IOSTAT=IERR) CHAR_list
        if (CHAR_list(1:10)=="CELL_TYPES") then
            read(CHAR_list,*) char_temp,cell_num
            exit
        end if
        if (IERR.NE.0) stop "CELL_TYPES not found"
    end do
    allocate(cell_type(cell_num))
    do i=1,cell_num
        read(vtk_unit,*) cell_type(i)
    end do
    DO
        READ(UNIT=vtk_unit,IOSTAT=IERR) CHAR_list
        IF (IERR.NE.0) EXIT
    end do
    !open(vtk_unit,file=trim(filename),status="old",position="append",action="write")
    write(vtk_unit,"(a)") ""
    ! e field
    ! HEDAER
    write(vtk_unit,"(a,i10)") "CELL_DATA",cell_num
    write(vtk_unit,"(a)") "VECTORS FIELD_ELECTRIC double"
    !write(89,"(a)") "LOOKUP_TABLE default"
    ! output
    K=0
    do i=1,cell_num
        if (cell_type(i)==1) then
            write(vtk_unit,"(3e16.8)") 0.0,0.0,0.0
        else if (cell_type(i)==3) then
            write(vtk_unit,"(3e16.8)") 0.0,0.0,0.0
        else if (cell_type(i)==10) then
            K=K+1
            write(vtk_unit,"(3e16.8)") (abs(cmod%e_field_global(K,j)),j=1,3)
        else
            stop "error"
        end if
    end do
    write(vtk_unit,"(a)") ""
    ! dbdt field
    !write(89,"(a,i10)") "CELL_DATA",cell_num
    write(vtk_unit,"(a)") "VECTORS FIELD_derivB double"
    K=0
    do i=1,cell_num
        if (cell_type(i)==1) then
            write(vtk_unit,"(3e16.8)") 0.0,0.0,0.0
        else if (cell_type(i)==3) then
            write(vtk_unit,"(3e16.8)") 0.0,0.0,0.0
        else if (cell_type(i)==10) then
            K=K+1
            write(vtk_unit,"(3e16.8)") (abs(cmod%d_field_global(K,j)),j=1,3)
        else
            stop "error"
        end if
    end do
    write(vtk_unit,"(a)") ""
    ! b field
    !write(89,"(a,i10)") "CELL_DATA",cell_num
    write(vtk_unit,"(a)") "VECTORS FIELD_B double"
    K=0
    do i=1,cell_num
        if (cell_type(i)==1) then
            write(vtk_unit,"(3e16.8)") 0.0,0.0,0.0
        else if (cell_type(i)==3) then
            write(vtk_unit,"(3e16.8)") 0.0,0.0,0.0
        else if (cell_type(i)==10) then
            K=K+1
            write(vtk_unit,"(3e16.8)") (abs(cmod%b_field_global(K,j)),j=1,3)
        else
            stop "error"
        end if
    end do
    IF (.NOT. COMMITQQ(vtk_unit)) stop 'vtk_unit writing failed, stop!'
    close(vtk_unit)
    end subroutine vtk_output

    subroutine mkl_to_csr (ai_mkl,n_csr,NNZ_CSR,a_csr,i_csr,j_csr)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    USE, INTRINSIC :: ISO_C_BINDING, only : C_PTR,c_f_pointer
    implicit none
    ! IN/OUT VARIABLES
    type(SPARSE_MATRIX_T),intent(in) :: ai_mkl
    integer,intent(out) :: n_csr,NNZ_CSR
    integer,allocatable,intent(out) :: i_csr(:),j_csr(:)
    real(8),allocatable,intent(out) :: a_csr(:)
    ! LOCAL VARIABLES
    integer :: info_1,info_2
    TYPE(C_PTR) :: isys_c_1,isys_c_2,jsys_c_1,vsys_c_1
    integer,pointer :: isys_1(:),isys_2(:),jsys_1(:)
    real(8),pointer :: vsys_1(:)
    ! process
    info_1 = mkl_sparse_d_export_csr (ai_mkl, info_2, n_csr, n_csr, isys_c_1, isys_c_2, jsys_c_1, vsys_c_1)
    if (info_1/=0) stop
    call c_f_pointer(isys_c_1, isys_1,[n_csr])
    call c_f_pointer(isys_c_2, isys_2,[n_csr])
    allocate(i_csr(n_csr+1))
    i_csr (1:n_csr) = isys_1
    i_csr (2:n_csr+1) = isys_2
    NNZ_CSR = i_csr(n_csr+1)-1
    call c_f_pointer(jsys_c_1, jsys_1,[NNZ_CSR])
    call c_f_pointer(vsys_c_1, vsys_1,[NNZ_CSR])
    allocate (j_csr(NNZ_CSR),a_csr(NNZ_CSR))
    j_csr = jsys_1
    a_csr = vsys_1
    end subroutine mkl_to_csr

    subroutine mkl_add_3 (num,av,ai_mkl,bv,bi_mkl,cv,ci_mkl,r_mkl)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    implicit none
    ! IN/OUT VARIABLES
    type(SPARSE_MATRIX_T),intent(in) :: ai_mkl,bi_mkl,ci_mkl
    integer,intent(in) :: num
    real(8),intent(in) :: av,bv,cv
    type(SPARSE_MATRIX_T),intent(out) :: r_mkl
    ! other variables
    integer :: info
    type(SPARSE_MATRIX_T) :: z_mkl,zt_mkl,tmp_mkl_1,tmp_mkl_2
    integer :: z_I(1),z_j(1)
    real(8) :: z_v(1)
    ! z_mkl
    z_i(1) = 1; z_j(1) = 1; z_v(1) = 0.d0
    info = mkl_sparse_d_create_coo (zt_mkl, SPARSE_INDEX_BASE_one, num, num, 1, z_i, z_j, z_v)
    if (info/=0) stop '2325'
    info = mkl_sparse_convert_csr (zt_mkl, SPARSE_OPERATION_NON_TRANSPOSE, z_mkl)
    if (info/=0) stop '2327'
    ! process
    info = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, Ai_mkl, av , z_mkl, tmp_mkl_1)
    if (info/=0) stop '2331'
    info = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, bi_mkl, bv , tmp_mkl_1, tmp_mkl_2)
    if (info/=0) stop '2333'
    info = mkl_sparse_d_add (SPARSE_OPERATION_NON_TRANSPOSE, ci_mkl, cv , tmp_mkl_2, r_mkl)
    if (info/=0) stop '2335'
    ! destroy
    info = mkl_sparse_destroy (z_mkl)
    info = mkl_sparse_destroy (zt_mkl)
    info = mkl_sparse_destroy (tmp_mkl_1)
    info = mkl_sparse_destroy (tmp_mkl_2)
    end subroutine mkl_add_3

    subroutine clear_model (cmod)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
    implicit none
    type(model),intent(inout) :: cmod
    integer :: info
    ! node
    cmod%nod_num=0
    if (allocated(cmod%nod_cor)) deallocate(cmod%nod_cor)
    ! source
    cmod%seg_num=0
    if (allocated(cmod%src_vex)) deallocate(cmod%src_vex)
    if (allocated(cmod%src_phy)) deallocate(cmod%src_phy)
    if (allocated(cmod%src_edg)) deallocate(cmod%src_edg)
    if (allocated(cmod%src_drc)) deallocate(cmod%src_drc)
    ! source info from the control file
    cmod%src_num=0
    cmod%src_sourcetype=0
    cmod%src_pulsetype=0
    cmod%src_pulsewidth=0.d0
    cmod%src_current=0.d0
    if (allocated(cmod%src_tp)) deallocate(cmod%src_tp)
    if (allocated(cmod%src_count_edges)) deallocate(cmod%src_count_edges)
    if (allocated(cmod%src_end_1)) deallocate(cmod%src_end_1)
    if (allocated(cmod%src_end_2)) deallocate(cmod%src_end_2)
    ! tet
    cmod%tet_num=0
    if (allocated(cmod%tet_vex)) deallocate(cmod%tet_vex)
    if (allocated(cmod%tet_phy)) deallocate(cmod%tet_phy)
    if (allocated(cmod%tet_drc)) deallocate(cmod%tet_drc)
    if (allocated(cmod%tet_vol)) deallocate(cmod%tet_vol)
    if (allocated(cmod%tet_edg)) deallocate(cmod%tet_edg)
    if (allocated(cmod%tet_len)) deallocate(cmod%tet_len)
    if (allocated(cmod%a)) deallocate(cmod%a)
    if (allocated(cmod%b)) deallocate(cmod%b)
    if (allocated(cmod%c)) deallocate(cmod%c)
    if (allocated(cmod%d)) deallocate(cmod%d)
    ! sur_edg
    if (allocated(cmod%sur_edg_list)) deallocate(cmod%sur_edg_list)
    if (allocated(cmod%sur_edg_index)) deallocate(cmod%sur_edg_index)
    ! edge
    cmod%edg_num=0
    if (allocated(cmod%edg_vex)) deallocate(cmod%edg_vex)
    if (allocated(cmod%edg_len)) deallocate(cmod%edg_len)
    ! receiver
    cmod%rec_num=0
    if (allocated(cmod%rec_cor)) deallocate(cmod%rec_cor)
    if (allocated(cmod%rec_tet)) deallocate(cmod%rec_tet)
    if (allocated(cmod%n_arr)) deallocate(cmod%n_arr)
    if (allocated(cmod%d_arr)) deallocate(cmod%d_arr)
    ! domain character
    cmod%iso=.true.
    cmod%dom_num=0
    if (allocated(cmod%dom_ind)) deallocate(cmod%dom_ind)
    if (allocated(cmod%sig_3d)) deallocate(cmod%sig_3d)
    if (allocated(cmod%per_3d)) deallocate(cmod%per_3d)
    if (allocated(cmod%mag_3d)) deallocate(cmod%mag_3d)
    if (allocated(cmod%sig_a3d)) deallocate(cmod%sig_a3d)
    if (allocated(cmod%per_a3d)) deallocate(cmod%per_a3d)
    if (allocated(cmod%mag_a3d)) deallocate(cmod%mag_a3d)
    ! system matrix
    info = mkl_sparse_destroy (cmod%a_mkl)
    info = mkl_sparse_destroy (cmod%b_mkl)
    info = mkl_sparse_destroy (cmod%c_mkl)
    ! surf_kind: 0 - PEC boundary ; 1 - 1st order ABC
    cmod%bound_cond=0
    ! time variables
    cmod%step_mthd=0
    cmod%step_dvsi=0
    cmod%step_dble=0
    cmod%step_dblesize=0
    cmod%step_init=0.d0
    cmod%time_maxm=0.d0
    ! field variables
    if (allocated(cmod%a_f)) deallocate(cmod%a_f)
    cmod%ind_f=0
    cmod%init_size_F=0
    ! output variables: ex ey ez dx dy dz bx by bz
    cmod%ind_f=0
    cmod%init_size_F=0
    cmod%o_f=.false.
    if (allocated(cmod%o_v)) deallocate(cmod%o_v)
    ! vtk output variables
    cmod%vtk_num = 0
    if (allocated(cmod%vtk_times))  deallocate(cmod%vtk_times)
    if (allocated(cmod%vtk_done))  deallocate(cmod%vtk_done)
    if (allocated(cmod%n_arr_global))  deallocate(cmod%n_arr_global)
    if (allocated(cmod%d_arr_global))  deallocate(cmod%d_arr_global)
    if (allocated(cmod%e_field_global))  deallocate(cmod%e_field_global)
    if (allocated(cmod%d_field_global))  deallocate(cmod%d_field_global)
    if (allocated(cmod%b_field_global))  deallocate(cmod%b_field_global)
    ! inversion variables
    cmod%deriv = .false.
    if (allocated(cmod%a_prime))  deallocate(cmod%a_prime)
    cmod%total_step_record = 0
    if (allocated(cmod%step_record))  deallocate(cmod%step_record)
    if (allocated(cmod%e_record))  deallocate(cmod%e_record)
    end subroutine clear_model

    end module fwd_modeling