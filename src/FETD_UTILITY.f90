	! included MKL headers
	include "mkl_service.f90"
	include "mkl_spblas.f90"
	include 'mkl_pardiso.f90'
	include 'mkl_sparse_handle.f90'
	include 'mkl_cluster_sparse_solver.f90'

	module utility
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	INTERFACE swap
	MODULE PROCEDURE swap_i,swap_d,masked_swap_i,masked_swap_d
	END INTERFACE
	INTERFACE assert_eq
	MODULE PROCEDURE assert_eq2,assert_eq3
	END INTERFACE
    contains
    
	real(8) function interp_lin (x1,y1,x2,y2,x0)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	real(8),intent(in) :: x1,y1,x2,y2,x0
	real(8) :: a,b
	a = (y1-y2) / (x1-x2)
	b = y1-a*x1
	interp_lin = a*x0+b
    end function interp_lin
    
	real(8) FUNCTION det2 (A)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	IMPLICIT NONE
	real(8), INTENT(IN)  :: A(:,:)
	det2 = A(1,1)*A(2,2) - A(1,2)*A(2,1)
	RETURN
    END FUNCTION det2
    
	real(8) FUNCTION det3 (A)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	IMPLICIT NONe
	real(8), INTENT(IN)  :: A(:,:)
	det3 = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3)  + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)
	RETURN
    eND FUNCTION det3
    
	real(8) FUNCTION det4 (A)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	IMPLICIT NONE
	real(8), INTENT(IN)  :: A(:,:)
	det4 =  A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))-A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+ &
		A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+ &
		A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
	RETURN
    END FUNCTION det4
    
	real(8) FUNCTION det_mkl_3 (A)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	real(8),intent(in) :: a(3,3)
	integer :: i,info,piv(4)
	call dgetrf(3,3,a,3,piv,info)
	det_mkl_3 = 0d0
	if (info/=0) return
	det_mkl_3 = 1d0
	do i=1,3
		if (piv(i).ne.i) then
			det_mkl_3 = -det_mkl_3 * a(i,i)
		else
			det_mkl_3 = det_mkl_3 * a(i,i)
		endif
	end do
    END FUNCTION det_mkl_3
    
	real(8) FUNCTION det_mkl_4 (A)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	real(8),intent(in) :: a(4,4)
	integer :: i,info,piv(4)
	call dgetrf(4,4,a,4,piv,info)
	det_mkl_4 = 0d0
	if (info/=0) return
	det_mkl_4 = 1d0
	do i=1,4
		if (piv(i).ne.i) then
			det_mkl_4 = -det_mkl_4 * a(i,i)
		else
			det_mkl_4 = det_mkl_4 * a(i,i)
		endif
	end do
    END FUNCTION det_mkl_4
    
	real(8) function tet_volume (x,y,z)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	real(8),intent(in) :: x(4),y(4),z(4)
	tet_volume = det_mkl_4(RESHAPE([1.d0,1.d0,1.d0,1.d0,x(1),x(2),x(3),x(4),y(1),y(2),y(3),y(4),z(1),z(2),z(3),z(4)],[4,4])) / 6.d0
    end function tet_volume
    
	real(8) function dist (p_1,p_2)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	real(8),intent(in) :: p_1(:),p_2(:)
	dist = sqrt( (p_1(1)-p_2(1))**2.d0 + (p_1(2)-p_2(2))**2.d0 + (p_1(3)-p_2(3))**2.d0 )
    end function dist
    
	real(8) function dirac (I,J)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	IMPLICIT NONE
	INTEGER,INTENT(IN) :: I,J
	dirac = 0.d0
	IF (I==J) dirac = 1.d0
    END function dirac
    
	function cross_p (a,b)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	real(8) :: cross_p(3),a(3),b(3)
	cross_p(1) = a(2) * b(3) - a(3) * b(2)
	cross_p(2) = a(3) * b(1) - a(1) * b(3)
	cross_p(3) = a(1) * b(2) - a(2) * b(1)
    end function cross_p
    
	SUBROUTINE swap_i(a,b)
	INTEGER :: a,b,dum
	dum=a;a=b;b=dum
    END SUBROUTINE swap_i
    
	SUBROUTINE swap_d(a,b)
	REAL(8) :: a,b,dum
	dum=a;a=b;b=dum
    END SUBROUTINE swap_d
    
	SUBROUTINE masked_swap_i(a,b,mask)
	integer :: a,b,swp
	LOGICAL, INTENT(IN) :: mask
	if (mask) then
		swp=a;a=b;b=swp
	end if
    END SUBROUTINE masked_swap_i
    
	SUBROUTINE masked_swap_d(a,b,mask)
	REAL(8) :: a,b,swp
	LOGICAL, INTENT(IN) :: mask
	if (mask) then
		swp=a;a=b;b=swp
	end if
	END SUBROUTINE masked_swap_d
	FUNCTION assert_eq2(n1,n2,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN) :: n1,n2
	INTEGER :: assert_eq2
	if (n1 == n2) then
		assert_eq2=n1
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', string
		STOP 'program terminated by assert_eq2'
	end if
	END FUNCTION assert_eq2
	FUNCTION assert_eq3(n1,n2,n3,string)
	CHARACTER(LEN=*), INTENT(IN) :: string
	INTEGER, INTENT(IN) :: n1,n2,n3
	INTEGER :: assert_eq3
	if (n1 == n2 .and. n2 == n3) then
		assert_eq3=n1
	else
		write (*,*) 'nrerror: an assert_eq failed with this tag:', string
		STOP 'program terminated by assert_eq3'
	end if
    END FUNCTION assert_eq3
    
	FUNCTION arth_i(first,increment,n)
	INTEGER, PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
	INTEGER, INTENT(IN) :: first,increment,n
	INTEGER :: arth_i(n)
	INTEGER :: k,k2,temp
	if (n > 0) arth_i(1)=first
	if (n <= NPAR_ARTH) then
		do k=2,n
			arth_i(k)=arth_i(k-1)+increment
		end do
	else
		do k=2,NPAR2_ARTH
			arth_i(k)=arth_i(k-1)+increment
		end do
		temp=increment*NPAR2_ARTH
		k=NPAR2_ARTH
		do
			if (k >= n) exit
			k2=k+k
			arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
			temp=temp+temp
			k=k2
		end do
	end if
    END FUNCTION arth_i
    
	SUBROUTINE indexx_i(iarr,index)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: iarr(:)
	INTEGER, INTENT(OUT) :: index(:)
	INTEGER, PARAMETER :: NN=15, NSTACK=50
	INTEGER :: a
	INTEGER :: n,k,i,j,indext,jstack,l,r
	INTEGER :: istack(NSTACK)
	n=assert_eq(size(index),size(iarr),'indexx_sp')
	index=arth_i(1,1,n)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				indext=index(j)
				a=iarr(indext)
				do i=j-1,l,-1
					if (iarr(index(i)) <= a) exit
					index(i+1)=index(i)
				end do
				index(i+1)=indext
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(index(k),index(l+1))
			call icomp_xchg(index(l),index(r))
			call icomp_xchg(index(l+1),index(r))
			call icomp_xchg(index(l),index(l+1))
			i=l+1
			j=r
			indext=index(l+1)
			a=iarr(indext)
			do
				do
					i=i+1
					if (iarr(index(i)) >= a) exit
				end do
				do
					j=j-1
					if (iarr(index(j)) <= a) exit
				end do
				if (j < i) exit
				call swap(index(i),index(j))
			end do
			index(l+1)=index(j)
			index(j)=indext
			jstack=jstack+2
			if (jstack > NSTACK) stop'indexx: NSTACK too small'
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	contains
	SUBROUTINE icomp_xchg(i,j)
	INTEGER :: i,j,swp
	if (iarr(j) < iarr(i)) then
		swp=i;i=j;j=swp
	end if
	END SUBROUTINE icomp_xchg
    END SUBROUTINE indexx_i
    
	SUBROUTINE sort_i(arr)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	IMPLICIT NONE
	integer, INTENT(INOUT) :: arr(:)
	INTEGER, PARAMETER :: NN=15, NSTACK=50
	integer :: a
	INTEGER :: n,k,i,j,jstack,l,r
	INTEGER :: istack(NSTACK)
	n=size(arr)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				a=arr(j)
				do i=j-1,l,-1
					if (arr(i) <= a) exit
					arr(i+1)=arr(i)
				end do
				arr(i+1)=a
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(arr(k),arr(l+1))
			call swap(arr(l),arr(r),arr(l)>arr(r))
			call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
			call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
			i=l+1
			j=r
			a=arr(l+1)
			do
				do
					i=i+1
					if (arr(i) >= a) exit
				end do
				do
					j=j-1
					if (arr(j) <= a) exit
				end do
				if (j < i) exit
				call swap(arr(i),arr(j))
			end do
			arr(l+1)=arr(j)
			arr(j)=a
			jstack=jstack+2
			if (jstack > NSTACK) stop'sort: NSTACK too small'
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
    END SUBROUTINE sort_i
    
	SUBROUTINE sort_ii(arr,slave)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	IMPLICIT NONE
	integer, INTENT(INOUT) :: arr(:),slave(:)
	INTEGER :: ndum
	INTEGER, DIMENSION(size(arr)) :: index
	integer :: ptr_1,ptr_2,n,i
	ndum=assert_eq(size(arr),size(slave),'sort2')
	! phase 1 sort
	call indexx_i(arr,index)
	arr=arr(index)
	slave=slave(index)
	! phase 2 sort
	n = size(arr)
	ptr_1 = 1
	do i=1,n-1
		if( arr(i)==arr(i+1) ) then
			ptr_2=i+1
			if (ptr_2==n) call sort_i(slave(ptr_1:ptr_2))
			cycle
		else
			ptr_2=i
			call sort_i(slave(ptr_1:ptr_2))
			ptr_1=i+1
			cycle
		end if
	end do
    END SUBROUTINE sort_ii
    
	subroutine unique_i(arr)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	integer, allocatable, INTENT(INOUT) :: arr(:)
	integer :: n,i
	logical :: mask(size(arr))
	integer :: arr_co(size(arr))
	n = size(arr)
	arr_co = arr
	call sort_i (arr_co)
	mask = .true.
	DO i = n,2,-1
		mask(i) = .NOT. (arr_co(i-1)==arr_co(i))
	END DO
	deallocate(arr)
	allocate(arr,source=pack(arr_co,mask))
    end subroutine unique_i
    
	subroutine unique_ii(arr,slave)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	integer, allocatable, INTENT(INOUT) :: arr(:),slave(:)
	integer :: n,i
	logical :: mask(size(arr))
	integer :: arr_co(size(arr)),slave_co(size(slave))
	n = size(arr)
	arr_co = arr
	slave_co = slave
	call sort_ii (arr_co,slave_co)
	mask = .true.
	DO i = n,2,-1
		mask(i) = .NOT. (arr_co(i-1)==arr_co(i).and.slave_co(i-1)==slave_co(i))
	END DO
	deallocate(arr,slave)
	allocate(arr,source=pack(arr_co,mask))
	allocate(slave,source=pack(slave_co,mask))
    end subroutine unique_ii
    
	SUBROUTINE sort_ir(arr,slave)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	IMPLICIT NONE
	integer, DIMENSION(:), INTENT(INOUT) :: arr
	real(8), DIMENSION(:), INTENT(INOUT) :: slave
	INTEGER :: ndum
	INTEGER, DIMENSION(size(arr)) :: index
	ndum=assert_eq(size(arr),size(slave),'sort2')
	call indexx_i(arr,index)
	arr=arr(index)
	slave=slave(index)
    END SUBROUTINE sort_ir
    
	SUBROUTINE sort_iir(arr,slave1,slave2)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	IMPLICIT NONE
	integer, DIMENSION(:), INTENT(INOUT) :: arr,slave1
	REAL(8), DIMENSION(:), INTENT(INOUT) :: slave2
	INTEGER :: ndum, ptr_1,ptr_2,n,i
	INTEGER, DIMENSION(size(arr)) :: index
	ndum=assert_eq(size(arr),size(slave1),size(slave2),'sort3')
	! phase 1 sort
	call indexx_i(arr,index)
	arr=arr(index)
	slave1=slave1(index)
	slave2=slave2(index)
	! phase 2 sort
	n=size(arr)
	ptr_1 = 1
	do i=1,n-1
		if( arr(i)==arr(i+1) ) then
			ptr_2=i+1
			if (ptr_2==n) call sort_ir(slave1(ptr_1:ptr_2),slave2(ptr_1:ptr_2))
			cycle
		else
			ptr_2=i
			call sort_ir(slave1(ptr_1:ptr_2),slave2(ptr_1:ptr_2))
			ptr_1=i+1
			cycle
		end if
	end do
    END SUBROUTINE sort_iir
    
	real(8) function inf_norm (a)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	include "mkl_vml.f90"
	! in/out variables
	real(8),intent(in),allocatable :: a(:)
	! local variablea
	real(8),allocatable :: tmp_arr(:)
	integer :: length
	! process
	length = size(a)
	allocate(tmp_arr(length))
	tmp_arr = abs(a)
	inf_norm = maxval(tmp_arr)
    end function inf_norm
    
	real(8) function one_norm (a)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	include "mkl_vml.f90"
	! in/out variables
	real(8),intent(in),allocatable :: a(:)
	! local variablea
	real(8),allocatable :: tmp_arr(:)
	integer :: length
	! process
	length = size(a)
	allocate(tmp_arr(length))
	tmp_arr = abs(a)
	one_norm = sum(tmp_arr)
    end function one_norm
    
	real(8) function two_norm (a)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	include "mkl_vml.f90"
	! in/out variables
	real(8),intent(in),allocatable :: a(:)
	! local variablea
	real(8),allocatable :: tmp_arr(:)
	integer :: length
	! process
	length = size(a)
	allocate(tmp_arr(length))
	call vdsqr ( length, a ,tmp_arr )
	two_norm = sqrt(sum(tmp_arr))
    end function two_norm
    
	real(8) function inf_error (a,b)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	! output: |a-b|(inf)/|b|(inf)
	implicit none
	include "mkl_vml.f90"
	! in/out variables
	real(8),intent(in),allocatable :: a(:),b(:)
	! local variablea
	real(8),allocatable :: tmp_arr(:)
	real(8) :: val_1,val_2
	integer :: length
	! process
	if (size(a)/=size(b)) stop "ERROR - array size do not consist, stop!"
	length = size(a)
	allocate(tmp_arr(length))
	call vdsub (length,a,b,tmp_arr)
	val_1 = inf_norm (tmp_arr)
	val_2 = inf_norm (b)
	inf_error = val_1 / val_2
    end function inf_error
    
	real(8) function one_error (a,b)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	! output: |a-b|/|b|
	implicit none
	include "mkl_vml.f90"
	! in/out variables
	real(8),intent(in),allocatable :: a(:),b(:)
	! local variablea
	real(8),allocatable :: tmp_arr(:)
	real(8) :: val_1,val_2
	integer :: length
	! process
	if (size(a)/=size(b)) stop "ERROR - array size do not consist, stop!"
	length = size(a)
	allocate(tmp_arr(length))
	call vdsub (length,a,b,tmp_arr)
	val_1 = one_norm (tmp_arr)
	val_2 = one_norm (b)
	one_error = val_1 / val_2
    end function one_error
    
	real(8) function two_error (a,b)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	! output: ||a-b||/||b||
	implicit none
	include "mkl_vml.f90"
	! in/out variables
	real(8),intent(in),allocatable :: a(:),b(:)
	! local variablea
	real(8),allocatable :: tmp_arr(:)
	real(8) :: val_1,val_2
	integer :: length
	! process
	if (size(a)/=size(b)) stop "ERROR - array size do not consist, stop!"
	length = size(a)
	allocate(tmp_arr(length))
	call vdsub (length,a,b,tmp_arr)
	val_1 = two_norm (tmp_arr)
	val_2 = two_norm (b)
	two_error = val_1 / val_2
    end function two_error
    
	integer function bSearch_I (a, value)
	integer, intent(in), target :: a(:)
	integer, intent(in)         :: value
	integer, pointer            :: p(:)
	integer                  :: mid, offset
	p => a
	bSearch_I = 0
	offset = 0
	do while (size(p) > 0)
		mid = size(p)/2 + 1
		if (p(mid) > value) then
			p => p(:mid-1)
		else if (p(mid) < value) then
			offset = offset + mid
			p => p(mid+1:)
		else
			bSearch_I = offset + mid
			return
		end if
	end do
    end function bSearch_I
    
	integer function bSearch_Ii (a, b, va, vb)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	integer, intent(in), target :: a(:),b(:)
	integer, intent(in)         :: va,vb
	integer, pointer            :: p(:)
	integer :: mid, offset, ind1, ind2, search2
	p => a
	bSearch_Ii = 0
	offset = 0
	do while (size(p) > 0)
		mid = size(p)/2 + 1
		if (p(mid) > va) then
			p => p(:mid-1)
		else if (p(mid) < va) then
			offset = offset + mid
			p => p(mid+1:)
		else
			bSearch_Ii = offset + mid
			ind1 = bSearch_Ii
			ind2 = bSearch_Ii
			do
				if (ind1>1) then
					if (a(ind1-1)==a(ind1)) then
						ind1 = ind1-1
					else
						exit
					end if
				else
					exit
				end if
			end do
			do
				if (ind2<size(a)) then
					if (a(ind2+1)==a(ind2)) then
						ind2 = ind2+1
					else
						exit
					end if
				else
					exit
				end if
			end do
			search2 = bSearch_I(b(ind1:ind2),vb)
			if (search2==0) then
				bSearch_Ii = 0
				return
			end if
			bSearch_Ii = ind1 + search2 - 1
			return
		end if
	end do
    end function bSearch_Ii
    
	subroutine parseCode( nLen, sLine, sCode, sValue, bComment )
	! parse a line read from a file into a code & value. Force the code to be all lowercase with no ending colon. Terminate the line at a '%' or '!' sign (these allow for user comments!)
	implicit none
	! Args
	integer, intent(in)   :: nLen
	character(nLen)       :: sLine
	character(nLen), intent(out) :: sCode, sValue
	logical, intent(out)    :: bComment
	! Local vars
	integer :: iFrom, iTo
	! Init returns
	bComment = .false.
	sCode = ' '
	sValue = ' '
	! Convert all tab characters to spaces
	forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '
	! Skip any beginning blanks
	do iFrom = 1,nLen
		if (sLine(iFrom:iFrom) .ne. ' ') exit
	enddo
	! If the first char is a comment char, then the whole line is a comment. Also, if the line is blank, consider it a comment.
	if (iFrom >= nLen) then !KWK may 2009 pulled this out in from since sometimes iFrom > nLen and this kills (iFrom:iFrom) below
		bComment = .true.
		return
	endif
	if(  sLine(iFrom:iFrom) == '%' .or. sLine(iFrom:iFrom) == '!' ) then
		bComment = .true.
		return
	endif
	! Pull off the code value. Cvt to lowercase as we go.
	iTo = index(sLine,':') - 1
	if (iTo < iFrom) then
		!write(*,*) 'Parsing Error: missing colon in line:',sLine
		return
	endif
	sCode = sLine(iFrom:iTo)
	call Lower(sCode)
	! Skip spaces after the colon
	do iFrom = iTo+2,nLen
		if (sLine(iFrom:iFrom) .ne. ' ') exit
	enddo
	! Get the rest, up to any comment
	sValue = sLine(iFrom:)
	iTo = len_trim(sValue)
	iFrom = index(sValue,'%')
	if (iFrom > 0 .and. iFrom < iTo) then
		sValue(iFrom:iTo) = ' '
	endif
	iFrom = index(sValue,'!')
	if (iFrom > 0 .and. iFrom < iTo) then
		sValue(iFrom:iTo) = ' '
	endif
	!call Lower(sValue)   ! No: Some values are filenames which are case-sensitive on UNIX!
    end subroutine parseCode
    
	subroutine parseLine( nLen, sLine, bComment )
	! Subroutine to check if the sLine is blank or a comment line, and if it isn't then any comment at the end of the line is blanked.
	! This is useful for reading in data tables that have user comment lines or comments at the end of the line, denoted by the ! and % symbols.
	! If the entire line is a comment, bComment = .true.
	implicit none
	! Args
	integer, intent(in)     :: nLen
	character(nLen)         :: sLine
	logical, intent(out)    :: bComment
	! Local vars
	integer :: iFrom, iTo
	! Init returns
	bComment = .false.
	! Convert all tab characters to spaces
	forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '
	! Skip any beginning blanks
	do iFrom = 1,nLen
		if (sLine(iFrom:iFrom) .ne. ' ') exit
	enddo
	! If the first char is a comment char, then the whole line is a comment.
	! DGM April 2008 Also, if the line is blank, consider it a comment.
	if( iFrom >= nLen .or. sLine(iFrom:iFrom) == '%' .or. sLine(iFrom:iFrom) == '!' ) then
		bComment = .true.
		return
	endif
	! Now trim off any comments at the end of the line
	iTo = len_trim(sLine)
	iFrom = index(sLine,'%')
	if (iFrom > 0 .and. iFrom < iTo) then
		sLine(iFrom:iTo) = ' '
	endif
	iFrom = index(sLine,'!')
	if (iFrom > 0 .and. iFrom < iTo) then
		sLine(iFrom:iTo) = ' '
	endif
	end subroutine parseLine
	subroutine Lower( s )
	! convert string to lower case
	character(*), intent(out)  :: s
	integer :: i
	do  i=1,len_trim(s)
		if  ( s(i:i) >= 'A' .and. s(i:i) <= 'Z' ) then
			s(i:i) = char(ichar(s(i:i)) + 32)
		endif
	enddo
    end subroutine Lower
    
	subroutine parseFields( nLen, sLine, nFields, sFields)
	! Routine to parse out mixed numeric and character fields from a line
	! This is useful for mixed format input tables, which are awkward for Fortran
	implicit none
	! Args
	integer, intent(in)     :: nLen
	character(nLen)         :: sLine
	integer, intent(in)     :: nFields
	character(nLen), intent(out) :: sFields(nFields)
	! Local vars
	integer :: iFrom, iTo, i
	! Convert all tab characters to spaces
	forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '
	iFrom = 1
	! Loop through the line and get the nFields:
	do i=1,nFields
		! Skip any beginning blanks
		do iFrom = iFrom,nLen
			if (sLine(iFrom:iFrom) .ne. ' ') exit
		enddo
		! Pull out nonblank character string:
		do iTo = iFrom,nLen
			if (sLine(iTo:iTo) .eq. ' ') exit
		enddo
		sFields(i) = trim(sLine(iFrom:iTo-1))
		iFrom = iTo
	enddo
    end subroutine parseFields
    
	subroutine  get_time_offset(timein,timeout)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	integer :: i,j,k,mjd,values(8)
	real(8) :: timein, timeout, fracday
	call date_and_time(values=values)
	i = values(1)
	j = values(2)
	k = values(3)
	mjd = -678927 + k + 1461*(i+(j-14)/12)/4 + 367*(j-2-12 * ((j-14)/12))/12 + (24002-12*i-j)/1200
	fracday = ((values(5)*60.d0 + values(6))*60.d0 + values(7) + values(8)/1000.d0 )/86400.d0
	timeout = mjd + fracday
	timeout =  timeout*86400.d0  - timein
	end subroutine  get_time_offset
	end module utility

	module pardiso_multi_real_mpi
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	! a module calling MKL Pardiso to solve symmetric/nonsymmetric real linear equations
	! ----------------------------------------------------------
	! public function protos for calling
	!pardiso_alloc (max_matrix)
	!pardiso_direct_1 (num_m,upper,n_in,nnz_in,ia,ja,va)
	!pardiso_direct_2 (num_m,b_in,x_out,ia,ja,va)
	!pardiso_direct_3 (num_m)
	!pardiso_dealloc
	! ----------------------------------------------------------
	use mkl_cluster_sparse_solver
	use mpi
	use mpi_parameters
	implicit none
	private
	! Internal solver memory pointer pt and iparm vector
	type(MKL_CLUSTER_SPARSE_SOLVER_HANDLE),allocatable :: pt(:,:)
	INTEGER, ALLOCATABLE :: iparm(:,:)
	! sparse matrix variables
	integer, allocatable ::	ia_par(:,:),ja_par(:,:)
	real(8), allocatable :: va_par(:,:)
	! temporal variables
	INTEGER :: idum(1)
	REAL(8) :: ddum(1)
	! other variables
	INTEGER :: maxfct=1, mnum=1, nrhs=1, msglvl=0, error, phase
	integer,allocatable :: mtype(:), n(:), nnz(:)
	! public subroutines
	public :: pardiso_alloc,pardiso_direct_1,pardiso_direct_2,pardiso_direct_3,pardiso_dealloc
	! module subroutine
    contains
    
	subroutine pardiso_alloc (max_matrix)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	! pre allocate variables
	implicit none
	! in/out variables
	integer,intent(in) :: max_matrix
	! process
	if (max_matrix<1) stop "maximum number of matrix error, stop!"
	maxfct = max_matrix
	! allocate internal solver memory pointer pt and iparm vector
	allocate ( pt(maxfct,64) )
	allocate ( iparm(maxfct,64) )
	! allocate other solver variables
	allocate(mtype(maxfct),n(maxfct),nnz(maxfct))
    end subroutine pardiso_alloc
    
	subroutine pardiso_direct_1 (num_m,symm,n_in,nnz_in,ia,ja,va,log)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	! set up initial variables and perform symbolic factorization
	IMPLICIT NONE
	! in/out variables
	integer,intent(in) :: num_m
	logical,intent(in) :: symm
	integer,intent(in) :: n_in,nnz_in
	INTEGER,ALLOCATABLE,intent(in) :: ia(:),ja(:)
	REAL(8),ALLOCATABLE,intent(in) :: va(:)
	logical,optional,intent(in) :: log
	integer :: i
	! local variables
	logical :: checker
	! check the input matrix and assign mtype_par
	if(rank==0) then
		call pardiso_checker (symm,n_in,nnz_in,ia,ja,checker)
		if (.not.checker) stop "matrix error!"
	end if
	if (symm) then
		mtype(num_m) = 2
	else
		mtype(num_m) = 11
	end if
	! read in n_par (row of A), nnz_par (non-zero of A) and matrix A
	n(num_m) = n_in
	nnz(num_m) = nnz_in
	! fill in pt_par and iparm_par with default values
	iparm(num_m,:) = 0
	pt(num_m,:)%dummy = 0
	! reordering and factorization
	phase = 12
	if (present(log)) then
		error = -1
	else
		call cluster_sparse_solver (pt(num_m,:),1,mnum,mtype(num_m),phase,n(num_m),va,ia,ja,idum,nrhs,iparm(num_m,:),msglvl,ddum,ddum,MPI_COMM_WORLD,error)
	end if
	!call cluster_sparse_solver (pt(num_m,:),1,mnum,mtype(num_m),phase,n(num_m),va_par(:,num_m),ia_par(:,num_m),ja_par(:,num_m),idum,nrhs,iparm(num_m,:),msglvl,ddum,ddum,MPI_COMM_WORLD,error)
	!CALL pardiso (pt(num_m,:),1,mnum,mtype(num_m),phase,n(num_m),va,ia,ja,idum,nrhs,iparm(num_m,:),msglvl,ddum,ddum,error)
	IF (error/=0 .and. mtype(num_m)==11) THEN
		if(rank==0) print*, 'ERROR - The following error was detected during analysis step of pardiso solver: ', error
		stop
	else if (error/=0 .and. mtype(num_m)==2) then
		! try indefinite matrix
		mtype(num_m) = -2
		call cluster_sparse_solver (pt(num_m,:),1,mnum,mtype(num_m),phase,n(num_m),va,ia,ja,idum,nrhs,iparm(num_m,:),msglvl,ddum,ddum,MPI_COMM_WORLD,error)
		!CALL pardiso (pt(num_m,:),1,mnum,mtype(num_m),phase,n(num_m),va,ia,ja,idum,nrhs,iparm(num_m,:),msglvl,ddum,ddum,error)
		! error occurs
		IF (error /= 0) THEN
			if(rank==0) print*, 'ERROR - The following error was detected during analysis step of pardiso solver: ', error
			stop
		else
			if(rank==0.and..not.present(log)) print*, 'WARNING - SYSTEM MATRIX IS INDEFINITE, CHANGE TO INDEFINITE MODE.'
		END IF
	END IF
    END subroutine pardiso_direct_1
    
	subroutine pardiso_direct_2 (num_m,b_in,x_out,ia,ja,va)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	! solve matrix, can be called repeatedly with different b_in
	implicit none
	! in/out variables
	integer,intent(in) :: num_m
	real(8),allocatable :: b_in(:), x_out(:)
	INTEGER,ALLOCATABLE,intent(in) :: ia(:),ja(:)
	REAL(8),ALLOCATABLE,intent(in) :: va(:)
	! only solving
	phase = 33
	call cluster_sparse_solver (pt(num_m,:),1,mnum,mtype(num_m),phase,n(num_m),va,ia,ja,idum,nrhs,iparm(num_m,:),msglvl,b_in,x_out,MPI_COMM_WORLD,error)
	!CALL pardiso (pt(num_m,:),1,mnum,mtype(num_m),phase,n(num_m),va,ia,ja,idum,nrhs,iparm(num_m,:),msglvl,b_in,x_out,error)
	IF (error/= 0) THEN
		if(rank==0) print*, 'ERROR - The following error was detected during solve step of pardiso solver: ', error
		stop
    END IF
    END subroutine pardiso_direct_2
    
	subroutine pardiso_direct_3 (num_m)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	! release internal memory
	implicit none
	! in/out variables
	integer,intent(in) :: num_m
	! process
	phase = -1
	call cluster_sparse_solver (pt(num_m,:),1,mnum,mtype(num_m),phase,n(num_m),ddum,idum,idum,idum,nrhs,iparm(num_m,:),msglvl,ddum,ddum,MPI_COMM_WORLD,error)
	!CALL pardiso (pt(num_m,:),1,mnum,mtype(num_m),phase,n(num_m),ddum,idum,idum,idum,nrhs,iparm(num_m,:),msglvl,ddum,ddum,error)
	IF (error /= 0) THEN
		if(rank==0) print*, 'ERR - The following error was detected during terminal step of pardiso solver: ',error
		stop
	END IF
    end subroutine pardiso_direct_3
    
	subroutine pardiso_dealloc
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	! decallocate all woring arrays of PARDISO after all calls to pardiso
	implicit none
	IF(ALLOCATED(pt)) DEALLOCATE(pt)
	IF(ALLOCATED(iparm)) DEALLOCATE(iparm)
	IF(ALLOCATED(mtype)) DEALLOCATE(mtype)
	IF(ALLOCATED(n)) DEALLOCATE(n)
	IF(ALLOCATED(nnz)) DEALLOCATE(nnz)
	if(rank==0) then
		IF(ALLOCATED(ia_par)) deallocate(ia_par,ja_par,va_par)
	end if
    END subroutine pardiso_dealloc
    
	subroutine pardiso_checker (sym,n_in,nnz_in,ia_in,ja_in,checker)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	USE MKL_SPARSE_HANDLE
	implicit none
	! variables
	logical,intent(in) :: sym
	integer,intent(in) :: n_in,nnz_in
	INTEGER,intent(in) :: ia_in(n_in+1),ja_in(nnz_in)
	logical,intent(out) :: checker
	TYPE(SPARSE_STRUCT) :: PT
	integer :: CHECK_RESULT_CODE
	! process
	CALL SPARSE_MATRIX_CHECKER_INIT(PT)
	PT % N = n_in
	PT % CSR_IA = LOC(ia_in)
	PT % CSR_JA = LOC(ja_in)
	PT % INDEXING         = MKL_ONE_BASED
	if (sym) then
		PT % MATRIX_STRUCTURE = MKL_UPPER_TRIANGULAR
	else
		PT % MATRIX_STRUCTURE = MKL_GENERAL_STRUCTURE
	end if
	PT % PRINT_STYLE      = MKL_C_STYLE
	PT % MESSAGE_LEVEL    = MKL_PRINT
	CHECK_RESULT_CODE = SPARSE_MATRIX_CHECKER(PT)
	IF (CHECK_RESULT_CODE/=MKL_SPARSE_CHECKER_SUCCESS) THEN
		WRITE(*,'("Matrix check details: (",i0,", ",i0,", ",i0,")")') PT % CHECK_RESULT
		WRITE(*,*)'matrix check failed, stop!'
		checker = .false.
	else
		checker = .true.
	end if
	end subroutine pardiso_checker
	end module pardiso_multi_real_mpi

	module source_mod
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	!source: ind source_type(1 electric, 2 magnetic) pulse_type (1 Turn-on step 2 Gaussian pulse 3 delta pulse) t0 t1 current
	!   51    1    1    1.0e-7    0.0    1.0
	!   52    1    2    5.0e-9    0.0    1.e-11
	!   53    1    3    1.0e-7    0.0    1.0
	! pulse_type 1: t0 - the width is t0
	! (0.35875-0.48829*cos(2.0*pi*t/t0)+0.14128*cos(4.0*pi*t/t0)-0.01168*cos(6.0*pi*t/t0))/(0.35875*t0)
	! pulse_type 2: t0 - the width is 2*t0
	! exp(-4*pi*(t-t0)^2/t0^2)
	! pulse_type 3: t0 - the width is t0
	! 1/abs(t0)/sqrt(pi)*exp(-((t*4*exp(1)-2*exp(1)*t0)/t0)^2)
	implicit none
	private
	!public :: time_const
	public :: time_const_p
    contains
    
	subroutine time_const_p (source_type,pulse_type,t,t0,cur,tmp)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	!logical,intent(in) :: deriv
	integer,intent(in) :: source_type,pulse_type
	real(8),intent(in) :: t,t0,cur
	real(8),intent(out) :: tmp
	if      (source_type==1 .and. pulse_type==1) then
		tmp = erfstepon (1,t,t0) * cur
		!else if (source_type==2 .and. pulse_type==1) then
		!tmp = erfstepon (1,t,t0) * cur
	else if (source_type==1 .and. pulse_type==2) then
		tmp = deltafunc (1,t,t0) * cur
		!else if (source_type==2 .and. pulse_type==2) then
		!tmp = gaussian(1,t,t0) * cur
	else if (source_type==1 .and. pulse_type==3) then
		tmp = -erfstepon (1,t,t0) * cur
		!else if (source_type==2 .and. pulse_type==3) then
		!tmp = -erfstepon (2,t,t0) * cur
		!else if (source_type==1 .and. pulse_type==4) then
		!tmp = gaussian(2,t,t0) * cur
		!else if (source_type==2 .and. pulse_type==4) then
		!tmp = gaussian(1,t,t0) * cur
	else
		stop 'no such kind of source'
	end if
    end subroutine time_const_p
    
	real(8) function erfstepon (deriv,t,t0)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	integer,intent(in) :: deriv
	real(8),intent(in) :: t,t0
	real(8),parameter :: pi = 3.141592653589793D0
	if (deriv==0) then
		stop "error - no representation for 0 order erf function!"
		!erfstepon = 1.d0
	else if (deriv==1) then
		erfstepon = (1.0D0/sqrt(3.141592653589793D0)*exp(1.0D0/t0**2*(t-t0/2.0D0)&
			&**2*(-2.5D+1))*5.0D0)/t0
	else if (deriv==2) then
		erfstepon = 1.0D0/t0**3*1.0D0/sqrt(3.141592653589793D0)*exp(1.0D0/t0**2*(&
			&t-t0/2.0D0)**2*(-2.5D+1))*(t*2.0D0-t0)*(-1.25D+2)
	else if (deriv==3) then
		erfstepon = 1.0D0/t0**3*1.0D0/sqrt(3.141592653589793D0)*exp(1.0D0/t0**2*(&
			&t-t0/2.0D0)**2*(-2.5D+1))*(-2.5D+2)+1.0D0/t0**5*1.0D0/sqrt(3.14159&
			&2653589793D0)*exp(1.0D0/t0**2*(t-t0/2.0D0)**2*(-2.5D+1))*(t*2.0D0-&
			&t0)**2*3.125D+3
	end if
    end function erfstepon
    
	real(8) function gaussian (deriv,t,t0)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	integer,intent(in) :: deriv
	real(8),intent(in) :: t,t0
	if (deriv==0) then
		gaussian = exp(1.0D0/t0**2*3.141592653589793D0*(t-t0/2.0D0)**2*(-1.6D+1)&
			&)
	else if (deriv==1) then
		gaussian = 1.0D0/t0**2*3.141592653589793D0*exp(1.0D0/t0**2*3.14159265358&
			&9793D0*(t-t0/2.0D0)**2*(-1.6D+1))*(t*2.0D0-t0)*(-1.6D+1)
	else if (deriv==2) then
		gaussian = 1.0D0/t0**2*3.141592653589793D0*exp(1.0D0/t0**2*3.14159265358&
			&9793D0*(t-t0/2.0D0)**2*(-1.6D+1))*(-3.2D+1)+1.0D0/t0**4*3.14159265&
			&3589793D0**2*exp(1.0D0/t0**2*3.141592653589793D0*(t-t0/2.0D0)**2*(&
			&-1.6D+1))*(t*2.0D0-t0)**2*2.56D+2
	else if (deriv==3) then
		gaussian = 1.0D0/t0**6*3.141592653589793D0**3*exp(1.0D0/t0**2*3.14159265&
			&3589793D0*(t-t0/2.0D0)**2*(-1.6D+1))*(t*2.0D0-t0)**3*(-4.096D+3)+1&
			&.0D0/t0**4*3.141592653589793D0**2*exp(1.0D0/t0**2*3.14159265358979&
			&3D0*(t-t0/2.0D0)**2*(-1.6D+1))*(t*2.0D0-t0)*5.12D+2+1.0D0/t0**4*3.&
			&141592653589793D0**2*exp(1.0D0/t0**2*3.141592653589793D0*(t-t0/2.0&
			&D0)**2*(-1.6D+1))*(t*8.0D0-t0*4.0D0)*2.56D+2
	else
		stop "gaussian error"
	end if
    end function gaussian
    
	real(8) function deltafunc (deriv,t,t0)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
	integer,intent(in) :: deriv
	real(8),intent(in) :: t,t0
	if (deriv==0) then
		deltafunc = exp(-1.0D0/t0**2*(t*5.436563656918091D0-t0*2.718281828459046D&
			&0)**2)
	else if (deriv==1) then
		deltafunc = -1.0D0/t0**2*exp(-1.0D0/t0**2*(t*5.436563656918091D0-t0*2.718&
			&281828459046D0)**2)*(t*5.911244879144521D+1-t0*2.955622439572261D+&
			&1)
	else if (deriv==2) then
		deltafunc = 1.0D0/t0**2*exp(-1.0D0/t0**2*(t*5.436563656918091D0-t0*2.7182&
			&81828459046D0)**2)*(-5.911244879144521D+1)+1.0D0/t0**4*exp(-1.0D0/&
			&t0**2*(t*5.436563656918091D0-t0*2.718281828459046D0)**2)*(t*5.9112&
			&44879144521D+1-t0*2.955622439572261D+1)**2
	else if (deriv==3) then
		deltafunc = 1.0D0/t0**4*exp(-1.0D0/t0**2*(t*5.436563656918091D0-t0*2.7182&
			&81828459046D0)**2)*(t*5.911244879144521D+1-t0*2.955622439572261D+1&
			&)*5.911244879144521D+1+1.0D0/t0**4*exp(-1.0D0/t0**2*(t*5.436563656&
			&918091D0-t0*2.718281828459046D0)**2)*(t*6.988563204242466D+3-t0*3.&
			&494281602121233D+3)-1.0D0/t0**6*exp(-1.0D0/t0**2*(t*5.436563656918&
			&091D0-t0*2.718281828459046D0)**2)*(t*5.911244879144521D+1-t0*2.955&
			&622439572261D+1)**3
	else
		stop "deltafunc error"
	end if
	end function deltafunc
	end module source_mod

	module system_lib
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	implicit none
    contains
    
	function itoa(i) result(res)
    ! TEMF3DT Copyright (C) 2021-2022 Xiaodong Yang
    ! This program comes with ABSOLUTELY NO WARRANTY.
    ! This is free software, and you are welcome to redistribute it under certain conditions.
	character(:),allocatable :: res
	integer,intent(in) :: i
	character(range(i)+2) :: tmp
	write(tmp,'(i0)') i
	res = trim(tmp)
    end function
	end module system_lib