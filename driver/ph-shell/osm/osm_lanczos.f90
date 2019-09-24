subroutine sm_lanczos
  use mpi_mapping
  use lanczos
  use m_sp
  integer,parameter :: seed = 86456
  integer:: i,k,cu_dim
  logical:: conv
  real*8:: cor,cor1
  

  if(num_mbasis<500)then
          call no_lanc
  else
   cor1=100000.d0

  !! intitial vector
  cu_dim=0
  call ini_startv
  do i=1,num_mbasis-1
    if(iam==0)write(*,*)i,'th arnoldi iteration ..'
  call mpv_ax
  call mpv_orth(i,conv)
    cu_dim=i
  if(conv)exit
  if(i==nstate)and_eig2=and_eig
  if(i>nstate+1)then
      cor=0.d0
 and_haml2=and_haml
  call  lapack_diag(and_haml2(1:cu_dim,1:cu_dim), &
    and_vec(1:cu_dim,1:cu_dim), and_lvec(1:cu_dim,1:cu_dim),and_eig(1:cu_dim), cu_dim )
     do k=1,nstate
       cor=abs(and_eig(k)-and_eig2(k))+cor
     enddo
     if((abs(cor-cor1)<1E-6) .or. (i==max_iter_lanc-1))exit
     and_eig2=and_eig
     cor1=cor
     endif

  enddo
  do i=1,nstate
    if(iam==0)write(*,*)i,'th state:',and_eig(i)
  enddo
  endif



  

end subroutine sm_lanczos

subroutine no_lanc
  use mpi_mapping
  use lanczos
  use m_sp
  implicit none
  COMPLEX*16, DIMENSION(num_mbasis,num_mbasis) :: lcvec, lcvecl
  COMPLEX*16, DIMENSION(num_mbasis) :: lceig
  integer:: i
  call lapack_diag(haml,lcvec,lcvecl,lceig,num_mbasis)
  do i=1,nstate
    if(iam==0)write(*,*)i,'th state:',lceig(i)
  enddo


end subroutine no_lanc


SUBROUTINE lapack_diag(h, cvec, cvecl,ceig, n )


  implicit none
  integer, intent(in) :: n
  complex*16, dimension(n,n), intent(in) :: h
  COMPLEX*16, DIMENSION(n,n), intent(out) :: cvec, cvecl
  COMPLEX*16, DIMENSION(n,n) ::  vl,vr
  COMPLEX*16, DIMENSION(n), intent(out) :: ceig
  DOUBLE PRECISION, DIMENSION(2*n) :: rwork
  COMPLEX*16, DIMENSION(10000) :: work1
  INTEGER :: lda, ldvl, ldvr, info, lwork
  CHARACTER*1 :: jobvl, jobvr
  complex*16 :: norm
  integer :: i,j

  jobvl = 'V' ;  jobvr = 'V';  lda = n
  ldvl = n;  ldvr = n;  lwork = 10000
  ceig = 0.; cvec = 0.; cvecl = 0.
!  write(*,*)h
  CALL zgeev( jobvl, jobvr, n, h, lda, ceig, cvecl, ldvl, cvec, ldvr, &
    work1, lwork, rwork, info )

  ! berggren normalization
  do i = 1, n
    norm = sum( cvec(:,i)*cvec(:,i) )
    cvec(:,i) = cvec(:,i)/sqrt(norm)
    norm = sum( cvecl(:,i)*cvec(:,i) )
    !write(*,*)norm
    if(real(norm) < 0.d0 )then
    cvecl(:,i)=-cvecl(:,i)/sqrt(-norm)
    else
    cvecl(:,i)=cvecl(:,i)/sqrt(norm)
    endif
  
     
  end do

do i=1,n
   do j=1,n
    norm=sum(cvecl(:,i)*cvecl(:,j))
   ! write(6,*)'ss',i,j, norm
enddo
enddo

  call hf_sort_cmplx2(ceig,cvec,cvecl, n)


end SUBROUTINE lapack_diag


SUBROUTINE hf_sort_cmplx2(a,b,c, n)
  IMPLICIT NONE

  INTEGER :: i, j, n
  complex*16, DIMENSION(n), INTENT(INOUT) :: a
  complex*16, DIMENSION(n,n), INTENT(INOUT) :: b,c
  complex*16 :: temp1, temp2
  complex*16, DIMENSION(n) :: temp3

  DO i = 1, n
     DO j = 1, n
        IF ( real( a(i) )  < real( a(j) ) ) THEN

           temp1 = a(i)
           a(i) = a(j)
           a(j) = temp1

           temp3(:) = b(:,i)
           b(:,i) = b(:,j)
           b(:,j) = temp3(:)


           temp3(:) = c(:,i)
           c(:,i) = c(:,j)
           c(:,j) = temp3(:)

        END IF
     END DO
  END DO

END SUBROUTINE hf_sort_cmplx2
