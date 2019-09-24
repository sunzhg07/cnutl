!--------------------------------------------------------------------                                                                                                                
FUNCTION calc_objective(work,mean,work_add,istart,istop,where_add)
  IMPLICIT NONE
  REAL :: calc_objective
  REAL :: work(1:istop)
  REAL :: mean, work_add
  INTEGER :: istart,istop,where_add
  INTEGER :: i

  calc_objective = 0.0
  DO i=istart,istop
     IF(i.NE.where_add)THEN
       calc_objective=calc_objective + ABS(work(i)-mean)
     ELSE
       calc_objective=calc_objective+ ABS(work(i) + work_add - mean)
     ENDIF
  ENDDO
end FUNCTION calc_objective

! sort routine 

subroutine sort_channel( c, b, n )
  implicit none 
  integer, intent(in) :: n
  integer, intent(inout) :: c(n), b(n)
  integer :: i,j, tmp1, tmp2


  do i = 1, n 
     do j = i+1, n 
        if ( c(j) < c(i) ) then
           tmp1 = c(j)
           c(j) = c(i)
           c(i) = tmp1
           
           tmp2 = b(j) 
           b(j) = b(i)
           b(i) = tmp2
           

        end if
     end do
  end do
  
end subroutine sort_channel


!     swaps values of 2 integers
SUBROUTINE SWAP_AB(a, b)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: a, b
  INTEGER :: c
  c = a; a = b; b = c
END SUBROUTINE SWAP_AB

!
!     Function to check triangular relations     
!
LOGICAL FUNCTION triag(i,j,k)
  IMPLICIT NONE
  INTEGER :: i, j, k
  triag = ((i-j-k)*(i-ABS(j-k)) > 0)

END FUNCTION triag
!
!      Function to calculate norm of g-mat    
!
real*8 FUNCTION dij(ja,jb)
  IMPLICIT NONE
  INTEGER :: ja, jb
  IF(ja == jb ) THEN
     dij=sqrt(2.d0)
  ELSE
     dij=1.0
  ENDIF
  
END FUNCTION dij

!
!     Function to calculate phase factors (-1)**(n)
!
INTEGER FUNCTION iph(n)
  IMPLICIT NONE
  INTEGER :: n
  iph=(-1)**n
END FUNCTION iph


FUNCTION delta(a,b)
  IMPLICIT NONE
  INTEGER :: a,b
  REAL*8  :: delta
  delta=0.0
  if(a == b)delta=1.0
end FUNCTION delta

SUBROUTINE zero(n,a)
  IMPLICIT NONE
  INTEGER :: n,i
  REAL*8  :: a(n)
  do i=1,n
    a(i)=0.0
  end do
end subroutine zero



