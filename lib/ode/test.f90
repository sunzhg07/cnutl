program main
implicit none
integer :: neq
real*8, allocatable :: y(:), ydot(:)
external deriv, deriv1
integer :: nstep,istep
real*8 :: x1,x2


! FOR TESTING DLSODE
 INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL,  LIW, LRW,&
         MF
 DOUBLE PRECISION  ATOL(2), RTOL
real*8,allocatable:: rwork(:)
integer,allocatable:: iwork(:)

! RK variables
Real*8                 :: eps
double precision ::  atol1


allocate(y(2),ydot(2))
nstep=50




        NEQ = 2
        ITOL = 2
        RTOL = 1.D-4
        ATOL(1) = 1.D-6
        ATOL(2) = 1.D-6
        ITASK = 1
        ISTATE = 1
        IOPT = 0

        MF = 10
        LRW = 20+16*neq
        LIW = 20
allocate(rwork(LRW))
allocate(iwork(LIW))

 

x1=0.d0
x2=0.d0
y=0.d0
ydot=0.d0
y(1)=sqrt(1.5d0)
y(2)=sqrt(3.d0)

do istep=1,nstep
x1=x2
x2=x2+0.1

!! solve ode with DLSODE, compile with dlsode.f

!    CALL DLSODE (deriv, NEQ, Y, x1, x2, ITOL, RTOL, ATOL, ITASK, &
!          ISTATE, IOPT, RWORK, LRW, IWORK, LIW, 'dum', MF)

!! solve ode with rk4 , compile with rk45.f
call rkf45(deriv1, neq, y,ydot, x1,x2,RTOL,ATOL1, ITASK)
write(*,'(I4, 2x, 3(F13.6,2x))')istep, x2, y(1),y(2)

enddo

end program main


subroutine deriv1(t,y,ydot)
implicit none

real*8,intent(in):: t
real*8,intent(in):: y(2)
real*8,intent(inout):: ydot(2)

real*8:: r,s,w,g,e,rd,sd,u,v
real*8:: aa(2,2),xx(2),bb(2)
e=0.5
g=-1e2
w=100

u=y(1)
v=y(2)


r=0.5*cos(t)
rd=-0.5*sin(t)
s=cos(w*t)
sd=-w*sin(w*t)

aa(1,1)=g
aa(1,2)=e
aa(2,1)=e
aa(2,2)=-1

xx(1)=(-1.d0+u**2.d0-r)/2.d0/u
xx(2)=(-2.d0+v**2.d0-s)/2.d0/v
bb(1)=rd/2.d0/u
bb(2)=sd/2.d0/sqrt(2.+s)

ydot=matmul(aa,xx)+bb


end subroutine deriv1

subroutine deriv(neq,t,y,ydot)
implicit none

integer,intent(in):: neq
real*8,intent(in):: t
real*8,intent(in):: y(neq)
real*8,intent(inout):: ydot(neq)

real*8:: r,s,w,g,e,rd,sd,u,v
real*8:: aa(2,2),xx(2),bb(2)
e=0.5
g=-1e2
w=100

u=y(1)
v=y(2)


r=0.5*cos(t)
rd=-0.5*sin(t)
s=cos(w*t)
sd=-w*sin(w*t)

aa(1,1)=g
aa(1,2)=e
aa(2,1)=e
aa(2,2)=-1

xx(1)=(-1.d0+u**2.d0-r)/2.d0/u
xx(2)=(-2.d0+v**2.d0-s)/2.d0/v
bb(1)=rd/2.d0/u
bb(2)=sd/2.d0/sqrt(2.+s)

ydot=matmul(aa,xx)+bb


end subroutine deriv


