module cplx
  implicit none
  type,public:: dc
    real*8:: a,b
    contains
       procedure :: disp =>display
  end type

  interface operator(+)
    module procedure add
  end interface
  interface operator(-)
    module procedure minus
  end interface
  interface operator(*)
    module procedure mult
  end interface
  interface operator(/)
    module procedure div
  end interface

  contains
    subroutine display(this)
      implicit none
      class(dc):: this
        write(*,*)this%a,this%b
    end subroutine display

    type(dc) function add(a,b)
      type(dc),intent(in)::a,b
      add%a=a%a+b%a
      add%b=a%b+b%b
      end function
      
    type(dc) function minus(a,b)
      type(dc),intent(in)::a,b
      minus%a=a%a-b%a
      minus%b=a%b-b%b
      end function

    type(dc) function mult(a,b)
      type(dc),intent(in)::a,b
      mult%a=a%a*b%a
      mult%b=a%b*b%b
      end function

    type(dc) function div(a,b)
      type(dc),intent(in)::a,b
      div%a=a%a/b%a
      div%b=a%b/b%b
      end function
end module cplx

program main
  use cplx
  implicit none
  type(dc):: a,b,c
  a%a=1
  a%b=2
  b%a=3
  b%b=4
  c=a+b
  write(*,*)'a+b='
  call c%disp
  c=a-b
  write(*,*)'a-b='
  call c%disp
  c=a*b
  write(*,*)'a*c='
  call c%disp
  c=a/b
  write(*,*)'a/c='
  call c%disp
  end program
