module m_const
  implicit none
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  integer, public :: jz,num_interaction

end module m_const


module m_sp
  implicit none
  integer, public :: np,nn,np_orb,nn_orb,tot_orbs,ntot
  integer ,public,dimension(:),pointer :: obs

  TYPE, PUBLIC :: single_particle_descript
     INTEGER :: total_orbits
     INTEGER, DIMENSION(:), POINTER :: nn, ll, jj, mm, itzp, nshell
  END TYPE single_particle_descript
  type(single_particle_descript) all_orbits
CONTAINS
  SUBROUTINE allocate_sp_array(this_array,n)
    type(single_particle_descript),intent(inout):: this_array
    integer, intent(inout):: n
    integer ::i
    IF (ASSOCIATED (this_array%mm) ) DEALLOCATE(this_array%mm)
    ALLOCATE(this_array%mm(n))
    IF (ASSOCIATED (this_array%nn) ) DEALLOCATE(this_array%nn)
    ALLOCATE(this_array%nn(n))
    IF (ASSOCIATED (this_array%ll) ) DEALLOCATE(this_array%ll)
    ALLOCATE(this_array%ll(n))
    IF (ASSOCIATED (this_array%jj) ) DEALLOCATE(this_array%jj)
    ALLOCATE(this_array%jj(n))
    IF (ASSOCIATED (this_array%itzp) ) DEALLOCATE(this_array%itzp)
    ALLOCATE(this_array%itzp(n))
    IF (ASSOCIATED (this_array%nshell) ) DEALLOCATE(this_array%nshell)
    ALLOCATE(this_array%nshell(n))
    !           blank all characters and zero all other values
    DO i= 1, n
       this_array%nn(i)=0
       this_array%mm(i)=0
       this_array%ll(i)=0
       this_array%jj(i)=0
       this_array%nshell(i)=0
       this_array%itzp(i)=0
    ENDDO

  END SUBROUTINE allocate_sp_array

  SUBROUTINE deallocate_sp_array(this_array)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    DEALLOCATE(this_array%nn) ; DEALLOCATE(this_array%ll)
    DEALLOCATE(this_array%jj) ;DEALLOCATE(this_array%itzp)
    DEALLOCATE(this_array%nshell)
    DEALLOCATE(this_array%mm)
  END SUBROUTINE deallocate_sp_array


end module m_sp


module c_funcs
  use iso_c_binding
  implicit none
interface
    subroutine qsort(array,elem_count,elem_size,compare) bind(C,name="qsort")
    use iso_c_binding
    implicit none

      type(c_ptr),value       :: array
      integer(c_size_t),value :: elem_count
      integer(c_size_t),value :: elem_size
      type(c_funptr),value    :: compare !int(*compare)(const void *, const void *)
    end subroutine qsort !standard C library qsort
  end interface
end module c_funcs

module m_int
  use c_funcs
  implicit none
  type, bind(c),public ::  mat
    integer(C_INT128_T) :: idx
    real(C_DOUBLE) :: val
  end type 
     interface operator(<=)
        module procedure LessOrEqual
     end interface operator(<=)


  integer:: n1b,n2b,n3b,no1b,no2b,no3b
  type(mat),allocatable,target  :: h1b(:),h2b(:),h3b(:)
  type(mat),allocatable,target  :: o1b(:),o2b(:),o3b(:)



  contains

        function compare(v1, v2) bind(C)
           integer(C_INT) compare
           type(C_PTR), value :: v1, v2
           type(mat), pointer :: p1, p2
  
           call C_F_POINTER(v1,p1)
           call C_F_POINTER(v2,p2)
           if(.NOT.(p1 <= p2)) then
              compare = 1
           else if(.NOT.(p2 <= p1)) then
              compare = -1
           else
              compare = 0
           end if

        end function compare

    function LessOrEqual(v1, v2)
       logical LessOrEqual
       type(mat), intent(in) :: v1, v2
    
       LessOrEqual = v1%idx <= v2%idx
    end function LessOrEqual


    subroutine ini_mat(m,k)
      integer, intent(in) :: m,k
      integer :: i
      integer*16 :: ia,ib,ic,id,ie,ik
      if(k==1)then
        write(6,*)'initializing onebd interaction',m
        allocate(h1b(m))
        n1b=m

        do i=1,m
       h1b(i)%idx=0
       h1b(i)%val=0.d0
     enddo
        do i=1,m
        read(8,*)ia,ib,h1b(i)%val
        ia=ishft(ia,7)
        h1b(i)%idx=ia+ib
        enddo
       close(8)
      else if(k==2)then
        write(6,*)'initializing twobd interaction',m
        allocate(h2b(m))
        n2b=m

        do i=1,m
       h2b(i)%idx=0
       h2b(i)%val=0.d0
     enddo
        do i=1,m
        read(9,*)ia,ib,ic,id,h2b(i)%val
        ia=ishft(ia,21)
        ib=ishft(ib,14)
        ic=ishft(ic,7)
        h2b(i)%idx=ia+ib+ic+id
!        write(*,*)ia+ib+ic+id
        enddo
       close(9)
      else if(k==3)then
        write(6,*)'initializing three-bd interaction',m
        allocate(h3b(m))
        n3b=m

        do i=1,m
       h3b(i)%idx=0
       h3b(i)%val=0.d0
     enddo
        do i=1,m
          read(10,*)ia,ib,ic,id,ie,ik,h3b(i)%val
        ia=ishft(ia,35)
        ib=ishft(ib,28)
        ic=ishft(ic,21)
        id=ishft(id,14)
        ie=ishft(ie,7)
        h3b(i)%idx=ia+ib+ic+id+ie+ik
        enddo
       close(10)

      elseif(k==4)then
        write(6,*)'initializing onebd interaction',m
        allocate(o1b(m))
        no1b=m

        do i=1,m
       o1b(i)%idx=0
       o1b(i)%val=0.d0
     enddo
        do i=1,m
        read(18,*)ia,ib,o1b(i)%val
        ia=ishft(ia,7)
        o1b(i)%idx=ia+ib
        enddo
       close(18)
      else if(k==5)then
        write(6,*)'initializing twobd interaction',m
        allocate(o2b(m))
        no2b=m

        do i=1,m
       o2b(i)%idx=0
       o2b(i)%val=0.d0
     enddo
        do i=1,m
        read(19,*)ia,ib,ic,id,o2b(i)%val
        ia=ishft(ia,21)
        ib=ishft(ib,14)
        ic=ishft(ic,7)
        o2b(i)%idx=ia+ib+ic+id
!        write(*,*)ia+ib+ic+id
        enddo
       close(19)
      else if(k==6)then
        write(6,*)'initializing three-bd interaction',m
        allocate(o3b(m))
        no3b=m

        do i=1,m
      o3b(i)%idx=0
      o3b(i)%val=0.d0
     enddo

        do i=1,m
          read(20,*)ia,ib,ic,id,ie,ik,o3b(i)%val
        ia=ishft(ia,35)
        ib=ishft(ib,28)
        ic=ishft(ic,21)
        id=ishft(id,14)
        ie=ishft(ie,7)
        o3b(i)%idx=ia+ib+ic+id+ie+ik
        enddo
       close(20)



      endif


    end  subroutine ini_mat





    subroutine fetch_mat(m,k,val)
      integer*16, intent(in) :: m
      integer,intent(in) :: k
      real*8, intent(inout):: val
      type(mat), pointer :: pt(:)
      integer :: low,high,mid,flag
      val=0.d0
      low=1
      flag=0
      select case(k)
      case(1)
        pt=>h1b
        high=n1b
      case(2)
        pt=>h2b
        high=n2b
      case(3)
        pt=>h3b
        high=n3b
      case(4)
        pt=>o1b
        high=no1b
      case(5)
        pt=>o2b
        high=no2b
      case(6)
        pt=>o3b
        high=no3b
      end select
      if(pt(1)%idx==m)then
        val=pt(1)%val
        return
      endif
      if(pt(high)%idx==m)then
        val=pt(high)%val
        return
      endif


      do
      mid=(low+high)/2
      if(mid==low)exit

      if(pt(mid)%idx==m)then
        val=pt(mid)%val
        exit
      elseif(pt(mid)%idx>m)then
        high=mid
      elseif(pt(mid)%idx<m)then
        low=mid
      endif
      
      enddo
      return
      



    end subroutine fetch_mat


    subroutine sort_mat(k)

           integer(c_size_t) :: n
           integer,intent(in):: k



           if(k==1)then
             n=n1b
             call qsort(c_loc(h1b(1)),n,sizeof(h1b(1)),c_funloc(compare))
           else if(k==2)then
             n=n2b
             call qsort(c_loc(h2b(1)),n,sizeof(h2b(1)),c_funloc(compare))
           else if(k==3)then
             n=n3b
             call qsort(c_loc(h3b(1)),n,sizeof(h3b(1)),c_funloc(compare))
           elseif(k==4)then
             n=no1b
             call qsort(c_loc(o1b(1)),n,sizeof(o1b(1)),c_funloc(compare))
           else if(k==5)then
             n=no2b
             call qsort(c_loc(o2b(1)),n,sizeof(o2b(1)),c_funloc(compare))
           else if(k==6)then
             n=no3b
             call qsort(c_loc(o3b(1)),n,sizeof(o3b(1)),c_funloc(compare))
           endif
   
  !  write(*,*)'after', n2b
  !  do n=1,n2b
  !    write(*,*)h2b(n)%idx,h2b(n)%val
  !  enddo

   !
    end subroutine sort_mat

    subroutine very_mat
      integer :: i,m
      integer*16 :: ia,ib,ic,id,ie,ik,idx
      real*8:: tmp, tmp1

      !open(unit=10,file=threebd_int)

      read(10,*)m
        do i=1,m
          read(10,*)ia,ib,ic,id,ie,ik,tmp
        ia=ishft(ia,35)
        ib=ishft(ib,28)
        ic=ishft(ic,21)
        id=ishft(id,14)
        ie=ishft(ie,7)
        idx=ia+ib+ic+id+ie+ik
        call fetch_mat(idx,3,tmp1)
        if(abs(tmp-tmp1)>1E-10)then
          write(*,*)ia,ib,ic,id,ie,ik,idx
          stop 'error'
        endif
        enddo


    end  subroutine very_mat



  end module m_int




  module m_file
  implicit none
      character(len=120):: onebd_int,twobd_int,threebd_int,&
                           onebd_opr,twobd_opr,threebd_opr,&
        & output,sp_file
  end module m_file

  module bit_op
  implicit none

    contains
subroutine decode(a)
  implicit none
  integer*16, intent(in):: a
  write(*,'(b12.12)')a
end subroutine decode

end  module bit_op


module combine
  use m_sp
  implicit none
  integer :: c(20)
  integer, dimension(:,:), pointer ::cn2,cn3,cm2,cn3t
  integer,public :: ncn2,ncm2,ncn3
  contains
    subroutine combin2( n, k, icnk )
       implicit none
  !
    integer i
    integer icnk
    integer k
    integer mn
    integer mx
    integer n
  !
    mn = min ( k, n-k )
  
    if ( mn < 0 ) then
  
      icnk = 0
  
    else if ( mn == 0 ) then
  
      icnk = 1
  
    else
  
      mx = max ( k, n-k )
      icnk = mx + 1
  
      do i = 2, mn
        icnk = ( icnk * ( mx + i ) ) / i
      end do
  
    end if
  
    return
  end subroutine combin2

  subroutine combin3(cn,n)
    integer,intent(inout),dimension(:,:),pointer:: cn
    integer,intent(in):: n
    integer :: i,k,m
       m=0
    do i=1,n-1
      do k=i+1,n
        m=m+1
        cn(1,m)=i
        cn(2,m)=k
      enddo
    enddo


    end subroutine combin3

  subroutine ini_combin
    integer :: n,k,num
    integer :: i,l,m
  ! number of C(N,2)
     n=nn+np
     k=2
     call combin2(n,k,num)
     ncn2=num
     allocate(cn2(2,num))
     cn2=0
     call combin3(cn2,n)



     num=0
     n=n-2
     call combin2(n,k,num)
     allocate(cm2(2,num))
     ncm2=num
     cm2=0
     call combin3(cm2,n)

     n=n+2
     m=0
     do i=1,n-2
     num=0
       call combin2(n-i,k,num)
       m=m+num
     enddo

     allocate(cn3(3,m))
     ncn3=m
     cn3=0
     m=0
     do i=1,n-2
     num=0
       call combin2(n-i,k,num)
       allocate(cn3t(2,num));cn3t=0
       call combin3(cn3t,n-i)
       cn3t=cn3t+i
       !write(*,*)cn3t(1,1),cn3t(2,1)

       do l=m+1,m+num
         cn3(1,l)=i
         cn3(2,l)=cn3t(1,l-m)
         cn3(3,l)=cn3t(2,l-m)
        enddo
       m=m+num
       deallocate(cn3t)
     enddo






    end subroutine ini_combin


  end module combine


module mbasis
    use bit_op
    use m_sp
  implicit none
  integer, public:: num_mbasis,np_mbasis,nn_mbasis
  integer*16,dimension(:),pointer:: mbsp,mbsn,idmaskn,idmaskp
  integer,dimension(:),pointer :: mzp,mzn,mz
  integer,dimension(:,:),pointer :: mbs
  complex*16, dimension(:,:),allocatable:: hm,jm,vec,vecl,om
  complex*16,dimension(:),allocatable :: em,oem
        
   contains 
        function ch_sign(m,nt,lt)
          integer*16,intent(in) ::m,nt,lt
          integer*16:: n,l
          integer :: ch_sign
          integer*16 :: sam_n,sam_p,sam
          n=nt; l=lt
          if(n>l)write(*,*)'faltal error in ch_sign',n,l
          !write(*,*)'nl',n,l
          if(l-n<=1)then
            ch_sign=0
          elseif(l<=nn_orb)then

  ! nn interaction
          sam_n=mbs(m,2)
          sam_n=mbsn(sam_n)
          sam=ibits(sam_n,n,l-n-1)
!  write(*,*)sam,sam_n,n,l
  ch_sign=popcnt(sam)
!  write(*,*)'num 1 in sam',ch_sign

 elseif(n>nn_orb)then
  ! pp interaction
  !write(*,*)'nl',n,l
          ch_sign=0
          sam_p=mbs(m,1)
          sam_p=mbsp(sam_p)
          n=n-nn_orb;
          l=l-nn_orb
  sam=ibits(sam_p,n,l-n-1)
  ch_sign=popcnt(sam)

        else
          ! np
          sam_n=mbs(m,2)
          sam_n=mbsn(sam_n)
          sam=ibits(sam_n,n,nn_orb-n)
          ch_sign=popcnt(sam)
          if(l==nn_orb+1)return
          sam_p=mbs(m,1)
          sam_p=mbsp(sam_p)
          sam=ibits(sam_p,0,l-nn_orb-1)
          ch_sign=popcnt(sam)+ch_sign
         endif
          return

      end function
end module mbasis
