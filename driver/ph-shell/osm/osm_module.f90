module parallel
  include'mpif.h'
  INTEGER, PUBLIC   :: ierror,iam,num_procs,master
end module parallel

module constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  integer,public::  jtot
  integer,public::  nu_parity
  integer,public:: number_interactions
  integer,public :: nstate
  integer, public :: max_iter_lanc
  REAL(DP) , PARAMETER, PUBLIC :: p_mass =938.926_dp
  REAL(DP), PARAMETER, PUBLIC :: hbarc = 197.326968_dp
  REAL(DP), PARAMETER, PUBLIC :: hb2ip = hbarc*hbarc/p_mass
  REAL(DP), PUBLIC, PARAMETER :: pi = 3.141592741012573_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_2 = 1.570796370506287_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_4 = 0.7853981852531433_dp
end module constants


module mpi_mapping
  integer,public:: n_submat  !! the threads are divided to nxn
  integer,public:: sub_row,sub_col
  integer,allocatable,public:: mapping(:,:)
  integer,allocatable,public:: mapping_1d(:,:)
end module mpi_mapping

module storage
  TYPE, PUBLIC :: block_storage
     REAL*8, DIMENSION(:), ALLOCATABLE :: val1
     REAL*8, DIMENSION(:,:), ALLOCATABLE :: val
     REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: val3
     complex*16, DIMENSION(:), ALLOCATABLE :: cval1
     complex*16, DIMENSION(:,:), ALLOCATABLE :: cval
     complex*16, DIMENSION(:,:,:), ALLOCATABLE :: cval3
  END TYPE block_storage
  
  TYPE, PUBLIC :: integer_storage
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: ival
     INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ival3
  END TYPE integer_storage



  real*8,dimension(:,:),allocatable:: vt1
  type(block_storage),dimension(:),allocatable:: vt2,vt3
  integer,dimension(:,:,:),allocatable:: ip2
  integer,dimension(:,:,:,:),allocatable:: ip3

  integer,dimension(:,:),allocatable:: gmat2_channels,gmat3_channels

end module storage

module m_sp
  implicit none
  integer, public :: np,nn,np_orb,nn_orb,tot_orbs,ntot
  integer ,public,dimension(:),pointer :: obs

  TYPE, PUBLIC :: single_particle_descript
     INTEGER :: total_orbits
     INTEGER, DIMENSION(:), POINTER :: nn, ll, jj, mm, itzp, nshell,occ,jord
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
    IF (ASSOCIATED (this_array%occ) ) DEALLOCATE(this_array%occ)
    ALLOCATE(this_array%occ(n))
    IF (ASSOCIATED (this_array%jord) ) DEALLOCATE(this_array%jord)
    ALLOCATE(this_array%jord(n))
    !           blank all characters and zero all other values
    DO i= 1, n
       this_array%nn(i)=0
       this_array%mm(i)=0
       this_array%ll(i)=0
       this_array%jj(i)=0
       this_array%nshell(i)=0
       this_array%itzp(i)=0
       this_array%occ(i)=0
       this_array%jord(i)=0
    ENDDO

  END SUBROUTINE allocate_sp_array

  SUBROUTINE deallocate_sp_array(this_array)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    DEALLOCATE(this_array%nn) ; DEALLOCATE(this_array%ll)
    DEALLOCATE(this_array%jj) ;DEALLOCATE(this_array%itzp)
    DEALLOCATE(this_array%nshell)
    DEALLOCATE(this_array%mm)
    DEALLOCATE(this_array%occ)
    DEALLOCATE(this_array%jord)
  END SUBROUTINE deallocate_sp_array


end module m_sp

module tiny_functions
  implicit none
  contains
  function ij_bits(m,nt,lt)
    !! return number of occ bits bettwen nt(not include) and lt(not include)
    implicit none
    integer*16,intent(in):: m
    integer,intent(in) ::nt,lt
    integer:: n,l
    integer :: ij_bits
    integer*16 ::sam
    n=nt; l=lt
    if(n>l)write(*,*)'faltal error in ij_bits',n,l

    if(n+1==l)then
      ij_bits=0
      return
    endif

    sam=ibits(m,n,l-n-1)
    ij_bits=popcnt(sam)
    return

  end function

    SUBROUTINE diagonalize_complex(h, eval, evect)
        COMPLEX*16, DIMENSION(:,:) :: h
        ! Variables  needed for dgeev
        CHARACTER*1 :: jobvl, jobvr
        INTEGER :: info, lwork, n, i, j
        COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: wr, work, eval
        COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: vr, vl, evect
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rwork

        eval = 0.0D0; evect = 0.0D0
        n = SIZE(eval)

        jobvl = 'N'; jobvr = 'V'
        lwork = 8*n
        ALLOCATE(wr(n), vr(n,n), vl(1,n), work(lwork))
        ALLOCATE(rwork(2*n))

        CALL ZGEEV(jobvl, jobvr, n, h, n, wr, vl, 1, vr, n, &
            work, lwork, rwork, info)

        ! Sort eigenvalues (lowest first)
        eval(1) = wr(1)
        evect(:,1) = vr(:,1)
        DO i = 2, n
            DO j = i-1, 0, -1
                IF (j == 0) THEN
                    eval(j+1) = wr(i)
                    evect(:,j+1) = vr(:,i)
                ELSE IF (REAL(wr(i)) < REAL(eval(j))) THEN
                    eval(j+1) = eval(j)
                    evect(:,j+1) = evect(:,j)
                ELSE
                    eval(j+1) = wr(i)
                    evect(:,j+1) = vr(:,i)
                    EXIT
                ENDIF
            ENDDO
        ENDDO

        DEALLOCATE(wr, vr, vl, work, rwork)
    END SUBROUTINE diagonalize_complex






  end module tiny_functions


module lanczos
  use constants
  use m_sp
  use tiny_functions
  use mpi_mapping
  use parallel
  implicit none
  integer*16,allocatable,public::mbsn(:),mbsp(:)
  integer,allocatable,public::mzn(:),mzp(:),iparn(:),iparp(:),mbsp_valid(:), mbsn_valid(:)
  integer,public:: np_mbasis,nn_mbasis,num_mbasis
  integer,allocatable,public::mbs_row(:,:),mbs_col(:,:)
  complex*16,dimension(:,:),allocatable:: haml,cvec,cvecl,vec,and_haml,and_haml2,and_vec,and_lvec
  complex*16,dimension(:),allocatable:: vi,vj,vk,and_eig,and_eig2

  
contains
  function ch_sign(m,n,nt,lt)
    !! m= neutron n=proton
    integer:: ch_sign
    integer*16,intent(in):: m,n
    integer,intent(in)::nt,lt
    integer:: l_nt,l_lt
    if(nt>lt)then
      write(*,*)'error in ch_sign'
      stop
    endif
    if(nt+1==lt)then
      ch_sign=0
      return
    endif

    if(lt<=nn_orb)then
      ch_sign=ij_bits(m,nt,lt)
      return
    endif
    if(nt>nn_orb)then
      l_nt=nt-nn_orb
      l_lt=lt-nn_orb

      ch_sign=ij_bits(n,l_nt,l_lt)
      return
    else
      l_lt=lt-nn_orb
    ch_sign=ij_bits(n,0,l_lt)+ij_bits(m,nt,nn_orb+1)
    return
    endif

    end function ch_sign

   subroutine ini_startv
     implicit none
      integer:: bra_min,bra_max
      real*8,allocatable::tmp(:)
      complex*16:: norm

      bra_min=mapping_1d(iam+1,1)
      bra_max=mapping_1d(iam+1,2)
      allocate(tmp(bra_min:bra_max))
      tmp=0.d0

     ! call random_number(tmp(bra_min:bra_max))
     ! write(*,*)tmp
    if(iam==0)tmp(1)=1.d0
      vec(bra_min:bra_max,1)=dcmplx(tmp)

      norm=0.0_dpc
      norm=dot_product(vec(bra_min:bra_max,1),vec(bra_min:bra_max,1))

    call mpi_allreduce(mpi_in_place,norm,1,mpi_double_complex,mpi_sum, &
      mpi_comm_world,ierror)

      norm=sqrt(norm)
      vec(bra_min:bra_max,1)=vec(bra_min:bra_max,1)/norm

      norm=dot_product(vec(bra_min:bra_max,1),vec(bra_min:bra_max,1))
     
    call mpi_allreduce(mpi_in_place,vec(:,1),num_mbasis,mpi_double_complex,mpi_sum, &
      mpi_comm_world,ierror)
    vi=vec(:,1)

   end subroutine ini_startv



    subroutine mpv_ax
      implicit none
      integer:: bra_min,bra_max,ket_min,ket_max,dk,dm,dn
      vj=0.0_dpc


   bra_min=mapping(sub_row,1)
   bra_max=mapping(sub_row,2)
   ket_min=mapping(sub_col,1)
   ket_max=mapping(sub_col,2)
   dm=bra_max-bra_min+1
   dn=1
   dk=ket_max-ket_min+1
  call zgemm('n','n',dm,dn,dk,1.d0,haml(bra_min:bra_max,ket_min:ket_max),&
    dm,vi(ket_min:ket_max),dk, 0.d0,vj(bra_min:bra_max),dm)
   call mpi_barrier(mpi_comm_world,ierror)
    call mpi_allreduce(mpi_in_place,vj,num_mbasis,mpi_double_complex,mpi_sum, &
      mpi_comm_world,ierror)

    end subroutine mpv_ax

    subroutine mpv_orth(k,conv)
      implicit none
      integer,intent(in)::k
      logical,intent(inout)::conv
      integer:: bra_min,bra_max
      integer:: dm,dn,dk,i
      complex*16:: norm
      conv=.FALSE.

      bra_min=mapping_1d(iam+1,1)
      bra_max=mapping_1d(iam+1,2)

      do i=1,k
      norm=dot_product(vj(bra_min:bra_max),vec(bra_min:bra_max,i))
    call mpi_allreduce(mpi_in_place,norm,1,mpi_double_complex,mpi_sum, &
      mpi_comm_world,ierror)
      and_haml(i,k)=norm
      vj=vj-norm*vec(:,i)
      enddo
      
      norm=dot_product(vj(bra_min:bra_max),vj(bra_min:bra_max))
    call mpi_allreduce(mpi_in_place,norm,1,mpi_double_complex,mpi_sum, &
      mpi_comm_world,ierror)
      and_haml(k+1,k)=sqrt(norm)
      if(abs(norm)<1E-16)conv=.TRUE.
     vj=vj/sqrt(norm)

     vec(:,k+1)=vj
     vi=vj
    end subroutine mpv_orth




end module lanczos

module file_list
  implicit none
      character(len=120):: onebd_int,twobd_int,threebd_int,&
                           onebd_opr,twobd_opr,threebd_opr,&
        & output,sp_file
end module file_list






module m_configurations
  use parallel
  use m_sp
  type m_configuration_descriptor
    integer:: number_confs
    integer, dimension(:),pointer:: config_ab
    integer, dimension(:),pointer:: config_abc
    end type m_configuration_descriptor

    contains

  subroutine number_vm_confs(ij,ipar,itz,struct,this)
    use constants
    use m_sp
    implicit none
    TYPE (m_configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, na, la, nb, lb, ja, jb, a, b, b_end, nconfs, &
         itza, itzb, struct,mb,ma
    LOGICAL triag

    

    nconfs=0
    DO a=1, all_orbits%total_orbits
       na = all_orbits%nn(a)
       la=all_orbits%ll(a)
       ja=all_orbits%jj(a)
       ma=all_orbits%mm(a)
       itza=all_orbits%itzp(a)
     
       
          b_end = all_orbits%total_orbits
       DO b=1, b_end
          nb = all_orbits%nn(b)
          lb=all_orbits%ll(b)
          itzb=all_orbits%itzp(b)
          jb=all_orbits%jj(b)
          mb=all_orbits%mm(b)
          
          if ( a == b ) cycle 
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          if ( ma + mb /= 2*ij ) cycle 
          if(a==b)cycle
          
       
          nconfs=nconfs+1
       ENDDO
    ENDDO
    this%number_confs=nconfs


  end subroutine number_vm_confs

  subroutine setup_vm_confs(ij,ipar,itz,struct,this)
    use constants
    use m_sp
    implicit none

    TYPE (m_configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, na, la, nb, lb, ja, jb, a, b, b_end, nconfs, &
         ma,mb,k1, k2, itza, itzb, struct
    LOGICAL triag
    nconfs=0
    DO a = 1, all_orbits%total_orbits
       la=all_orbits%ll(a)
       ja=all_orbits%jj(a)
       na = all_orbits%nn(a)
       itza=all_orbits%itzp(a)
       ma=all_orbits%mm(a)
       
          b_end = all_orbits%total_orbits
       DO b = 1, b_end 
          lb=all_orbits%ll(b)
          jb=all_orbits%jj(b)
          nb = all_orbits%nn(b)
          itzb=all_orbits%itzp(b)
          mb=all_orbits%mm(b)
          
          if ( a == b ) cycle 
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          if ( ma + mb /= 2*ij ) cycle 
          if(a==b)cycle
          
          
          
          nconfs=nconfs+1
          k2=nconfs*2
          k1=k2-1
          this%config_ab(k1)=b
          this%config_ab(k2)=a
       ENDDO
    ENDDO
    IF ( nconfs /= this%number_confs ) THEN
       WRITE(6,*) ' Error in configuration allocation ' ; STOP
    ENDIF

  end subroutine setup_vm_confs


  SUBROUTINE number_vm3_confs(jj,ipar,itz,struct,this)
    use constants
    USE m_sp
    
    IMPLICIT NONE
    TYPE (m_configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: jj, ipar, itz, la, lb, lc, ja, jb, jc, na, nb, nc, &
         a, b, c, nconfs, itza, itzb, itzc, j_ab, &
         j_ab_min, j_ab_max, j_bc, j_bc_min, j_bc_max, b_end, c_end, struct, count
    INTEGER :: orb_low, orb_high, j_ac, j_ac_min, j_ac_max, count_ac, count_bc, ma,mb,mc
    DOUBLE PRECISION ::  wf_norm, factor
    LOGICAL triag
    CHARACTER(LEN = 100 ) :: coup_order
    nconfs=0

       orb_low = 1; orb_high = all_orbits%total_orbits
    
    
    DO a = orb_low, orb_high
              
       la=all_orbits%ll(a)
       ja=all_orbits%jj(a)
       na=all_orbits%nn(a)
       ma=all_orbits%mm(a)
       itza=all_orbits%itzp(a)
       
       DO b = orb_low, orb_high 
          
          lb=all_orbits%ll(b)
          jb=all_orbits%jj(b)
          nb=all_orbits%nn(b)
          mb=all_orbits%mm(b)
          itzb=all_orbits%itzp(b)
          
          DO c = orb_low, orb_high             
             
             lc=all_orbits%ll(c)
             jc=all_orbits%jj(c)
             nc=all_orbits%nn(c)
             mc=all_orbits%mm(c)
             itzc=all_orbits%itzp(c)
                          
             ! check on isospin and parity of three particle state
             IF ( ma + mb + mc /= jj ) CYCLE
             IF ( itzc+itza+itzb /= itz ) CYCLE
             IF ( ipar /= MOD ( la+lb+lc, 2) ) CYCLE
             if(a==b)cycle
             if(b==c)cycle
             if(a==c)cycle
             
             
             nconfs=nconfs+1
             
             
          ENDDO
       ENDDO
    ENDDO
    
    this%number_confs=nconfs

  END SUBROUTINE number_vm3_confs


  SUBROUTINE setup_vm3_confs(jj,ipar,itz,struct,this)
    use constants
    USE m_sp
    
    IMPLICIT NONE
    TYPE (m_configuration_descriptor), INTENT(INOUT) :: this
    DOUBLE PRECISION :: wf_norm
    INTEGER :: jj, ipar, itz, la, lb, lc, ja, jb, jc, na, nb, nc, &
         a, b, c, nconfs, k1, k2, k3, k4, itza, itzb, itzc, struct
    INTEGER :: orb_low, orb_high, j_ac, j_ac_min, j_ac_max,  ma,mb,mc
    LOGICAL triag
    nconfs=0
    
    
    !WRITE(6,*) ' Three-body configurations for J, TZ, and parity: ',jj,itz,ipar
    
       orb_low = 1; orb_high = all_orbits%total_orbits
    
    DO a = orb_low, orb_high
              
       la=all_orbits%ll(a)
       ja=all_orbits%jj(a)
       na=all_orbits%nn(a)
       ma=all_orbits%mm(a)
       itza=all_orbits%itzp(a)
       
       DO b = orb_low, orb_high 
          
          lb=all_orbits%ll(b)
          jb=all_orbits%jj(b)
          nb=all_orbits%nn(b)
          mb=all_orbits%mm(b)
          itzb=all_orbits%itzp(b)
          
          DO c = orb_low, orb_high             
             
             lc=all_orbits%ll(c)
             jc=all_orbits%jj(c)
             nc=all_orbits%nn(c)
             mc=all_orbits%mm(c)
             itzc=all_orbits%itzp(c)
             
             ! check on isospin and parity of three particle state
             IF ( itzc+itza+itzb /= itz ) CYCLE
             IF ( ipar /= MOD ( la+lb+lc, 2) ) CYCLE
             IF ( ma + mb + mc /= jj ) CYCLE

             
             if(a==b)cycle
             if(b==c)cycle
             if(a==c)cycle
             nconfs=nconfs+1
                   
             k3=nconfs*3
             k2=nconfs*3-1
             k1=nconfs*3-2
             
             this%config_ab(k1)=a
             this%config_ab(k2)=b
             this%config_ab(k3)=c
                          
          END DO
       END DO
    END DO

    IF ( nconfs /= this%number_confs ) THEN
       WRITE(6,*) ' Error in configuration allocation ' ; STOP
    ENDIF
  END SUBROUTINE setup_vm3_confs



end module m_configurations

module truncations
  use  m_sp
  use lanczos
  implicit none
  integer:: number_spj_orbs
  integer,allocatable::occ_lim(:,:)
  integer*16,allocatable:: occ_mask(:),obtype(:)
  logical:: trunc
  contains
    subroutine setup_occ_mask
      implicit none
      integer:: i,j,bit_ahead,obsize,itz
      integer*16:: one,bit_mask
      integer:: it_tmp,ob_start
      bit_ahead=0

      one=1
      do i=1,number_spj_orbs
        obsize=0
        bit_mask=0
        bit_ahead=0
        do j=1,all_orbits%total_orbits
          if(all_orbits%jord(j)/=i)then
            bit_ahead=bit_ahead+1
          else
            ob_start=j
            obsize=all_orbits%jj(j)+1
            it_tmp=all_orbits%itzp(j)
            obtype(i)=it_tmp
            exit
          endif
        enddo
          bit_mask=2**obsize-1

        if(it_tmp==-1)then
          occ_mask(i)=ishft(bit_mask,bit_ahead-nn_orb)
        else
          occ_mask(i)=ishft(bit_mask,bit_ahead)
        endif
        write(*,*)obsize,ob_start,bit_ahead,popcnt(occ_mask(i)),occ_mask(i),it_tmp


      enddo
      write(*,*)'number jj orbs:', number_spj_orbs
      write(*,*)occ_mask

    end subroutine setup_occ_mask

    logical function check_occ(a,b)
      integer,intent(in):: a,b
      integer*16:: tp_long
      integer:: i,n
      check_occ=.true.
      do i=1,number_spj_orbs
        if(obtype(i)==-1)then
          tp_long=mbsp(a)
        else
          tp_long=mbsn(b)
        endif
          n=popcnt(iand(tp_long,occ_mask(i)))
          if(n<occ_lim(1,i).or. n> occ_lim(2,i))then       
            check_occ=.FALSE.
            exit
          endif
      enddo
      return
      end function check_occ


end module




