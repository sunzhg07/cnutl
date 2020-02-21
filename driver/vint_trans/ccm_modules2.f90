

MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  ! min and max isospin projection
  INTEGER, PUBLIC :: itzmin, itzmax
  ! min and max total two-body angular momentum in lab frame
  INTEGER, PUBLIC :: j_lab_min, j_lab_max
  ! min and max total two-body angular momentum in Rel-CoM frame
  INTEGER, PUBLIC :: jmin, jmax, n_nucleus, z_nucleus
  INTEGER, PUBLIC :: occ_protons, occ_neutrons
  real*8, public :: com_beta,com_switch, j2_beta, eom_energy
  REAL*8, public :: hcom_value,occnum3_cut,e0,o0,total_energy, eccsd
  INTEGER, PUBLIC :: mass_nucleus, tot_orbs, below_ef, above_ef,nlmax_cut, & 
       e3max_cut, l3max_cut, j3_cut, j2_max
  REAL(DP), PUBLIC ::  oscl, hbar_omega, efermi
  REAL(DP) , PARAMETER, PUBLIC :: p_mass = 938.926_dp
  REAL(DP) , PARAMETER, PUBLIC :: g_A = 1.29d0
  REAL(DP) , PARAMETER, PUBLIC :: fpi = 92.4d0
  REAL(DP), PARAMETER, PUBLIC :: theta_rot = 0.0_dp! 0.125_dp
  REAL(DP), PARAMETER, PUBLIC :: hbarc = 197.326968_dp
  REAL(DP), PARAMETER, PUBLIC :: hb2ip = hbarc*hbarc/p_mass
  REAL(DP), PUBLIC, PARAMETER :: pi = 3.141592741012573_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_2 = 1.570796370506287_dp
  REAL(DP), PUBLIC, PARAMETER :: pi_4 = 0.7853981852531433_dp
  LOGICAL :: switch_density, switch_triples
  CHARACTER(LEN=100) :: operator_type, gs_type
  CHARACTER(LEN=40) :: clem_j
  integer, public :: stime=0 ! whether show the time of running
  integer,public:: run_unit
  
END MODULE constants



MODULE REORDER 

  INTEGER, ALLOCATABLE, public :: INDX(:), indx_inv(:)

END MODULE REORDER




MODULE single_particle_orbits
  TYPE, PUBLIC :: single_particle_descript
     INTEGER :: total_orbits
     INTEGER, DIMENSION(:), POINTER :: nn, ll, jj, itzp, nshell
     CHARACTER (LEN=10), DIMENSION(:), POINTER :: orbit_status, orb_type, model_space
     REAL*8, DIMENSION(:), POINTER :: e, occ_numbers
  END TYPE single_particle_descript
  TYPE (single_particle_descript), PUBLIC :: all_orbit, neutron_data, &
       proton_data
  
CONTAINS
  SUBROUTINE allocate_sp_array(this_array,n)
    TYPE (single_particle_descript), INTENT(INOUT) :: this_array
    INTEGER , INTENT(IN) :: n
    integer :: i 

    IF (ASSOCIATED (this_array%nn) ) DEALLOCATE(this_array%nn)
    ALLOCATE(this_array%nn(n))
    IF (ASSOCIATED (this_array%ll) ) DEALLOCATE(this_array%ll)
    ALLOCATE(this_array%ll(n))
    IF (ASSOCIATED (this_array%jj) ) DEALLOCATE(this_array%jj)
    ALLOCATE(this_array%jj(n))
    IF (ASSOCIATED (this_array%itzp) ) DEALLOCATE(this_array%itzp)
    ALLOCATE(this_array%itzp(n))
    IF (ASSOCIATED (this_array%e) ) DEALLOCATE(this_array%e)
    ALLOCATE(this_array%e(n))
    IF (ASSOCIATED (this_array%occ_numbers) ) DEALLOCATE(this_array%occ_numbers)
    ALLOCATE(this_array%occ_numbers(n))
    IF (ASSOCIATED (this_array%nshell) ) DEALLOCATE(this_array%nshell)
    ALLOCATE(this_array%nshell(n))
    IF (ASSOCIATED (this_array%orbit_status) ) DEALLOCATE(this_array%orbit_status)
    ALLOCATE(this_array%orbit_status(n))
    IF (ASSOCIATED (this_array%model_space) ) DEALLOCATE(this_array%model_space)
    ALLOCATE(this_array%model_space(n))
    !           blank all characters and zero all other values
    DO i= 1, n
       this_array%model_space(i)= ' '
       this_array%orbit_status(i)= ' '
       this_array%e(i)=0.
       this_array%occ_numbers(i)=0.
       this_array%nn(i)=0
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
    DEALLOCATE(this_array%e) ;DEALLOCATE(this_array%nshell)
    DEALLOCATE(this_array%orbit_status); DEALLOCATE(this_array%model_space)
    DEALLOCATE(this_array%occ_numbers)
  END SUBROUTINE deallocate_sp_array




END MODULE single_particle_orbits





MODULE gmat_storage
  
  
  TYPE, PUBLIC :: block_storage
     COMPLEX*16, DIMENSION(:), ALLOCATABLE :: val1
     REAL*8, DIMENSION(:,:), ALLOCATABLE :: val
     REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: val3
     COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: cval
     REAL*8, DIMENSION(:), ALLOCATABLE  :: angmom_factor 
     INTEGER(8), DIMENSION(:), ALLOCATABLE :: nconf_high, nconf_low
     INTEGER, DIMENSION(:), ALLOCATABLE :: r3bra, r3ket, r3channel, vbra, vchannel, vket, r2bra, r2ket, r2channel, j2, r1a, r1i
  END TYPE block_storage
  
  TYPE, PUBLIC :: integer_storage
     INTEGER :: number_confs
     INTEGER, DIMENSION(:), ALLOCATABLE :: ival1
     INTEGER, DIMENSION(:,:), ALLOCATABLE :: ival
     INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: ival3
     INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: ival4
  END TYPE integer_storage

  type(block_storage),dimension(:),allocatable :: gmatrix
    
  TYPE (integer_storage), ALLOCATABLE, PUBLIC :: lookup_pq_configs(:),pq_configs(:)
  integer,dimension(:,:),allocatable:: quantum_numbers_ch3
  integer,dimension(:,:,:,:),allocatable:: locate_gmatchannel
  real*8,dimension(:,:),allocatable:: tkin

END MODULE gmat_storage



!     module which defines configuration type, general structure
!     either for lab frame case or relcm system
MODULE configurations
  USE single_particle_orbits
  
  TYPE configuration_descriptor        
     INTEGER :: number_confs, model_space_confs
     INTEGER, DIMENSION(:), pointer :: config_ab, config_abc
     INTEGER, DIMENSION(:), pointer :: config_J_ipar_tz
     INTEGER, DIMENSION(:), pointer :: channel_rabij
     INTEGER, DIMENSION(:), pointer :: channel_r2vbb
     INTEGER, DIMENSION(:), pointer :: channel_t3
     INTEGER, DIMENSION(:), pointer :: channel_rabcijk
     INTEGER, DIMENSION(:), pointer :: channel_rabcijk_titus
     INTEGER, DIMENSION(:), pointer :: channel_2vbb_3p3h
     INTEGER, DIMENSION(:), pointer :: channel_rabi, channel_rabcij, channel_rabijk
     INTEGER, DIMENSION(:), pointer :: channel_rkcj
     INTEGER, DIMENSION(:), pointer :: channel_rcki
     INTEGER, DIMENSION(:), pointer :: channel_rbcd
     INTEGER, DIMENSION(:), pointer :: channel_raic
     INTEGER, DIMENSION(:), pointer :: channel_rmij
     INTEGER, DIMENSION(:), pointer :: channel_rabij_cross
     INTEGER, DIMENSION(:), pointer :: channel_rabi_cross
     INTEGER, DIMENSION(:), pointer :: channel_rijb
     INTEGER, DIMENSION(:), pointer :: channel_rijb_cross
     INTEGER, DIMENSION(:), pointer :: channel_r3_reduced
     INTEGER, DIMENSION(:), pointer :: channel_l3_reduced
     
  END TYPE configuration_descriptor
  
  TYPE (configuration_descriptor), PUBLIC :: channels
  TYPE (configuration_descriptor), PUBLIC :: ph_channels
  
CONTAINS

  
  !
  !     setting up all configurations for given J, Tz and parity for
  !     the single particle model space
  !
  SUBROUTINE number_sp_confs(ij,ipar,itz,this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, a, nconfs

    nconfs=0
    DO a=1, all_orbit%total_orbits
       IF ( MOD(all_orbit%ll(a),2) /= ipar ) CYCLE
       IF ( all_orbit%jj(a) /= ij) CYCLE
       IF ( all_orbit%itzp(a) /= itz ) CYCLE
       nconfs=nconfs+1
    ENDDO
    this%number_confs=nconfs

  END SUBROUTINE number_sp_confs
  !
  !     Allocates space for the sp configurations and sets them up
  !
  SUBROUTINE setup_sp_configurations(ij,ipar,itz,this)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: a, ij, ipar, itz, nconfs
    nconfs=0
    DO a = 1, all_orbit%total_orbits
       IF ( MOD(all_orbit%ll(a),2) /= ipar ) CYCLE
       IF ( all_orbit%jj(a) /= ij) CYCLE
       IF ( all_orbit%itzp(a) /= itz ) CYCLE
       nconfs=nconfs+1
       this%config_ab(nconfs)=a
    ENDDO
    IF ( nconfs /= this%number_confs ) THEN
       WRITE(6,*) ' Error in configuration allocation ' ; STOP
    ENDIF

  END SUBROUTINE setup_sp_configurations



  !                       
  !     setting up all configurations for given J, Tz and parity for
  !     the g-matrix model space

  SUBROUTINE number_gmatrix_confs(ij,ipar,itz,struct,this)
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, na, la, nb, lb, ja, jb, a, b, b_end, nconfs, &
         itza, itzb, struct
    LOGICAL triag

    

    nconfs=0
    DO a=1, all_orbit%total_orbits
       na = all_orbit%nn(a)
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       itza=all_orbit%itzp(a)
     
       if ( struct == 1 .and. all_orbit%orbit_status(a) /= 'hole' ) cycle 
       if ( struct == 2 .and. all_orbit%orbit_status(a) /= 'particle' ) cycle 
       if ( struct == 3 .and. all_orbit%orbit_status(a) /= 'particle' ) cycle 
       IF ( itz /= 0 .and. struct == 4 ) THEN
          b_end = a
       ELSE
          b_end = all_orbit%total_orbits
       ENDIF
       DO b=1, b_end
          nb = all_orbit%nn(b)
          lb=all_orbit%ll(b)
          itzb=all_orbit%itzp(b)
          jb=all_orbit%jj(b)
          
          IF ( itzb > itza .and. struct == 4 ) CYCLE
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
          !if ( all_orbit%orbit_status(b) == 'particle' .and. & 
          !     all_orbit%orbit_status(a) == 'hole' ) cycle 
!!$          if ( abs(all_orbit%occ_numbers(a) ) < occnum3_cut ) cycle 
!!$          if ( abs(all_orbit%occ_numbers(b) ) < occnum3_cut ) cycle 
          
          
          if ( struct == 1 .and. all_orbit%orbit_status(b) /= 'hole' ) cycle 
          if ( struct == 2 .and. all_orbit%orbit_status(b) /= 'hole' ) cycle 
          if ( struct == 3 .and. all_orbit%orbit_status(b) /= 'particle' ) cycle 
          nconfs=nconfs+1
       ENDDO
    ENDDO
    this%number_confs=nconfs

  END SUBROUTINE number_gmatrix_confs
  !
  !
  !
  SUBROUTINE setup_gmatrix_configurations(ij,ipar,itz,struct,this)
    USE constants
    USE single_particle_orbits
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: ij, ipar, itz, na, la, nb, lb, ja, jb, a, b, b_end, nconfs, &
         k1, k2, itza, itzb, struct
    LOGICAL triag
    nconfs=0


    DO a = 1, all_orbit%total_orbits
       la=all_orbit%ll(a)
       ja=all_orbit%jj(a)
       na = all_orbit%nn(a)
       itza=all_orbit%itzp(a)
       
       if ( struct == 1 .and. all_orbit%orbit_status(a) /= 'hole' ) cycle 
       if ( struct == 2 .and. all_orbit%orbit_status(a) /= 'particle' ) cycle 
       if ( struct == 3 .and. all_orbit%orbit_status(a) /= 'particle' ) cycle 
       
       
       IF ( itz /= 0 .and. struct == 4  ) THEN
          b_end = a
       ELSE
          b_end = all_orbit%total_orbits
       ENDIF
       DO b = 1, b_end 
          lb=all_orbit%ll(b)
          jb=all_orbit%jj(b)
          nb = all_orbit%nn(b)
          itzb=all_orbit%itzp(b)
          
          
          IF ( itzb > itza .and. struct == 4 ) CYCLE
          IF (itzb+itza /= itz*2 ) CYCLE
          IF ( (-1)**(la+lb+ipar) < 0) CYCLE
          IF ( triag( 2*ij, ja, jb) ) CYCLE
          IF ( (a == b) .AND. (MOD(ij,2) /= 0) ) CYCLE
!!$          if ( abs(all_orbit%occ_numbers(a) ) < occnum3_cut ) cycle 
!!$          if ( abs(all_orbit%occ_numbers(b) ) < occnum3_cut ) cycle 

          
          if ( struct == 1 .and. all_orbit%orbit_status(b) /= 'hole' ) cycle 
          if ( struct == 2 .and. all_orbit%orbit_status(b) /= 'hole' ) cycle 
          if ( struct == 3 .and. all_orbit%orbit_status(b) /= 'particle' ) cycle 
          !if ( all_orbit%orbit_status(b) == 'particle' .and. & 
          !     all_orbit%orbit_status(a) == 'hole' ) cycle 
          
          
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
    CALL sort_configs(this,nconfs)

!!$    IF( struct == 0 ) then
!!$       SELECT CASE (itz )
!!$       CASE ( -1)
!!$          if ( iam == 0 ) WRITE(6,*) nconfs,' proton-proton configurations for J',ij
!!$       CASE (1)
!!$          if ( iam == 0 ) WRITE(6,*) nconfs,' neutron-neutron configurations for J ',ij
!!$       CASE (0)
!!$          if ( iam == 0 ) WRITE(6,*) nconfs,' proton-neutron configurations for J ',ij
!!$       END SELECT
!!$    end IF

  END SUBROUTINE setup_gmatrix_configurations



  SUBROUTINE sort_configs(this,n_confs)
    IMPLICIT NONE
    TYPE (configuration_descriptor), INTENT(INOUT) :: this
    INTEGER :: i, j, tempa, tempb, n_confs

    DO i = 1, n_confs
       DO j =i, n_confs
          IF( (this%config_ab(2*i-1)*100 +this%config_ab(2*i)) >  &
               (this%config_ab(2*j-1)*100 +this%config_ab(2*j)) ) THEN 
             tempa=this%config_ab(2*i-1)
             tempb=this%config_ab(2*i)
             this%config_ab(2*i-1)=this%config_ab(2*j-1)
             this%config_ab(2*i)=this%config_ab(2*j)
             this%config_ab(2*j-1)=tempa
             this%config_ab(2*j)=tempb                             
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE sort_configs
  

  
END MODULE configurations

module transe
  integer, public:: amass
 CHARACTER (LEN=100):: sp_file,kin_file,vint_file,snt_file

end module transe
