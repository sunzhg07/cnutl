!             Spherical Coupled Cluster Project
!
!             Author:   Thomas Papenbrock and Gaute Hagen 
!             ADDRESS:  Physics Division, Oak Ridge National Laboratory
!             E-MAIL:   hageng@ornl.gov, tpapenbr@utk.edu
!             LANGUAGE: F90/F95, MPI
!             LAST UPGRADE : September 2007
!
PROGRAM ccsd
  USE ang_mom_functions
  use parallel
  use m_scheme_orbits
  use cc_config_module
  
  implicit none
  real*8  ::  startwtime , endwtime, res, cd, c3,c4
  character(len=100) :: argument, filename
  !TYPE(cc_config) :: config
  
!!$  real*8, allocatable :: sendbuf(:),recvbuf(:)
!!$  integer :: sendcount
!!$  integer, allocatable :: recvcount(:),rdispl(:),nconf_low(:), nconf_high(:)
!!$  integer :: nconfs_tot, number_mtxel_iam,diff,bra_confs,bra_min,bra_max,n,nconfs,i
  
  

  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world,iam,ierror)
  call mpi_comm_size(mpi_comm_world,num_procs,ierror)
  
  master=0
  
  call read_cmd_arguments(argument)
  startwtime = MPI_WTIME()
  
  filename = argument
  call config%read_config_file(filename)
  CALL config%distribute(iam, mpi_comm_world)
  
  CALL commons_to_angmom
  
  if ( iam == 0 ) write(6,*) 'hf-format', config%hf_format
  if ( config%hf_format == 'standard' ) then 
     CALL read_nuclear_sp_data
     call setup_fockmtx_and_interaction
  elseif ( config%hf_format == 'oslo' ) then 
     CALL read_nuclear_sp_data_imsrg
     call setup_fockmtx_and_interaction_imsrg
  end if
  
  call setup_particle_hole
  if ( .not. config%precalc_groundstate ) call ccsd_iter
  if ( config%precalc_groundstate ) call ccsd_stored
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Total execution time for CCSD(T) code', endwtime - startwtime
  
  call mpi_finalize(ierror)

END PROGRAM ccsd

subroutine mytimer(tm)
  implicit none
  real(4):: tm
  integer:: val(8)

  val = 0
  call date_and_time(values=val)
  tm = 0.
  tm = real(val(6))*60.d0 + real(val(7)) + real(val(8))*10.d-3
end subroutine mytimer

SUBROUTINE read_cmd_arguments(argument)
  CHARACTER(len=100) :: argument
  INTEGER :: numargs
  
  numargs = command_argument_count()
  IF (numargs < 1) STOP 'Need configfile as first argument!'
  CALL get_command_argument(1, argument)
  
END SUBROUTINE read_cmd_arguments


!
! read in the single-particle data
! and set up the single-particle orbitals; 
! allocate sp space.
!
SUBROUTINE read_nuclear_sp_data
  USE PARALLEL
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use reorder 
  use gmat_storage, only : scale
  use eom_module, only : Nmax, Nucleus,Interaction
  use cc_config_module
  
  IMPLICIT NONE
  CHARACTER(LEN=100) :: input1, input2
  INTEGER :: ntot_n, ntot_p, junk1, junk2, junk3, junk4, t_z, nn, tza, tzc, pp
  REAL*8 :: re_e1, im_e1, re_hcom, im_hcom
  INTEGER :: i,j, a,b,c, n1, l1, j1, t1, nlmax, number_p, number_n, number_orbs
  real*8 :: delta 
  integer:: ii, iph
  integer, allocatable, dimension(:):: orbtype
  
  
  !
  ! Number of occupied j-shells 
  !
  occ_protons = config%occ_protons
  occ_neutrons = config%occ_neutrons
  mass_nucleus = config%mass_nucleus
  com_beta = config%com_beta
  com_switch = config%com_beta
  j2_beta = config%com_switch
  input1 = config%hf_input_file_orbits
  input2 = config%hf_input_file_onebody
  
  !     Read neutron and proton single-particle data
  !read(5,*);read(5,*) input1
  OPEN(1, FILE = input1) 
  !read(5,*);read(5,*) input2
  OPEN(2, FILE = input2)
  
  if ( iam == 0 ) write(6,*) 'sp-energy file', input1
  if ( iam == 0 ) write(6,*) 'kinetic-energy file', input2
  
  !
  ! Count number of j-scheme orbitals
  !
  read(1,*) 
  ntot_n = 0
  ntot_p = 0
1 READ(1,*,end=2) junk1, junk2, junk3, junk4, t_z
  
  if ( t_z == -1 ) ntot_p = ntot_p + 1
  if ( t_z ==  1 ) ntot_n = ntot_n + 1
  
  GOTO 1
2 CONTINUE
  REWIND 1
  
  proton_data%total_orbits = ntot_p
  neutron_data%total_orbits = ntot_n
  
  
  
  if(iam == 0)write(6,*)'neutron_total orbits ',neutron_data%total_orbits
  if(iam == 0)write(6,*)'proton_total orbits  ', proton_data%total_orbits
  if(iam == 0)write(6,*)'occupation p,n',occ_protons, occ_neutrons
  !if(iam == 0)write(6,*)'center of mass switch, beta',com_switch,com_beta, j2_beta
  
  !     Setup all possible orbit information
  all_orbit%total_orbits=neutron_data%total_orbits +    &
       proton_data%total_orbits
  CALL allocate_sp_array(neutron_data,neutron_data%total_orbits)
  CALL allocate_sp_array(proton_data,proton_data%total_orbits)
  CALL allocate_sp_array(all_orbit,all_orbit%total_orbits)
  
  !
  !
  !
  below_ef = occ_protons + occ_neutrons
  tot_orbs = all_orbit%total_orbits
  
  !     Read neutron single-particle data
  number_n = 0 
  number_p = 0
  nlmax = 0
  read(1,*) hbar_omega
  
  
  allocate(orbtype(all_orbit%total_orbits), indx(all_orbit%total_orbits))
  allocate( indx_inv( all_orbit%total_orbits) ) 
  
  do a = 1, all_orbit%total_orbits 
     
     read(1,*) i, n1, l1, j1, t1, re_e1, im_e1 
     
     !read(1,'(1X,I3,1X,I3,1X,I3,1X,I3,1X,I3,5X,2(E16.10,1X))') &
     ! i, n1, l1, j1, t1, re_e1, im_e1
          
     if ( a /= i ) then
        if ( iam == 0 ) write(6,*) 'wrong orbit number'
        stop
     end if
     
     !
     !     Neutrons are in the internal structure always even numbers
     !     Protons are in the internal structure always odd numbers
     !
     if ( t1 == 1 ) then
        number_n = number_n + 1
        if ( number_n <= occ_neutrons ) then
           neutron_data%orbit_status(number_n) = 'hole'
           all_orbit%orbit_status(i) = 'hole'
        else
           neutron_data%orbit_status(number_n) = 'particle'
           all_orbit%orbit_status(i) = 'particle'
        end if
     elseif ( t1 == -1 ) then
        number_p = number_p + 1
        if ( number_p <= occ_protons ) then
           proton_data%orbit_status(number_p) = 'hole'
           all_orbit%orbit_status(i) = 'hole'
        else
           proton_data%orbit_status(number_p) = 'particle'
           all_orbit%orbit_status(i) = 'particle'
        end if
     end if
     
     if(all_orbit%orbit_status(i).eq.'hole')then
        orbtype(i)=0
     else
        orbtype(i)=1
     end if
       
  end do
  
  !write(*,*) orbtype
  !call index(all_orbit%total_orbits,orbtype,indx_inv) ! indexing in ascending order
  
  !write(*,*) indx
  do i=1, all_orbit%total_orbits
     
     indx_inv(i) = i
     indx(i) = i 
     !indx(indx_inv(i))=i
  end do
  

  !do i = 1, all_orbit%total_orbits
  !   write(6,*) i, indx(i), indx_inv(indx(i))
  !end do
  
  !stop
  
  !write(*,*) indx

  deallocate(orbtype)
  number_n=0
  number_p=0
  rewind(1)
  read(1,*) hbar_omega

  do a = 1, all_orbit%total_orbits 
     
     !read(1,'(1X,I3,1X,I3,1X,I3,1X,I3,1X,I3,5X,2(E16.10,1X))') &
     !     i, n1, l1, j1, t1, re_e1, im_e1
     read(1,*) i, n1, l1, j1, t1, re_e1, im_e1
          
     if ( a /= i ) then
        if ( iam == 0 ) write(6,*) 'wrong orbit number'
        stop
     end if

     !ii = a 
     ii=indx(i)

     
     !
     !     Neutrons are in the internal structure always even numbers
     !     Protons are in the internal structure always odd numbers
     !
     if ( t1 == 1 ) then
        number_n = number_n + 1
        neutron_data%nn(number_n)=n1
        neutron_data%ll(number_n)=l1
        neutron_data%jj(number_n)=j1
        neutron_data%itzp(number_n)=t1
        neutron_data%e(number_n)= cmplx(re_e1,im_e1) 
        if ( number_n <= occ_neutrons ) then
           neutron_data%orbit_status(number_n) = 'hole'
           all_orbit%orbit_status(ii) = 'hole'
        else
           neutron_data%orbit_status(number_n) = 'particle'
           all_orbit%orbit_status(ii) = 'particle'
        end if
     elseif ( t1 == -1 ) then
        number_p = number_p + 1
        proton_data%nn(number_p)=n1 
        proton_data%ll(number_p)=l1 
        proton_data%jj(number_p)=j1
        proton_data%itzp(number_p)=t1 
        proton_data%e(number_p)= cmplx(re_e1,im_e1) 
        if ( number_p <= occ_protons ) then
           proton_data%orbit_status(number_p) = 'hole'
           all_orbit%orbit_status(ii) = 'hole'
        else
           proton_data%orbit_status(number_p) = 'particle'
           all_orbit%orbit_status(ii) = 'particle'
        end if
     end if
     
     all_orbit%nn(ii)=n1
     all_orbit%ll(ii)=l1
     all_orbit%jj(ii)=j1
     all_orbit%e(ii)= re_e1 
     all_orbit%itzp(ii)=t1
        
  end do


  ! set shell model modelspace
  open(33,file='model.dat')

  DO i=1, all_orbit%total_orbits
            all_orbit%model_space(i)='outside'
        !@!    if(all_orbit%itzp(i)==-1)cycle
        !@!    if(all_orbit%nn(i)*2+all_orbit%ll(i)==1 .and. all_orbit%orbit_status(i)=='particle' )then
        !@!    all_orbit%model_space(i)='inside'
        !@!  endif
        read(33,*)junk1,junk1,junk1,junk1,junk1,a
        if(a==1)all_orbit%model_space(i)='inside'
  ENDDO
  close(33)
  

  
  DO i=1, all_orbit%total_orbits
     if(iam == 0)WRITE(6,'(5i5,e16.8,2x,a,2x,a)')  &
          i, all_orbit%nn(i), all_orbit%ll(i), all_orbit%jj(i), &
          all_orbit%itzp(i), real(all_orbit%e(i)),   all_orbit%orbit_status(i),all_orbit%model_space(i)
     ! set initial occupations
     all_orbit%occ_numbers(i) = 1.d0
  ENDDO
  
  
  nn = 0
  pp = 0 
  if(iam == 0)write(6,*)'Occupied j-orbits' 
  efermi = 0.d0
  DO i=1,all_orbit%total_orbits
     if ( all_orbit%orbit_status(i) /= 'hole' ) cycle 
     if(iam == 0)WRITE(6,'(4i5,e16.8,2x,a)')  &
          all_orbit%nn(i), all_orbit%ll(i), all_orbit%jj(i), &
          all_orbit%itzp(i), real(all_orbit%e(i)),   all_orbit%orbit_status(i)
     
     if ( all_orbit%itzp(i) == -1 ) pp = pp + all_orbit%jj(i) + 1 
     if ( all_orbit%itzp(i) ==  1 ) nn = nn + all_orbit%jj(i) + 1  
     if ( all_orbit%nn(i)*2 + all_orbit%ll(i) > efermi ) efermi =  all_orbit%nn(i)*2 + all_orbit%ll(i) 
  ENDDO
  
  if ( iam == 0 ) write(6,*) 'Fermi level: 2n+l', efermi
  if ( iam == 0 ) write(6,*) 'Number of neutrons', nn
  if ( iam == 0 ) write(6,*) 'Number of protons', pp



  if(iam == 0)write(6,*)'Unoccupied j-orbits' 
  !DO i= occ_protons + occ_neutrons +1, all_orbit%total_orbits
  DO i= 1, all_orbit%total_orbits
     if ( all_orbit%orbit_status(i) /= 'particle' ) cycle 
     if(iam == 0)WRITE(6,'(4i5,e16.8,3x,2a)')  &
          all_orbit%nn(i), all_orbit%ll(i), all_orbit%jj(i), &
          all_orbit%itzp(i), real(all_orbit%e(i)),   all_orbit%orbit_status(i),all_orbit%model_space(i)
  ENDDO
  
  

  !
  ! set up kinetic energy matrix t2
  !
  !
  ! set up kinetic energy matrix t2
  !
  allocate( u_hf( all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( tkin( all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( hcom_1p( all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( theta_1b( all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( theta_1bdag( all_orbit%total_orbits, all_orbit%total_orbits) )
  
  theta_1bdag =0.d0
  theta_1b = 0.d0 
  re_e1 = 0.d0
  im_e1 = 0.d0
  re_hcom = 0.d0
  im_hcom = 0.d0
  u_hf = 0.d0 
  tkin = 0.d0
  hcom_1p = 0.d0
  do a = 1, all_orbit%total_orbits
     do c = 1, all_orbit%total_orbits
        
        !read(2,'(1X,I3,1X,I3,5X,4(E16.10,1X))') i,j,re_e1, im_e1, re_hcom, im_hcom
        read(2,*) i,j,re_e1, im_e1, re_hcom, im_hcom
        if ( a /= i .or. c /= j  ) then
           if ( iam == 0 ) write(6,*) 'wrong orbit number'
           stop
        end if
        
        i = indx(i)
        j = indx(j) 

        tza = all_orbit%itzp(i)
        tzc = all_orbit%itzp(j)
        
        !
        !  Make it a spherical tensor of rank 0 
        !       
        u_hf(i,j) = im_e1 
        
        tkin(i,j) = re_e1 
        theta_1b(i,j) = re_hcom*g_A/sqrt(6.*pi) !+ im_hcom 
        theta_1bdag(i,j) = im_hcom 
     end do
  end do
  
  do a = 1, all_orbit%total_orbits
     if ( all_orbit%itzp(a) /= 1 )cycle 
     do c = 1, all_orbit%total_orbits
        if ( all_orbit%itzp(c) /= -1 )cycle 
        theta_1bdag(a,c) = theta_1bdag(c,a)*iph( (all_orbit%jj(a)-all_orbit%jj(c))/2 )* & 
             sqrt(all_orbit%jj(c)+1.)/sqrt(all_orbit%jj(a)+1.)
     end do
  end do
  
  if ( iam == 0 ) write(6,*) 'Charge Z, and Mass number A', pp, mass_nucleus
  
  j2_max = -10 
  do a = 1, all_orbit%total_orbits
     do b = 1, all_orbit%total_orbits
        if ( (all_orbit%jj(a)+all_orbit%jj(b))/2 > j2_max ) j2_max = (all_orbit%jj(a)+all_orbit%jj(b))/2
     end do
  end do
  if ( iam == 0 ) write(6,*) 'Jab-max', j2_max 
  
  
END SUBROUTINE read_nuclear_sp_data


!
! read in the single-particle data
! and set up the single-particle orbitals; 
! allocate sp space.
!
SUBROUTINE read_nuclear_sp_data_imsrg
  USE PARALLEL
  USE single_particle_orbits
  USE constants
  use one_body_operators
  use reorder 
  use gmat_storage, only : scale
  use eom_module, only : Nmax, Nucleus,Interaction
  use cc_config_module

  IMPLICIT NONE
  CHARACTER(LEN=100) :: input1, input2
  INTEGER :: ntot_n, ntot_p, junk1, junk2, junk3, junk4, t_z, nn, tza, tzc, pp
  REAL*8 :: re_e1, im_e1, re_hcom, im_hcom
  INTEGER :: i,j, a,c, n1, l1, j1, t1, nlmax, number_p, number_n, number_orbs
  real*8 :: delta 
  integer:: ii, iph
  integer, allocatable, dimension(:):: orbtype
  
  !
  ! Number of occupied j-shells 
  !
  occ_protons = config%occ_protons
  occ_neutrons = config%occ_neutrons
  mass_nucleus = config%mass_nucleus
  com_beta = config%com_beta
  com_switch = config%com_beta
  j2_beta = config%com_switch

!!$  input1 = config%hf_input_file_orbits
!!$  input2 = config%hf_input_file_onebody
!!$
!!$  !     Read neutron and proton single-particle data
!!$  !read(5,*);read(5,*) input1
!!$  OPEN(1, FILE = input1)
!!$  !read(5,*);read(5,*) input2
!!$  OPEN(2, FILE = input2)
!!$  
!!$  if ( iam == 0 ) write(6,*) 'sp-energy file', input1
!!$  if ( iam == 0 ) write(6,*) 'kinetic-energy file', input2

  
  !
  ! Number of occupied j-shells 
  !
  !READ(5,*) ; READ(5,*) occ_protons, occ_neutrons
  !
  ! 
  !
  !READ(5,*) ; READ(5,*) com_beta,com_switch, j2_beta 
  !READ(5,*) ; READ(5,*) mass_nucleus
  hbar_omega = 14.d0
  
  call setup_sp_data_imsrg
  
  if(iam == 0)write(6,*)'neutron_total orbits ',neutron_data%total_orbits
  if(iam == 0)write(6,*)'proton_total orbits  ', proton_data%total_orbits
  if(iam == 0)write(6,*)'occupation p,n',occ_protons, occ_neutrons
  if(iam == 0)write(6,*)'center of mass switch, beta',com_switch,com_beta, j2_beta
  
  allocate( indx(all_orbit%total_orbits))
  allocate( indx_inv( all_orbit%total_orbits) ) 
  do i=1, all_orbit%total_orbits
     
     indx_inv(i) = i
     indx(i) = i 
  end do
  
  !
  !
  !
  below_ef = occ_protons + occ_neutrons
  tot_orbs = all_orbit%total_orbits
  
  number_n = 0 
  number_p = 0 
  do a = 1, all_orbit%total_orbits 
     t1 = all_orbit%itzp(a)

     !
     !     Neutrons are in the internal structure always even numbers
     !     Protons are in the internal structure always odd numbers
     !
     if ( t1 == 1 ) then
        number_n = number_n + 1
        if ( number_n <= occ_neutrons ) then
           neutron_data%orbit_status(number_n) = 'hole'
           all_orbit%orbit_status(a) = 'hole'
        else
           neutron_data%orbit_status(number_n) = 'particle'
           all_orbit%orbit_status(a) = 'particle'
        end if
     elseif ( t1 == -1 ) then
        number_p = number_p + 1
        if ( number_p <= occ_protons ) then
           proton_data%orbit_status(number_p) = 'hole'
           all_orbit%orbit_status(a) = 'hole'
        else
           proton_data%orbit_status(number_p) = 'particle'
           all_orbit%orbit_status(a) = 'particle'
        end if
     end if
          
  end do

  
  DO i=1, all_orbit%total_orbits
     if(iam == 0)WRITE(6,'(5i5,e16.8,2x,a)')  &
          i, all_orbit%nn(i), all_orbit%ll(i), all_orbit%jj(i), &
          all_orbit%itzp(i), real(all_orbit%e(i)),   all_orbit%orbit_status(i)
  ENDDO
    
  nn = 0
  pp = 0 
  if(iam == 0)write(6,*)'Occupied j-orbits' 
  DO i=1,all_orbit%total_orbits
     if ( all_orbit%orbit_status(i) /= 'hole' ) cycle 
     if(iam == 0)WRITE(6,'(4i5,e16.8,2x,a)')  &
          all_orbit%nn(i), all_orbit%ll(i), all_orbit%jj(i), &
          all_orbit%itzp(i), real(all_orbit%e(i)),   all_orbit%orbit_status(i)
     
     if ( all_orbit%itzp(i) == -1 ) pp = pp + all_orbit%jj(i) + 1 
     if ( all_orbit%itzp(i) ==  1 ) nn = nn + all_orbit%jj(i) + 1  

  ENDDO
  
  if ( iam == 0 ) write(6,*) 'Number of neutrons', nn
  if ( iam == 0 ) write(6,*) 'Number of protons', pp



  if(iam == 0)write(6,*)'Unoccupied j-orbits' 
  !DO i= occ_protons + occ_neutrons +1, all_orbit%total_orbits
  DO i= 1, all_orbit%total_orbits
     if ( all_orbit%orbit_status(i) /= 'particle' ) cycle 
     if(iam == 0)WRITE(6,'(4i5,e16.8,2x,a)')  &
          all_orbit%nn(i), all_orbit%ll(i), all_orbit%jj(i), &
          all_orbit%itzp(i), real(all_orbit%e(i)),   all_orbit%orbit_status(i)
  ENDDO
  
  
  
  !
  ! set up kinetic energy matrix t2
  !
  !
  ! set up kinetic energy matrix t2
  !
  allocate( u_hf( all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( tkin( all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( hcom_1p( all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( theta_1b( all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( theta_1bdag( all_orbit%total_orbits, all_orbit%total_orbits) )
  
  theta_1bdag =0.d0
  theta_1b = 0.d0 
  re_e1 = 0.d0
  im_e1 = 0.d0
  re_hcom = 0.d0
  im_hcom = 0.d0
  u_hf = 0.d0 
  tkin = 0.d0
  hcom_1p = 0.d0
  
END SUBROUTINE read_nuclear_sp_data_imsrg

!
!
! 
subroutine setup_fockmtx_and_interaction 

  use parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  use gmat_storage
  use one_body_operators
  use reorder 
  use HF_mod
  use gmat_storage
  use cc_config_module

  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmat_bra_configs 
  TYPE (configuration_descriptor) :: gmat_ket_configs 
  TYPE (configuration_descriptor) :: gmatrix_configs 
  INTEGER :: isospin_z, p_parity, ang_mom, pq_confs, ket_confs, bra_confs
  INTEGER :: i,j,ia,ib,ic,id, number_channels, ii, bra_conf, ket_conf, test_chan
  REAL*8 :: gmat,  angmom_fact, sum
  CHARACTER(LEN=100) :: g_matrix_file1, g_matrix_file2, channel_file, u_hf_file, fock_mtx_file 
  INTEGER :: j1, h, ja, jc, jh, jab_max, jcd_max, jab_min, jcd_min, jtmin, jtmax, jtot
  INTEGER :: aa, bb,cc,dd, hh, hh_number, hp_number, pp_number, k1, k2, k3, i1
  REAL*4 :: TM1, TM2
  INTEGER :: channel,channel1, channel2, channel3, channel4, channel5, channel6, & 
       bra, ket, nconfs, hf_switch, amount_bra, amount_ket,amount, iph
  REAL*8 :: gmat1, gmat2, gmat3,gmat4, phase_ab, phase_cd, dij, fnorm, memory 
  INTEGER :: total_work, work_per_proc, curr_work, curr_proc, curr_channel
  integer :: processor_work(num_procs), mtx_dim, temp, tot_work, work_pr_proc, local_ch_id
  integer :: row_dim, column_dim, ket_dim, bra_dim, my_channel_low, my_channel_high 
  INTEGER :: bra_min, bra_max, ket_min, ket_max
  integer, allocatable :: quantum_numbers_ch(:,:), quantum_numbers_ch3(:,:)
  integer(8), allocatable :: channel_map(:,:)
  LOGICAL :: has_ch_been_added
  REAL*8, allocatable :: work(:)
  integer :: vec, vec_sum, count, ierr, ichan, reclen, a, b
  integer(8):: istart, iend, jlong, numrec2
  logical :: precalc_fockmtx
  
  
  !
  ! read in j-coupled gmatrix elements
  !if(iam ==0)
  !read(5,*) ; read(5,*) g_matrix_file1
  !read(5,*) g_matrix_file2
  g_matrix_file1 = config%hf_input_file_twobody_pppp
  g_matrix_file2 = config%hf_input_file_twobody_rest
  channel_file = config%hf_input_file_channels    
  u_hf_file = config%hf_input_file_coefficients
  fock_mtx_file = config%hf_input_file_fock_matrix
  precalc_fockmtx = config%precalc_fockmtx
  !read(5,*);read(5,*) u_hf_file, fock_mtx_file 
  !read(5,*);read(5,*) precalc_fockmtx
  
  
  ! For complex files 
  inquire(IOLENGTH=reclen) channel, ia, ib, ic, id, gmat1, gmat2
  OPEN(UNIT=71, FILE=g_matrix_file1, ACCESS='direct', RECL=reclen, CONVERT='LITTLE_ENDIAN')
  OPEN(UNIT=72, FILE=g_matrix_file2, ACCESS='direct', RECL=reclen, CONVERT='LITTLE_ENDIAN')
  
  ! Let root read channels from file and broadcast 
  !read(5,*) ; read(5,*) channel_file

  open(UNIT=73, FILE=channel_file)
  i=0
11 read(73,*,end=12) ichan, istart, iend
  i=i+1
  goto 11
12 continue
  rewind(73)
  i=i-1 !last record is number of records for gmat file number 2
  !write(6,*) 'NUMBER OF CHANNELS READ FROM FILE: ', i
  allocate(channel_map(2,i))
  do j=1, i
     read(73,*) ichan, istart, iend
     channel_map(1,ichan)=istart
     channel_map(2,ichan)=iend
  end do
  read(73,*) i, numrec2, numrec2

  
  u_hf = 0.d0 
  !
  ! read and setup HF transformation matrix 
  !
  nlmax_cut = 1000
  !read(5,*);read(5,*) u_hf_file, fock_mtx_file 
  !read(5,*);read(5,*) precalc_fockmtx
  OPEN(UNIT=11, FILE=u_hf_file )
  OPEN(UNIT=12, FILE=fock_mtx_file)

  if ( iam == 0 ) write(6,*) u_hf_file, fock_mtx_file
  u_hf = 0.d0
  do i = 1, tot_orbs
     do j = 1, tot_orbs
        !read(11,2002) a,b, gmat
        read(11,*) a,b, gmat
        
        a = indx(a)
        b = indx(b) 
        
        u_hf(a,b) = gmat
     end do
  end do
  
2002 FORMAT(1X,I3,1X,I3,5X,(E16.10,1X))


  if(iam==0)write(6,*) g_matrix_file1
  if(iam==0)write(6,*) g_matrix_file2
  
  ! read in j-coupled gmatrix elements
  call setup_proc_mappings_pp
  
  number_channels = size( check_my_channel ) 
  allocate( quantum_numbers_ch3( number_channels, 3 ))
  
  !
  ! Count total number of channels
  !
  number_channels = 0
  !     loop over isospin projection
  DO isospin_z=-1,1 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60 

           
           !CALL  number_gmatrix_confs&
           !     (ang_mom,p_parity,isospin_z,1,gmat_bra_configs)
           !IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
          
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,3,gmat_ket_configs)
           IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
           
           number_channels = number_channels + 1 
           quantum_numbers_ch3(number_channels,1) = ang_mom
           quantum_numbers_ch3(number_channels,2) = p_parity
           quantum_numbers_ch3(number_channels,3) = isospin_z
           

        end DO
     end DO
  end DO
  
  
  number_channels = 0
  !     loop over isospin projection
  DO isospin_z=-1,1 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60 
           
           
                      
           !
           ! call number of bra_confs
           !
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,0,gmat_bra_configs)
           IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
           number_channels = number_channels + 1 
        end DO
     end DO
  end DO
  allocate( quantum_numbers_ch( number_channels, 3 ))
  ALLOCATE( lookup_ab_configs(1:3,number_channels))
  
  allocate( locate_gmatchannel(0:7,-1:1, 0:1, 0:60) ) 
  locate_gmatchannel = 0
  channels%number_confs = number_channels
  ALLOCATE( channels%config_J_ipar_tz(3*number_channels) )
  number_channels = 0
  !     loop over isospin projection
  DO isospin_z=-1,1 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60 

           
                      
           !
           ! call number of bra_confs
           !
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,0,gmat_bra_configs)
           IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
           nconfs = gmat_bra_configs%number_confs
           ALLOCATE(gmat_bra_configs%config_ab(nconfs + nconfs) )
           CALL setup_gmatrix_configurations &
                (ang_mom,p_parity,isospin_z,0,gmat_bra_configs)
           


           number_channels = number_channels + 1 
           k1 = number_channels*3 
           k2 = number_channels*3 - 1
           k3 = number_channels*3 - 2
           channels%config_J_ipar_tz(k1) = ang_mom
           channels%config_J_ipar_tz(k2) = p_parity
           channels%config_J_ipar_tz(k3) = isospin_z
           
           quantum_numbers_ch(number_channels,1) = ang_mom
           quantum_numbers_ch(number_channels,2) = p_parity
           quantum_numbers_ch(number_channels,3) = isospin_z
           locate_gmatchannel(0,isospin_z, p_parity, ang_mom) = number_channels
        end DO
     end DO
  end DO
  
  number_channels = 0
  !     loop over isospin projection
  DO isospin_z=-1,1 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60 

           
           !
           ! call number of bra_confs
           !
           !CALL number_gmatrix_confs&
           !     (ang_mom,p_parity,isospin_z,1,gmat_ket_configs)
           !IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
           
           
           !
           ! call number of bra_confs
           !
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,3,gmat_bra_configs)
           IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
           nconfs = gmat_bra_configs%number_confs
           ALLOCATE(gmat_bra_configs%config_ab(nconfs + nconfs) )
           CALL setup_gmatrix_configurations &
                (ang_mom,p_parity,isospin_z,3,gmat_bra_configs)
           
           number_channels = number_channels + 1 
           locate_gmatchannel(6,isospin_z, p_parity, ang_mom) = number_channels
           
           allocate( pq_conf( number_channels)%ival( tot_orbs, tot_orbs) )
           allocate( rs_conf( number_channels)%ival( tot_orbs, tot_orbs) )
           
           
           pq_conf(number_channels)%ival = 0
           rs_conf(number_channels)%ival = 0
           
           if ( check_my_channel(number_channels) == 0 ) cycle
           
           local_ch_id = number_channels
           !if ( iam ==0 ) write(6,'(a20,2x,i3)') 'Column index start:', mapping(iam+1,local_ch_id,2)
           !if ( iam ==0 ) write(6,'(a20,2x,i3)') 'Column index end:', mapping(iam+1,local_ch_id,3)
           !if ( iam ==0 ) write(6,'(a20,2x,i3)') 'Row index start:', mapping(iam+1,local_ch_id,4)
           !if ( iam ==0 ) write(6,'(a20,2x,i3)') 'Row index end:', mapping(iam+1,local_ch_id,5)
           
           ket_min = mapping(iam+1,local_ch_id,2)
           ket_max = mapping(iam+1,local_ch_id,3)
           bra_min = mapping(iam+1,local_ch_id,4)
           bra_max = mapping(iam+1,local_ch_id,5)
!!$
!!$
!!$
!!$           bra_min = 1; bra_max = nconfs
!!$           ket_min = 1; ket_max = nconfs

           if ( ket_max == bra_max .and. ket_min == bra_min ) fully_stored_channel( local_ch_id ) = 1 
           
           DO i=bra_min, bra_max
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              
              pq_conf( number_channels)%ival(ia,ib) = i 
           end DO
           
           DO i=ket_min, ket_max
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              rs_conf( number_channels)%ival(ia,ib) = i 
              
           end DO
           
           
           allocate( gmat_mtx_pppp%depot(number_channels)%val(bra_min:bra_max, ket_min:ket_max) )
           gmat_mtx_pppp%depot(number_channels)%val = 0.d0
           !write(6,*) number_channels, bra_min, bra_max, size(gmat_mtx_pppp%depot(number_channels)%val,1 ) 

        end DO
     end DO
  end DO

  
  
  call mpi_barrier(mpi_comm_world,ierror)
  
  do i = 1, number_channels 
     if ( check_my_channel(i) == 0 ) cycle
     if ( fully_stored_channel(i) == 0 ) then
        !write(6,*) 'iam', iam, 'channel', i 
     end if
  end do
  
  memory = 0.d0 
  do ii = 1, 3
     !
     ! Count number of channels 
     !
     do channel = 1, channels%number_confs 
        ang_mom   = channels%config_J_ipar_tz(channel*3)
        p_parity  = channels%config_J_ipar_tz(channel*3-1)
        isospin_z = channels%config_J_ipar_tz(channel*3-2)
        
        if ( ii == 1 ) then
           bra_conf = 1 
        elseif ( ii == 2 ) then
           bra_conf = 2 
        elseif ( ii == 3 ) then
           bra_conf = 3 
        end if
              
        !
        ! call number of bra_confs
        !
        CALL  number_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        
        nconfs = gmat_bra_configs%number_confs
        ALLOCATE( lookup_ab_configs(ii,channel)%ival(2,nconfs) )
        
        IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
        ALLOCATE(gmat_bra_configs%config_ab(2*nconfs) )
        CALL setup_gmatrix_configurations &
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        
        memory = memory + dble(ii*channel*2*nconfs)*4./1.e9 
        lookup_ab_configs(ii,channel)%ival = 0
        DO i=1,gmat_bra_configs%number_confs
           ia= gmat_bra_configs%config_ab(i*2-1)
           ib= gmat_bra_configs%config_ab(i*2)
           
           lookup_ab_configs(ii,channel)%ival(1,i) = ia
           lookup_ab_configs(ii,channel)%ival(2,i) = ib
        end DO
                
     end DO
  end do
  if ( iam == 0 ) write(6,*) 'Memory for lookup ab configs', memory, 'GByte'
  
  
  test_chan = 0
  do ii = 1, 5
     !
     ! Count number of channels 
     !
     number_channels = 0
     do i1 = 1, channels%number_confs 
        ang_mom   = channels%config_J_ipar_tz(i1*3)
        p_parity  = channels%config_J_ipar_tz(i1*3-1)
        isospin_z = channels%config_J_ipar_tz(i1*3-2)
                           
        if ( ii == 1 ) then
           bra_conf = 1 
           ket_conf = 1
        elseif ( ii == 2 ) then
           bra_conf = 1 
           ket_conf = 2 
        elseif ( ii == 3 ) then
           bra_conf = 1 
           ket_conf = 3 
        elseif ( ii == 4 ) then
           bra_conf = 2 
           ket_conf = 2 
        elseif ( ii == 5 ) then
           bra_conf = 2 
           ket_conf = 3 
        elseif ( ii == 6 ) then
           bra_conf = 3
           ket_conf = 3 
        end if
              
        !
        ! call number of bra_confs
        !
        CALL  number_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
        !
        ! call number of ket_confs
        !
        CALL  number_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
        IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
        
        
        number_channels = number_channels + 1 
        locate_gmatchannel(ii,isospin_z, p_parity, ang_mom) = number_channels
        
              
              
     end DO
     
     if ( ii == 1 ) allocate( gmat_mtx_hhhh%depot(1:number_channels) )
     if ( ii == 2 ) allocate( gmat_mtx_hhhp%depot(1:number_channels) )
     if ( ii == 3 ) allocate( gmat_mtx_hhpp%depot(1:number_channels) )
     if ( ii == 4 ) allocate( gmat_mtx_hphp%depot(1:number_channels) )
     if ( ii == 5 ) allocate( gmat_mtx_hppp%depot(1:number_channels) )
     

     if ( ii == 1 ) allocate( j_hhhh(1:number_channels) )
     if ( ii == 2 ) allocate( j_hhhp(1:number_channels) )
     if ( ii == 3 ) allocate( j_hhpp(1:number_channels) )
     if ( ii == 4 ) allocate( j_hphp(1:number_channels) )
     if ( ii == 5 ) allocate( j_hppp(1:number_channels) )
     
     
     
     if ( ii == 1 ) allocate( hh_hhhh(1:number_channels) )
     if ( ii == 2 ) allocate( hh_hhhp(1:number_channels),hp_hhhp(1:number_channels) )
     if ( ii == 3 ) allocate( hh_hhpp(1:number_channels),pp_hhpp(1:number_channels) )
     if ( ii == 4 ) allocate( hp_hphp(1:number_channels) )
     if ( ii == 5 ) allocate( hp_hppp(1:number_channels),pp_hppp(1:number_channels) )
               
  end do
  
  
  do i=1, size(hh_hhhh)
     allocate( hh_hhhh(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     hh_hhhh(i)%ival = 0
  end do
  do i=1, size(hh_hhhp)
     allocate( hp_hhhp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     allocate( hh_hhhp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     hh_hhhp(i)%ival = 0
     hp_hhhp(i)%ival = 0
  end do
  
  do i=1, size(hh_hhpp)
     allocate( hh_hhpp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     allocate( pp_hhpp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     hh_hhpp(i)%ival = 0
     pp_hhpp(i)%ival = 0
  end do
  do i=1, size(hp_hphp)
     allocate( hp_hphp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     hp_hphp(i)%ival = 0
  end do
  do i=1, size(hp_hppp)
     allocate( hp_hppp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     allocate( pp_hppp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     hp_hppp(i)%ival = 0
     pp_hppp(i)%ival = 0
  end do
!  do i=1, size(pp_pppp)
!     allocate( pp_pppp(i)%ival(below_ef+1:tot_orbs, below_ef+1:tot_orbs) )
!     pp_pppp(i)%ival = 0
!  end do
  
  call mytimer(tm1)

  
  if ( iam ==0 ) write(6,*) 'setting up interaction'
  do ii = 1, 5
     number_channels = 0
     
     do i1 = 1, channels%number_confs 
        ang_mom   = channels%config_J_ipar_tz(i1*3)
        p_parity  = channels%config_J_ipar_tz(i1*3-1)
        isospin_z = channels%config_J_ipar_tz(i1*3-2)
     
        if ( ii == 1 ) then
           bra_conf = 1 
           ket_conf = 1
        elseif ( ii == 2 ) then
           bra_conf = 1 
           ket_conf = 2 
        elseif ( ii == 3 ) then
           bra_conf = 1 
           ket_conf = 3 
        elseif ( ii == 4 ) then
           bra_conf = 2 
           ket_conf = 2 
        elseif ( ii == 5 ) then
           bra_conf = 2 
           ket_conf = 3 
        elseif ( ii == 6 ) then
           bra_conf = 3
           ket_conf = 3 
        end if
              
        !
        ! call number of bra_confs
        !
        CALL  number_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
        !
        ! call number of ket_confs
        !
        CALL  number_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
        IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
              
        bra_confs=gmat_bra_configs%number_confs
        ket_confs=gmat_ket_configs%number_confs
              
        ALLOCATE(gmat_bra_configs%config_ab(bra_confs+bra_confs) )
        ALLOCATE(gmat_ket_configs%config_ab(ket_confs+ket_confs) )
              
        CALL setup_gmatrix_configurations &
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        
        CALL setup_gmatrix_configurations &
             (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
        
        number_channels = number_channels + 1 
     

              
        
        if ( ii == 1 ) allocate( gmat_mtx_hhhh%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 2 ) allocate( gmat_mtx_hhhp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 3 ) allocate( gmat_mtx_hhpp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 4 ) allocate( gmat_mtx_hphp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        !if ( ii == 5 ) allocate( gmat_mtx_hppp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        !if ( ii == 6 ) allocate( gmat_mtx_pppp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 5 ) then
           amount = 0
           DO i=1,gmat_ket_configs%number_confs
              ia= gmat_ket_configs%config_ab(i*2-1)
              ib= gmat_ket_configs%config_ab(i*2)
                    
              if ( ia > ib ) cycle
              amount = amount + 1
                                  
           end DO
           allocate( gmat_mtx_hppp%depot(number_channels)%val(1:bra_confs, 1:amount) ) 
        end if
        
        if ( ii == 1 ) gmat_mtx_hhhh%depot(number_channels)%val = 0.d0 
        if ( ii == 2 ) gmat_mtx_hhhp%depot(number_channels)%val = 0.d0
        if ( ii == 3 ) gmat_mtx_hhpp%depot(number_channels)%val = 0.d0
        if ( ii == 4 ) gmat_mtx_hphp%depot(number_channels)%val = 0.d0
        if ( ii == 5 ) gmat_mtx_hppp%depot(number_channels)%val = 0.d0
        !if ( ii == 6 ) gmat_mtx_pppp%depot(number_channels)%val = 0.d0
        
        if ( ii == 1 ) j_hhhh(number_channels) = ang_mom 
        if ( ii == 2 ) j_hhhp(number_channels) = ang_mom
        if ( ii == 3 ) j_hhpp(number_channels) = ang_mom
        if ( ii == 4 ) j_hphp(number_channels) = ang_mom
        if ( ii == 5 ) j_hppp(number_channels) = ang_mom
        !if ( ii == 6 ) j_pppp(number_channels) = ang_mom
              
        
        if  ( ii == 1 ) then
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
                    
              hh_hhhh(number_channels)%ival(ia,ib) = i 
              
           end DO
           
        elseif  ( ii == 2 ) then
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              
              hh_hhhp(number_channels)%ival(ia,ib) = i 
                    
           end DO
                 
           DO i=1,gmat_ket_configs%number_confs
              ia= gmat_ket_configs%config_ab(i*2-1)
              ib= gmat_ket_configs%config_ab(i*2)
                    
              hp_hhhp(number_channels)%ival(ia,ib) = i 
                    
           end DO
        elseif( ii  == 3 ) then 
               
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              
              hh_hhpp(number_channels)%ival(ia,ib) = i 
                    
           end DO
                 
           DO i=1,gmat_ket_configs%number_confs
              ia= gmat_ket_configs%config_ab(i*2-1)
              ib= gmat_ket_configs%config_ab(i*2)
                    
              pp_hhpp(number_channels)%ival(ia,ib) = i 
                    
           end DO
                 
        elseif ( ii == 4 ) then
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
                    
              hp_hphp(number_channels)%ival(ia,ib) = i 
              
           end DO
              
        elseif( ii  == 5 ) then 
           
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              
              hp_hppp(number_channels)%ival(ia,ib) = i 
                    
           end DO
           
           amount = 0
           DO i=1,gmat_ket_configs%number_confs
              ia= gmat_ket_configs%config_ab(i*2-1)
              ib= gmat_ket_configs%config_ab(i*2)
                    
              if ( ia > ib ) cycle
              amount = amount + 1
              pp_hppp(number_channels)%ival(ia,ib) = amount 
                    
           end DO
                            
!        elseif( ii == 6 ) then
!           DO i=1,gmat_bra_configs%number_confs
!              ia= gmat_bra_configs%config_ab(i*2-1)
!              ib= gmat_bra_configs%config_ab(i*2)
!                    
!              pp_pppp(number_channels)%ival(ia,ib) = i 
!                    
!           end DO
        end if
              
              
                            
                    
                    
     end DO
  end DO


  count = 0
  do i = 1, size( gmat_mtx_hppp%depot(: ))
     bra_conf = size( gmat_mtx_hppp%depot(i)%val,1 )
     ket_conf = size( gmat_mtx_hppp%depot(i)%val,2 )
     count = count + bra_conf*ket_conf
  end do
  if ( iam == 0 ) write(6,*) 'Memory of hppp matrix', real(count)*8.d0/1.e9, 'Gb'
  
  !call mpi_barrier(mpi_comm_world,ierror)
  
  if(iam ==0) write(6,*) 'read in g-matrix'

  !
  ! Read pppp block
  ! 

  amount = 0
  
  do i = 1, size( channel_map(1,:) ) 
     
     channel  = mapping( iam+1, i, 1)
     if ( channel == 0 ) cycle
     istart = channel_map(1,channel)
     iend   = channel_map(2,channel)     
!  DO WHILE (.TRUE.)
     !read in            t_z,     pi,     J, a, b, c, d, G
!5    READ(71,*,end=6) isospin_z, p_parity, ang_mom, aa, bb, cc, dd, gmat1 !, gmat3 !, gmat3 , gmat4
!5    READ(71,*,end=6) channel, aa, bb, cc, dd, gmat1 !, gmat3 !, gmat3 , gmat4
     do jlong=istart, iend
        !write(*,*) i, channel, jlong
        READ(71,rec=jlong) channel, aa, bb, cc, dd, gmat1, gmat2  
        

        !write(6,('5I3,1x,g')) channel, aa, bb, cc, dd, gmat1 
        gmat = gmat1 ! - gmat3*com_switch/mass_nucleus
        ang_mom = ang_mom/2 
        
        ang_mom = quantum_numbers_ch3(channel,1) 
        p_parity = quantum_numbers_ch3(channel,2) 
        isospin_z = quantum_numbers_ch3(channel,3) 
                
        if ( isospin_z /= 0 ) then
           fnorm = dij(aa,bb)*dij(cc,dd)
        elseif ( isospin_z == 0 ) then
           fnorm =1.d0
        end if
        gmat = gmat * fnorm
        
        !
        ! Map orbitals to new orbitals
        !
        ia = indx(aa)
        ib = indx(bb)
        ic = indx(cc)
        id = indx(dd)
        
        if ( 2*all_orbit%nn(ia) + all_orbit%ll(ia) > nlmax_cut ) gmat = 0.d0 
        if ( 2*all_orbit%nn(ib) + all_orbit%ll(ib) > nlmax_cut ) gmat = 0.d0 
        if ( 2*all_orbit%nn(ic) + all_orbit%ll(ic) > nlmax_cut ) gmat = 0.d0 
        if ( 2*all_orbit%nn(id) + all_orbit%ll(id) > nlmax_cut ) gmat = 0.d0 
        
        !
        ! Define phases from interchange of particles
        !
        phase_ab = iph( (all_orbit%jj(ia) + all_orbit%jj(ib))/2 + ang_mom +1 )
        phase_cd = iph( (all_orbit%jj(ic) + all_orbit%jj(id))/2 + ang_mom +1 )
     
     
     !
     ! Find channel
     !
        channel6 = locate_gmatchannel(6,isospin_z, p_parity, ang_mom) 
     !
     ! check if channel6 is my channel6
     !
        if ( channel6 == 0 ) cycle
        if ( check_my_channel(channel6) == 0 ) cycle
     
        
     !
     ! 1 
     !
        bra= pq_conf(channel6)%ival(ia,ib)
        ket= rs_conf(channel6)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = gmat 
        end if
        bra= pq_conf(channel6)%ival(ic,id)
        ket= rs_conf(channel6)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = gmat 
        end if
     !
     ! (ab)
     !
        bra= pq_conf(channel6)%ival(ib,ia)
        ket= rs_conf(channel6)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = phase_ab*gmat 
        end if
        bra= pq_conf(channel6)%ival(ic,id)
        ket= rs_conf(channel6)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = phase_ab*gmat 
        end if
     !
     ! (cd)
     !
        bra= pq_conf(channel6)%ival(ia,ib)
        ket= rs_conf(channel6)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = phase_cd*gmat 
        end if
        bra= pq_conf(channel6)%ival(id,ic)
        ket= rs_conf(channel6)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = phase_cd*gmat 
        end if
     !
     ! (ab)(cd)
     !
        bra= pq_conf(channel6)%ival(ib,ia)
        ket= rs_conf(channel6)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        bra= pq_conf(channel6)%ival(id,ic)
        ket= rs_conf(channel6)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        
        amount = amount + 1
     end do
     
  ENDDO
6 CONTINUE
  


!  DO WHILE (.TRUE.)
  do jlong=1, numrec2
     !read in            t_z,     pi,     J, a, b, c, d, G
!7    READ(72,*,end=8) isospin_z, p_parity, ang_mom, aa, bb, cc, dd, gmat1 !, gmat3 !, gmat3 , gmat4
!7    READ(72,*,end=8) channel, aa, bb, cc, dd, gmat1 !, gmat3 !, gmat3 , gmat4
     !write(6,*) jlong 
     READ(72,rec=jlong) channel, aa, bb, cc, dd, gmat1, gmat2 
     
     gmat = gmat1 ! - gmat3*com_switch/mass_nucleus
     ang_mom = ang_mom/2 
     
     
     ang_mom = quantum_numbers_ch(channel,1) 
     p_parity = quantum_numbers_ch(channel,2) 
     isospin_z = quantum_numbers_ch(channel,3) 
     
     if ( isospin_z /= 0 ) then
        fnorm = dij(aa,bb)*dij(cc,dd)
     elseif ( isospin_z == 0 ) then
        fnorm =1.d0
     end if
     gmat = gmat * fnorm
     
     !
     ! Map orbitals to new orbitals
     !
     ia = indx(aa)
     ib = indx(bb)
     ic = indx(cc)
     id = indx(dd)

     if ( 2*all_orbit%nn(ia) + all_orbit%ll(ia) > nlmax_cut ) gmat = 0.d0 
     if ( 2*all_orbit%nn(ib) + all_orbit%ll(ib) > nlmax_cut ) gmat = 0.d0 
     if ( 2*all_orbit%nn(ic) + all_orbit%ll(ic) > nlmax_cut ) gmat = 0.d0 
     if ( 2*all_orbit%nn(id) + all_orbit%ll(id) > nlmax_cut ) gmat = 0.d0 
     
     !
     ! Define phases from interchange of particles
     !
     phase_ab = iph( (all_orbit%jj(ia) + all_orbit%jj(ib))/2 + ang_mom +1 )
     phase_cd = iph( (all_orbit%jj(ic) + all_orbit%jj(id))/2 + ang_mom +1 )
     
     
     !
     ! Find channel
     !
     channel1 = locate_gmatchannel(1,isospin_z, p_parity, ang_mom) 
     channel2 = locate_gmatchannel(2,isospin_z, p_parity, ang_mom) 
     channel3 = locate_gmatchannel(3,isospin_z, p_parity, ang_mom) 
     channel4 = locate_gmatchannel(4,isospin_z, p_parity, ang_mom) 
     channel5 = locate_gmatchannel(5,isospin_z, p_parity, ang_mom) 
     
     
     if ( channel1 /= 0 ) then
        
        !
        ! 1 
        !
        bra= hh_hhhh(channel1)%ival(ia,ib)
        ket= hh_hhhh(channel1)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = gmat 
        end if
        bra= hh_hhhh(channel1)%ival(ic,id)
        ket= hh_hhhh(channel1)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = gmat 
        end if
        !
        ! (ab)
        !
        bra= hh_hhhh(channel1)%ival(ib,ia)
        ket= hh_hhhh(channel1)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = phase_ab*gmat 
        end if
        bra= hh_hhhh(channel1)%ival(ic,id)
        ket= hh_hhhh(channel1)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = phase_ab*gmat 
        end if
        !
        ! (cd)
        !
        bra= hh_hhhh(channel1)%ival(ia,ib)
        ket= hh_hhhh(channel1)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = phase_cd*gmat 
        end if
        bra= hh_hhhh(channel1)%ival(id,ic)
        ket= hh_hhhh(channel1)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = phase_cd*gmat 
        end if
        !
        ! (ab)(cd)
        !
        bra= hh_hhhh(channel1)%ival(ib,ia)
        ket= hh_hhhh(channel1)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        bra= hh_hhhh(channel1)%ival(id,ic)
        ket= hh_hhhh(channel1)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
     end if
     
     if ( channel2 /= 0 ) then
        
        !
        ! 1 
        !
        bra= hh_hhhp(channel2)%ival(ia,ib)
        ket= hp_hhhp(channel2)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = gmat 
        end if
        bra= hh_hhhp(channel2)%ival(ic,id)
        ket= hp_hhhp(channel2)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = gmat 
        end if
        !
        ! (ab)
        !
        bra= hh_hhhp(channel2)%ival(ib,ia)
        ket= hp_hhhp(channel2)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = phase_ab*gmat 
        end if
        bra= hh_hhhp(channel2)%ival(ic,id)
        ket= hp_hhhp(channel2)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = phase_ab*gmat 
        end if
        !
        ! (cd)
        !
        bra= hh_hhhp(channel2)%ival(ia,ib)
        ket= hp_hhhp(channel2)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = phase_cd*gmat 
        end if
        bra= hh_hhhp(channel2)%ival(id,ic)
        ket= hp_hhhp(channel2)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = phase_cd*gmat 
        end if
        !
        ! (ab)(cd)
        !
        bra= hh_hhhp(channel2)%ival(ib,ia)
        ket= hp_hhhp(channel2)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        bra= hh_hhhp(channel2)%ival(id,ic)
        ket= hp_hhhp(channel2)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
     end if
     
     
     if ( channel3 /= 0 ) then
        
        !
        ! 1 
        !
        bra= hh_hhpp(channel3)%ival(ia,ib)
        ket= pp_hhpp(channel3)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = gmat 
        end if
        bra= hh_hhpp(channel3)%ival(ic,id)
        ket= pp_hhpp(channel3)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = gmat 
        end if
        !
        ! (ab)
        !
        bra= hh_hhpp(channel3)%ival(ib,ia)
        ket= pp_hhpp(channel3)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = phase_ab*gmat 
        end if
        bra= hh_hhpp(channel3)%ival(ic,id)
        ket= pp_hhpp(channel3)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = phase_ab*gmat 
        end if
        !
        ! (cd)
        !
        bra= hh_hhpp(channel3)%ival(ia,ib)
        ket= pp_hhpp(channel3)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = phase_cd*gmat 
        end if
        bra= hh_hhpp(channel3)%ival(id,ic)
        ket= pp_hhpp(channel3)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = phase_cd*gmat 
        end if
        !
        ! (ab)(cd)
        !
        bra= hh_hhpp(channel3)%ival(ib,ia)
        ket= pp_hhpp(channel3)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        bra= hh_hhpp(channel3)%ival(id,ic)
        ket= pp_hhpp(channel3)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
     end if
     
     if ( channel4 /= 0 ) then
        
        !
        ! 1 
        !
        bra= hp_hphp(channel4)%ival(ia,ib)
        ket= hp_hphp(channel4)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = gmat 
        end if
        bra= hp_hphp(channel4)%ival(ic,id)
        ket= hp_hphp(channel4)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = gmat 
        end if
        !
        ! (ab)
        !
        bra= hp_hphp(channel4)%ival(ib,ia)
        ket= hp_hphp(channel4)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = phase_ab*gmat 
        end if
        bra= hp_hphp(channel4)%ival(ic,id)
        ket= hp_hphp(channel4)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = phase_ab*gmat 
        end if
        !
        ! (cd)
        !
        bra= hp_hphp(channel4)%ival(ia,ib)
        ket= hp_hphp(channel4)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = phase_cd*gmat 
        end if
        bra= hp_hphp(channel4)%ival(id,ic)
        ket= hp_hphp(channel4)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = phase_cd*gmat 
        end if
        !
        ! (ab)(cd)
        !
        bra= hp_hphp(channel4)%ival(ib,ia)
        ket= hp_hphp(channel4)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        bra= hp_hphp(channel4)%ival(id,ic)
        ket= hp_hphp(channel4)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
     end if
     
     if ( channel5 /= 0 ) then
        
        !
        ! 1 
        !
        bra= hp_hppp(channel5)%ival(ia,ib)
        ket= pp_hppp(channel5)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = gmat 
        end if
        bra= hp_hppp(channel5)%ival(ic,id)
        ket= pp_hppp(channel5)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = gmat 
        end if
        !
        ! (ab)
        !
        bra= hp_hppp(channel5)%ival(ib,ia)
        ket= pp_hppp(channel5)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = phase_ab*gmat 
        end if
        bra= hp_hppp(channel5)%ival(ic,id)
        ket= pp_hppp(channel5)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = phase_ab*gmat 
        end if
        !
        ! (cd)
        !
        bra= hp_hppp(channel5)%ival(ia,ib)
        ket= pp_hppp(channel5)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = phase_cd*gmat 
        end if
        bra= hp_hppp(channel5)%ival(id,ic)
        ket= pp_hppp(channel5)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = phase_cd*gmat 
        end if
        !
        ! (ab)(cd)
        !
        bra= hp_hppp(channel5)%ival(ib,ia)
        ket= pp_hppp(channel5)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        bra= hp_hppp(channel5)%ival(id,ic)
        ket= pp_hppp(channel5)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
     end if
               
     amount = amount + 1
       
  ENDDO
8 CONTINUE

  
  
  deallocate( quantum_numbers_ch )
  deallocate( quantum_numbers_ch3 )
  
  !call mpi_barrier(mpi_comm_world,ierror)
  if ( iam ==0 )  write(6,*) 'Amount of G-matrix elements', amount
  
  !allocate( u_hf(tot_orbs, tot_orbs) )
  !u_hf = 0.d0
  !do i = 1, tot_orbs 
  !   u_hf(i,i) = 1.d0
  !end do
  
!  if ( iam == 0 ) write(6,*) 'HF-solver' 
!  READ(5,*);READ(5,*) hf_switch
!  if ( hf_switch == 1 ) call HF_solver 
  
  
  !if ( iam == 0 ) write(6,*) 'checking transformation' 
  !deallocate( rho )
  !call mpi_barrier(mpi_comm_world,ierror)
  
  !call HF_solver
  !return
  !stop
  
  
  
  call mytimer(tm2)
  if ( iam ==0 ) write(6,*) 'setting up interaction done!'
  



  if ( iam ==0 ) write(6,*) 'Time setting up interaction in pp-coupling', tm2-tm1, 'seconds'
  !
  ! calculate vacuum expectation value
  !
  e0 = 0.d0 
  do i = 1, all_orbit%total_orbits 
     if ( all_orbit%orbit_status(i) /= 'hole' ) cycle 
     j1 = all_orbit%jj(i)
     e0 = e0 + (j1+1.d0)*tkin(i,i) 
  end do
  !write(6,*) 'e0 kinetic only', e0
  
  number_channels = 0
  !     loop over isospin projection
  DO isospin_z=-1,1 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60 
           
           !
           ! call number of hh_confs 
           !
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,1,gmatrix_configs)
           IF (gmatrix_configs%number_confs <= 0 ) CYCLE
                                    
           number_channels = number_channels + 1 
           DO i=1,gmatrix_configs%number_confs
              
              gmat = gmat_mtx_hhhh%depot(number_channels)%val(i,i)
              e0 = e0 + 0.5d0* gmat*(2.d0*ang_mom + 1.d0) 
                            
           end DO
        end DO
     end do
  end do
  
  if ( iam == 0 ) write(6,*) 'Vacuum energy e0', e0
  
  ! 
  ! Setup Fock-matrix in j-scheme
  !
  allocate( fock_mtx(all_orbit%total_orbits ,all_orbit%total_orbits )) 

  fock_mtx = 0.d0
  !
  ! Add kinetic energy to Fock matrix
  !
  fock_mtx = fock_mtx + tkin 
       
  


  !
  ! Here we set of one-body Fock-matrix
  !
  do ia = 1, all_orbit%total_orbits
     do ic = 1, all_orbit%total_orbits
        do h = 1, all_orbit%total_orbits
           if ( all_orbit%orbit_status(h) /= 'hole') cycle
        
           aa = indx_inv(ia) 
           cc = indx_inv(ic) 
           hh = indx_inv(h) 
           
           ja = all_orbit%jj(ia)
           jc = all_orbit%jj(ic)
           jh = all_orbit%jj(h)
           
           jab_min = ABS(ja-jh); jab_max = (ja+jh)
           jcd_min = ABS(jc-jh); jcd_max = (jc+jh)
           jtmin = MAX( jcd_min, jab_min); jtmax = MIN(jab_max,jcd_max)
           IF ( jtmin > jtmax) CYCLE
           if ( all_orbit%jj(ia) /= all_orbit%jj(ic) ) cycle
           if ( all_orbit%ll(ia) /= all_orbit%ll(ic) ) cycle
           if ( all_orbit%itzp(ia) /= all_orbit%itzp(ic) ) cycle
           
           isospin_z =  (all_orbit%itzp(ia) + all_orbit%itzp(h))/2
           p_parity  = ( 1 - (-1)**( all_orbit%ll(ia) + all_orbit%ll(h) ))/2
           
           DO jtot = jtmin, jtmax, 2
              
              if ( ia > below_ef .and. ic > below_ef ) then
                 channel = locate_gmatchannel(4,isospin_z, p_parity, jtot/2 ) 
                 if ( channel == 0 ) cycle

                 bra = hp_hphp(channel)%ival(h,ia)
                 ket = hp_hphp(channel)%ival(h,ic)
                 if ( bra * ket == 0 ) cycle
                 gmat = gmat_mtx_hphp%depot(channel)%val(bra,ket) 
              elseif ( ia <= below_ef .and. ic <= below_ef ) then
                 channel = locate_gmatchannel(1,isospin_z, p_parity, jtot/2 ) 
                 if ( channel == 0 ) cycle

                 bra = hh_hhhh(channel)%ival(h,ia)
                 ket = hh_hhhh(channel)%ival(h,ic)
                 if ( bra * ket == 0 ) cycle
                 gmat = gmat_mtx_hhhh%depot(channel)%val(bra,ket) 
              elseif ( ia > below_ef .and. ic <= below_ef ) then
                 channel = locate_gmatchannel(2,isospin_z, p_parity, jtot/2 ) 
                 if ( channel == 0 ) cycle

                 bra = hh_hhhp(channel)%ival(h,ic)
                 ket = hp_hhhp(channel)%ival(h,ia)
                 if ( bra * ket == 0 ) cycle
                 gmat = gmat_mtx_hhhp%depot(channel)%val(bra,ket) 
              elseif ( ia <= below_ef .and. ic > below_ef ) then
                 channel = locate_gmatchannel(2,isospin_z, p_parity, jtot/2 ) 
                 if ( channel == 0 ) cycle

                 bra = hh_hhhp(channel)%ival(h,ia)
                 ket = hp_hhhp(channel)%ival(h,ic)
                 if ( bra * ket == 0 ) cycle
                 gmat = gmat_mtx_hhhp%depot(channel)%val(bra,ket) 
              
              end if
              
              angmom_fact = (jtot+1.d0)/( ja + 1.d0) 
              !CALL jcoupled_gmatrix_element_new (aa,hh,cc,hh,jtot/2,gmat)
              fock_mtx(ia,ic) = fock_mtx(ia,ic) + angmom_fact * gmat
           end DO
        end do
     end do
  end do
  
  sum = 0.d0
  do i = 1, below_ef
     sum = sum + 0.5d0*( tkin(i,i) + fock_mtx(i,i) ) * ( all_orbit%jj(i) + 1.d0) 
  end do
  
  
!!$  do ia = 1, all_orbit%total_orbits
!!$     write(6,*) 'fock matrix', all_orbit%ll(ia),all_orbit%jj(ia), all_orbit%itzp(ia),fock_mtx(ia,ia)
!!$  end do
!!$  stop
  
  !precalc_fockmtx = .TRUE.
  if ( precalc_fockmtx ) then 
     fock_mtx = 0.d0 
     read(12,*) e0 
     do ia = 1, all_orbit%total_orbits
        do ic = 1, all_orbit%total_orbits
           
           read(12,*) aa,cc, gmat
           
           if ( aa /= ia ) then
              write(6,*) 'error', aa, ia 
              stop
           end if
           if ( cc /= ic ) then
              write(6,*) 'error', cc, ic 
              stop
           end if

           if ( 2*all_orbit%nn(ia) + all_orbit%ll(ia) > nlmax_cut ) gmat = 0.d0 
           if ( 2*all_orbit%nn(ic) + all_orbit%ll(ic) > nlmax_cut ) gmat = 0.d0 
           fock_mtx(ia,ic) = gmat 
           
        ENDDO
     ENDDO
     if ( iam == 0 ) write(6,*) 'Vacuum energy e0', e0
  end if
  
  if ( iam ==0 ) write(6,*) sum
  !deallocate ( work ) 

end subroutine setup_fockmtx_and_interaction

!
!
! 
subroutine setup_fockmtx_and_interaction_imsrg 
  
  use parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  use gmat_storage
  use one_body_operators
  use reorder 
  use HF_mod
  use gmat_storage
  use cc_config_module
  
  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmat_bra_configs 
  TYPE (configuration_descriptor) :: gmat_ket_configs 
  TYPE (configuration_descriptor) :: gmatrix_configs 
  INTEGER :: isospin_z, p_parity, ang_mom, pq_confs, ket_confs, bra_confs
  INTEGER :: i,j,ia,ib,ic,id, number_channels, ii, bra_conf, ket_conf, test_chan
  REAL*8 :: gmat,  angmom_fact, sum
  CHARACTER(LEN=100) :: g_matrix_file1, g_matrix_file2, channel_file, u_hf_file, fock_mtx_file 
  INTEGER :: j1, h, ja, jc, jh, jab_max, jcd_max, jab_min, jcd_min, jtmin, jtmax, jtot
  INTEGER :: aa, bb,cc,dd, hh, hh_number, hp_number, pp_number, k1, k2, k3, i1
  REAL*4 :: TM1, TM2
  INTEGER :: channel,channel1, channel2, channel3, channel4, channel5, channel6, & 
       bra, ket, nconfs, hf_switch, amount_bra, amount_ket,amount, iph
  REAL*8 :: gmat1, gmat2, gmat3,gmat4, phase_ab, phase_cd, dij, fnorm, memory 
  INTEGER :: total_work, work_per_proc, curr_work, curr_proc, curr_channel
  integer :: processor_work(num_procs), mtx_dim, temp, tot_work, work_pr_proc, local_ch_id
  integer :: row_dim, column_dim, ket_dim, bra_dim, my_channel_low, my_channel_high 
  INTEGER :: bra_min, bra_max, ket_min, ket_max
  integer, allocatable :: quantum_numbers_ch(:,:), quantum_numbers_ch3(:,:)
  integer(8), allocatable :: channel_map(:,:)
  LOGICAL :: has_ch_been_added
  REAL*8, allocatable :: work(:)
  integer :: vec, vec_sum, count, ierr, ichan, reclen, a, b
  integer(8):: istart, iend, jlong, numrec2
  logical :: precalc_fockmtx

  g_matrix_file1 = config%hf_input_file_twobody_pppp
  g_matrix_file2 = config%hf_input_file_twobody_rest
  channel_file = config%hf_input_file_channels    
  u_hf_file = config%hf_input_file_coefficients
  fock_mtx_file = config%hf_input_file_fock_matrix
  precalc_fockmtx = config%precalc_fockmtx
  
  
  !
  ! read in j-coupled gmatrix elements
  !if(iam ==0)
  !read(5,*) ; read(5,*) g_matrix_file1
  !read(5,*) g_matrix_file2
  ! Let root read channels from file and broadcast 
  !read(5,*) ; read(5,*) channel_file
  
  
  OPEN(UNIT=71, FILE=g_matrix_file1, STATUS='OLD',form='formatted')
  u_hf = 0.d0 
  !
  ! read and setup HF transformation matrix 
  !
  nlmax_cut = 1000
  !read(5,*);read(5,*) u_hf_file, fock_mtx_file 
  !read(5,*);read(5,*) precalc_fockmtx
  
  
  OPEN(UNIT=11, FILE=u_hf_file )
  OPEN(UNIT=12, FILE=fock_mtx_file)
  
  if(iam==0)write(6,*) g_matrix_file1
  if(iam==0)write(6,*) g_matrix_file2
  
  ! read in j-coupled gmatrix elements
  call setup_proc_mappings_pp
  
  number_channels = size( check_my_channel ) 
  allocate( quantum_numbers_ch3( number_channels, 3 ))
  
  !
  ! Count total number of channels
  !
  number_channels = 0
  !     loop over isospin projection
  DO isospin_z=-1,1 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60 

           
           !CALL  number_gmatrix_confs&
           !     (ang_mom,p_parity,isospin_z,1,gmat_bra_configs)
           !IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
          
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,3,gmat_ket_configs)
           IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
           
           number_channels = number_channels + 1 
           quantum_numbers_ch3(number_channels,1) = ang_mom
           quantum_numbers_ch3(number_channels,2) = p_parity
           quantum_numbers_ch3(number_channels,3) = isospin_z
           

        end DO
     end DO
  end DO
  
  
  number_channels = 0
  !     loop over isospin projection
  DO isospin_z=-1,1 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60 
           
           
                      
           !
           ! call number of bra_confs
           !
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,0,gmat_bra_configs)
           IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
           number_channels = number_channels + 1 
        end DO
     end DO
  end DO
  allocate( quantum_numbers_ch( number_channels, 3 ))
  ALLOCATE( lookup_ab_configs(1:3,number_channels))
  
  allocate( locate_gmatchannel(0:6,-1:1, 0:1, 0:60) ) 
  locate_gmatchannel = 0
  channels%number_confs = number_channels
  ALLOCATE( channels%config_J_ipar_tz(3*number_channels) )
  number_channels = 0
  !     loop over isospin projection
  DO isospin_z=-1,1 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60 

           
                      
           !
           ! call number of bra_confs
           !
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,0,gmat_bra_configs)
           IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
           nconfs = gmat_bra_configs%number_confs
           ALLOCATE(gmat_bra_configs%config_ab(nconfs + nconfs) )
           CALL setup_gmatrix_configurations &
                (ang_mom,p_parity,isospin_z,0,gmat_bra_configs)
           


           number_channels = number_channels + 1 
           k1 = number_channels*3 
           k2 = number_channels*3 - 1
           k3 = number_channels*3 - 2
           channels%config_J_ipar_tz(k1) = ang_mom
           channels%config_J_ipar_tz(k2) = p_parity
           channels%config_J_ipar_tz(k3) = isospin_z
           
           quantum_numbers_ch(number_channels,1) = ang_mom
           quantum_numbers_ch(number_channels,2) = p_parity
           quantum_numbers_ch(number_channels,3) = isospin_z
           locate_gmatchannel(0,isospin_z, p_parity, ang_mom) = number_channels
        end DO
     end DO
  end DO
  
  number_channels = 0
  !     loop over isospin projection
  DO isospin_z=-1,1 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60 

           
           !
           ! call number of bra_confs
           !
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,3,gmat_bra_configs)
           IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
           nconfs = gmat_bra_configs%number_confs
           ALLOCATE(gmat_bra_configs%config_ab(nconfs + nconfs) )
           CALL setup_gmatrix_configurations &
                (ang_mom,p_parity,isospin_z,3,gmat_bra_configs)
           
           number_channels = number_channels + 1 
           locate_gmatchannel(6,isospin_z, p_parity, ang_mom) = number_channels
           
           allocate( pq_conf( number_channels)%ival( tot_orbs, tot_orbs) )
           allocate( rs_conf( number_channels)%ival( tot_orbs, tot_orbs) )
           
           
           pq_conf(number_channels)%ival = 0
           rs_conf(number_channels)%ival = 0
           
           if ( check_my_channel(number_channels) == 0 ) cycle
           
           local_ch_id = number_channels
           
           ket_min = mapping(iam+1,local_ch_id,2)
           ket_max = mapping(iam+1,local_ch_id,3)
           bra_min = mapping(iam+1,local_ch_id,4)
           bra_max = mapping(iam+1,local_ch_id,5)

           if ( ket_max == bra_max .and. ket_min == bra_min ) fully_stored_channel( local_ch_id ) = 1 
           
           DO i=bra_min, bra_max
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              
              pq_conf( number_channels)%ival(ia,ib) = i 
           end DO
           
           DO i=ket_min, ket_max
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              rs_conf( number_channels)%ival(ia,ib) = i 
              
           end DO
           
           
           allocate( gmat_mtx_pppp%depot(number_channels)%val(bra_min:bra_max, ket_min:ket_max) )
           gmat_mtx_pppp%depot(number_channels)%val = 0.d0
           
           
        end DO
     end DO
  end DO

  
  
  
  memory = 0.d0 
  do ii = 1, 3
     !
     ! Count number of channels 
     !
     do channel = 1, channels%number_confs 
        ang_mom   = channels%config_J_ipar_tz(channel*3)
        p_parity  = channels%config_J_ipar_tz(channel*3-1)
        isospin_z = channels%config_J_ipar_tz(channel*3-2)
        
        if ( ii == 1 ) then
           bra_conf = 1 
        elseif ( ii == 2 ) then
           bra_conf = 2 
        elseif ( ii == 3 ) then
           bra_conf = 3 
        end if
              
        !
        ! call number of bra_confs
        !
        CALL  number_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        
        nconfs = gmat_bra_configs%number_confs
        ALLOCATE( lookup_ab_configs(ii,channel)%ival(2,nconfs) )
        
        IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
        ALLOCATE(gmat_bra_configs%config_ab(2*nconfs) )
        CALL setup_gmatrix_configurations &
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        
        memory = memory + dble(ii*channel*2*nconfs)*4./1.e9 
        lookup_ab_configs(ii,channel)%ival = 0
        DO i=1,gmat_bra_configs%number_confs
           ia= gmat_bra_configs%config_ab(i*2-1)
           ib= gmat_bra_configs%config_ab(i*2)
           
           lookup_ab_configs(ii,channel)%ival(1,i) = ia
           lookup_ab_configs(ii,channel)%ival(2,i) = ib
        end DO
                
     end DO
  end do
  if ( iam == 0 ) write(6,*) 'Memory for lookup ab configs', memory, 'GByte'
  
  
  test_chan = 0
  do ii = 1, 5
     !
     ! Count number of channels 
     !
     number_channels = 0
     do i1 = 1, channels%number_confs 
        ang_mom   = channels%config_J_ipar_tz(i1*3)
        p_parity  = channels%config_J_ipar_tz(i1*3-1)
        isospin_z = channels%config_J_ipar_tz(i1*3-2)
                           
        if ( ii == 1 ) then
           bra_conf = 1 
           ket_conf = 1
        elseif ( ii == 2 ) then
           bra_conf = 1 
           ket_conf = 2 
        elseif ( ii == 3 ) then
           bra_conf = 1 
           ket_conf = 3 
        elseif ( ii == 4 ) then
           bra_conf = 2 
           ket_conf = 2 
        elseif ( ii == 5 ) then
           bra_conf = 2 
           ket_conf = 3 
        elseif ( ii == 6 ) then
           bra_conf = 3
           ket_conf = 3 
        end if
              
        !
        ! call number of bra_confs
        !
        CALL  number_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
        !
        ! call number of ket_confs
        !
        CALL  number_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
        IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
        
        
        number_channels = number_channels + 1 
        locate_gmatchannel(ii,isospin_z, p_parity, ang_mom) = number_channels
        
              
              
     end DO
     
     if ( ii == 1 ) allocate( gmat_mtx_hhhh%depot(1:number_channels) )
     if ( ii == 2 ) allocate( gmat_mtx_hhhp%depot(1:number_channels) )
     if ( ii == 3 ) allocate( gmat_mtx_hhpp%depot(1:number_channels) )
     if ( ii == 4 ) allocate( gmat_mtx_hphp%depot(1:number_channels) )
     if ( ii == 5 ) allocate( gmat_mtx_hppp%depot(1:number_channels) )
     

     if ( ii == 1 ) allocate( j_hhhh(1:number_channels) )
     if ( ii == 2 ) allocate( j_hhhp(1:number_channels) )
     if ( ii == 3 ) allocate( j_hhpp(1:number_channels) )
     if ( ii == 4 ) allocate( j_hphp(1:number_channels) )
     if ( ii == 5 ) allocate( j_hppp(1:number_channels) )
     
     
     
     if ( ii == 1 ) allocate( hh_hhhh(1:number_channels) )
     if ( ii == 2 ) allocate( hh_hhhp(1:number_channels),hp_hhhp(1:number_channels) )
     if ( ii == 3 ) allocate( hh_hhpp(1:number_channels),pp_hhpp(1:number_channels) )
     if ( ii == 4 ) allocate( hp_hphp(1:number_channels) )
     if ( ii == 5 ) allocate( hp_hppp(1:number_channels),pp_hppp(1:number_channels) )
               
  end do
  
  
  do i=1, size(hh_hhhh)
     allocate( hh_hhhh(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     hh_hhhh(i)%ival = 0
  end do
  do i=1, size(hh_hhhp)
     allocate( hp_hhhp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     allocate( hh_hhhp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     hh_hhhp(i)%ival = 0
     hp_hhhp(i)%ival = 0
  end do
  
  do i=1, size(hh_hhpp)
     allocate( hh_hhpp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     allocate( pp_hhpp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     hh_hhpp(i)%ival = 0
     pp_hhpp(i)%ival = 0
  end do
  do i=1, size(hp_hphp)
     allocate( hp_hphp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     hp_hphp(i)%ival = 0
  end do
  do i=1, size(hp_hppp)
     allocate( hp_hppp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     allocate( pp_hppp(i)%ival(1:tot_orbs, 1:tot_orbs) ) 
     hp_hppp(i)%ival = 0
     pp_hppp(i)%ival = 0
  end do
!  do i=1, size(pp_pppp)
!     allocate( pp_pppp(i)%ival(below_ef+1:tot_orbs, below_ef+1:tot_orbs) )
!     pp_pppp(i)%ival = 0
!  end do
  
  call mytimer(tm1)

  
  if ( iam ==0 ) write(6,*) 'setting up interaction'
  do ii = 1, 5
     number_channels = 0
     
     do i1 = 1, channels%number_confs 
        ang_mom   = channels%config_J_ipar_tz(i1*3)
        p_parity  = channels%config_J_ipar_tz(i1*3-1)
        isospin_z = channels%config_J_ipar_tz(i1*3-2)
     
        if ( ii == 1 ) then
           bra_conf = 1 
           ket_conf = 1
        elseif ( ii == 2 ) then
           bra_conf = 1 
           ket_conf = 2 
        elseif ( ii == 3 ) then
           bra_conf = 1 
           ket_conf = 3 
        elseif ( ii == 4 ) then
           bra_conf = 2 
           ket_conf = 2 
        elseif ( ii == 5 ) then
           bra_conf = 2 
           ket_conf = 3 
        elseif ( ii == 6 ) then
           bra_conf = 3
           ket_conf = 3 
        end if
              
        !
        ! call number of bra_confs
        !
        CALL  number_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
        !
        ! call number of ket_confs
        !
        CALL  number_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
        IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
              
        bra_confs=gmat_bra_configs%number_confs
        ket_confs=gmat_ket_configs%number_confs
              
        ALLOCATE(gmat_bra_configs%config_ab(bra_confs+bra_confs) )
        ALLOCATE(gmat_ket_configs%config_ab(ket_confs+ket_confs) )
              
        CALL setup_gmatrix_configurations &
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        
        CALL setup_gmatrix_configurations &
             (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
        
        number_channels = number_channels + 1 
     

              
        
        if ( ii == 1 ) allocate( gmat_mtx_hhhh%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 2 ) allocate( gmat_mtx_hhhp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 3 ) allocate( gmat_mtx_hhpp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 4 ) allocate( gmat_mtx_hphp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        !if ( ii == 5 ) allocate( gmat_mtx_hppp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        !if ( ii == 6 ) allocate( gmat_mtx_pppp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 5 ) then
           amount = 0
           DO i=1,gmat_ket_configs%number_confs
              ia= gmat_ket_configs%config_ab(i*2-1)
              ib= gmat_ket_configs%config_ab(i*2)
                    
              if ( ia > ib ) cycle
              amount = amount + 1
                                  
           end DO
           allocate( gmat_mtx_hppp%depot(number_channels)%val(1:bra_confs, 1:amount) ) 
        end if
        
        if ( ii == 1 ) gmat_mtx_hhhh%depot(number_channels)%val = 0.d0 
        if ( ii == 2 ) gmat_mtx_hhhp%depot(number_channels)%val = 0.d0
        if ( ii == 3 ) gmat_mtx_hhpp%depot(number_channels)%val = 0.d0
        if ( ii == 4 ) gmat_mtx_hphp%depot(number_channels)%val = 0.d0
        if ( ii == 5 ) gmat_mtx_hppp%depot(number_channels)%val = 0.d0
        !if ( ii == 6 ) gmat_mtx_pppp%depot(number_channels)%val = 0.d0
        
        if ( ii == 1 ) j_hhhh(number_channels) = ang_mom 
        if ( ii == 2 ) j_hhhp(number_channels) = ang_mom
        if ( ii == 3 ) j_hhpp(number_channels) = ang_mom
        if ( ii == 4 ) j_hphp(number_channels) = ang_mom
        if ( ii == 5 ) j_hppp(number_channels) = ang_mom
        !if ( ii == 6 ) j_pppp(number_channels) = ang_mom
              
        
        if  ( ii == 1 ) then
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
                    
              hh_hhhh(number_channels)%ival(ia,ib) = i 
              
           end DO
           
        elseif  ( ii == 2 ) then
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              
              hh_hhhp(number_channels)%ival(ia,ib) = i 
                    
           end DO
                 
           DO i=1,gmat_ket_configs%number_confs
              ia= gmat_ket_configs%config_ab(i*2-1)
              ib= gmat_ket_configs%config_ab(i*2)
                    
              hp_hhhp(number_channels)%ival(ia,ib) = i 
                    
           end DO
        elseif( ii  == 3 ) then 
               
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              
              hh_hhpp(number_channels)%ival(ia,ib) = i 
                    
           end DO
                 
           DO i=1,gmat_ket_configs%number_confs
              ia= gmat_ket_configs%config_ab(i*2-1)
              ib= gmat_ket_configs%config_ab(i*2)
                    
              pp_hhpp(number_channels)%ival(ia,ib) = i 
                    
           end DO
                 
        elseif ( ii == 4 ) then
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
                    
              hp_hphp(number_channels)%ival(ia,ib) = i 
              
           end DO
              
        elseif( ii  == 5 ) then 
           
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              
              hp_hppp(number_channels)%ival(ia,ib) = i 
                    
           end DO
           
           amount = 0
           DO i=1,gmat_ket_configs%number_confs
              ia= gmat_ket_configs%config_ab(i*2-1)
              ib= gmat_ket_configs%config_ab(i*2)
                    
              if ( ia > ib ) cycle
              amount = amount + 1
              pp_hppp(number_channels)%ival(ia,ib) = amount 
                    
           end DO
                            
!        elseif( ii == 6 ) then
!           DO i=1,gmat_bra_configs%number_confs
!              ia= gmat_bra_configs%config_ab(i*2-1)
!              ib= gmat_bra_configs%config_ab(i*2)
!                    
!              pp_pppp(number_channels)%ival(ia,ib) = i 
!                    
!           end DO
        end if
              
              
                            
                    
                    
     end DO
  end DO


  count = 0
  do i = 1, size( gmat_mtx_hppp%depot(: ))
     bra_conf = size( gmat_mtx_hppp%depot(i)%val,1 )
     ket_conf = size( gmat_mtx_hppp%depot(i)%val,2 )
     count = count + bra_conf*ket_conf
  end do
  if ( iam == 0 ) write(6,*) 'Memory of hppp matrix', real(count)*8.d0/1.e9, 'Gb'
  
  !call mpi_barrier(mpi_comm_world,ierror)
  
  if(iam ==0) write(6,*) 'read in g-matrix'
  
  !
  ! Read pppp block
  ! 
  read(71,*) 
  amount = 0
  DO WHILE (.TRUE.)
     !read in            t_z,     pi,     J, a, b, c, d, G
5    READ(71,*,end=6) isospin_z, p_parity, ang_mom, aa, bb, cc, dd, gmat1 
     
     aa = aa + 1
     bb = bb + 1
     cc = cc + 1
     dd = dd + 1
     
     gmat = gmat1 
     ang_mom = ang_mom/2 
     
     if ( isospin_z /= 0 ) then
        fnorm = dij(aa,bb)*dij(cc,dd)
     elseif ( isospin_z == 0 ) then
        fnorm =1.d0
     end if
     gmat = gmat * fnorm
     
     !
     ! Map orbitals to new orbitals
     !
     ia = aa
     ib = bb
     ic = cc
     id = dd
     
     !
     ! Define phases from interchange of particles
     !
     phase_ab = iph( (all_orbit%jj(ia) + all_orbit%jj(ib))/2 + ang_mom +1 )
     phase_cd = iph( (all_orbit%jj(ic) + all_orbit%jj(id))/2 + ang_mom +1 )
     
     
     !
     ! Find channel
     !
     channel1 = locate_gmatchannel(1,isospin_z, p_parity, ang_mom) 
     channel2 = locate_gmatchannel(2,isospin_z, p_parity, ang_mom) 
     channel3 = locate_gmatchannel(3,isospin_z, p_parity, ang_mom) 
     channel4 = locate_gmatchannel(4,isospin_z, p_parity, ang_mom) 
     channel5 = locate_gmatchannel(5,isospin_z, p_parity, ang_mom) 
     channel6 = locate_gmatchannel(6,isospin_z, p_parity, ang_mom) 
     
     if ( channel1 /= 0 ) then
        
        !
        ! 1 
        !
        bra= hh_hhhh(channel1)%ival(ia,ib)
        ket= hh_hhhh(channel1)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = gmat 
        end if
        bra= hh_hhhh(channel1)%ival(ic,id)
        ket= hh_hhhh(channel1)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = gmat 
        end if
        !
        ! (ab)
        !
        bra= hh_hhhh(channel1)%ival(ib,ia)
        ket= hh_hhhh(channel1)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = phase_ab*gmat 
        end if
        bra= hh_hhhh(channel1)%ival(ic,id)
        ket= hh_hhhh(channel1)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = phase_ab*gmat 
        end if
        !
        ! (cd)
        !
        bra= hh_hhhh(channel1)%ival(ia,ib)
        ket= hh_hhhh(channel1)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = phase_cd*gmat 
        end if
        bra= hh_hhhh(channel1)%ival(id,ic)
        ket= hh_hhhh(channel1)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = phase_cd*gmat 
        end if
        !
        ! (ab)(cd)
        !
        bra= hh_hhhh(channel1)%ival(ib,ia)
        ket= hh_hhhh(channel1)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        bra= hh_hhhh(channel1)%ival(id,ic)
        ket= hh_hhhh(channel1)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhh%depot(channel1)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
     end if
     
     if ( channel2 /= 0 ) then
        
        !
        ! 1 
        !
        bra= hh_hhhp(channel2)%ival(ia,ib)
        ket= hp_hhhp(channel2)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = gmat 
        end if
        bra= hh_hhhp(channel2)%ival(ic,id)
        ket= hp_hhhp(channel2)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = gmat 
        end if
        !
        ! (ab)
        !
        bra= hh_hhhp(channel2)%ival(ib,ia)
        ket= hp_hhhp(channel2)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = phase_ab*gmat 
        end if
        bra= hh_hhhp(channel2)%ival(ic,id)
        ket= hp_hhhp(channel2)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = phase_ab*gmat 
        end if
        !
        ! (cd)
        !
        bra= hh_hhhp(channel2)%ival(ia,ib)
        ket= hp_hhhp(channel2)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = phase_cd*gmat 
        end if
        bra= hh_hhhp(channel2)%ival(id,ic)
        ket= hp_hhhp(channel2)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = phase_cd*gmat 
        end if
        !
        ! (ab)(cd)
        !
        bra= hh_hhhp(channel2)%ival(ib,ia)
        ket= hp_hhhp(channel2)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        bra= hh_hhhp(channel2)%ival(id,ic)
        ket= hp_hhhp(channel2)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhhp%depot(channel2)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
     end if
     
     
     if ( channel3 /= 0 ) then
        
        !
        ! 1 
        !
        bra= hh_hhpp(channel3)%ival(ia,ib)
        ket= pp_hhpp(channel3)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = gmat 
        end if
        bra= hh_hhpp(channel3)%ival(ic,id)
        ket= pp_hhpp(channel3)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = gmat 
        end if
        !
        ! (ab)
        !
        bra= hh_hhpp(channel3)%ival(ib,ia)
        ket= pp_hhpp(channel3)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = phase_ab*gmat 
        end if
        bra= hh_hhpp(channel3)%ival(ic,id)
        ket= pp_hhpp(channel3)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = phase_ab*gmat 
        end if
        !
        ! (cd)
        !
        bra= hh_hhpp(channel3)%ival(ia,ib)
        ket= pp_hhpp(channel3)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = phase_cd*gmat 
        end if
        bra= hh_hhpp(channel3)%ival(id,ic)
        ket= pp_hhpp(channel3)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = phase_cd*gmat 
        end if
        !
        ! (ab)(cd)
        !
        bra= hh_hhpp(channel3)%ival(ib,ia)
        ket= pp_hhpp(channel3)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        bra= hh_hhpp(channel3)%ival(id,ic)
        ket= pp_hhpp(channel3)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hhpp%depot(channel3)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
     end if
     
     if ( channel4 /= 0 ) then
        
        !
        ! 1 
        !
        bra= hp_hphp(channel4)%ival(ia,ib)
        ket= hp_hphp(channel4)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = gmat 
        end if
        bra= hp_hphp(channel4)%ival(ic,id)
        ket= hp_hphp(channel4)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = gmat 
        end if
        !
        ! (ab)
        !
        bra= hp_hphp(channel4)%ival(ib,ia)
        ket= hp_hphp(channel4)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = phase_ab*gmat 
        end if
        bra= hp_hphp(channel4)%ival(ic,id)
        ket= hp_hphp(channel4)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = phase_ab*gmat 
        end if
        !
        ! (cd)
        !
        bra= hp_hphp(channel4)%ival(ia,ib)
        ket= hp_hphp(channel4)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = phase_cd*gmat 
        end if
        bra= hp_hphp(channel4)%ival(id,ic)
        ket= hp_hphp(channel4)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = phase_cd*gmat 
        end if
        !
        ! (ab)(cd)
        !
        bra= hp_hphp(channel4)%ival(ib,ia)
        ket= hp_hphp(channel4)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        bra= hp_hphp(channel4)%ival(id,ic)
        ket= hp_hphp(channel4)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hphp%depot(channel4)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
     end if
     
     if ( channel5 /= 0 ) then
        
        !
        ! 1 
        !
        bra= hp_hppp(channel5)%ival(ia,ib)
        ket= pp_hppp(channel5)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = gmat 
        end if
        bra= hp_hppp(channel5)%ival(ic,id)
        ket= pp_hppp(channel5)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = gmat 
        end if
        !
        ! (ab)
        !
        bra= hp_hppp(channel5)%ival(ib,ia)
        ket= pp_hppp(channel5)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = phase_ab*gmat 
        end if
        bra= hp_hppp(channel5)%ival(ic,id)
        ket= pp_hppp(channel5)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = phase_ab*gmat 
        end if
        !
        ! (cd)
        !
        bra= hp_hppp(channel5)%ival(ia,ib)
        ket= pp_hppp(channel5)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = phase_cd*gmat 
        end if
        bra= hp_hppp(channel5)%ival(id,ic)
        ket= pp_hppp(channel5)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = phase_cd*gmat 
        end if
        !
        ! (ab)(cd)
        !
        bra= hp_hppp(channel5)%ival(ib,ia)
        ket= pp_hppp(channel5)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        bra= hp_hppp(channel5)%ival(id,ic)
        ket= pp_hppp(channel5)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_hppp%depot(channel5)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
     end if
     
     !
     ! check if channel6 is my channel6
     !
     if ( channel6 /= 0 ) then 
        
        if ( check_my_channel(channel6) == 0 ) cycle
        !
        ! 1 
        !
        bra= pq_conf(channel6)%ival(ia,ib)
        ket= rs_conf(channel6)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = gmat 
        end if
        bra= pq_conf(channel6)%ival(ic,id)
        ket= rs_conf(channel6)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = gmat 
        end if
        !
        ! (ab)
        !
        bra= pq_conf(channel6)%ival(ib,ia)
        ket= rs_conf(channel6)%ival(ic,id)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = phase_ab*gmat 
        end if
        bra= pq_conf(channel6)%ival(ic,id)
        ket= rs_conf(channel6)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = phase_ab*gmat 
        end if
        !
        ! (cd)
        !
        bra= pq_conf(channel6)%ival(ia,ib)
        ket= rs_conf(channel6)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = phase_cd*gmat 
        end if
        bra= pq_conf(channel6)%ival(id,ic)
        ket= rs_conf(channel6)%ival(ia,ib)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = phase_cd*gmat 
        end if
        !
        ! (ab)(cd)
        !
        bra= pq_conf(channel6)%ival(ib,ia)
        ket= rs_conf(channel6)%ival(id,ic)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        bra= pq_conf(channel6)%ival(id,ic)
        ket= rs_conf(channel6)%ival(ib,ia)
        if ( bra*ket /= 0 ) then
           gmat_mtx_pppp%depot(channel6)%val(bra,ket) = phase_ab*phase_cd*gmat 
        end if
        
     end if
     amount = amount + 1
     
  ENDDO
6 CONTINUE

  
  
  deallocate( quantum_numbers_ch )
  deallocate( quantum_numbers_ch3 )
  
  !call mpi_barrier(mpi_comm_world,ierror)
  if ( iam ==0 )  write(6,*) 'Amount of G-matrix elements', amount
  
  call mytimer(tm2)
  if ( iam ==0 ) write(6,*) 'setting up interaction done!'
  


  
  if ( iam ==0 ) write(6,*) 'Time setting up interaction in pp-coupling', tm2-tm1, 'seconds'
  
  ! 
  ! Setup Fock-matrix in j-scheme
  !
  allocate( fock_mtx(all_orbit%total_orbits ,all_orbit%total_orbits )) 
  fock_mtx = 0.d0
  
  read(12,*)
  fock_mtx = 0.d0 
  DO WHILE (.TRUE.)
     !read in            t_z,     pi,     J, a, b, c, d, G
7    READ(12,*,end=8) aa,cc,gmat
     
     aa = aa + 1
     cc = cc + 1
        
     fock_mtx(aa,cc) = gmat 
     
  ENDDO
8 CONTINUE
  
  e0 = 0 
  do ia = 1, below_ef 
     e0 = e0 + fock_mtx(ia,ia)
  end do
  if ( iam == 0 ) write(6,*) 'Vacuum energy e0', e0
  !e0 = 45.367

  !deallocate ( work ) 
  
end subroutine setup_fockmtx_and_interaction_imsrg


!
! Setup two-body interaction in particle-hole coupled scheme
!
subroutine setup_particle_hole
  USE ang_mom_functions
  use parallel
  USE single_particle_orbits
  USE configurations
  USE constants
  use gmat_storage
  use one_body_operators
  use reorder 
  use HF_mod

  

  IMPLICIT NONE
  TYPE (configuration_descriptor) :: gmat_bra_configs 
  TYPE (configuration_descriptor) :: gmat_ket_configs 
  TYPE (configuration_descriptor) :: gmatrix_configs 
  INTEGER :: isospin_z, p_parity, ang_mom, number_channels
  INTEGER :: ph_channel, channel, channel1, channel2, nconfs, a, b,c,d, i1, j1, i,j, jt, j_min,j_max
  INTEGER :: ji,jj, ja,jb, aa, bb,cc,dd, bra, ket, iph, jbra_min, jket_min, jbra_max ,jket_max
  LOGICAL :: triag
  REAL*8 :: sixj, fnorm, ang_mom_factor, sum1, phase, angmom_fact, gmat, cross, dij
  INTEGER :: h1, p1, h2, p2, jtot, count
  REAL*8 :: sg, de, val1, val2, fact, ans1, ans2, gmat2, memory
  INTEGER :: li, itzi, lj, itzj, la, itza, lb, itzb, jc, jd
  INTEGER :: k1, k2, k3, ii, bra_conf, ket_conf, bra_confs, ket_confs, ia, ib,ic,id
  INTEGER :: bra1,ket1,bra2,ket2, pp_channel, itz_pp, ipar_pp
  real*4 :: tm1, tm2
  INTEGER :: total_work, work_per_proc, curr_work, curr_proc, curr_channel
  integer :: processor_work(num_procs), mtx_dim, mtx_dim_work,temp, tot_work, work_pr_proc, local_ch_id
  integer :: row_dim, column_dim, ket_dim, bra_dim, my_channel_low, my_channel_high 
  INTEGER :: bra_min, bra_max, ket_min, ket_max
  LOGICAL :: has_ch_been_added
  REAL*8, allocatable :: work(:)
  integer :: phase_cd, vec, vec_sum
  integer, allocatable :: ket_dim_channel(:)

  
  
  call setup_proc_mappings_ph
  number_channels = 0
  !     loop over isospin projection
  DO isospin_z=-1,1 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60
           
              
           
           !
           ! call number of bra_confs
           !
           CALL  number_ph_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,0,gmat_bra_configs)
           IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
           number_channels = number_channels + 1 
           
           

        end DO
     end DO
  end DO
  if ( iam ==0 ) write(6,*) 'number of  cross coupled channels', number_channels 
  
  ALLOCATE( lookup_ph_configs(1:4,number_channels))

  ph_channels%number_confs = number_channels
  ALLOCATE( ph_channels%config_J_ipar_tz(3*number_channels) )
  number_channels = 0
  !     loop over isospin projection
  DO isospin_z=-1,1 
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60 

           
                      
           !
           ! call number of bra_confs
           !
           CALL  number_ph_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,0,gmat_bra_configs)
           IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
           
           
           number_channels = number_channels + 1 
           k1 = number_channels*3 
           k2 = number_channels*3 - 1
           k3 = number_channels*3 - 2
           ph_channels%config_J_ipar_tz(k1) = ang_mom
           ph_channels%config_J_ipar_tz(k2) = p_parity
           ph_channels%config_J_ipar_tz(k3) = isospin_z
           

        end DO
     end DO
  end DO
  
  memory = 0.d0 
  do ii = 1, 4 
     !
     ! Count number of channels 
     !
     number_channels = 0
     do i1 = 1, ph_channels%number_confs 
        ang_mom   = ph_channels%config_J_ipar_tz(i1*3)
        p_parity  = ph_channels%config_J_ipar_tz(i1*3-1)
        isospin_z = ph_channels%config_J_ipar_tz(i1*3-2)


        if ( ii == 1 ) then
           bra_conf = 1 
        elseif ( ii == 2 ) then
           bra_conf = 2 
        elseif ( ii == 3 ) then
           bra_conf = 3 
        elseif ( ii == 4 ) then
           bra_conf = 4
        end if
        
        !
        ! call number of bra_confs
        !
        CALL  number_ph_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
        
        bra_confs=gmat_bra_configs%number_confs
        
        ALLOCATE( lookup_ph_configs(ii,i1)%ival(2,bra_confs) )
        ALLOCATE(gmat_bra_configs%config_ab(bra_confs+bra_confs) )
        CALL setup_ph_gmatrix_configurations &
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        
        memory = memory + dble(ii*i1*2*bra_confs)*4./1.e9 
        lookup_ph_configs(ii,i1)%ival = 0
        
        DO i=1,gmat_bra_configs%number_confs
           ia= gmat_bra_configs%config_ab(i*2-1)
           ib= gmat_bra_configs%config_ab(i*2)
           
           lookup_ph_configs(ii,i1)%ival(1,i) = ia
           lookup_ph_configs(ii,i1)%ival(2,i) = ib
        end DO
        
        deallocate( gmat_bra_configs%config_ab ) 
     end do
  end do
  
  if ( iam == 0 ) write(6,*) 'Memory for lookup ph configs', memory, 'GByte'
 

  ALLOCATE( locate_crosscoupled_channel(1:11,-1:1, 0:1, 0:60) ) 
  
  locate_crosscoupled_channel = 0 
  
  do ii = 1, 11
     !
     ! Count number of channels 
     !
     number_channels = 0
     do i1 = 1, ph_channels%number_confs 
        ang_mom   = ph_channels%config_J_ipar_tz(i1*3)
        p_parity  = ph_channels%config_J_ipar_tz(i1*3-1)
        isospin_z = ph_channels%config_J_ipar_tz(i1*3-2)
        
        if ( ii == 1 ) then
           bra_conf = 3 
           ket_conf = 2
        elseif ( ii == 2 ) then
           bra_conf = 2 
           ket_conf = 2 
        elseif ( ii == 3 ) then
           bra_conf = 2 
           ket_conf = 1 
        elseif ( ii == 4 ) then
           bra_conf = 2 
           ket_conf = 3 
        elseif ( ii == 5 ) then
           bra_conf = 2 
           ket_conf = 4 
        elseif ( ii == 6 ) then
           bra_conf = 3 
           ket_conf = 4 
        elseif ( ii == 7 ) then
           bra_conf = 3 
           ket_conf = 1 
        elseif ( ii == 8 ) then
           bra_conf = 4 
           ket_conf = 4 
        elseif ( ii == 9 ) then
           bra_conf = 1 
           ket_conf = 4 
        elseif ( ii == 10 ) then
           bra_conf = 3 
           ket_conf = 3 
        elseif ( ii == 11 ) then
           bra_conf = 1 
           ket_conf = 1 
        end if
        
        !
        ! call number of bra_confs
        !
        CALL  number_ph_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
        !
        ! call number of ket_confs
        !
        CALL  number_ph_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
        IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
        
        
        
        number_channels = number_channels + 1 
        locate_crosscoupled_channel(ii,isospin_z, p_parity, ang_mom) = number_channels
              
              
     end DO
     
     if ( ii == 1 ) allocate( gmat_hp_phhp%depot(1:number_channels) )
     if ( ii == 1 ) allocate( t2_cross_ccm%depot(1:number_channels) )
     if ( ii == 1 ) allocate( t2_cross_ccm_eqn%depot(1:number_channels) )
     
     
     if ( ii == 2 ) allocate( gmat_hp_hphp%depot(1:number_channels) )
     if ( ii == 3 ) allocate( gmat_hp_hphh%depot(1:number_channels) )
     if ( ii == 4 ) allocate( gmat_hp_hpph%depot(1:number_channels) )
     !if ( ii == 5 ) allocate( gmat_hp_hppp%depot(1:number_channels) )
          
     
     if ( ii == 1 ) allocate( ph_phhp_cross(1:number_channels), hp_phhp_cross(1:number_channels))
     if ( ii == 2 ) allocate( hp_hphp_cross(1:number_channels) )
     if ( ii == 3 ) allocate( hp_hphh_cross(1:number_channels), hh_hphh_cross(1:number_channels) )
     if ( ii == 4 ) allocate( hp_hpph_cross(1:number_channels), ph_hpph_cross(1:number_channels))
     !if ( ii == 5 ) allocate( hp_hppp_cross(1:number_channels), pp_hppp_cross(1:number_channels))
     if ( ii == 5 ) allocate( hp_hppp_cross_full(1:number_channels), pp_hppp_cross_full(1:number_channels))
     if(ii==10) allocate(ph_phph_cross(1:number_channels))
     if(ii==11) allocate(hh_hhhh_cross(1:number_channels))
     if ( ii == 7 ) allocate( hh_hhph_cross(1:number_channels), ph_hhph_cross(1:number_channels))

     
  end do
  
  

  if ( iam ==0 ) write(6,*) 'Number of hhpp channels', size( gmat_mtx_hhpp%depot(:) )
  if ( iam ==0 ) write(6,*) 'Number of cross-coupled hhpp channels', size( gmat_hp_phhp%depot(:) )
  
  do i=1, size(ph_phhp_cross)
     allocate( ph_phhp_cross(i)%ival(below_ef+1:tot_orbs, 1:below_ef) )
     allocate( hp_phhp_cross(i)%ival(1:below_ef, below_ef+1:tot_orbs) )
     ph_phhp_cross(i)%ival = 0
     hp_phhp_cross(i)%ival = 0
  end do
  
  do i=1, size(hp_hphp_cross)
     allocate( hp_hphp_cross(i)%ival(1:below_ef, below_ef+1:tot_orbs) )
     hp_hphp_cross(i)%ival = 0
  end do
  do i=1, size(hp_hphh_cross)
     allocate( hp_hphh_cross(i)%ival(1:below_ef,below_ef+1:tot_orbs) )
     allocate( hh_hphh_cross(i)%ival(1:below_ef, 1:below_ef) )
     hp_hphh_cross(i)%ival = 0
     hh_hphh_cross(i)%ival = 0
  end do
  
  do i=1, size(hp_hpph_cross)
     allocate( ph_hpph_cross(i)%ival(below_ef+1:tot_orbs, 1:below_ef) )
     allocate( hp_hpph_cross(i)%ival(1:below_ef, below_ef+1:tot_orbs) )
     ph_hpph_cross(i)%ival = 0
     hp_hpph_cross(i)%ival = 0
  end do

  do i=1, size(hp_hppp_cross)
     allocate( hp_hppp_cross(i)%ival(1:below_ef, below_ef+1:tot_orbs) )
     allocate( pp_hppp_cross(i)%ival(below_ef+1:tot_orbs, below_ef+1:tot_orbs) )
     pp_hppp_cross(i)%ival = 0
     hp_hppp_cross(i)%ival = 0
  end do


  do i=1, size(hp_hppp_cross_full)
     allocate( hp_hppp_cross_full(i)%ival(1:below_ef, below_ef+1:tot_orbs) )
     allocate( pp_hppp_cross_full(i)%ival(below_ef+1:tot_orbs, below_ef+1:tot_orbs) )
     pp_hppp_cross_full(i)%ival = 0
     hp_hppp_cross_full(i)%ival = 0
  end do



  do i=1, size(hh_hhph_cross)
     allocate( hh_hhph_cross(i)%ival(1:below_ef, 1:below_ef) )
     allocate( ph_hhph_cross(i)%ival(below_ef+1:tot_orbs, 1:below_ef) )
     hh_hhph_cross(i)%ival = 0
     ph_hhph_cross(i)%ival = 0
  end do

  do i=1,size(hh_hhhh_cross)
    allocate(hh_hhhh_cross(i)%ival(1:below_ef,1:below_ef))
    hh_hhhh_cross(i)%ival=0
  enddo

  
  call mytimer(tm1)
  do ii = 1, 4
     number_channels = 0
     
     do i1 = 1, ph_channels%number_confs 
        ang_mom   = ph_channels%config_J_ipar_tz(i1*3)
        p_parity  = ph_channels%config_J_ipar_tz(i1*3-1)
        isospin_z = ph_channels%config_J_ipar_tz(i1*3-2)
        
        if ( ii == 1 ) then
           bra_conf = 3 
           ket_conf = 2
        elseif ( ii == 2 ) then
           bra_conf = 2 
           ket_conf = 2 
        elseif ( ii == 3 ) then
           bra_conf = 2 
           ket_conf = 1 
        elseif ( ii == 4 ) then
           bra_conf = 2 
           ket_conf = 3
        elseif ( ii == 5 ) then
           bra_conf = 2 
           ket_conf = 4
        end if
              
        !
        ! call number of bra_confs
        !
        CALL  number_ph_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
        !
        ! call number of ket_confs
        !
        CALL  number_ph_gmatrix_confs&
             (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
        IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
              
        bra_confs=gmat_bra_configs%number_confs
        ket_confs=gmat_ket_configs%number_confs
              
        ALLOCATE(gmat_bra_configs%config_ab(bra_confs+bra_confs) )
        ALLOCATE(gmat_ket_configs%config_ab(ket_confs+ket_confs) )
              
        CALL setup_ph_gmatrix_configurations &
             (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        
        CALL setup_ph_gmatrix_configurations &
             (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
        
        number_channels = number_channels + 1 
        

              
        
        if ( ii == 1 ) allocate( gmat_hp_phhp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 1 ) allocate( t2_cross_ccm%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 1 ) allocate( t2_cross_ccm_eqn%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
     
        if ( ii == 2 ) allocate( gmat_hp_hphp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 3 ) allocate( gmat_hp_hphh%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 4 ) allocate( gmat_hp_hpph%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        if ( ii == 5 ) allocate( gmat_hp_hppp%depot(number_channels)%val(1:bra_confs, 1:ket_confs) ) 
        
        if ( ii == 1 ) gmat_hp_phhp%depot(number_channels)%val = 0.d0 
        if ( ii == 1 ) t2_cross_ccm_eqn%depot(number_channels)%val = 0.d0
        if ( ii == 1 ) t2_cross_ccm%depot(number_channels)%val = 0.d0
        
        
        if ( ii == 2 ) gmat_hp_hphp%depot(number_channels)%val = 0.d0
        if ( ii == 3 ) gmat_hp_hphh%depot(number_channels)%val = 0.d0
        if ( ii == 4 ) gmat_hp_hpph%depot(number_channels)%val = 0.d0 
        if ( ii == 5 ) gmat_hp_hppp%depot(number_channels)%val = 0.d0 
        
        
        if  ( ii == 1 ) then
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)

              ph_phhp_cross(number_channels)%ival(ia,ib) = i 
              
           end DO
          
           DO i=1,gmat_ket_configs%number_confs
              ia= gmat_ket_configs%config_ab(i*2-1)
              ib= gmat_ket_configs%config_ab(i*2)
                    
              hp_phhp_cross(number_channels)%ival(ia,ib) = i 
              
           end DO
          
 
        elseif  ( ii == 2 ) then
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              
              hp_hphp_cross(number_channels)%ival(ia,ib) = i 
                    
           end DO
        
        elseif( ii  == 3 ) then 
               
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              
              hp_hphh_cross(number_channels)%ival(ia,ib) = i 
                    
           end DO
                 
           DO i=1,gmat_ket_configs%number_confs
              ia= gmat_ket_configs%config_ab(i*2-1)
              ib= gmat_ket_configs%config_ab(i*2)
                    
              hh_hphh_cross(number_channels)%ival(ia,ib) = i 
                    
           end DO
           
        elseif  ( ii == 4 ) then
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              
              hp_hpph_cross(number_channels)%ival(ia,ib) = i 
              
           end DO
           
           DO i=1,gmat_ket_configs%number_confs
              ia= gmat_ket_configs%config_ab(i*2-1)
              ib= gmat_ket_configs%config_ab(i*2)
                    
              ph_hpph_cross(number_channels)%ival(ia,ib) = i 
              
           end DO
           
        elseif  ( ii == 5 ) then
           DO i=1,gmat_bra_configs%number_confs
              ia= gmat_bra_configs%config_ab(i*2-1)
              ib= gmat_bra_configs%config_ab(i*2)
              
              hp_hppp_cross(number_channels)%ival(ia,ib) = i 
              
           end DO
           
           DO i=1,gmat_ket_configs%number_confs
              ia= gmat_ket_configs%config_ab(i*2-1)
              ib= gmat_ket_configs%config_ab(i*2)
                    
              pp_hppp_cross(number_channels)%ival(ia,ib) = i 
              
           end DO
        
        end if
        
        
        !
        ! Setup cross-coupled matrix elements gmat_ph_hhpp%depot(:)%val(:,:)
        ! 
        ! Coupling direction: ingoing -> outgoing 
        ! <ij||ab> (ai->J) and (bj->J) if one changes coupling direction 
        ! to outgoing -> ingoing one gets a phase (-1)**(ji+jj+ja+jb - 2J ) 
        ! 
        
        DO i=1,gmat_bra_configs%number_confs
           a= gmat_bra_configs%config_ab(i*2-1)
           c= gmat_bra_configs%config_ab(i*2)
           
           ja = all_orbit%jj(a)
           jc = all_orbit%jj(c)

           aa = indx_inv(a) 
           cc = indx_inv(c) 
                    
           DO j=1,gmat_ket_configs%number_confs
              d= gmat_ket_configs%config_ab(j+j-1)
              b= gmat_ket_configs%config_ab(j+j)
              
              jd = all_orbit%jj(d)
              jb = all_orbit%jj(b)
              
              dd = indx_inv(d) 
              bb = indx_inv(b) 
              
              
                 

              
              
              
              itz_pp = (all_orbit%itzp(a) + all_orbit%itzp(b))/2
              ipar_pp = ( 1 - (-1)**( all_orbit%ll(a) + all_orbit%ll(b) ))/2
              j_min=max( ABS((all_orbit%jj(a)-all_orbit%jj(b))/2), ABS((all_orbit%jj(c)-all_orbit%jj(d))/2) )
              j_max=min( (all_orbit%jj(a)+all_orbit%jj(b))/2, (all_orbit%jj(c)+all_orbit%jj(d))/2 ) 
              

              !fnorm= sqrt(2.*ang_mom+1.)*iph((ji+jb)/2+ang_mom) 
              fnorm=iph((ja+jd)/2+ang_mom) 
              
              IF(j_min > j_max) CYCLE
              cross=0.
              DO jt=j_min,j_max
                 
                 IF(( c == d ).AND.(MOD(jt,2)/=0)) CYCLE
                 IF(( a == b ).AND.(MOD(jt,2)/=0)) CYCLE
                 IF(triag(all_orbit%jj(c),all_orbit%jj(d),2*jt)) CYCLE
                 IF(triag(all_orbit%jj(a),all_orbit%jj(b),2*jt)) CYCLE
              
                 ang_mom_factor=sjs(jc,ja,2*ang_mom, jb,jd,2*jt)* &
                      (2.*jt+1.)*fnorm*iph(jt)
                    gmat2=0.d0

                 if ( ii == 1 ) then
                    pp_channel = locate_gmatchannel(3,itz_pp, ipar_pp, jt ) 
                    bra = hh_hhpp(pp_channel)%ival(c,d) 
                    ket = pp_hhpp(pp_channel)%ival(a,b) 
                    gmat2 = gmat_mtx_hhpp%depot(pp_channel)%val(bra,ket) 
 
                 elseif ( ii == 2 ) then 
                    pp_channel = locate_gmatchannel(4,itz_pp, ipar_pp, jt ) 
                    
                    
                    bra = hp_hphp(pp_channel)%ival(a,b)
                    ket = hp_hphp(pp_channel)%ival(d,c)
                    
                    phase = (-1.d0)**( (jd + jc)/2 + jt + 1 )
                    gmat2 = phase * gmat_mtx_hphp%depot(pp_channel)%val(bra,ket) 
                    
                 elseif ( ii == 3 ) then 
                    pp_channel = locate_gmatchannel(2,itz_pp, ipar_pp, jt ) 
                    
                    bra = hh_hhhp(pp_channel)%ival(a,b)
                    ket = hp_hhhp(pp_channel)%ival(d,c) 
                    phase = (-1.d0)**( (jd + jc)/2 + jt + 1 )
                    gmat2 = phase * gmat_mtx_hhhp%depot(pp_channel)%val(bra,ket) 
                                        
                 elseif ( ii == 4 ) then 
                    pp_channel = locate_gmatchannel(3,itz_pp, ipar_pp, jt ) 
                    bra = hh_hhpp(pp_channel)%ival(a,b)
                    ket = pp_hhpp(pp_channel)%ival(c,d) 
                    gmat2 = gmat_mtx_hhpp%depot(pp_channel)%val(bra,ket) 
                    

                 elseif ( ii == 5 ) then
                    pp_channel = locate_gmatchannel(5,itz_pp, ipar_pp, jt ) 
                    bra = hp_hppp(pp_channel)%ival(a,b)
                    ket = pp_hppp(pp_channel)%ival(c,d) 
                    
                    phase_cd = 1.d0
                    if ( c > d ) then
                       ket = pp_hppp(pp_channel)%ival(d,c)
                       phase_cd = iph( (all_orbit%jj(c) + all_orbit%jj(d))/2 + jt +1 )
                    end if
                    
                    gmat2 = phase_cd * gmat_mtx_hppp%depot(pp_channel)%val(bra,ket) 
                    
                    
                 end if

                 !CALL jcoupled_gmatrix_element_new (aa,bb,cc,dd,jt,gmat)
                 cross = cross + ang_mom_factor*gmat2
                 !if ( abs( gmat - gmat2) > 1.e-6  ) write(6,*) gmat , gmat2

              end DO
              
              if ( ii == 1 ) gmat_hp_phhp%depot(number_channels)%val(i,j) = cross
              if ( ii == 2 ) gmat_hp_hphp%depot(number_channels)%val(i,j) = cross
              if ( ii == 3 ) gmat_hp_hphh%depot(number_channels)%val(i,j) = cross
              if ( ii == 4 ) gmat_hp_hpph%depot(number_channels)%val(i,j) = cross
              if ( ii == 5 ) gmat_hp_hppp%depot(number_channels)%val(i,j) = cross
              
           end DO
        end DO
     end DO
          
  end do




  number_channels = 0
  do i1 = 1, ph_channels%number_confs 
     ang_mom   = ph_channels%config_J_ipar_tz(i1*3)
     p_parity  = ph_channels%config_J_ipar_tz(i1*3-1)
     isospin_z = ph_channels%config_J_ipar_tz(i1*3-2)
     
     bra_conf = 2 
     ket_conf = 4
     
     !
     ! call number of bra_confs
     !
     CALL  number_ph_gmatrix_confs&
          (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
     IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
     !
     ! call number of ket_confs
     !
     CALL  number_ph_gmatrix_confs&
          (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
     IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
     
     bra_confs=gmat_bra_configs%number_confs
     ket_confs=gmat_ket_configs%number_confs
     
     ALLOCATE(gmat_bra_configs%config_ab(bra_confs+bra_confs) )
     ALLOCATE(gmat_ket_configs%config_ab(ket_confs+ket_confs) )
              
     CALL setup_ph_gmatrix_configurations &
          (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        
     CALL setup_ph_gmatrix_configurations &
          (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
        number_channels=number_channels+1
             

     DO i=1,gmat_bra_configs%number_confs
        ia= gmat_bra_configs%config_ab(i*2-1)
        ib= gmat_bra_configs%config_ab(i*2)
        
        hp_hppp_cross_full(number_channels)%ival(ia,ib) = i 
              
     end DO
           
     DO i=1,gmat_ket_configs%number_confs
        ia= gmat_ket_configs%config_ab(i*2-1)
        ib= gmat_ket_configs%config_ab(i*2)
                    
        pp_hppp_cross_full(number_channels)%ival(ia,ib) = i 
              
     end DO
  end DO
    


  
  number_channels = 0
  do i1 = 1, ph_channels%number_confs 
     ang_mom   = ph_channels%config_J_ipar_tz(i1*3)
     p_parity  = ph_channels%config_J_ipar_tz(i1*3-1)
     isospin_z = ph_channels%config_J_ipar_tz(i1*3-2)
     
     
     number_channels = locate_crosscoupled_channel(5,isospin_z, p_parity, ang_mom)
     
     if ( number_channels == 0 ) cycle
     local_ch_id = number_channels
     if ( check_my_channel_hp_hppp(local_ch_id) == 0 ) cycle
     
     bra_conf = 2 
     ket_conf = 4
     
     !
     ! call number of bra_confs
     !
     CALL  number_ph_gmatrix_confs&
          (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
     IF (gmat_bra_configs%number_confs <= 0 ) CYCLE
     !
     ! call number of ket_confs
     !
     CALL  number_ph_gmatrix_confs&
          (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
     IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
     
     bra_confs=gmat_bra_configs%number_confs
     ket_confs=gmat_ket_configs%number_confs
     
     ALLOCATE(gmat_bra_configs%config_ab(bra_confs+bra_confs) )
     ALLOCATE(gmat_ket_configs%config_ab(ket_confs+ket_confs) )
              
     CALL setup_ph_gmatrix_configurations &
          (ang_mom,p_parity,isospin_z,bra_conf,gmat_bra_configs)
        
     CALL setup_ph_gmatrix_configurations &
          (ang_mom,p_parity,isospin_z,ket_conf,gmat_ket_configs)
             
     !
     ! ket side if fully stored on each proc
     !
     ket_min = mapping_hp_hppp(iam+1,local_ch_id,2)
     ket_max = mapping_hp_hppp(iam+1,local_ch_id,3)
     !
     ! bra side is distributed 
     !
     bra_min = mapping_hp_hppp(iam+1,local_ch_id,4)
     bra_max = mapping_hp_hppp(iam+1,local_ch_id,5)
     
     
     allocate( gmat_hp_hppp%depot(number_channels)%val(bra_min:bra_max, ket_min:ket_max) )
     gmat_hp_hppp%depot(number_channels)%val = 0.d0 
        



     DO i=bra_min, bra_max !1,gmat_bra_configs%number_confs
        ia= gmat_bra_configs%config_ab(i*2-1)
        ib= gmat_bra_configs%config_ab(i*2)
        
        hp_hppp_cross(number_channels)%ival(ia,ib) = i 
              
     end DO
           
     DO i=ket_min, ket_max !1,gmat_ket_configs%number_confs
        ia= gmat_ket_configs%config_ab(i*2-1)
        ib= gmat_ket_configs%config_ab(i*2)
                    
        pp_hppp_cross(number_channels)%ival(ia,ib) = i 
              
     end DO
     
     !
     ! Setup cross-coupled matrix elements gmat_ph_hhpp%depot(:)%val(:,:)
     ! 
     ! Coupling direction: ingoing -> outgoing 
     ! <ij||ab> (ai->J) and (bj->J) if one changes coupling direction 
     ! to outgoing -> ingoing one gets a phase (-1)**(ji+jj+ja+jb - 2J ) 
     ! 
     
     DO i=bra_min, bra_max !1,gmat_bra_configs%number_confs
        a= gmat_bra_configs%config_ab(i*2-1)
        c= gmat_bra_configs%config_ab(i*2)
        
        ja = all_orbit%jj(a)
        jc = all_orbit%jj(c)

        aa = indx_inv(a) 
        cc = indx_inv(c) 
                    
        DO j=ket_min, ket_max !1,gmat_ket_configs%number_confs
           d= gmat_ket_configs%config_ab(j+j-1)
           b= gmat_ket_configs%config_ab(j+j)
           
           jd = all_orbit%jj(d)
           jb = all_orbit%jj(b)
           
           dd = indx_inv(d) 
           bb = indx_inv(b) 
           
           itz_pp = (all_orbit%itzp(a) + all_orbit%itzp(b))/2
           ipar_pp = ( 1 - (-1)**( all_orbit%ll(a) + all_orbit%ll(b) ))/2
           j_min=max( ABS((all_orbit%jj(a)-all_orbit%jj(b))/2), ABS((all_orbit%jj(c)-all_orbit%jj(d))/2) )
           j_max=min( (all_orbit%jj(a)+all_orbit%jj(b))/2, (all_orbit%jj(c)+all_orbit%jj(d))/2 ) 
              

           !fnorm= sqrt(2.*ang_mom+1.)*iph((ji+jb)/2+ang_mom) 
           fnorm=iph((ja+jd)/2+ang_mom) 
              
           IF(j_min > j_max) CYCLE
           cross=0.
           DO jt=j_min,j_max
                 
              IF(( c == d ).AND.(MOD(jt,2)/=0)) CYCLE
              IF(( a == b ).AND.(MOD(jt,2)/=0)) CYCLE
              IF(triag(all_orbit%jj(c),all_orbit%jj(d),2*jt)) CYCLE
              IF(triag(all_orbit%jj(a),all_orbit%jj(b),2*jt)) CYCLE
              
              ang_mom_factor=sjs(jc,ja,2*ang_mom, jb,jd,2*jt)* &
                   (2.*jt+1.)*fnorm*iph(jt)
              
              pp_channel = locate_gmatchannel(5,itz_pp, ipar_pp, jt ) 
              bra = hp_hppp(pp_channel)%ival(a,b)
              ket = pp_hppp(pp_channel)%ival(c,d) 
                    
              phase_cd = 1.d0
              if ( c > d ) then
                 ket = pp_hppp(pp_channel)%ival(d,c)
                 phase_cd = iph( (all_orbit%jj(c) + all_orbit%jj(d))/2 + jt +1 )
              end if
              
              gmat2 = phase_cd * gmat_mtx_hppp%depot(pp_channel)%val(bra,ket) 
                    
              cross = cross + ang_mom_factor*gmat2 
           end DO
           gmat_hp_hppp%depot(number_channels)%val(i,j) = cross
              
        end DO
     end DO
  end DO
    

  call mytimer(tm2)

  call mpi_barrier(mpi_comm_world,ierror)
  if ( iam ==0 ) write(6,*) 'Time setting up interaction in cross-coupling', tm2-tm1, 'seconds'
  !stop
  
  
  
!  deallocate( indx, indx_inv )

end subroutine setup_particle_hole


!*****7**-*********-*********-*********-*********-*********-*********-72
      SUBROUTINE index(n,arr,indx)
! integer array arr is indexed in ascending order
      INTEGER :: n,indx(*),M,NSTACK
      integer :: arr(*)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      integer :: a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) write(*,*)'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END subroutine index
!*****7**-*********-*********-*********-*********-*********-*********-72

  

! jobvl, jobvr, np, temp1, lda, omega2_eig, vl, ldvl, eigen_vec, ldvr, &
!       work1, lwork, rwork, info
SUBROUTINE lapack_diag(h, cvec, ceig, n )
  

  implicit none
  integer, intent(in) :: n
  complex*16, dimension(n,n), intent(in) :: h
  COMPLEX*16, DIMENSION(n,n), intent(out) :: cvec
  COMPLEX*16, DIMENSION(n,n) ::  vl
  COMPLEX*16, DIMENSION(n), intent(out) :: ceig
  DOUBLE PRECISION, DIMENSION(2*n) :: rwork
  COMPLEX*16, DIMENSION(10000) :: work1
  INTEGER :: lda, ldb, ldvl, ldvr, info, lwork, ilo , ihi
  CHARACTER*1 :: jobvl, jobvr
  DOUBLE PRECISION, DIMENSION(n) :: scale, rconde, rcondv
  complex*16 :: norm
  integer :: i, j

  jobvl = 'N' ;  jobvr = 'V';  lda = n
  ldvl = 1;  ldvr = n;  lwork = 10000
  ceig = 0.; cvec = 0.
  CALL zgeev( jobvl, jobvr, n, h, lda, ceig, vl, ldvl, cvec, ldvr, &
       work1, lwork, rwork, info )
  
  ! berggren normalization
  do i = 1, n
     norm = sum( cvec(:,i)*cvec(:,i) )
     cvec(:,i) = cvec(:,i)/sqrt(norm)
     !write(6,*) i, norm

  end do
  
  do i = 1, n
     do j = 1, n
!        if ( i /=j .and. abs(sum( (cvec(:,i))*cvec(:,j) ) ) > 1.e-5 ) then
!           write(6,*)  i,j, sum( (cvec(:,i))*cvec(:,j) ) 
!           !stop
!        end if
     end do
  end do
  
  call hf_sort_cmplx2(ceig,cvec, n)

end SUBROUTINE lapack_diag

!
! eigenvalue sort
! sort cmplx vector real(a(1))<real(a(2)) < ... < real(a(n))
!
SUBROUTINE hf_sort_cmplx2(a,b, n)
  IMPLICIT NONE
  
  INTEGER :: i, j, n
  complex*16, DIMENSION(n), INTENT(INOUT) :: a
  complex*16, DIMENSION(n,n), INTENT(INOUT) :: b
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
           


        END IF
     END DO
  END DO

END SUBROUTINE hf_sort_cmplx2

subroutine setup_dipole_operator
  USE partial_waves
  USE relcm_gmatrix
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  USE constants
  use ang_mom_functions
  use one_body_operators
  
  implicit none
  !
  integer :: iph, a, na, la,ja,tza,c, nc,lc,jc,tzc, i, nn
  double precision :: ang_fac
  double precision :: osc_a, osc_c, vsum,vsum1, osc_c_deriv, factor
  double precision, allocatable, dimension(:) :: rr ,wrr
  double precision, allocatable :: hol(:,:), temp_kin(:,:), temp_kin2(:,:)
  DOUBLE PRECISION ::  sum_rel, cx(0:200), xp
  real*8 :: zlab, delta, re_hcom 
  integer :: ph, pp, j, lambda
  
  allocate( reduced_r1(all_orbit%total_orbits, all_orbit%total_orbits) )
  
  lambda = 1
  CALL commons_to_angmom  
  
  oscl=hbarc/SQRT(p_mass*hbar_omega)

  nn = 600
  allocate( rr(nn), wrr(nn) )
  call gauss_legendre(0.d0, 25.d0, rr, wrr, nn)
  allocate( hol(all_orbit%total_orbits, nn) )
  
  hol = 0.d0
  do a = 1, all_orbit%total_orbits
     
     na=all_orbit%nn(a) 
     la=all_orbit%ll(a)
     ja=all_orbit%jj(a)
     tza = all_orbit%itzp(a)

     ph = 1.d0 
     factor = 0.5D0*((na+la+2)*LOG(2.D0)+fac(na)-dfac(2*na+2*la+1)-0.5D0*LOG(pi))
     factor = EXP(factor)
     sum_rel=0.
     DO i=1,nn
        
        zlab= rr(i)/oscl
        
        CALL laguerre_general( na, la+0.5D0, zlab*zlab, cx )
        xp = cx(na)*exp(-zlab*zlab*0.5d0)*(zlab**la)
        hol(a,i) = xp*factor*(1.d0/oscl**(1.5D0))
        sum_rel = sum_rel + hol(a, i)**2 * rr(i)**2*wrr(i)
        
     end DO
     if ( iam == 0 ) write(6,*) 'norm', sum_rel

  end do
  
  reduced_r1 = 0.d0
  do a = 1, all_orbit%total_orbits
     
     na=all_orbit%nn(a) 
     la=all_orbit%ll(a)
     ja=all_orbit%jj(a)
     tza = all_orbit%itzp(a)
     
     do c = 1, all_orbit%total_orbits
        nc=all_orbit%nn(c)
        lc=all_orbit%ll(c)
        jc=all_orbit%jj(c)
        tzc = all_orbit%itzp(c)
        

        if ( tza /= tzc ) cycle
        
        ang_fac = iph( (jc + 1)/2 + la + lc )* dsqrt(jc+1.d0) * & 
             dsqrt( 2.*la+1.d0)*dsqrt( 2.*lc+1.d0) * tjs(2*lc,2,2*la,0,0,0) * sjs( 2*lc, 1, jc, ja, 2, 2*la)

        vsum = 0.d0
        do i = 1, nn
           
           
           osc_a = hol(a,i)
           osc_c = hol(c,i)
           
           vsum = vsum + wrr(i)*rr(i)**2 * osc_a * rr(i)**(lambda) *osc_c
     
        end do
        
        if ( lambda == 1 ) reduced_r1(a,c) = ang_fac * vsum
        if ( lambda == 2 ) reduced_r1(a,c) = ang_fac * vsum/dsqrt(4.d0*pi) ! presumably the quadrupole
        
        
     end do
  end do
  
  nn = 0
  pp = 0 
  DO i=1,below_ef

     if ( all_orbit%itzp(i) == -1 ) pp = pp + all_orbit%jj(i) + 1 
     if ( all_orbit%itzp(i) ==  1 ) nn = nn + all_orbit%jj(i) + 1  

  ENDDO
  
  if ( iam == 0 ) write(6,*) 'Number of neutrons', nn
  if ( iam == 0 ) write(6,*) 'Number of protons', pp
  

  do i = 1, all_orbit%total_orbits
     do j = 1, all_orbit%total_orbits
        tza = all_orbit%itzp(i)
        tzc = all_orbit%itzp(j)
        re_hcom = reduced_r1(i,j)
        if ( lambda == 1 ) theta_1b(i,j) = re_hcom * delta(tza,-1) * delta(tza,tzc) - delta(tza,tzc) * pp * re_hcom /mass_nucleus
        if ( lambda == 2 ) theta_1b(i,j) = re_hcom * delta(tza,tzc) 
     end do
  end do
  
  allocate( temp_kin(all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( temp_kin2(all_orbit%total_orbits, all_orbit%total_orbits) )
  
  
  temp_kin2 = 0.d0
  temp_kin2 = theta_1b
  temp_kin = 0.d0
  temp_kin = matmul( theta_1b,(u_hf) )
  theta_1b = 0.d0
  theta_1b = matmul( transpose(u_hf), temp_kin )
  
  
!  do i = 1, all_orbit%total_orbits
!     do j = 1, all_orbit%total_orbits
        
        !if ( abs(theta_1b(i,j) ) < 1.e-10 ) cycle 
        !write(6,'(a,1x,4(I3,1x),g)') 'check dipole', all_orbit%jj(i), all_orbit%ll(i), all_orbit%jj(j), all_orbit%ll(j), real(theta_1b(i,j))
!     end do
!  end do
  
  
  deallocate( temp_kin, temp_kin2,reduced_r1 )
  
end subroutine setup_dipole_operator



!============================================
subroutine setup_dipole_no_CMcorr_operator
  USE partial_waves
  USE relcm_gmatrix
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  USE constants
  use ang_mom_functions
  use one_body_operators
  
  implicit none
  !double precision :: reduced_r1(all_orbit%total_orbits, all_orbit%total_orbits) 
  integer :: iph, a, na, la,ja,tza,c, nc,lc,jc,tzc, i, nn
  double precision :: ang_fac
  double precision :: osc_a, osc_c, vsum,vsum1, osc_c_deriv, factor
  double precision, allocatable, dimension(:) :: rr ,wrr
  double precision, allocatable :: hol(:,:), temp_kin(:,:), temp_kin2(:,:)
  DOUBLE PRECISION ::  sum_rel, cx(0:200), xp
  real*8 :: zlab, delta, re_hcom 
  integer :: ph, pp, j, lambda
  
  allocate( reduced_r1(all_orbit%total_orbits, all_orbit%total_orbits) )
  lambda = 1
  CALL commons_to_angmom  
  
  oscl=hbarc/SQRT(p_mass*hbar_omega)

  nn = 600
  allocate( rr(nn), wrr(nn) )
  call gauss_legendre(0.d0, 25.d0, rr, wrr, nn)
  allocate( hol(all_orbit%total_orbits, nn) )
  
  hol = 0.d0
  do a = 1, all_orbit%total_orbits
     
     na=all_orbit%nn(a) 
     la=all_orbit%ll(a)
     ja=all_orbit%jj(a)
     tza = all_orbit%itzp(a)

     ph = 1.d0 
     factor = 0.5D0*((na+la+2)*LOG(2.D0)+fac(na)-dfac(2*na+2*la+1)-0.5D0*LOG(pi))
     factor = EXP(factor)
     sum_rel=0.
     DO i=1,nn
        
        zlab= rr(i)/oscl
        
        CALL laguerre_general( na, la+0.5D0, zlab*zlab, cx )
        xp = cx(na)*exp(-zlab*zlab*0.5d0)*(zlab**la)
        hol(a,i) = xp*factor*(1.d0/oscl**(1.5D0))
        sum_rel = sum_rel + hol(a, i)**2 * rr(i)**2*wrr(i)
        
     end DO
     if ( iam == 0 ) write(6,*) 'norm', sum_rel

  end do
  
  reduced_r1 = 0.d0
  do a = 1, all_orbit%total_orbits
     
     na=all_orbit%nn(a) 
     la=all_orbit%ll(a)
     ja=all_orbit%jj(a)
     tza = all_orbit%itzp(a)
     
     do c = 1, all_orbit%total_orbits
        nc=all_orbit%nn(c)
        lc=all_orbit%ll(c)
        jc=all_orbit%jj(c)
        tzc = all_orbit%itzp(c)
        
        if ( tza /= tzc ) cycle
 
!        if ( tza /= -1 ) cycle 

        ang_fac = iph( (jc + 1)/2 + 1 ) * dsqrt(jc+1.d0) *  & 
             dsqrt( 2.*la+1.d0)*dsqrt( 2.*lc+1.d0) * tjs(2*la,2*lambda,2*lc,0,0,0) * sjs( 2*la, ja, 1, jc, 2*lc, 2)*0.5d0*real(tzc)
        

        vsum = 0.d0
        do i = 1, nn
           
           
           osc_a = hol(a,i)
           osc_c = hol(c,i)
           
           vsum = vsum + wrr(i)*rr(i)**2 * osc_a * rr(i)**(lambda) *osc_c
     
        end do
        
        reduced_r1(a,c) = ang_fac * vsum
        
        
     end do
  end do
  
  nn = 0
  pp = 0 
  DO i=1,below_ef

     if ( all_orbit%itzp(i) == -1 ) pp = pp + all_orbit%jj(i) + 1 
     if ( all_orbit%itzp(i) ==  1 ) nn = nn + all_orbit%jj(i) + 1  

  ENDDO
  
  if ( iam == 0 ) write(6,*) 'Number of neutrons', nn
  if ( iam == 0 ) write(6,*) 'Number of protons', pp
  

  do i = 1, all_orbit%total_orbits
     do j = 1, all_orbit%total_orbits
        tza = all_orbit%itzp(i)
        tzc = all_orbit%itzp(j)
        re_hcom = reduced_r1(i,j)
!        theta_1b(i,j) = re_hcom * delta(tza,-1) * delta(tza,tzc) 
        theta_1b(i,j) = re_hcom * delta(tza,tzc) 
     end do
  end do
  
  allocate( temp_kin(all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( temp_kin2(all_orbit%total_orbits, all_orbit%total_orbits) )
  
  
  temp_kin2 = 0.d0
  temp_kin2 = theta_1b
  temp_kin = 0.d0
  temp_kin = matmul( theta_1b,(u_hf) )
  theta_1b = 0.d0
  theta_1b = matmul( transpose(u_hf), temp_kin )
    
  deallocate( temp_kin, temp_kin2, reduced_r1 )
  
end subroutine setup_dipole_no_CMcorr_operator
!====================================

subroutine setup_CoulombMultipole_operator(lambda, curr, q_rel) 
  USE partial_waves
  USE relcm_gmatrix
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  USE constants
  use ang_mom_functions
  use one_body_operators
  
  implicit none
  !double precision :: reduced_r1(all_orbit%total_orbits, all_orbit%total_orbits) 
  integer :: iph, a, na, la,ja,tza,c, nc,lc,jc,tzc, i, nn, np
  double precision :: ang_fac
  double precision :: osc_a, osc_c, vsum,vsum1, osc_c_deriv, factor
  double precision, allocatable, dimension(:) :: rr ,wrr, sph_bessel
  double precision, allocatable :: hol(:,:), temp_kin(:,:), temp_kin2(:,:)
  DOUBLE PRECISION ::  sum_rel, cx(0:200), xp
  real*8 :: zlab, q_trans, delta, re_hcom, y,sj,sy,sjp,syp 
  integer :: ph, pp, j
  real*8, intent(in) :: q_rel 
  integer, intent(in) :: lambda 
  real*8 :: iso_t
  character(1), intent(in) :: curr
  
  ! convert momentum transfer in fm-1
  q_trans = q_rel / hbarc 
  nn = size(rmesh)
  
  allocate( reduced_r1(all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( sph_bessel(nn) )

  !
  ! compute the spherical bessel functions
  !
  sph_bessel = 0.d0 
  do i = 1, nn
     
     y = dble(q_trans*rmesh(i))
     call sphbes(lambda,y,sj,sy,sjp,syp)
     sph_bessel(i) = sj
     
   
  END DO
  
  reduced_r1 = 0.d0
  do a = 1, all_orbit%total_orbits
     
     na=all_orbit%nn(a) 
     la=all_orbit%ll(a)
     ja=all_orbit%jj(a)
     tza = all_orbit%itzp(a)
     
     do c = 1, all_orbit%total_orbits
        nc=all_orbit%nn(c)
        lc=all_orbit%ll(c)
        jc=all_orbit%jj(c)
        tzc = all_orbit%itzp(c)


        if(curr.eq.'V')  iso_t=real(tzc) ! isovector
        if(curr.eq.'S')  iso_t=1.d0 ! isoscalar
        
        if ( tza /= tzc ) cycle
        
        ang_fac =  iph( (jc + 1)/2 + lambda ) * dsqrt(jc+1.d0) * (2.*lambda+1.d0) *  & 
             dsqrt( 2.*la+1.d0)*dsqrt( 2.*lc+1.d0) * tjs(2*lc,2*lambda,2*la,0,0,0) * sjs( 2*lc, 1, jc, ja, 2*lambda, 2*la)* &
             0.5d0*iso_t
        
        vsum = 0.d0
        do i = 1, nn
           
           osc_a = hol_wfunc(a,i)
           osc_c = hol_wfunc(c,i)
           
           vsum = vsum + wrmesh(i)*rmesh(i)**2 * osc_a * sph_bessel(i) * osc_c
           
        end do

       reduced_r1(a,c) = ang_fac * vsum
       
     end do
  end do

  
 
  ! write reduced matrix elements to file
  theta_1b = 0.d0
  do i = 1, all_orbit%total_orbits
     do j = 1, all_orbit%total_orbits
        tza = all_orbit%itzp(i)
        tzc = all_orbit%itzp(j)
        re_hcom = reduced_r1(i,j)
        
        theta_1b(i,j) = re_hcom * delta(tza,tzc) 

!!$        if ( iam == 0 .and. abs(theta_1b(i,j) ) > 1.e-10 ) & 
!!$             write(333,'(1x,8(I3,1x),E)') & 
!!$             all_orbit%nn(i), all_orbit%ll(i), all_orbit%jj(i), all_orbit%itzp(i), & 
!!$             all_orbit%nn(j), all_orbit%ll(j), all_orbit%jj(j), all_orbit%itzp(j), & 
!!$             theta_1b(i,j)


     end do
  end do
  
  allocate( temp_kin(all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( temp_kin2(all_orbit%total_orbits, all_orbit%total_orbits) )
  
  !
  ! transform theta_1b to Hartree-Fock basis
  !
  
  temp_kin2 = 0.d0
  temp_kin2 = theta_1b
  temp_kin = 0.d0
  temp_kin = matmul( theta_1b,(u_hf) )
  theta_1b = 0.d0
  theta_1b = matmul( transpose(u_hf), temp_kin )
  
  
  deallocate( sph_bessel )
  deallocate( temp_kin, temp_kin2,reduced_r1 )
  
end subroutine setup_CoulombMultipole_operator


subroutine setup_gt_operator
  USE partial_waves
  USE relcm_gmatrix
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  USE constants
  use ang_mom_functions
  use one_body_operators
  
  implicit none
  !double precision :: reduced_r1(all_orbit%total_orbits, all_orbit%total_orbits) 
  integer :: iph, a, na, la,ja,tza,c, nc,lc,jc,tzc, i, nn
  double precision :: ang_fac
  double precision :: osc_a, osc_c, vsum,vsum1, osc_c_deriv, factor
  double precision, allocatable, dimension(:) :: rr ,wrr
  double precision, allocatable :: hol(:,:), temp_kin(:,:), temp_kin2(:,:)
  DOUBLE PRECISION ::  sum_rel, cx(0:200), xp
  real*8 :: zlab, delta, re_hcom, im_hcom
  integer :: ph, pp, j


 
  CALL commons_to_angmom  
  
  allocate( reduced_r1(all_orbit%total_orbits, all_orbit%total_orbits) )
  oscl=hbarc/SQRT(p_mass*hbar_omega)

  nn = 600
  allocate( rr(nn), wrr(nn) )
  call gauss_legendre(0.d0, 25.d0, rr, wrr, nn)
  allocate( hol(all_orbit%total_orbits, nn) )
  
  hol = 0.d0
  do a = 1, all_orbit%total_orbits
     
     na=all_orbit%nn(a) 
     la=all_orbit%ll(a)
     ja=all_orbit%jj(a)
     tza = all_orbit%itzp(a)

     ph = 1.d0 
     factor = 0.5D0*((na+la+2)*LOG(2.D0)+fac(na)-dfac(2*na+2*la+1)-0.5D0*LOG(pi))
     factor = EXP(factor)
     sum_rel=0.
     DO i=1,nn
        
        zlab= rr(i)/oscl
        
        CALL laguerre_general( na, la+0.5D0, zlab*zlab, cx )
        xp = cx(na)*exp(-zlab*zlab*0.5d0)*(zlab**la)
        hol(a,i) = xp*factor*(1.d0/oscl**(1.5D0))
        sum_rel = sum_rel + hol(a, i)**2 * rr(i)**2*wrr(i)
        
     end DO
     if ( iam == 0 ) write(6,*) 'norm', sum_rel

  end do
  
  reduced_r1 = 0.d0
  do a = 1, all_orbit%total_orbits
     
     na=all_orbit%nn(a) 
     la=all_orbit%ll(a)
     ja=all_orbit%jj(a)
     tza = all_orbit%itzp(a)
     

     do c = 1, all_orbit%total_orbits
        nc=all_orbit%nn(c)
        lc=all_orbit%ll(c)
        jc=all_orbit%jj(c)
        tzc = all_orbit%itzp(c)
        
        
        !
        ! check: la + lc + 1 = even
        !
        if ( na /= nc ) cycle 
        if ( la /= lc ) cycle 
        
        ang_fac = dsqrt(6.d0)*iph( (ja + 1)/2 + la + 1 )* dsqrt(jc+1.d0) * sjs(1, ja, 2*lc, jc, 1, 2)
        
        ! sqrt(2j+1) from the summing over projection in GT strength 
        reduced_r1(a,c) = ang_fac*dsqrt(3.d0)
        
     end do
  end do
  
  do i = 1, all_orbit%total_orbits
     do j = 1, all_orbit%total_orbits
        tza = all_orbit%itzp(i)
        tzc = all_orbit%itzp(j)
        re_hcom = reduced_r1(i,j)
        
        !write(6,*) dsqrt(all_orbit%jj(i)+1.d0) * reduced_r1(i,j), dsqrt(all_orbit%jj(j)+1.d0) * reduced_r1(j,i) 
        theta_1b(i,j) = re_hcom !*g_A/sqrt(6.*pi) 
        theta_1bdag(i,j) = im_hcom
                
        if ( abs(theta_1b(i,j)) > 1.e-10 .and. iam == 0 )  write(6,'(a,1x,4(I3,1x),g20.10)') 'ESPM for GT:', & 
             all_orbit%jj(i), all_orbit%ll(i), & 
             all_orbit%jj(j), all_orbit%ll(j), re_hcom * sqrt(all_orbit%jj(i)+1.)/sqrt(3.)
     end do
  end do
  !stop

!!$  
!!$  theta_1b = 0.d0
!!$  do i = 1, all_orbit%total_orbits
!!$     theta_1b(i,i) = 1.d0
!!$  end do
  
!!$  do a = 1, all_orbit%total_orbits
!!$     if ( all_orbit%itzp(a) /= 1 )cycle 
!!$     do c = 1, all_orbit%total_orbits
!!$        if ( all_orbit%itzp(c) /= -1 )cycle 
!!$        theta_1bdag(a,c) = theta_1bdag(c,a)*iph( (all_orbit%jj(a)-all_orbit%jj(c))/2 )*sqrt(all_orbit%jj(c)+1.)/sqrt(all_orbit%jj(a)+1.)
!!$     end do
!!$  end do
!!$  
!!$  theta_1b = theta_1b + theta_1bdag 
    
  allocate( temp_kin(all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( temp_kin2(all_orbit%total_orbits, all_orbit%total_orbits) )
  
  
  temp_kin2 = 0.d0
  temp_kin2 = theta_1b
  temp_kin = 0.d0
  temp_kin = matmul( theta_1b,(u_hf) )
  theta_1b = 0.d0
  theta_1b = matmul( transpose(u_hf), temp_kin )
  
  
!  do i = 1, all_orbit%total_orbits
!     do j = 1, all_orbit%total_orbits
        
        !if ( abs(theta_1b(i,j) ) < 1.e-10 ) cycle 
        !write(6,*) all_orbit%jj(i), all_orbit%ll(i), all_orbit%jj(j), all_orbit%ll(j), theta_1b(i,j)
!     end do
!  end do
  
  deallocate( temp_kin, temp_kin2 )
  deallocate( rr, wrr )
  deallocate( hol, reduced_r1 )

end subroutine setup_gt_operator


subroutine setup_quadrupole_operator
  USE partial_waves
  USE relcm_gmatrix
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  USE constants
  use ang_mom_functions
  use one_body_operators
  
  implicit none
  integer :: lambda
  !double precision :: reduced_r1(all_orbit%total_orbits, all_orbit%total_orbits) 
  integer :: iph, a, na, la,ja,tza,c, nc,lc,jc,tzc, i, nn
  double precision :: ang_fac
  double precision :: osc_a, osc_c, vsum,vsum1, osc_c_deriv, factor
  double precision, allocatable, dimension(:) :: rr ,wrr
  double precision, allocatable :: hol(:,:), temp_kin(:,:), temp_kin2(:,:)
  DOUBLE PRECISION ::  sum_rel, cx(0:200), xp
  real*8 :: zlab, delta, re_hcom, eff_ch
  integer :: ph, pp, j

  lambda = 2
  CALL commons_to_angmom  
  
  allocate( reduced_r1(all_orbit%total_orbits, all_orbit%total_orbits) )
  oscl=hbarc/SQRT(p_mass*hbar_omega)
  write(6,*) 'Oscillator length', oscl
  
  nn = 600
  allocate( rr(nn), wrr(nn) )
  call gauss_legendre(0.d0, 25.d0, rr, wrr, nn)
  allocate( hol(all_orbit%total_orbits, nn) )
  
  hol = 0.d0
  do a = 1, all_orbit%total_orbits
     
     na=all_orbit%nn(a) 
     la=all_orbit%ll(a)
     ja=all_orbit%jj(a)
     tza = all_orbit%itzp(a)

     ph = 1.d0 
     factor = 0.5D0*((na+la+2)*LOG(2.D0)+fac(na)-dfac(2*na+2*la+1)-0.5D0*LOG(pi))
     factor = EXP(factor)
     sum_rel=0.
     DO i=1,nn
        
        zlab= rr(i)/oscl
        
        CALL laguerre_general( na, la+0.5D0, zlab*zlab, cx )
        xp = cx(na)*exp(-zlab*zlab*0.5d0)*(zlab**la)
        hol(a,i) = xp*factor*(1.d0/oscl**(1.5D0))
        sum_rel = sum_rel + hol(a, i)**2 * rr(i)**2*wrr(i)
        
     end DO
     if ( iam == 0 ) write(6,*) 'norm', sum_rel

  end do
  
  reduced_r1 = 0.d0
  do a = 1, all_orbit%total_orbits

     na=all_orbit%nn(a)
     la=all_orbit%ll(a)
     ja=all_orbit%jj(a)
     tza = all_orbit%itzp(a)

     do c = 1, all_orbit%total_orbits
        nc=all_orbit%nn(c)
        lc=all_orbit%ll(c)
        jc=all_orbit%jj(c)
        tzc = all_orbit%itzp(c)
        

        ! 
        ! check: la + lc + 1 = even  
        !                  
        if ( tza /= tzc ) cycle
        
!!$        ang_fac = iph( (jc + 1)/2 + la + lc ) * dsqrt(jc+1.d0) * dsqrt( 2.*lambda + 1.d0) * &
!!$             dsqrt( 2.*la+1.d0)*dsqrt( 2.*lc+1.d0) * tjs(2*lc,2*lambda,2*la,0,0,0) * sjs( 2*lc, 1, jc, ja, 2*lambda, 2*la)
        
        ang_fac = 0.5d0*iph( (jc - 1)/2 + lambda) * ( 1.d0 + iph(la+lc+lambda) ) * dsqrt(jc+1.d0) * dsqrt( 2.*lambda + 1.d0) * &
             tjs(ja,jc,2*lambda,1,-1,0)
        
        vsum = 0.d0
        do i = 1, nn
           
           osc_a = hol(a,i)
           osc_c = hol(c,i)

           vsum = vsum + wrr(i)*rr(i)**2 * osc_a * rr(i)**(lambda) * osc_c

        end do

        reduced_r1(a,c) = ang_fac * vsum/dsqrt(4.d0*pi)
     end do
  end do
  
  eff_ch = 0.d0
  do i = 1, all_orbit%total_orbits
     do j = 1, all_orbit%total_orbits
        tza = all_orbit%itzp(i)
        tzc = all_orbit%itzp(j)
        re_hcom = reduced_r1(i,j)
        if ( lambda == 2 ) theta_1b(i,j) = re_hcom * ((1.d0-eff_ch)*delta(tza,-1) + eff_ch*delta(tza,1))
        !if ( lambda == 2 ) theta_1b(i,j) = re_hcom * delta(tza,-1) 
     end do
  end do
  
  
  allocate( temp_kin(all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( temp_kin2(all_orbit%total_orbits, all_orbit%total_orbits) )
  
  
  temp_kin2 = 0.d0
  temp_kin2 = theta_1b
  temp_kin = 0.d0
  temp_kin = matmul( theta_1b,(u_hf) )
  theta_1b = 0.d0
  theta_1b = matmul( transpose(u_hf), temp_kin )
  
  
!!$  do i = 1, all_orbit%total_orbits
!!$     do j = 1, all_orbit%total_orbits
!!$        
!!$        if ( abs(theta_1b(i,j) ) < 1.e-10 ) cycle 
!!$        write(6,'(a,1x,4(I3,1x),g20.10)') 'check dipole', & 
!!$             all_orbit%jj(i), all_orbit%ll(i), all_orbit%jj(j), all_orbit%ll(j), real(theta_1b(i,j))
!!$     end do
!!$  end do
!!$  
  
  deallocate( temp_kin, temp_kin2,reduced_r1 )
  
end subroutine setup_quadrupole_operator

subroutine setup_M1_operator
  USE partial_waves
  USE relcm_gmatrix
  USE single_particle_orbits
  USE wave_functions
  USE configurations
  USE constants
  use ang_mom_functions
  use one_body_operators
  
  implicit none
  !double precision :: reduced_r1(all_orbit%total_orbits, all_orbit%total_orbits) 
  integer :: iph, a, na, la,ja,tza,c, nc,lc,jc,tzc, i, nn
  double precision :: ang_fac
  double precision :: osc_a, osc_c, vsum,vsum1, osc_c_deriv, factor
  double precision, allocatable, dimension(:) :: rr ,wrr
  double precision, allocatable :: hol(:,:), temp_kin(:,:), temp_kin2(:,:)
  DOUBLE PRECISION ::  sum_rel, cx(0:200), xp
  real*8 :: g_p, g_n, zlab, delta, re_hcom, im_hcom, sum1
  integer :: ph, pp, j, lambda
  
  lambda = 1

  allocate( reduced_r1(all_orbit%total_orbits, all_orbit%total_orbits) )
  
  g_p =  5.586d0
  g_n = -3.826d0
  
  nn = 400
  call setup_hol_wfunc(nn)
  CALL commons_to_angmom  
  
!!$  oscl=hbarc/SQRT(p_mass*hbar_omega)
!!$  nn = 600
!!$  allocate( rr(nn), wrr(nn) )
!!$  call gauss_legendre(0.d0, 25.d0, rr, wrr, nn)
!!$  allocate( hol(all_orbit%total_orbits, nn) )
!!$  
!!$  hol = 0.d0
!!$  do a = 1, all_orbit%total_orbits
!!$     
!!$     na=all_orbit%nn(a) 
!!$     la=all_orbit%ll(a)
!!$     ja=all_orbit%jj(a)
!!$     tza = all_orbit%itzp(a)
!!$
!!$     ph = 1.d0 
!!$     factor = 0.5D0*((na+la+2)*LOG(2.D0)+fac(na)-dfac(2*na+2*la+1)-0.5D0*LOG(pi))
!!$     factor = EXP(factor)
!!$     sum_rel=0.
!!$     DO i=1,nn
!!$        
!!$        zlab= rr(i)/oscl
!!$        
!!$        CALL laguerre_general( na, la+0.5D0, zlab*zlab, cx )
!!$        xp = cx(na)*exp(-zlab*zlab*0.5d0)*(zlab**la)
!!$        hol(a,i) = xp*factor*(1.d0/oscl**(1.5D0))
!!$        sum_rel = sum_rel + hol(a, i)**2 * rr(i)**2*wrr(i)
!!$        
!!$     end DO
!!$     if ( iam == 0 ) write(6,*) 'norm', sum_rel
!!$
!!$  end do
  
      
  reduced_r1 = 0.d0

  do a = 1, all_orbit%total_orbits
     
     na=all_orbit%nn(a) 
     la=all_orbit%ll(a)
     ja=all_orbit%jj(a)
     tza = all_orbit%itzp(a)
     

     do c = 1, all_orbit%total_orbits
        nc=all_orbit%nn(c)
        lc=all_orbit%ll(c)
        jc=all_orbit%jj(c)
        tzc = all_orbit%itzp(c)
        

        if ( na /= nc ) cycle 
        if ( la /= lc ) cycle        
        if ( tza /= tzc ) cycle 
 
        ang_fac = 0.d0
        if     ( tza == -1) then     ! proton
           ang_fac = delta(tza,-1)* dsqrt(jc+1.d0)* &
                (iph( (jc + 1)/2 + lc + 1 ) * dsqrt(2.*lc+1.d0)*dsqrt(lc*(lc+1.d0)) * sjs(2*lc, 1, jc, ja, 2, 2*lc) + & 
                g_p*dsqrt(1.5d0)*iph( (ja + 1)/2 + lc + 1 ) * sjs(1, 2*lc, jc, ja, 2, 1))
        elseif ( tza ==  1) then     ! neutron
           ang_fac = delta(tza,1)* dsqrt(jc+1.d0) * iph( (ja + 1)/2 + lc + 1 ) * g_n * dsqrt(1.5d0)     * sjs(1, 2*lc, jc, ja, 2, 1) 
        endif
        

!!$        
!!$        if     ( tza == -1) then     ! proton
!!$           ang_fac = dsqrt(jc+1.d0)* &
!!$                ((1-2*mod((jc+2*lc)/2,2)) * dsqrt(2.*lc+1.d0)*dsqrt(lc*(lc+1.d0)) * sjs(2*lc, 1, jc, ja, 2, 2*lc) + & 
!!$                g_p*dsqrt(1.5d0)*(1-2*mod((ja+2*lc)/2,2))* sjs(1, 2*lc, jc, ja, 2, 1))
!!$        elseif ( tza == +1) then     ! neutron
!!$           ang_fac = dsqrt(jc+1.d0) *(1-2*mod((ja+2*lc)/2,2)) * g_n * dsqrt(1.5d0)   * sjs(1, 2*lc, jc, ja, 2, 1) 
!!$        endif
        
!!$        vsum = 0.d0
!!$        do i = 1, nn
!!$           
!!$           osc_a = hol_wfunc(a,i)
!!$           osc_c = hol_wfunc(c,i)
!!$           
!!$           vsum = vsum + wrmesh(i)*rmesh(i)**2 * osc_a * osc_c
!!$           
!!$        end do
        reduced_r1(a,c) = ang_fac 
           
        !if ( abs(reduced_r1(a,c)) > 1.e-10 ) write(66,*) a,c, reduced_r1(a,c), vsum
        
     end do
  end do
  
  
  do i = 1, all_orbit%total_orbits
     do j = 1, all_orbit%total_orbits
        tza = all_orbit%itzp(i)
        tzc = all_orbit%itzp(j)
        re_hcom = reduced_r1(i,j)
        theta_1b(i,j) = re_hcom
     end do
  end do
  
  
  
  allocate( temp_kin(all_orbit%total_orbits, all_orbit%total_orbits) )
  allocate( temp_kin2(all_orbit%total_orbits, all_orbit%total_orbits) )
  

  
  temp_kin2 = 0.d0
  temp_kin2 = theta_1b
  temp_kin = 0.d0
  temp_kin = matmul( theta_1b,(u_hf) )
  theta_1b = 0.d0
  theta_1b = matmul( transpose(u_hf), temp_kin )
  
  
  deallocate( temp_kin, temp_kin2, reduced_r1 )
  
end subroutine setup_M1_operator

!SONIA
!OTHER QUADRATURE
      SUBROUTINE gauleg2(x1,x2,x,w,n,np)
      INTEGER n
      real*8 x1,x2,x(np),w(np)
      real*8 EPS
      PARAMETER (EPS=1.d-14)
      INTEGER i,j,m
      real*8 p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software .)2">`!.


module setup_orbitals_imsrg

  integer :: lab_lmax, lab_nmax, nlmax

end module setup_orbitals_imsrg
!
!     Reads in and allocates sp data for the nuclear case
!     Read comments below in order to properly understand the input
! 
SUBROUTINE setup_sp_data_imsrg
  USE single_particle_orbits  
  USE constants
  USE setup_orbitals_imsrg
  use cc_config_module
  use parallel 
  
  IMPLICIT NONE
  INTEGER :: i, n1, l1, j1, t1, shell1, nnn, nn1, ll1, lll, j, l12, max_space
  CHARACTER(LEN= 10) space, model
  INTEGER :: model_space, states, neutron_hole, proton_hole 
  INTEGER, ALLOCATABLE, DIMENSION(:) ::  n_big,l_big,j_big,tzp_big, shell_big
  CHARACTER(LEN= 10), ALLOCATABLE, DIMENSION(:) :: space_big, model_big
  INTEGER :: a,na,la,ja,tza, c,nc,lc,jc,tzc, neutron_tot_orbs, proton_tot_orbs
  INTEGER :: pp, nn, n , l_ang, j_ang, itz, nmax_proton, np
  REAL(DP) :: ekin

  ! 2n+ l 
  max_space = config%twobody_nmax
  lab_lmax = max_space
  lab_nmax = max_space/2
  nlmax = max_space
  
  proton_hole = 1 
  neutron_hole = 1 
  CALL count_single_particle_orbits(states)
  neutron_data%total_orbits = states
  proton_data%total_orbits = states
  ALLOCATE( n_big(2*states), l_big(2*states), j_big(2*states), tzp_big(2*states), shell_big(2*states))
  ALLOCATE( space_big(2*states), model_big(2*states))
  !     Setup all possible orbit information
  all_orbit%total_orbits=neutron_data%total_orbits+ proton_data%total_orbits
  !     Allocate space for all single-particle data
  CALL allocate_sp_array(neutron_data,neutron_data%total_orbits) 
  CALL allocate_sp_array(proton_data,proton_data%total_orbits) 
  CALL allocate_sp_array(all_orbit,all_orbit%total_orbits) 
  !     Read neutron single-particle data
  CALL setup_orbitals(model_space, states, &
       n_big,l_big,j_big,tzp_big, shell_big, space_big, model_big,&
       proton_hole,neutron_hole)
  DO i=1, neutron_data%total_orbits 
     neutron_data%nn(i)=n_big(2*i)
     neutron_data%ll(i)=l_big(2*i); neutron_data%jj(i)=j_big(2*i)
     neutron_data%itzp(i)=tzp_big(2*i)
     neutron_data%nshell(i)=shell_big(2*i)
     neutron_data%orbit_status(i)=space_big(2*i) 
     neutron_data%model_space(i)=model_big(2*i)
     neutron_data%e(i)=hbar_omega*(shell_big(2*i)+1.5)
     !     Neutrons are in the internal structure always even numbers 
     all_orbit%nn(i*2)=neutron_data%nn(i)
     all_orbit%ll(i*2)=neutron_data%ll(i)
     all_orbit%jj(i*2)=neutron_data%jj(i)
     all_orbit%nshell(i*2)=neutron_data%nshell(i)
     all_orbit%e(i*2)=neutron_data%e(i)
     all_orbit%itzp(i*2)=neutron_data%itzp(i)
     all_orbit%orbit_status(i*2)=neutron_data%orbit_status(i)
     all_orbit%model_space(i*2)=neutron_data%model_space(i)
  ENDDO
  DO i=1,  proton_data%total_orbits 
     proton_data%nn(i)=n_big(2*i-1)
     proton_data%ll(i)=l_big(2*i-1)
     proton_data%jj(i)=j_big(2*i-1)
     proton_data%itzp(i)=tzp_big(2*i-1)
     proton_data%nshell(i)=shell_big(2*i-1)
     proton_data%orbit_status(i)=space_big(2*i-1) 
     proton_data%model_space(i)=model_big(2*i-1)
     proton_data%e(i)=hbar_omega*(shell_big(2*i-1)+1.5)
     !     protons are in the internal structure always odd numbers
     all_orbit%nn(i*2-1)=proton_data%nn(i)
     all_orbit%ll(i*2-1)=proton_data%ll(i)
     all_orbit%jj(i*2-1)=proton_data%jj(i)
     all_orbit%nshell(i*2-1)=proton_data%nshell(i)
     all_orbit%e(i*2-1)=proton_data%e(i)
     all_orbit%itzp(i*2-1)=proton_data%itzp(i)
     all_orbit%orbit_status(i*2-1)=proton_data%orbit_status(i)
     all_orbit%model_space(i*2-1)=proton_data%model_space(i)
  ENDDO

  
  
  
  DO i=1, all_orbit%total_orbits
     if ( iam == 0 ) WRITE(6,'(7HNumber:,6(I3,2X),2X,E12.6)') i, all_orbit%nn(i), all_orbit%ll(i), &
          all_orbit%jj(i), &
          all_orbit%nshell(i), all_orbit%itzp(i), all_orbit%e(i) 
  ENDDO
  DEALLOCATE( space_big, model_big); DEALLOCATE( n_big, l_big, j_big, tzp_big, shell_big)
  
  
  nlmax=0 ! find max possible 2*n+l for large Hilbert space and model space
  DO j=1, all_orbit%total_orbits
     ll1=all_orbit%ll(j)
     nn1=all_orbit%nn(j)
     DO i=j,all_orbit%total_orbits
        lll=all_orbit%ll(i)
        nnn=all_orbit%nn(i)
        l12=ll1+lll+nn1+nn1+nnn+nnn
        IF ( l12 > nlmax) nlmax = l12
      
     ENDDO
  ENDDO
  if ( iam == 0 ) WRITE(6,*) ' Max value of 2*n + l+ cm 2*N +L for large space:', nlmax
 
  


END SUBROUTINE setup_sp_data_imsrg
!
!    setup of single-particle basis
!
SUBROUTINE setup_orbitals(model_space, norb, nn, ll, jj, itzp, nshell, space, model, &
     proton_hole,neutron_hole)
  USE single_particle_orbits
  USE constants
  use setup_orbitals_imsrg
  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: norb, model_space, proton_hole, neutron_hole 
  INTEGER, DIMENSION(norb) :: na,la,ja
  INTEGER, DIMENSION(2*norb), INTENT(INOUT) ::  nn,ll,jj,itzp, nshell
  CHARACTER(LEN= 10), DIMENSION(2*norb), INTENT(INOUT)::   space, model
  INTEGER :: icount, j_max, j_min, j, n, l
  
  
  icount = 0
  DO n = 0, lab_nmax
     DO l = 0, lab_lmax
        IF ( n+n+l > nlmax) CYCLE
        
        j_min = l+l-1; j_max = l+l+1
        IF (j_min < 0) j_min = 1
        DO j = j_max, j_min, -2
           icount = icount + 1
           ja(icount)=j
           la(icount)=l
           na(icount)=n
        ENDDO
     ENDDO
  ENDDO
  IF ( icount /= norb) THEN
     WRITE(6,*) 'Error in allocation of single-particle states'
     STOP 'setup_single_particle_orbits'
  ENDIF
  ! now sort the sp orbits with increasing 2n+l and largest j value first
  CALL sort_spbasis_imsrg(icount,na,la,ja)
  DO j = 1, icount
     
     jj(j*2)=ja(j)                     ! neutrons are
     ll(j*2)=la(j)                     ! always even numbers
     nn(j*2)=na(j)                     ! protons are odd
     nshell(j*2)=2*na(j)+la(j)
     
!!$     IF ( nshell(j*2) <= model_space-1) THEN
!!$        model(j*2)='inside'
!!$     ELSE
!!$        model(j*2)= 'inside' !'outside'
!!$     ENDIF
!!$     IF ( j <= neutron_hole) THEN
!!$        space(j*2)='hole'
!!$     ELSEIF ( (j > neutron_hole).AND.( nshell(j*2) <= model_space-1)) THEN
!!$        space(j*2)='particle'
!!$     ELSE
!!$        space(j*2)='passive'
!!$     ENDIF
     jj(j*2-1)=ja(j)                    
     ll(j*2-1)=la(j)
     nn(j*2-1)=na(j)
     nshell(j*2-1)=2*na(j)+la(j)
!!$     IF ( nshell(j*2-1) <= model_space-1) THEN
!!$        model(j*2-1)='inside'
!!$     ELSE
!!$        model(j*2-1)='inside' !'outside'
!!$     ENDIF
!!$     IF ( j <= proton_hole) THEN
!!$        space(j*2-1)='hole'
!!$     ELSEIF ( (j > proton_hole).AND.( nshell(j*2-1) <= model_space-1)) THEN
!!$        space(j*2-1)='particle'
!!$     ELSE
!!$        space(j*2-1)='passive'
!!$     ENDIF
     itzp(j*2-1)=-1                   ! protons, nucl. phys. def.
     itzp(j*2)=1                      ! Neutrons, nucl. phys. def.
  ENDDO

END SUBROUTINE setup_orbitals
!
!   Find maximum number of single-particle states constrained by n, l and j
!
SUBROUTINE count_single_particle_orbits(icount)
  USE single_particle_orbits
  USE constants
  USE setup_orbitals_imsrg
  
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: icount
  INTEGER :: n,l, j, j_max, j_min

  icount = 0
  DO n = 0, lab_nmax
     DO l = 0, lab_lmax
        IF ( n+n+l > nlmax) CYCLE
        j_min = l+l-1; j_max = l+l+1
        IF (j_min < 0) j_min = 1
        DO j = j_max, j_min, -2
           icount = icount + 1
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE count_single_particle_orbits

