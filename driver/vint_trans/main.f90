program main
 use  constants
 use  REORDER
 use single_particle_orbits
 use gmat_storage
 use transe
 implicit none

 open(unit=10,file='input.dat')
 read(10,*)all_orbit%total_orbits
 read(10,*)sp_file
 read(10,*)kin_file
 read(10,*)vint_file
 read(10,*)snt_file
 read(10,*)amass

     CALL read_nuclear_sp_data
     call read_tbme
     call generate_int


end program main



SUBROUTINE read_nuclear_sp_data
  USE single_particle_orbits
  USE constants
  use reorder 
  use gmat_storage
  use transe
  
  IMPLICIT NONE
  INTEGER :: i, a,b, n1, l1,j1, t1
  integer:: aa,bb
  real*8:: gmat
  
  
  OPEN(1, FILE = sp_file) 
  !read(5,*);read(5,*) input2
  OPEN(2, FILE = kin_file)
  
   write(6,*) 'sp-energy file', sp_file
   write(6,*) 'kinetic-energy file', kin_file
  
  CALL allocate_sp_array(all_orbit,all_orbit%total_orbits)
  allocate(indx_inv(all_orbit%total_orbits))
  indx_inv=0
  
  
  do a = 1, all_orbit%total_orbits 
     
     read(1,*) n1, l1, j1, t1
     
     
           all_orbit%orbit_status(a) = 'particle'
           all_orbit%model_space(a) = 'inside'
           all_orbit%nn(a)=n1
           all_orbit%ll(a)=l1
           all_orbit%jj(a)=j1
           all_orbit%itzp(a)=t1
       
  end do

  a=0
  b=0
  do i=1, all_orbit%total_orbits
    if(all_orbit%itzp(i)==-1)then
            a=a+1
    else
            b=b+1
    endif
  end do
  proton_data%total_orbits=a
  neutron_data%total_orbits=b

  
  a=0
  b=proton_data%total_orbits
  do i=1, all_orbit%total_orbits
    if(all_orbit%itzp(i)==-1)then
            a=a+1
            indx_inv(i)=a
    else
            b=b+1
            indx_inv(i)=b
    endif
     
  end do
  

  do i = 1, all_orbit%total_orbits
     write(6,*) i, indx_inv(i)
  end do
  
  allocate( tkin( all_orbit%total_orbits, all_orbit%total_orbits) )
  tkin=0.d0
  
  do a=1,all_orbit%total_orbits
  do b=1,all_orbit%total_orbits
  read(2,*)aa,bb,gmat
  tkin(aa,bb)=gmat*(1.-1./amass)
  enddo
  enddo
  
END SUBROUTINE read_nuclear_sp_data

subroutine read_tbme
 use  constants
 use  REORDER
 use single_particle_orbits
 use gmat_storage
 use configurations
 use transe
 implicit none
  TYPE (configuration_descriptor) :: gmat_ket_configs 
  integer:: number_channels,channel
  integer:: a,b,c,d,ket_dim
  real*8:: gmat,p2,hc
  integer:: ang_mom,p_parity,isospin_z
  integer:: aa,bb,cc,dd,bra,ket
  real:: phase_ab,phase_cd
  integer:: iph




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
                (ang_mom,p_parity,isospin_z,4,gmat_ket_configs)
           IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
           
           number_channels = number_channels + 1 
           

        end DO
     end DO
  end DO
allocate(quantum_numbers_ch3(number_channels,3))
quantum_numbers_ch3=0
allocate(gmatrix(number_channels))
allocate(lookup_pq_configs(number_channels))
allocate(pq_configs(number_channels))

  allocate( locate_gmatchannel(0:7,-1:1, 0:1, 0:60) ) 
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
                (ang_mom,p_parity,isospin_z,4,gmat_ket_configs)
           IF (gmat_ket_configs%number_confs <= 0 ) CYCLE
           ket_dim=gmat_ket_configs%number_confs
           allocate(gmat_ket_configs%config_ab(2*ket_dim))
           call setup_gmatrix_configurations&
                (ang_mom,p_parity,isospin_z,4,gmat_ket_configs)
           
           number_channels = number_channels + 1 
           quantum_numbers_ch3(number_channels,1) = ang_mom
           quantum_numbers_ch3(number_channels,2) = p_parity
           quantum_numbers_ch3(number_channels,3) = isospin_z
           allocate(gmatrix(number_channels)%val(ket_dim,ket_dim))
           gmatrix(number_channels)%val=0.d0
           allocate(pq_configs(number_channels)%ival(all_orbit%total_orbits,all_orbit%total_orbits))
           allocate(lookup_pq_configs(number_channels)%ival(2,ket_dim))
           pq_configs(number_channels)%ival=0
           lookup_pq_configs(number_channels)%ival=0


           do ket=1,ket_dim
           a=gmat_ket_configs%config_ab(2*ket-1)
           b=gmat_ket_configs%config_ab(2*ket)
           pq_configs(number_channels)%ival(a,b)=ket
           lookup_pq_configs(number_channels)%ival(1,ket)=a
           lookup_pq_configs(number_channels)%ival(2,ket)=b
           enddo
           locate_gmatchannel(0,isospin_z, p_parity, ang_mom) = number_channels
           deallocate(gmat_ket_configs%config_ab)
        end DO
     end DO
  end DO

  open(unit=3,file=vint_file)

5 READ(3,*,end=6) isospin_z,p_parity,ang_mom,aa,bb,cc,dd,gmat,p2,hc
   
           channel=locate_gmatchannel(0,isospin_z, p_parity, ang_mom/2)
           bra=pq_configs(channel)%ival(aa,bb)
           ket=pq_configs(channel)%ival(cc,dd)
           if(bra*ket==0)write(*,*)isospin_z,p_parity,ang_mom,aa,bb,cc,dd,gmat,p2,hc

           gmatrix(channel)%val(bra,ket)=gmat-p2/dble(amass)
  
  GOTO 5
6 CONTINUE
  REWIND 1

  number_channels=size(gmatrix)
  do channel=1,number_channels
           ang_mom   = quantum_numbers_ch3(number_channels,1) 
           p_parity  = quantum_numbers_ch3(number_channels,2) 
           isospin_z = quantum_numbers_ch3(number_channels,3) 
    do bra=1,size(lookup_pq_configs(channel)%ival,2)
      a=lookup_pq_configs(channel)%ival(1,bra)
      b=lookup_pq_configs(channel)%ival(2,bra)
      phase_ab=iph((all_orbit%jj(a)+all_orbit%jj(b))/2+ang_mom+1)
    do ket=1,size(lookup_pq_configs(channel)%ival,2)
      c=lookup_pq_configs(channel)%ival(1,ket)
      d=lookup_pq_configs(channel)%ival(2,ket)
      phase_cd=iph((all_orbit%jj(c)+all_orbit%jj(d))/2+ang_mom+1)
      if(bra==ket)cycle
      gmat=0.d0
      if(b==d .and. a/=c)then
              gmat= gmat+tkin(a,c)
      endif
      if(a==c .and. b/=d)then
              gmat= gmat+tkin(b,d)
      endif
      if(b==c .and. a/=d)then
              gmat= gmat+tkin(a,d)*phase_cd
      endif
      if(a==d .and. b/=c)then
              gmat= gmat+tkin(b,c)*phase_ab
      endif
       gmatrix(channel)%val(bra,ket)=gmatrix(channel)%val(bra,ket)+gmat


  enddo

  enddo
  enddo
end subroutine read_tbme



subroutine generate_int
 use  constants
 use  REORDER
 use single_particle_orbits
 use gmat_storage
 use configurations
 use transe
 implicit none

  integer:: channel
  integer:: a,b,c,d,ket_dim,i
  integer:: aa,bb,cc,dd,bra,ket
  integer:: ang_mom
  open(unit=20,file=snt_file)
  write(20,*)'! shell model interaction from this following files:'
  write(20,*) '!', sp_file
  write(20,*) '!', kin_file
  write(20,*) '!', vint_file
  write(20,*) '! model space'
  write(20,'(2(I4,1x))')proton_data%total_orbits,neutron_data%total_orbits 
  do i=1,all_orbit%total_orbits
      if(all_orbit%itzp(i)==1)cycle
      write(20,'(6(I4,1x))')i,indx_inv(i),&
              all_orbit%nn(i),&
              all_orbit%ll(i),&
              all_orbit%jj(i),&
              all_orbit%itzp(i)
      end do
  do i=1,all_orbit%total_orbits
      if(all_orbit%itzp(i)==-1)cycle
      write(20,'(6(I4,1x))')i,indx_inv(i),&
              all_orbit%nn(i),&
              all_orbit%ll(i),&
              all_orbit%jj(i),&
              all_orbit%itzp(i)
      end do


  write(20,*) '! '
  write(20,*) '! '
  write(20,*) '! '
  write(20,*)all_orbit%total_orbits,'0'
  do i=1,all_orbit%total_orbits
  write(20,'(2(I4,1x),F13.6)')indx_inv(i),indx_inv(i),tkin(i,i)
  enddo

  write(20,*) '! '

  write(20,*) 'number 1 1 0'
  do channel=1,size(gmatrix)
  ang_mom=quantum_numbers_ch3(channel,1)
  ket_dim=size(gmatrix(channel)%val,1)
  do bra=1,ket_dim
  a=lookup_pq_configs(channel)%ival(1,bra)
  b=lookup_pq_configs(channel)%ival(2,bra)
  aa=indx_inv(a)
  bb=indx_inv(b)
  do ket=1,ket_dim
  c=lookup_pq_configs(channel)%ival(1,ket)
  d=lookup_pq_configs(channel)%ival(2,ket)
  cc=indx_inv(c)
  dd=indx_inv(d)
  write(20,'(5(I4,1x),F13.6)')aa,bb,cc,dd,ang_mom,gmatrix(channel)%val(bra,ket)
  enddo
  enddo
  enddo


end subroutine generate_int
