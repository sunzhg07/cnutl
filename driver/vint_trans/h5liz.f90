program main
  USE single_particle_orbits
  USE configurations
  USE constants
  use gmat_storage
  implicit none
  integer:: a
  CHARACTER(LEN=100) :: vint_name,vout_name
  CALL read_nuclear_sp_data
  call generate_configurations
  do a=1,16
  read(1,*)vint_name
  write(vout_name,'(a,a)')trim(vint_name),".h5"
  write(*,*)vint_name,vout_name
  call read_twobody_matrix(vint_name)
  call write_twobody_matrix(vout_name)
  enddo

end program main


SUBROUTINE read_nuclear_sp_data
  USE single_particle_orbits
  USE constants
  use gmat_storage
  use reorder
  
  IMPLICIT NONE
  CHARACTER(LEN=100) :: input1, input2
  REAL*8 :: re_e1, im_e1
  INTEGER :: i,j, n1, l1, j1, t1
  integer:: io_status
  
  
  OPEN(unit=1, FILE = 'input.dat') 
  read(1,*)input1
  read(1,*)input2
  
  write(6,*) 'sp-energy file', input1
  write(6,*) 'kinetic-energy file', input2
  open(unit=2,file=input1)

  read(2,*)hbar_omega
  j=0
  do
   read(2,*,IOSTAT=io_status)i,n1,l1,j1,t1,re_e1,im_e1
   if(io_status>0)then
           write(*,*)'error in read single particle file'
   elseif(io_status<0)then
           !! EOF
           exit
    else
           j=j+1

   end if
  enddo

  all_orbit%total_orbits=j
  CALL allocate_sp_array(all_orbit,all_orbit%total_orbits)
  rewind(2)
  read(2,*)hbar_omega
  do i=1,all_orbit%total_orbits
   read(2,*,IOSTAT=io_status)j,n1,l1,j1,t1,re_e1,im_e1
   all_orbit%nn(i)=n1
   all_orbit%ll(i)=l1
   all_orbit%jj(i)=j1
   all_orbit%itzp(i)=t1
   all_orbit%nshell(i)=l1+2*n1
   all_orbit%e(i)=re_e1
   all_orbit%orbit_status='particle'
   all_orbit%model_space='inside'
  enddo

end subroutine read_nuclear_sp_data

subroutine generate_configurations
  USE single_particle_orbits
  USE configurations
  USE constants
  use gmat_storage
  implicit none
  TYPE (configuration_descriptor) :: bra_configs
  integer:: isospin_z, p_parity, ang_mom
  integer:: ndim,ndim_tot
  integer:: number_channels
  integer:: isob(3),tt



  isob(1)=-1
  isob(2)=0
  isob(3)=1


!@!  number_channels = 0
!@!  !     loop over isospin projection
!@!  do tt=1,3
!@!   isospin_z=isob(tt)
!@!     !     loop over parity values, here positive parity is 0, negative 1
!@!     DO p_parity=0,1           
!@!        !     loop over angular momenta
!@!        DO ang_mom=0, 60 
!@!
!@!           
!@!          
!@!           CALL  number_gmatrix_confs&
!@!                (ang_mom,p_parity,isospin_z,4,bra_configs)
!@!           IF (bra_configs%number_confs <= 0 ) CYCLE
!@!           
!@!           number_channels = number_channels + 1 
!@!        end DO
!@!     end DO
!@!  end DO
!@!
!@!  allocate(quantum_numbers(number_channels,3))
!@!  allocate(vpot1(number_channels))
!@!  allocate(vpot2(number_channels))
!@!  allocate(vpot3(number_channels))
ndim_tot=0

 number_channels=0
 do tt=1,3
  isospin_z=isob(tt)
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60 

           
          
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,4,bra_configs)
           IF (bra_configs%number_confs <= 0 ) CYCLE
           ndim=bra_configs%number_confs
           
           number_channels = number_channels + 1 
           !quantum_numbers(number_channels,1) = ang_mom
           !quantum_numbers(number_channels,2) = p_parity
           !quantum_numbers(number_channels,3) = isospin_z
           ndim_tot=ndim_tot+ndim*(ndim+1)/2
!@!           allocate(vpot1(number_channels)%val(ndim,ndim))
!@!           vpot1(number_channels)%val=0.d0
!@!           allocate(vpot2(number_channels)%val(ndim,ndim))
!@!           vpot2(number_channels)%val=0.d0
!@!           allocate(vpot3(number_channels)%val(ndim,ndim))
!@!           vpot3(number_channels)%val=0.d0
        end DO
     end DO
  end DO
  allocate(vpot1(ndim_tot));vpot1=0.d0
  allocate(vpot2(ndim_tot));vpot2=0.d0
  allocate(vpot3(ndim_tot));vpot3=0.d0

end subroutine generate_configurations

subroutine read_twobody_matrix(filename)
  USE single_particle_orbits
  USE configurations
  USE constants
  use gmat_storage
    implicit none
  TYPE (configuration_descriptor) :: bra_configs
  integer:: isospin_z, p_parity, ang_mom
  integer:: iostatus
  integer:: tz,ipar,jtot,a,b,c,d
  real*8 :: gmat,p2,hc
  integer:: bra,ket,aa,bb,cc,dd
  integer:: ndim,ndim_tot
  integer:: number_channels
  CHARACTER(LEN=100) :: filename
  integer:: isob(3),tt


  
      OPEN (UNIT=12,FILE=filename, STATUS='OLD', &
           ACCESS='STREAM', FORM='unformatted' )


  vpot1=0.d0
  vpot2=0.d0
  vpot3=0.d0



  isob(1)=-1
  isob(2)=0
  isob(3)=1

  
ndim_tot=0

 number_channels=0
 do tt=1,3
  isospin_z=isob(tt)
     !     loop over parity values, here positive parity is 0, negative 1
     DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=0, 60 


           
          
           CALL  number_gmatrix_confs&
                (ang_mom,p_parity,isospin_z,4,bra_configs)
           IF (bra_configs%number_confs <= 0 ) CYCLE
           ndim=bra_configs%number_confs
           allocate(bra_configs%config_ab(2*ndim))
           CALL  setup_gmatrix_configurations&
                (ang_mom,p_parity,isospin_z,4,bra_configs)

        number_channels=number_channels+1




         do bra=1,ndim
         aa=bra_configs%config_ab(2*bra-1)
         bb=bra_configs%config_ab(2*bra)
         do ket=bra,ndim
         cc=bra_configs%config_ab(2*ket-1)
         dd=bra_configs%config_ab(2*ket)

        READ(12,IOSTAT=iostatus) tz, ipar, jtot, a, b, c, d, gmat, p2, hc
        ndim_tot=ndim_tot+1
        vpot1(ndim_tot)=gmat 
        vpot2(ndim_tot)=p2
        vpot3(ndim_tot)=hc 
!        write(71,*)'sdsds',ndim
       write(71,'(11(i4,x))')tz,ipar,jtot,a,b,c,d,aa,bb,cc,dd
      !  write(72,'(7(i4,x),3(F23.6,1x))')isospin_z,p_parity,ang_mom*2,aa, bb, cc, dd,gmat,p2,hc
         !@! if(2*ang_mom/=jtot)write(*,*)'error'
         !@! if(ipar/=p_parity)write(*,*)'error'
         !@! if(isospin_z/=tz)write(*,*)'error'
         !@! if(aa/=a)write(*,*)'error'
         !@! if(bb/=b)write(*,*)'error'
         !@! if(cc/=c)write(*,*)'error'
         !@! if(dd/=d)write(*,*)'error'
              enddo
              enddo
              deallocate(bra_configs%config_ab)


           
        end DO
     end DO
  end DO
  close(12)
    
end subroutine read_twobody_matrix


subroutine write_twobody_matrix(filename)
  USE single_particle_orbits
  USE configurations
  USE constants
  use gmat_storage
  USE hdf5
  USE h5lt
  use hdf5_integer_matrix_wrapper
  use hdf5_double_matrix_wrapper
    implicit none
  integer:: ndim
  CHARACTER(LEN=100) :: filename

  INTEGER(HID_T) :: file_id
  INTEGER :: hdferr
  INTEGER, DIMENSION(:,:), POINTER :: dummy_int
  real*8,allocatable:: tvec(:,:)
  INTEGER :: error

      ndim=size(vpot1)
      allocate(tvec(ndim,1))

        call h5open_f(error)
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, hdferr)
        allocate(dummy_int(1,1))
        dummy_int=(ndim)
        error=write_integer_matrix(file_id, dummy_int, "number_mtx_elements")
        tvec(:,1) = vpot1
        error = write_double_matrix(file_id, tvec, "gmat")
        tvec(:,1) = vpot2
        error = write_double_matrix(file_id, tvec, "p2")
        tvec(:,1) = vpot3
        error = write_double_matrix(file_id, tvec, "hc")
        CALL h5fclose_f(file_id, hdferr)
        deallocate(tvec)
    
end subroutine write_twobody_matrix
