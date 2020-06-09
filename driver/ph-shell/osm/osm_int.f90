subroutine ini_interactions
  use constants
  use storage
  use m_sp
  use m_configurations
  use file_list
  implicit none
  
  integer:: itz,p_parity,ang_mom
  integer:: number_channels,channel
  type(m_configuration_descriptor)::bra_configs
  integer:: nmtx,a,b,c,d,e,f,i,j,k,bra,ket,ii
  integer:: nconfs
  real*8:: gmat,gmat1
  !! one body interaction
  write(*,*)tot_orbs
  allocate(vt1(tot_orbs,tot_orbs))
  vt1=0.d0
    open(unit=20,file=onebd_int)
    read(20,*)nconfs
    do i=1,nconfs
      read(20,*) a, b, gmat
      vt1(a,b)=gmat
    enddo
    close(20)
  

  !! construct channels 2b
  number_channels=0
  nconfs=0
    DO itz=-1,1 
      !     loop over parity values, here positive parity is 0, negative 1
      DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=-101,101
          call number_vm_confs(ang_mom,p_parity,itz,3,bra_configs)
          if(bra_configs%number_confs==0)cycle
          number_channels=number_channels+1
        enddo
      enddo
    enddo

    allocate(vt2(number_channels))


    allocate(ip2(tot_orbs,tot_orbs,2))
    ip2=0

  number_channels=0
    DO itz=-1,1 
      !     loop over parity values, here positive parity is 0, negative 1
      DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=-101,101
          call number_vm_confs(ang_mom,p_parity,itz,3,bra_configs)
          if(bra_configs%number_confs==0)cycle
          number_channels=number_channels+1
          nconfs=bra_configs%number_confs
          allocate(bra_configs%config_ab(2*nconfs))
          call setup_vm_confs(ang_mom,p_parity,itz,3,bra_configs)
          allocate(vt2(number_channels)%val(nconfs,nconfs))
          vt2(number_channels)%val=0.d0

          do  bra=1,nconfs
            a=bra_configs%config_ab(2*bra-1)
            b=bra_configs%config_ab(2*bra)
            ip2(a,b,1)=number_channels
            ip2(a,b,2)=bra
          enddo
          deallocate(bra_configs%config_ab)

        enddo
      enddo
   
    enddo

    open(unit=21,file=twobd_int)

    read(21,*)nmtx
    do  ii=1,nmtx
      read(21,*)a,b,c,d,gmat
      if(a*b<=0)cycle
      if(c*d<=0)cycle
      channel=ip2(a,b,1)
      if(channel==0)then
        cycle
      endif
      bra=ip2(a,b,2)
      ket=ip2(c,d,2)
      if(bra*ket==0)then
        cycle
      endif
      vt2(channel)%val(bra,ket)=gmat
    enddo
    read(21,*)nmtx
    do  ii=1,nmtx
      read(21,*)a,b,c,d,gmat
      if(a*b<=0)cycle
      if(c*d<=0)cycle
      channel=ip2(a,b,1)
      if(channel==0)then
        cycle
      endif
      bra=ip2(a,b,2)
      ket=ip2(c,d,2)
      if(bra*ket==0)then
        cycle
      endif
      vt2(channel)%val(bra,ket)=gmat
    enddo
    read(21,*)nmtx
    do  ii=1,nmtx
      read(21,*)a,b,c,d,gmat
      if(a*b<=0)cycle
      if(c*d<=0)cycle
      channel=ip2(a,b,1)
      if(channel==0)then
        cycle
      endif
      bra=ip2(a,b,2)
      ket=ip2(c,d,2)
      if(bra*ket==0)then
        cycle
      endif
      vt2(channel)%val(bra,ket)=gmat
    enddo
    close(21)



    if(number_interactions==3)then
  number_channels=0
    DO itz=-3,3,2 
      !     loop over parity values, here positive parity is 0, negative 1
      DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=-101,101
          call number_vm3_confs(ang_mom,p_parity,itz,3,bra_configs)
          if(bra_configs%number_confs==0)cycle
          number_channels=number_channels+1
        enddo
      enddo
    enddo

    allocate(vt3(number_channels))
    allocate(ip3(tot_orbs,tot_orbs,tot_orbs,2))
    ip3=0

  number_channels=0
    DO itz=-3,3,2
      !     loop over parity values, here positive parity is 0, negative 1
      DO p_parity=0,1           
        !     loop over angular momenta
        DO ang_mom=-101,101
          call number_vm3_confs(ang_mom,p_parity,itz,3,bra_configs)
          if(bra_configs%number_confs==0)cycle
          number_channels=number_channels+1
          nconfs=bra_configs%number_confs
          allocate(bra_configs%config_ab(3*nconfs))
          call setup_vm3_confs(ang_mom,p_parity,itz,3,bra_configs)
          allocate(vt3(number_channels)%val(nconfs,nconfs))
          vt3(number_channels)%val=0.d0

          do  bra=1,nconfs
            a=bra_configs%config_ab(3*bra-2)
            b=bra_configs%config_ab(3*bra-1)
            c=bra_configs%config_ab(3*bra)
            ip3(a,b,c,1)=number_channels
            ip3(a,b,c,2)=bra
          enddo
          deallocate(bra_configs%config_ab)

        enddo
      enddo
   
    enddo


    open(unit=22,file=threebd_int)

    read(22,*)nmtx
    do  ii=1,nmtx
      read(22,*)a,b,c,i,j,k,gmat
      channel=ip3(a,b,c,1)
      if(channel==0)then
        cycle
      endif
      bra=ip3(a,b,c,2)
      ket=ip3(i,j,k,2)
      if(bra*ket==0)cycle
      vt3(channel)%val(bra,ket)=gmat
    enddo
  endif


end subroutine ini_interactions
