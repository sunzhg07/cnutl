subroutine ini_files
  use parallel
  use constants
  use storage
  use m_sp
  use file_list
  use truncations
  implicit none
  integer :: i,n
  open(file='sm.dat',unit=5)
  read(5,*);read(5,*)np,nn
  read(5,*);
  read(5,*)number_interactions
  read(5,*)
  read(5,*)jtot,nu_parity
  read(5,*)
  read(5,*)onebd_int
  read(5,*)
  read(5,*)twobd_int
  read(5,*)
  read(5,*)threebd_int
  read(5,*)
  read(5,*)onebd_opr
  read(5,*)
  read(5,*)twobd_opr
  read(5,*)
  read(5,*)threebd_opr
  read(5,*)
  read(5,*)output
!  open(file=output,unit=6)
  read(5,*)
  read(5,*)sp_file
  read(5,*)
  read(5,*)nstate
  read(5,*)
  read(5,*)max_iter_lanc
  read(5,*)
  read(5,*)trunc
  if(trunc)then
  read(5,*)
    read(5,*)number_spj_orbs
    allocate(occ_lim(2,number_spj_orbs))
    allocate(occ_mask(number_spj_orbs))
    allocate(obtype(number_spj_orbs))
    occ_mask=0
    occ_lim=0
    obtype=0
    do i=1,number_spj_orbs
      read(5,*)occ_lim(1,i),occ_lim(2,i)
    enddo
  endif


  open(unit=7,file=sp_file)
  read(7,*)tot_orbs
  call allocate_sp_array(all_orbits,tot_orbs)

  np_orb=0
  nn_orb=0
  all_orbits%total_orbits=tot_orbs
  do i=1,tot_orbs
    read(7,*)all_orbits%nn(i),&
      & all_orbits%nshell(i),&
      & all_orbits%ll(i),&
      & all_orbits%jj(i),&
      & all_orbits%mm(i),&
      & all_orbits%itzp(i),&
       all_orbits%jord(i),&
       all_orbits%occ(i)
    if(all_orbits%itzp(i)==-1)then
      np_orb=np_orb+1
    elseif(all_orbits%itzp(i)==1)then
      nn_orb=nn_orb+1
    else
      write(*,*)'error in single particle itzp'
    endif

  enddo
  close(7)

  call setup_occ_mask
  allocate(obs(tot_orbs))
  obs=0
  n=0
  do i=1,all_orbits%total_orbits
    if(all_orbits%itzp(i)==1)then
      n=n+1
      obs(n)=i
    endif
  enddo
  do i=1,all_orbits%total_orbits
    if(all_orbits%itzp(i)==-1)then
      n=n+1
      obs(n)=i
    endif
  enddo

  if(iam==0)write(*,*)'SM calculation starting..'
  if(iam==0)write(*,*)'# of protons and neutrons ',np,nn
  if(iam==0)write(*,*)'# interaction files ',number_interactions
  if(iam==0)write(*,*)'jtot = ',jtot
  if(iam==0)write(*,*)'interaction files:'
  if(iam==0)write(*,*)onebd_int
  if(iam==0)write(*,*)twobd_int
  if(iam==0)write(*,*)threebd_int
  if(iam==0)write(*,*)'operator files:'
  if(iam==0)write(*,*)onebd_opr
  if(iam==0)write(*,*)twobd_opr
  if(iam==0)write(*,*)threebd_opr
  if(iam==0)write(*,*)'single partitcle file'
  if(iam==0)write(*,*)sp_file
  if(iam==0)write(*,*)'# of converge lanczos vectors'
  if(iam==0)write(*,*)nstate
  if(iam==0)write(*,*)'# of max lanczos iterations'
  if(iam==0)write(*,*)max_iter_lanc

end subroutine ini_files
