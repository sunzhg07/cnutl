subroutine sm_basis
  use parallel
  use mpi_mapping
  use lanczos
  use constants
  use m_sp
  implicit none
  integer*16 :: low, high,one,tp
  integer*16:: u,ur,nh,no
  integer :: i,j,a,b,tmp,tmp1
  integer:: bra,ket,bra_start,bra_end,ket_start,ket_end
  integer:: bra_max,bra_min,ket_min,ket_max
  integer:: nv,nvec,nvres

  one=1
  if(np==0)then 
    write(*,*) 'number protons 0'
    np_mbasis=1
    allocate(mbsp(1))
    allocate(mzp(1))
    allocate(iparp(1))
    mbsp(1)=0
    mzp(1)=0
    iparp(1)=0
  else
    tp=0
    low=0
    do i=1,np
      one=1
      tp=ishft(one,i-1)
      low=low+tp
    enddo
    high=ishft(low,np_orb-np)
    !write(*,*)low,high
    u=low
    np_mbasis=1
    do while((u>=low) .and.( u<high))
      ur=and(u , -u)
      nh=u+ur
      no=xor(u,nh)
      no=no/ur
      no=ishft(no,-2)
      u=or(nh,no)
      np_mbasis=np_mbasis+1
      !   call decode(u)
    enddo
    !write(*,*)'proton basis',np_mbasis
    !call decode(low)
    !call decode(high)

    allocate(mbsp(np_mbasis))
    allocate(mzp(np_mbasis)); mzp=0
    allocate(iparp(np_mbasis)); iparp=0
    np_mbasis=1
    tp=0
    low=0
    do i=1,np
      one=1
      tp=ishft(one,i-1)
      low=low+tp
    enddo
    high=ishft(low,np_orb-np)
    mbsp(1)=low
    u=low
    !write(*,*)low,high

    do while((u>=low) .and.( u<high))
      ur=and(u,-u)
      nh=u+ur
      no=xor(u,nh)
      no=no/ur
      no=ishft(no,-2)
      u=or(nh, no)
      np_mbasis=np_mbasis+1
      mbsp(np_mbasis)=u
    enddo
    do i=1,np_mbasis
      tp=mbsp(i)
      tmp=0
      tmp1=0
      do j=1,np_orb
        if(and(tp,one) ==one)then
          tmp=tmp+all_orbits%mm(j+nn_orb)
          tmp1=tmp1+all_orbits%ll(j+nn_orb)
        endif
        tp=ishft(tp,-1)
      enddo
      mzp(i)=tmp
      iparp(i)=tmp1
    enddo


  endif

  if(iam==0)write(*,*)'number proton basis', np_mbasis

  tp=0
  low=0
  if(nn==0)then 
    write(*,*) 'number protons 0'
    nn_mbasis=1
    allocate(mbsn(1))
    allocate(mzn(1)); mzn=0
    allocate(iparn(1)); iparn=0
    mbsn(1)=0
  else
    do i=1,nn
      one=1
      tp=ishft(one,i-1)
      low=low+tp
    enddo
    high=ishft(low,nn_orb-nn)

    u=low
    nn_mbasis=1
    do while((u>=low) .and.( u<high))
      ur=and(u,-u)
      nh=u+ur
      no=xor(u,nh)
      no=no/ur
      no=ishft(no,-2)
      u=or(nh, no)
      nn_mbasis=nn_mbasis+1
      ! call decode(u)
    enddo
    !call decode(low)
    !call decode(high)

    allocate(mbsn(nn_mbasis))
    allocate(mzn(nn_mbasis)); mzn=0
    allocate(iparn(nn_mbasis)); iparn=0

    tp=0
    low=0
    do i=1,nn
      one=1
      tp=ishft(one,i-1)
      low=low+tp
    enddo
    high=ishft(low,nn_orb-nn)

    u=low
    nn_mbasis=1
    mbsn(1)=low
    do while((u>=low) .and.( u<high))
      ur=and(u,-u)
      nh=u+ur
      no=xor(u,nh)
      no=no/ur
      no=ishft(no,-2)
      u=or(nh,no)
      nn_mbasis=nn_mbasis+1
      mbsn(nn_mbasis)=u

    enddo

    do i=1,nn_mbasis
      tp=mbsn(i)
      tmp=0
      tmp1=0
      do j=1,nn_orb
        if(and(tp,one) ==one)then
          tmp=tmp+all_orbits%mm(j)
          tmp1=tmp1+all_orbits%ll(j)
        endif
        tp=ishft(tp,-1)
      enddo
      mzn(i)=tmp
      iparn(i)=tmp1
    enddo


  endif
  if(iam==0)write(*,*)'number neutron states: ', nn_mbasis






  num_mbasis=0
  do a=1,np_mbasis
    do b=1,nn_mbasis
      if(mzp(a)+mzn(b)==jtot .and. mod(iparp(a)+iparn(b),2)==nu_parity)then
        num_mbasis=num_mbasis+1
      endif
    enddo
  enddo

  if(iam==0)write(*,*)'total dimension of hamiltonian', num_mbasis
  if(iam==0)write(*,*)'distribute matrix on',n_submat**2,'procs'



  !! mapping matrix
  do bra=1,n_submat
    do ket=1,n_submat
      if(((bra-1)*n_submat+ket-1)==iam)then
        sub_row=bra
        sub_col=ket
      endif
    enddo
  enddo

  nvec=num_mbasis/n_submat 
  nvres=num_mbasis-nvec*n_submat

  allocate(mapping(1:n_submat,2))
  nv=0
  do a=1,n_submat
    if(a<=nvres)then
    bra_start=nv+1
    bra_end=nv+nvec+1
    mapping(a,1)=bra_start
    mapping(a,2)=bra_end
    nv=nv+nvec+1
  else
    bra_start=nv+1
    bra_end=nv+nvec
    mapping(a,1)=bra_start
    mapping(a,2)=bra_end
    nv=nv+nvec
  endif
enddo
if(iam==0)write(*,*)nv,num_mbasis





   bra_min=mapping(sub_row,1)
   bra_max=mapping(sub_row,2)
  allocate(mbs_row(bra_min:bra_max,2));mbs_row=0
   bra_min=mapping(sub_col,1)
   bra_max=mapping(sub_col,2)
  allocate(mbs_col(bra_min:bra_max,2));mbs_col=0

   bra_min=mapping(sub_row,1)
   bra_max=mapping(sub_row,2)
   num_mbasis=0

  do a=1,np_mbasis
    do b=1,nn_mbasis
      if(mzp(a)+mzn(b)==jtot .and. mod(iparp(a)+iparn(b),2)==nu_parity)then
        num_mbasis=num_mbasis+1
        if(num_mbasis>bra_max)cycle
        if(num_mbasis<bra_min)cycle

        mbs_row(num_mbasis,1)=a
        mbs_row(num_mbasis,2)=b
      endif
    enddo
  enddo

   ket_min=mapping(sub_col,1)
   ket_max=mapping(sub_col,2)

   num_mbasis=0
  do a=1,np_mbasis
    do b=1,nn_mbasis
      if(mzp(a)+mzn(b)==jtot .and. mod(iparp(a)+iparn(b),2)==nu_parity)then
        num_mbasis=num_mbasis+1
        if(num_mbasis>ket_max)cycle
        if(num_mbasis<ket_min)cycle

        mbs_col(num_mbasis,1)=a
        mbs_col(num_mbasis,2)=b
      endif
    enddo
  enddo

  allocate(haml(bra_min:bra_max,ket_min:ket_max))
  allocate(and_haml(max_iter_lanc,max_iter_lanc))
  allocate(and_haml2(max_iter_lanc,max_iter_lanc))
  allocate(and_eig(max_iter_lanc))
  allocate(and_eig2(max_iter_lanc))
  allocate(and_vec(max_iter_lanc,max_iter_lanc))
  allocate(and_lvec(max_iter_lanc,max_iter_lanc))
!@!  allocate(cvec(bra_min:bra_max,ket_min:ket_max))
!@!  allocate(cvecl(bra_min:bra_max,ket_min:ket_max))
  allocate(vj(1:num_mbasis))
  allocate(vi(1:num_mbasis))
!@!  allocate(vk(1:num_mbasis))
  haml=0.0_dpc
  and_haml=0.0_dpc
!@!  cvec=0.0_dpc
  vi=0.0_dpc
  vj=0.0_dpc
!@!  vk=0.0_dpc
  allocate(vec(num_mbasis,max_iter_lanc))
  vec=0.d0
if(iam==0)then
  write(6,*)'Size of Hamiltonian',dble(num_mbasis*num_mbasis)*4./1.e9
  write(6,*)'Size of Vectors',dble((max_iter_lanc+2)*num_mbasis)*4./1.e9

endif



  allocate(mapping_1d(1:num_procs,2))


  nvec=num_mbasis/num_procs
  nvres=num_mbasis-nvec*num_procs

  nv=0
  do a=1,num_procs
    if(a<=nvres)then
    bra_start=nv+1
    bra_end=nv+nvec+1
    mapping_1d(a,1)=bra_start
    mapping_1d(a,2)=bra_end
    nv=nv+nvec+1
  else
    bra_start=nv+1
    bra_end=nv+nvec
    mapping_1d(a,1)=bra_start
    mapping_1d(a,2)=bra_end
    nv=nv+nvec
  endif
enddo



end subroutine sm_basis
