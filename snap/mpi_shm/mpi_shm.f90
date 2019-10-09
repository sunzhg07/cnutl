program shm_test
  use mpi
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR,C_F_POINTER
  implicit none
  integer,parameter:: dp=kind(1.d0)
  integer:: win
  integer:: hostcomm,rootcomm
  integer:: nproc,iam,h_nproc,h_iam,r_nproc,r_iam,nroot
  INTEGER(KIND=MPI_ADDRESS_KIND) :: wsize
  integer:: disp,ierror
  integer:: color,key
  integer:: i,j,disp_unit
  type(c_ptr):: blk_a
  real(dp), POINTER :: a(:)
  real(dp),allocatable,dimension(:)::b
  ! all threads in whole word
  call mpi_init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD,iam,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

!! split threads based on nodes, all threads belog to same nodes build a new communicator,and have a internal rank
  call MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,hostcomm,ierror)
  call MPI_Comm_rank(hostcomm,h_iam,ierror)
  call MPI_COMM_SIZE(hostcomm,h_nproc,ierror)
  color=0
  if(h_iam==0)color=1
  key=iam

  !! all threads have intenal rank 0(first threads of each node) form a new communacator
  call MPI_Comm_split(MPI_COMM_WORLD,color,key,rootcomm,ierror)
  call MPI_Comm_rank(rootcomm,r_iam,ierror)
  call MPI_COMM_SIZE(rootcomm,nroot,ierror)
   if(color==1)then
  write(*,*)'root_order',r_iam,'world_order',iam,'node_order',h_iam
endif

allocate(b(10))
b=0.d0



   !! first thread on each node allocate a share memory
  if(h_iam==0)then
    wsize=int(10,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
  else
    wsize=0_MPI_ADDRESS_KIND
  endif
  disp_unit=1
   !! other threads on the node have access to the share memory of thread 0
    CALL MPI_Win_allocate_shared(wsize, disp_unit, MPI_INFO_NULL, hostcomm, blk_a, win, ierror)
    if(h_iam/=0)then
      call MPI_Win_shared_query(win,0,wsize,disp_unit,blk_a,ierror)
    endif

    CALL C_F_POINTER(blk_a, a,[10])


    if(h_iam==0)then
      do i=1,10
        a(i)=1.d0/i
      enddo
    endif

    CALL MPI_WIN_FENCE(0,win,ierror)
    if(h_iam==10)then
      write(*,*)a(1:10)
    endif

    CALL MPI_WIN_FENCE(0,win,ierror)

    !@!if(h_iam==0 .and.r_iam==0)then
    !@!  do i=1,10
    !@!    a(i)=2.d0/i
    !@!  enddo
    !@!endif
    call mpi_barrier(hostcomm,ierror)
    call mpi_barrier(MPI_COMM_WORLD,ierror)
    b=a
    call mpi_barrier(hostcomm,ierror)
    call mpi_barrier(MPI_COMM_WORLD,ierror)

    !! operation of 
    !call mpi_bcast(a,size(a),mpi_real8,0,rootcomm,ierror)
    call mpi_allreduce(mpi_in_place, b,size(b),mpi_real8,mpi_sum,rootcomm,ierror)
    !call mpi_reduce(mpi_in_place,a,size(a),mpi_real8,mpi_sum,0,rootcomm,ierror)

    call mpi_barrier(rootcomm,ierror)
    call mpi_barrier(hostcomm,ierror)
    call mpi_barrier(MPI_COMM_WORLD,ierror)

    CALL MPI_WIN_FENCE(0,win,ierror)


    if(h_iam==1)then
      write(*,*)b(1:10)
    endif

    CALL MPI_WIN_FENCE(0,win,ierror)


  call MPI_FINALIZE(ierror)




end program shm_test
