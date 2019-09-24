program osm_main
  use parallel
  use constants
  use storage
  use m_sp
  use mpi_mapping


  implicit none
  real*8  ::  startwtime , endwtime

  call mpi_init(ierror)
  call mpi_comm_rank(mpi_comm_world,iam,ierror)
  call mpi_comm_size(mpi_comm_world,num_procs,ierror)
  startwtime = MPI_WTIME()
  master=0
  n_submat=int(sqrt(real(num_procs)))
  if(iam==0) write(*,*)n_submat,num_procs
  if(n_submat**2 /= num_procs)then
    if(iam==0)write(*,*)'number of procs must be int**2'
   call mpi_finalize
   stop
  endif


  call ini_files

  call ini_interactions

  call mpi_barrier(mpi_comm_world,ierror)

  call sm_basis

  call mpi_barrier(mpi_comm_world,ierror)
  call sm_matrix

  if(iam==0) write(*,*)'Starting arnoldi iteration:'

  call sm_lanczos

  call mpi_barrier(mpi_comm_world,ierror)
  endwtime = MPI_WTIME()
  if ( iam == 0 ) write(6,*) 'Total execution time:', endwtime - startwtime
  call mpi_finalize(ierror)


end program osm_main
