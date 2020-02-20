program sharedmemtest
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR, C_F_POINTER
  use mpi
  implicit none
  integer, parameter :: dp = selected_real_kind(14,200)
  integer :: win,win2,hostcomm,hostrank
  INTEGER(KIND=MPI_ADDRESS_KIND) :: windowsize
  INTEGER :: disp_unit,my_rank,ierr,total
  TYPE(C_PTR) :: baseptr,baseptr2
  real(dp), POINTER :: matrix_elementsy(:,:,:,:)
  integer,allocatable :: arrayshape(:)

  call MPI_INIT( ierr )

  call MPI_COMM_RANK(MPI_COMM_WORLD,MY_RANK,IERR)  !GET THE RANK OF ONE PROCESS                                                                                                                                                                                                
  call MPI_COMM_SIZE(MPI_COMM_WORLD,Total,IERR)  !GET THE TOTAL PROCESSES OF THE COMM                                                                                                                                                                                          
  CALL MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, hostcomm,ierr)
  CALL MPI_Comm_rank(hostcomm, hostrank,ierr)

  allocate(arrayshape(4))
  arrayshape=(/ 10,10,10,10 /)
  if (hostrank == 0) then
     windowsize = int(10**4,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND !*8 for double ! Put the actual data size here                                                                                                                                                                
  else
     windowsize = 0_MPI_ADDRESS_KIND
  end if

  disp_unit = 1
  CALL MPI_Win_allocate_shared(windowsize, disp_unit, MPI_INFO_NULL, hostcomm, baseptr, win, ierr)    

  ! Obtain the location of the memory segment                                                                                                                                                                                                                                  
  if (hostrank /= 0) then
     CALL MPI_Win_shared_query(win, 0, windowsize, disp_unit, baseptr, ierr)     
  end if

  ! baseptr can now be associated with a Fortran pointer                                                                                                                                                                                                                       
  ! and thus used to access the shared data                                                                                                                                                                                                                                    
  CALL C_F_POINTER(baseptr, matrix_elementsy,arrayshape)

  !!! your code here!                                                                                                                                                                                                                                                          
  !!! sample below                                                                                                                                                                                                                                                             

  if (hostrank == 0) then
     matrix_elementsy=0.0_dp
     matrix_elementsy(1,2,3,4)=1.0_dp
  end if
  CALL MPI_WIN_FENCE(0, win, ierr)

  print *,"my_rank=",my_rank,matrix_elementsy(1,2,3,4),matrix_elementsy(1,2,3,5)

  !!! end sample code                                                                                                                                                                                                                                                          

  call MPI_WIN_FENCE(0, win, ierr) 
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
  call MPI_Win_free(win,ierr)     
  call MPI_FINALIZE(IERR)

  end program
