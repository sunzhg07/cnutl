module data_storage
  use mpi
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR,C_F_POINTER

  type, public :: block_data
    type(c_ptr):: blk_a
    real(dp),pointer:: a(:)
  end type block_data
end module data_storage
