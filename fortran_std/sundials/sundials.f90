module ode_mod
    use, intrinsic :: iso_c_binding
    use fcvode_mod                    ! Fortran interface to CVODE
    use fsundials_context_mod         ! Fortran interface to SUNContext
    use fsundials_nvector_mod         ! Fortran interface to generic N_Vector
    use fnvector_serial_mod           ! Fortran interface to serial N_Vector
    use fsundials_nonlinearsolver_mod ! Fortran interface to generic SUNNonlinearSolver
    use fsunnonlinsol_fixedpoint_mod  ! Fortran interface to fixed point SUNNonlinearSolver
    implicit none

    integer(c_long):: neq 
    type(c_ptr)    :: cvode_mem  ! CVODE memory
    type(c_ptr)    :: ctx        ! SUNDIALS simulation context
    type(N_Vector), pointer           :: sunvec_y,sunvec_f ! sundials vector
    type(SUNNonlinearSolver), pointer :: sunnls   ! sundials fixed-point nonlinear solver
    integer(c_int) :: ierr       ! error flag from C functions
    real(c_double) :: rtol, atol ! relative and absolute tolerance
    real(c_double) :: tstart     ! initial time

contains
    !! the right hand side function
    integer(c_int) function RhsFn(tn, sunvec_y, sunvec_f, user_data) &
            result(ierr) bind(C,name='RhsFn')
        implicit none

        ! calling variables
        real(c_double), value :: tn        ! current time
        type(N_Vector)        :: sunvec_y  ! solution N_Vector
        type(N_Vector)        :: sunvec_f  ! rhs N_Vector
        type(c_ptr), value    :: user_data ! user-defined data

        ! pointers to data in SUNDIALS vectors
        real(c_double), pointer :: yvec(:)
        real(c_double), pointer :: fvec(:)

        !! this interface is essential to keep our derivs external
        interface
            subroutine derivs(neq,yvec,fvec,x) 
                use, intrinsic :: iso_c_binding
                real(c_double),dimension(:),intent(inout):: yvec,fvec
                integer(c_long),intent(in):: neq
                real(c_double),intent(in) :: x
            end subroutine derivs
        end interface

        !======= Internals ============

        ! get data arrays from SUNDIALS vectors
        yvec => FN_VGetArrayPointer(sunvec_y)
        fvec => FN_VGetArrayPointer(sunvec_f)
        ! the actual calcualtion of the derivs is through the external subroutine

        call derivs(neq, yvec, fvec,tn)
        ierr = 0
        return
    end function RhsFn

    subroutine cvode_init(yvec, neqs)
        implicit none
    real(c_double),dimension(:),intent(inout) :: yvec(:)
    integer(c_long):: neqs

        ierr = FSUNContext_Create(c_null_ptr, ctx)


        ! create SUNDIALS N_Vector
        sunvec_y => FN_VMake_Serial(neqs, yvec, ctx)
        if (.not. associated(sunvec_y)) then
            print *, 'ERROR: sunvec = NULL'
            stop 1
        end if

        ! create CVode memory
        cvode_mem = FCVodeCreate(CV_ADAMS, ctx)
        if (.not. c_associated(cvode_mem)) then
            print *, 'ERROR: cvode_mem = NULL'
            stop 1
        end if

        ! initialize CVode
        ierr = FCVodeInit(cvode_mem, c_funloc(RhsFn), tstart, sunvec_y)
        if (ierr /= 0) then
            print *, 'Error in FCVodeInit, ierr = ', ierr, '; halting'
            stop 1
        end if

        ! set relative and absolute tolerances
        rtol = 1.0d-6
        atol = 1.0d-10

        ierr = FCVodeSStolerances(cvode_mem, rtol, atol)
        if (ierr /= 0) then
            print *, 'Error in FCVodeSStolerances, ierr = ', ierr, '; halting'
            stop 1
        end if

        ! create fixed point nonlinear solver object
        sunnls => FSUNNonlinSol_FixedPoint(sunvec_y, 0, ctx)

        if (.not. associated(sunnls)) then
            print *,'ERROR: sunnls = NULL'
            stop 1
        end if

        ! attache nonlinear solver object to CVode
        ierr = FCVodeSetNonlinearSolver(cvode_mem, sunnls)
        if (ierr /= 0) then
            print *, 'Error in FCVodeSetNonlinearSolver, ierr = ', ierr, '; halting'
            stop 1
        end if
    end subroutine cvode_init



    subroutine cvode_finalize
        ! clean up
        call FCVodeFree(cvode_mem)
        ierr = FSUNNonLinSolFree(sunnls)
        call FN_VDestroy(sunvec_y)
        ierr = FSUNContext_Free(ctx)
    end subroutine cvode_finalize

end module ode_mod


subroutine derivs(neq,yvec,fvec,x) 
    use, intrinsic :: iso_c_binding
    implicit none
    real(c_double),dimension(:),intent(inout):: yvec,fvec
    integer(c_long),intent(in):: neq
    real(c_double),intent(in) :: x

    fvec(1) = yvec(1)*yvec(2) + cos(x)-0.5*sin(2*x)
    fvec(2) = yvec(1)*yvec(1) +yvec(2)*yvec(2) - (1.+sin(x))

end subroutine derivs

program main
    use ode_mod                       ! ODE functions
    implicit none

    ! local variables
    real(c_double) :: tend       ! final time
    real(c_double) :: dtout      ! output time interval
    real(c_double) :: tout       ! output time
    real(c_double) :: tcur(1)    ! current time
    integer(c_int) :: nout       ! number of outputs



    ! solution vector, neq is set in the ode_mod module
    real(c_double),allocatable :: yvec(:)
    integer:: outstep

    neq=2
    allocate(yvec(neq))
    yvec=0.d0
    nout=25


    call cvode_init(yvec,neq)



        tstart=0.d0
        tout=0.d0
        dtout  = 0.08d0
        tend=2.0

    do outstep = 1,nout

        ! call CVode
        tout = min(tout + dtout, tend)

        ierr = FCVode(cvode_mem, tout, sunvec_y, tcur, CV_NORMAL)
        if (ierr /= 0) then
            print *, 'Error in FCVODE, ierr = ', ierr, '; halting'
            stop 1
        endif
        ! output current solution
        write(*,'(3(F13.6,1x))') tcur, yvec(1),yvec(2)

    enddo

    call cvode_finalize

end program main


