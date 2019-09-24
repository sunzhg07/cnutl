MODULE arnoldi
    USE parallel
    use mbasis
    IMPLICIT NONE


    DOUBLE PRECISION, EXTERNAL :: DDOT
    DOUBLE PRECISION :: zero = 1.0D-4
    DOUBLE PRECISION :: tol = 1.0D-6
    INTEGER :: maxit = 150, debug
CONTAINS
    SUBROUTINE solve_eigensystem_arnoldi_complex(eval, evect, start_vector, mydebug)
        COMPLEX*16, DIMENSION(:) :: eval
        COMPLEX*16, DIMENSION(:,:) :: evect
        COMPLEX*16, DIMENSION(:) :: start_vector
        INTEGER, OPTIONAL :: mydebug

        COMPLEX*16 :: local_buf
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: work_vector, local_vec
        COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: arnoldi, h
        
        
        DOUBLE PRECISION :: memory
        INTEGER*8 :: n_eval, n_evect
        INTEGER :: i, k, n_converged, proc, start, stopv, sizev

        DOUBLE PRECISION :: norm, t_start, t_start2, t_stop, t_stop2, t_start3, &
                            t_stop3, mvp_tsum

        t_start3 = MPI_WTIME()
        debug = 0
        IF (PRESENT(mydebug)) debug = mydebug

        evect = 0.0D0; eval = 0.0D0; mvp_tsum = 0.0D0
        n_eval = SIZE(eval); n_evect = SIZE(start_vector)

        ! Normalize start vector
        norm = SUM(CONJG(start_vector)*start_vector)
        start_vector(:) = start_vector(:)/sqrt(norm)
        IF (ABS(norm) < 1.0E-8) THEN
            IF ( iam == master) WRITE(*,*) 'solve_eigensystem_arnoldi: Zero start_vector!'
            STOP
        ENDIF
        
        ! Setup variables for parallel computation.
        CALL setup_parallel_lookup(tmin, n_evect)
        memory = 16.*tmin%my_size*(2*maxit+2) + 16.*n_evect*2
        IF (debug > 1) WRITE(*,*) 'Memory arnoldi: ', iam, memory/1024D0/1024D0


        ALLOCATE(work_vector(n_evect), arnoldi(tmin%my_size,0:maxit))
        ALLOCATE(local_vec(tmin%my_size), h(0:maxit,0:maxit))

        arnoldi = 0.0D0; h = 0.0D0
        ! All the work goes here
        arnoldi(:,0) = start_vector(tmin%my_start:tmin%my_stop)
        DO i = 1, maxit
            t_start2 = MPI_WTIME()
            !Calculate the new Arnoldi vector
            t_start = MPI_WTIME()
            CALL mvp_sub(start_vector, work_vector, debug)
            t_stop = MPI_WTIME()
            IF (iam == master .AND. debug > 0) WRITE(*,*) 'Time mvp iteration: ', i, t_stop-t_start
            mvp_tsum = mvp_tsum + (t_stop-t_start)

            ! Orthogonalize the work vector and setup the Hessenberg matrix
            t_start = MPI_WTIME()
            local_vec = work_vector(tmin%my_start:tmin%my_stop)
            DO k = 0,i-1
                local_buf = SUM(CONJG(arnoldi(:,k))*local_vec)

                ! Everybody needs h(k,i-1) for next step
                CALL MPI_ALLREDUCE(local_buf, h(k,i-1), 1, MPI_DOUBLE_COMPLEX, MPI_SUM, &
                        mpi_comm_world, ierror)

                local_vec(:) = local_vec(:) - h(k,i-1)*arnoldi(:,k)
            ENDDO

            ! Update work_vector and broadcast to all nodes
            t_start = MPI_WTIME()
            IF (iam == master) THEN
                work_vector(tmin%my_start:tmin%my_stop) = local_vec
                DO proc = 1, num_procs - 1
                    start = tmin%all_starts(proc+1)
                    stopv = tmin%all_stops(proc+1)
                    sizev = tmin%all_sizes(proc+1)
                    CALL MPI_RECV(work_vector(start:stopv), sizev, MPI_DOUBLE_COMPLEX, proc, &
                        500, mpi_comm_world, mystatus, ierror)
                ENDDO
            ELSE
                CALL MPI_SEND(local_vec, tmin%my_size, MPI_DOUBLE_COMPLEX, master, &
                        500, mpi_comm_world,ierror)
            ENDIF
            CALL MPI_BCAST(work_vector, n_evect, MPI_DOUBLE_COMPLEX, master, mpi_comm_world, ierror)

            norm = SUM(CONJG(work_vector)*work_vector)
            h(i,i-1) = sqrt(norm)

            ! Normalize the workvector to get the i'th arnoldi vector and new
            ! start vector.
            start_vector(:) = work_vector(:)/sqrt(norm)
            arnoldi(:,i) = start_vector(tmin%my_start:tmin%my_stop)
            t_stop = MPI_WTIME()
            IF (iam == master .AND. debug > 0) WRITE(*,*) 'Time update of work vector: ', i, t_stop-t_start

            ! Compute the Ritz values and vectors
            t_start = MPI_WTIME()
            CALL compute_ritz_pairs_complex(h, arnoldi(:,0:i-1), eval, evect, n_converged)
            t_stop = MPI_WTIME()
            IF (iam == master .AND. debug > 0) WRITE(*,*) 'Time calculation of ritz pairs: ', i, t_stop-t_start

            t_stop2 = MPI_WTIME()
            IF (iam == master .AND. debug > 0) WRITE(*,*) 'Time for one arnoldi iteration: ', i, t_stop2-t_start2
            IF (n_converged == n_eval) EXIT
        ENDDO

        CALL MPI_BARRIER(mpi_comm_world, ierror)
        CALL takedown_parallel_lookup(tmin)
        t_stop3 = MPI_WTIME()
        IF (iam == master ) WRITE(*,*) 'Time for complete arnoldi process: ', i, t_stop3-t_start3
        IF (iam == master ) WRITE(*,*) 'Time in matrix vector product: ', i, mvp_tsum

    END SUBROUTINE solve_eigensystem_arnoldi_complex


    SUBROUTINE compute_ritz_pairs_complex(hessen, qmat, ev, evec, n_converged)
        COMPLEX*16, DIMENSION(:,:) :: hessen, qmat, evec
        COMPLEX*16, DIMENSION(:), INTENT(INOUT) :: ev
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: ritz_eval
        COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: ritz_evect, temp_mat
        INTEGER :: k, i, n_converged, j,  start, stopv, sizev, proc
        INTEGER, ALLOCATABLE :: picked(:)
        DOUBLE PRECISION :: diff, t_stop, t_start

        evec = 0.0D0

        i = SIZE(qmat,2)
        ALLOCATE(picked(SIZE(ev)))
        picked = 0

        ALLOCATE(ritz_eval(i), ritz_evect(i,i), temp_mat(i,i))
        ritz_eval = 0.0D0; ritz_evect = 0.0D0; temp_mat = hessen(1:i,1:i)

        CALL diagonalize_complex(temp_mat, ritz_eval, ritz_evect)
        DEALLOCATE(temp_mat)

        IF (iam == master .AND. debug > 1) WRITE(*,*) 'Ritz values after iteration: ', i
        IF (iam == master .AND. debug > 1) &
                    WRITE(*,'(A4,A20,A20,A20)') 'i', 'REAL', 'IMAG', 'diff'

        k = 1
        n_converged = 0
        DO j = 1, i
            IF (ABS(ritz_eval(j)) > zero .AND. k <= SIZE(ev)) THEN
                diff = ABS(ritz_eval(j) - ev(k))
                picked(k) = j
                ev(k) = ritz_eval(j)
                IF (diff < tol) n_converged = n_converged + 1
                k = k+1
                IF (iam == master .AND. debug > 1) &
                    WRITE(*,'(I4,3E20.8)') j, REAL(ritz_eval(j)), AIMAG(ritz_eval(j)), diff
            ELSE IF (iam == master .AND. debug > 1 .AND. j < 10) THEN
                WRITE(*,'(I4,2E20.8)') j, REAL(ritz_eval(j)), AIMAG(ritz_eval(j))
            ENDIF
        ENDDO

        t_start = MPI_WTIME()
        ALLOCATE(temp_mat(tmin%my_size, i))
        CALL ZGEMM('N', 'N', tmin%my_size, i, i, DCMPLX(1.0D0,0.0D0), qmat, tmin%my_size, &
                ritz_evect, i, DCMPLX(0.0D0, 0.0D0), temp_mat, tmin%my_size)

        ! Send all ritz_vectors to master.
        DO j = 1, k-1
            IF ( j>i) CYCLE
            IF (iam == master) THEN
                evec(tmin%my_start:tmin%my_stop, j) = temp_mat(:,picked(j))
                DO proc = 1, num_procs - 1
                    start = tmin%all_starts(proc+1)
                    stopv = tmin%all_stops(proc+1)
                    sizev = tmin%all_sizes(proc+1)
                    CALL MPI_RECV(evec(start:stopv, j), sizev, MPI_DOUBLE_COMPLEX, proc, &
                        500, mpi_comm_world, mystatus, ierror)
                ENDDO
            ELSE
                CALL MPI_SEND(temp_mat(:,picked(j)), tmin%my_size, MPI_DOUBLE_COMPLEX, master, &
                        500, mpi_comm_world,ierror)
            ENDIF
        ENDDO

        DEALLOCATE(ritz_eval, ritz_evect, picked, temp_mat)
        t_stop = MPI_WTIME()
        IF (iam == master .AND. debug > 0 ) WRITE(*,*) 'Time for calulating ritz_vectors: ', &
                    i, t_stop-t_start
                
    END SUBROUTINE compute_ritz_pairs_complex



    SUBROUTINE diagonalize_complex(h, eval, evect)
        COMPLEX*16, DIMENSION(:,:) :: h
        ! Variables  needed for dgeev
        CHARACTER*1 :: jobvl, jobvr
        INTEGER :: info, lwork, n, i, j
        COMPLEX*16, ALLOCATABLE, DIMENSION(:) :: wr, work, eval
        COMPLEX*16, ALLOCATABLE, DIMENSION(:,:) :: vr, vl, evect
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: rwork

        eval = 0.0D0; evect = 0.0D0
        n = SIZE(eval)

        jobvl = 'N'; jobvr = 'V'
        lwork = 8*n
        ALLOCATE(wr(n), vr(n,n), vl(1,n), work(lwork))
        ALLOCATE(rwork(2*n))

        CALL ZGEEV(jobvl, jobvr, n, h, n, wr, vl, 1, vr, n, &
            work, lwork, rwork, info)

        ! Sort eigenvalues (lowest first)
        eval(1) = wr(1)
        evect(:,1) = vr(:,1)
        DO i = 2, n
            DO j = i-1, 0, -1
                IF (j == 0) THEN
                    eval(j+1) = wr(i)
                    evect(:,j+1) = vr(:,i)
                ELSE IF (REAL(wr(i)) < REAL(eval(j))) THEN
                    eval(j+1) = eval(j)
                    evect(:,j+1) = evect(:,j)
                ELSE
                    eval(j+1) = wr(i)
                    evect(:,j+1) = vr(:,i)
                    EXIT
                ENDIF
            ENDDO
        ENDDO

        DEALLOCATE(wr, vr, vl, work, rwork)
    END SUBROUTINE diagonalize_complex

END MODULE arnoldi
