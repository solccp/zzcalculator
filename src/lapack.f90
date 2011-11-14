module LAPACK

    implicit none

 private
 public :: syev, my_dsygvd, my_dsyevd, my_dsygv


contains

SUBROUTINE SYEV(JOBZ, UPLO, N, P, LDP, EVALS)
        INTEGER, INTENT(IN) :: N, LDP
        REAL(KIND=8), INTENT(INOUT), DIMENSION( LDP, N ) :: P
        REAL(KIND=8), INTENT(OUT), DIMENSION(N) :: EVALS
        REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: WORK

        CHARACTER, INTENT(IN) :: JOBZ, UPLO

        INTEGER :: INFO,LWORK

        LWORK=-1
        ALLOCATE(WORK(1))
        CALL DSYEV(JOBZ, UPLO, N, P, LDP, EVALS, WORK, LWORK, INFO)
        LWORK=INT(WORK(1))
        DEALLOCATE(WORK)
        ALLOCATE(WORK(LWORK))
        CALL DSYEV(JOBZ, UPLO, N, P, LDP, EVALS, WORK, LWORK, INFO)
        DEALLOCATE(WORK)

END SUBROUTINE

subroutine my_dsyevd(N, A, LDA, W, INFO)
    integer, intent(in) :: N, LDA
    double precision, intent(inout) :: A(:,:)
    double precision, intent(out) :: W(:)
    integer, intent(out) :: INFO
    
    
    double precision, allocatable :: work(:)
    integer, allocatable :: iwork(:)
    double precision :: temp(1)
    integer :: liwork, lwork

    call dsyevd(1, 'V', 'L', N, A, LDA, W, temp, -1, liwork, -1, INFO)
    lwork=floor(temp(1))
    allocate(work(lwork))
    allocate(iwork(liwork))
    
    call dsyevd(1, 'V', 'L', N, A, LDA, W, work, lwork, iwork, liwork, INFO)

    deallocate(work, iwork)



end subroutine

subroutine my_dsygvd(N, A, LDA, B, LDB, W, INFO)
    
    integer, intent(in) :: N, LDA, LDB
    double precision, intent(inout) :: A(:,:), B(:,:)
    double precision, intent(out) :: W(:)
    integer, intent(out) :: INFO
    
    
    double precision, allocatable :: work(:)
    integer, allocatable :: iwork(:)
    double precision :: temp(1)
    integer :: liwork, lwork

    call dsygvd(1, 'V', 'L', N, A, LDA, B, LDB, W, temp, -1, liwork, -1, INFO)
    lwork=floor(temp(1))
    allocate(work(lwork))
    allocate(iwork(liwork))
    
    call dsygvd(1, 'V', 'L', N, A, LDA, B, LDB, W, work, lwork, iwork, liwork, INFO)


    deallocate(work, iwork)

end subroutine
subroutine my_dsygv(N, A, LDA, B, LDB, W, INFO)
    
    integer, intent(in) :: N, LDA, LDB
    double precision, intent(inout) :: A(:,:), B(:,:)
    double precision, intent(out) :: W(:)
    integer, intent(out) :: INFO
    
    
    double precision, allocatable :: work(:)
    double precision :: temp(1)
    integer :: lwork

    call dsygv(1, 'V', 'L', N, A, LDA, B, LDB, W, temp, -1, INFO)
    lwork=floor(temp(1))
    allocate(work(lwork))
    
    call dsygv(1, 'V', 'L', N, A, LDA, B, LDB, W, work, lwork, INFO)


    deallocate(work)

end subroutine



end module
