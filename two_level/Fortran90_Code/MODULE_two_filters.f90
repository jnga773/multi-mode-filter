! This module contains the subroutines used in any of the two-filter
! programs [g2_cross_RK4.f90].

! This module contains the following subroutines:
! - SquareMatrixInverse: Calculates the inverse of an input matrix.
!
! - SquareMatrixZeroEigenvalue: Calculates the eigenvector of a matrix
!                               corresponding to the zero-valued eigenvalue.
!
! - MatrixInverseSS: Uses SquareMatrixInverse to return the steady state array
!                    without having to type out the multiplication.
!
! - SteadyStateMoments: Calculates the steady states of the various operator
!                       moment equations for the atom-filter coupled system.
!
! - G2_InitialConditions: Calculates the initial conditions for the second-
!                         order cross correlation function for two filters.
!
! - G2_InitialValues : Calculates only the initial value of the second-
!                      order correlation function.
!
! - G2_CalculateRK4: Calculates the time evolution of the second-order
!                    cross correlation function for two filters using
!                    Runge-Kutta 4th Order.

! This file must be added to the compilation command when compiling any of the 
! single-filter programs. Eg, with Intel oneAPI:
!     (LINUX): ifort -qmkl -O3 -o [NAME] ./[FILENAME].f90 ./MODULE_two_filters.f90
!   (WINDOWS): ifort /Qmkl /O3 /o [NAME] ./[FILENAME].f90 ./MODULE_two_filters.f90

MODULE TWO_FILTER_SUBROUTINES

CONTAINS

!==============================================================================!
!                          LAPACK MATRIX SUBROUTINES                           !
!==============================================================================!
! Subroutine to calculate the inverse of a matrix USING LAPACK LIBRARY
SUBROUTINE SquareMatrixInverse(N_in, MatrixInv_out)
  ! Import the MKL Library LAPACK95, from which the eigenvalue/eigenvector and
  ! matrix inversion subroutines come from.
  ! MUST INCLUDE THE -mkl OR /Qmkl FLAG AS A COMPILER OPTION IF USING INTEL.
  ! Otherwise you'll have to link it yourself and I don't know how to do that :)

  USE LAPACK95

  ! The subroutines used from LAPACK are:
  ! - zGETRF - Calculates LU-factorisation of a complexmatrix so it can be
  !            inverted by...,
  ! - zGETRI - Calculates the inverse of a complex matrix.

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! Dimension of matrix (N_in x N_in)
  INTEGER, INTENT(IN)                                   :: N_in

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Inverted matrix to be output
  COMPLEX(KIND=8), DIMENSION(N_in, N_in), INTENT(INOUT) :: MatrixInv_out

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
  ! Work space dimension
  INTEGER, PARAMETER                                    :: LWMAX = 300
  INTEGER                                               :: LWORK
  ! Work space array
  COMPLEX(KIND=8), DIMENSION(LWMAX)                     :: WORK
  REAL(KIND=8), DIMENSION(2*N_in)                       :: RWORK
  ! LU-factorisation array
  INTEGER, DIMENSION(N_in)                              :: IPIV
  ! Info IO
  INTEGER                                               :: INFO

  ! Perform LU-factorization of matrix
  CALL zGETRF(N_in, N_in, MatrixInv_out, N_in, IPIV, INFO)
  IF (INFO .NE. 0) THEN
    PRINT*, "zGETRF M failed :( INFO = ", INFO
    STOP
  END IF

  ! Set LWORK to N_in
  LWORK = N_in

  CALL zGETRI(N_in, MatrixInv_out, N_in, IPIV, WORK, LWORK, INFO)

  ! End of subroutine
END SUBROUTINE SquareMatrixInverse

! Subroutine to calculate the inverse of a matrix USING LAPACK LIBRARY
SUBROUTINE SquareMatrixZeroEigenvalue(N_in, Matrix_in, SS_out)
  ! Import the MKL Library LAPACK95, from which the eigenvalue/eigenvector and
  ! matrix inversion subroutines come from.
  ! MUST INCLUDE THE -mkl OR /Qmkl FLAG AS A COMPILER OPTION IF USING INTEL.
  ! Otherwise you'll have to link it yourself and I don't know how to do that :)

  USE LAPACK95

  ! The subroutines used from LAPACK are:
  ! - zGEEV  - Calculates the eigenvalues and eigevectors of a complex matrix.

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! Dimension of matrix (N_in x N_in)
  INTEGER, INTENT(IN)                                :: N_in
  ! Matrix to calculate eigenvalues/vectors from
  COMPLEX(KIND=8), DIMENSION(N_in, N_in), INTENT(IN) :: Matrix_in

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Steady state vector out
  COMPLEX(KIND=8), DIMENSION(N_in), INTENT(OUT)      :: SS_out

  !---------------------!
  !     OTHER STUFF     !
  !---------------------!
  INTEGER                                            :: j, x

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
  ! Work space dimension
  INTEGER, PARAMETER                                 :: LWMAX = 300
  INTEGER                                            :: LWORK
  ! Work space array
  COMPLEX(KIND=8), DIMENSION(LWMAX)                  :: WORK
  REAL(KIND=8), DIMENSION(2*N_in)                    :: RWORK
  ! Eigenvalues of matrix M
  COMPLEX(KIND=8), DIMENSION(N_in)                   :: eigval
  ! S, and S^{-1} matrix for diagonalising eigenvectors
  COMPLEX(KIND=8), DIMENSION(N_in, N_in)             :: S, Sinv
  ! Info IO
  INTEGER                                            :: INFO

  ! Calculate eigenvalues and eigenvectors (Optimal LWORK = 264)
  LWORK = 264
  CALL zGEEV('N', 'V', N_in, Matrix_in, N_in, eigval, S, N_in, S, N_in, WORK, LWORK, RWORK, INFO)
  ! Check convergence
  IF (INFO .GT. 0) THEN
     PRINT*, "zGEEV failed ON eigenvalues/vectors of Mat_OG"
     STOP
  END IF

  SS_out = 0.0d0
  ! Cycle through eigenvalues and, for the eigenvalue that = 0, use that
  ! eigenvector as the steady state
  DO x = 1, N_in
    IF (ABS(REAL(eigval(x))) .LT. 1D-10 .AND. ABS(REAL(eigval(x))) .LT. 1D-10) THEN
      ! Save steady state eigenvector
      DO j = 1, N_in
        SS_out(j) = S(j, x)
      END DO
    END IF
  END DO

  ! End of subroutine
END SUBROUTINE SquareMatrixZeroEigenvalue

!==============================================================================!
!                            STEADY STATE SUBROUTINE                           !
!==============================================================================!
! Subroutine to calculate steady state coupled moments using SquareMatrixInverse
SUBROUTINE MatrixInverseSS(N_in, Matrix_in, Bvec_in, SS_out)
  ! Calculates the steady state of a system of coupled equations,
  !   d/dt < A > = M < A > + B,
  ! by multiplying the non-homogeneous vector with the inverse matrix:
  !   < A >_{ss} = -M^{-1} B.

  ! Parameters
  ! ----------
  ! N_in : integer
  !   The dimension of Matrix_in (N_in x N_in)
  ! Matrix_in : complex matrix, dimension(N_in, N_in)
  !   The matrix which we will invert.
  ! Bvec_in : complex array, dimension(N_in)
  !   The non-homogeneous vector for the system.

  ! Output
  ! ------
  ! SS_out : complex array, dimension(N_in)
  !   Output eigenvector for zero eigenvalue

  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!

  IMPLICIT NONE

  !-------------------------!
  !     INPUT ARGUMENTS     !
  !-------------------------!
  ! Dimension of matrix (N_in x N_in)
  INTEGER, INTENT(IN)                                   :: N_in
  ! Input evolution matrix
  COMPLEX(KIND=8), DIMENSION(N_in, N_in), INTENT(IN)    :: Matrix_in
  ! Non-homogeneous vecotr
  COMPLEX(KIND=8), DIMENSION(N_in), INTENT(IN)          :: Bvec_in

  !--------------------------!
  !     OUTPUT ARGUMENTS     !
  !--------------------------!
  ! Inverted matrix to be output
  COMPLEX(KIND=8), DIMENSION(N_in), INTENT(INOUT)       :: SS_out

  !--------------------------!
  !     SUBROUTINE STUFF     !
  !--------------------------!
  ! Input evolution matrix
  COMPLEX(KIND=8), DIMENSION(N_in, N_in)                :: MatrixInverse

  !============================================================================!
  !                         END OF VARIABLE DELCARATION                        !
  !============================================================================!
  ! Set the matrix to be inverted
  MatrixInverse = Matrix_in

  ! Invert the matrix
  CALL SquareMatrixInverse(N_IN, MatrixInverse)

  ! Calculate steady states
  SS_out = 0.0d0
  SS_out = -MATMUL(MatrixInverse, Bvec_in)

END SUBROUTINE MatrixInverseSS

! Subroutine to calculate the steady states for the atom-filter operator moments
SUBROUTINE SteadyStateMoments(gamma_in, Omega_in, &
                            & epsilon_in, N_in, phase_in, &
                            & w0a_in, kappaa_in, dwa_in, &
                            & w0b_in, kappab_in, dwb_in, &
                            & photon_out, sigma_out, &
                            & f1_out, f1sig_out, &
                            & f2_out, f2sig_out, &
                            & f3_out, f3sig_out, &
                            & f4_out)

  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!

  IMPLICIT NONE

  !---------------!
  !     INPUT     !
  !---------------!
  ! Atomic decay rate
  REAL(KIND=8), INTENT(IN)                  :: gamma_in
  ! Driving amplitude
  REAL(KIND=8), INTENT(IN)                  :: Omega_in

  ! Filter parameter stuff
  ! Percentage of fluorecence aimed at cavity
  REAL(KIND=8), INTENT(IN)                  :: epsilon_in
  ! Number of mode either side of w0, 2N + 1 total mode
  INTEGER, INTENT(IN)                       :: N_in
  ! Phase modulation of mode coupling
  REAL(KIND=8), INTENT(IN)                  :: phase_in
  ! Central mode frequency of the filter cavity, with N mode frequencies either side
  REAL(KIND=8), INTENT(IN)                  :: w0a_in, w0b_in
  ! Cavity linewidth/transmission of cavity mode
  REAL(KIND=8), INTENT(IN)                  :: kappaa_in, kappab_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                  :: dwa_in, dwb_in

  !------------------------------------!
  !     MOMENT EQUATION ARRAY STUFF    !
  !------------------------------------!
  ! Dimension of M matrix
  INTEGER, PARAMETER                        :: N_mat = 3
  ! M matrix (filled as transpose)
  COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)  :: Mat, Mat_OG
  ! Non-homogeneous vector
  COMPLEX(KIND=8), DIMENSION(N_mat)         :: B_vec, B_OG

  ! Integer indices for sigma operators
  INTEGER, PARAMETER                        :: sm = 1, sp = 2, sz = 3
  ! Integer indices for: f, f^{\dagger}, f^{2}, f^{\dagger}^{2} ... etc
  INTEGER, PARAMETER                        :: f = 1, ft = 3, g = 2, gt = 4
  INTEGER, PARAMETER                        :: fg = 1, ftgt = 2, ftf = 3, gtg = 4, ftg = 5, gtf = 6
  INTEGER, PARAMETER                        :: ftfg = 1, gtgf = 2, gtftf = 3, ftgtg = 4

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Steady state photon number
  REAL(KIND=8), DIMENSION(2), INTENT(OUT)                                                 :: photon_out

  ! Steady state arrays
  ! First-order moments: Atomic equations (< \sigma >)
  COMPLEX(KIND=8), DIMENSION(N_mat), INTENT(OUT)                                          :: sigma_out
  ! First-order moments: Cavity (< a >, < f^{\dagger} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 4), INTENT(OUT)                                  :: f1_out
  ! Second-order moments: Cavity and atom (< a \sigma >, < f^{\dagger} \sigma >
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 4, N_mat), INTENT(OUT)                           :: f1sig_out
  ! Second-order moments: Cavity (< f^{\dagger} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 6), INTENT(OUT)                      :: f2_out
  ! Third-order moments: Cavity and atom (< f^{2} \sigma >, < f^{\dagger 2} \sigma >, < f^{\dagger} a \sigma >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 6, N_mat), INTENT(OUT)               :: f2sig_out
  ! Third-order moments: Cavity (< f^{2} f^{\dagger} >, < f^{\dagger 2} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, 4), INTENT(OUT)          :: f3_out
  ! Fourth-order moments: Cavity and atom ( < f^{\dagger} f^{2} \sigma >, < f^{\dagger 2} a \sigma >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, 4, N_mat), INTENT(OUT)   :: f3sig_out
  ! Fourth-order moments: Cavity (< f^{\dagger 2} f^{2} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, -N_in:N_in), INTENT(OUT) :: f4_out

  !----------------------------!
  !     OTHER USEFUL STUFF     !
  !----------------------------!
  ! Integer counter
  INTEGER                                   :: j, k, l, m, x
  ! Imaginary i
  COMPLEX(KIND=8), PARAMETER                :: i = CMPLX(0.0d0, 1.0d0, 8)
  ! pi
  REAL(KIND=8), PARAMETER                   :: pi = 3.1415926535897932384d0
  ! List of Delta values
  REAL(KIND=8), DIMENSION(-N_in:N_in, 2)    :: wl
  ! List of mode dependent cascade coupling values
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2) :: gkl
  ! Blackman window coefficient
  REAL(KIND=8)                              :: blackman
  ! Parameters
  REAL(KIND=8), DIMENSION(2)                :: kappa, w0, dw
  ! Temporary values
  COMPLEX(KIND=8)                           :: moment_out, moment_out2

  !============================================================================!
  !               DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS               !
  !============================================================================!
  kappa(f) = kappaa_in; kappa(g) = kappab_in
  w0(f) = w0a_in; w0(g) = w0b_in
  dw(f) = dwa_in; dw(g) = dwb_in

  !------------------------!
  !     BLOCH MATRIX M     !
  !------------------------!
  Mat_OG = 0.0d0
  ! Row 1: d/dt |g><e|
  Mat_OG(1, 1) = -0.5d0 * gamma_in
  Mat_OG(1, 2) = 0.0d0
  Mat_OG(1, 3) = i * 0.5d0 * Omega_in
  ! Row 2: d/dt |e><g|
  Mat_OG(2, 1) = 0.0d0
  Mat_OG(2, 2) = -0.5d0 * gamma_in
  Mat_OG(2, 3) = -i * 0.5d0 * Omega_in
  ! Row 3: d/dt |e><e| - |g><g|
  Mat_OG(3, 1) = i * Omega_in
  Mat_OG(3, 2) = -i * Omega_in
  Mat_OG(3, 3) = -gamma_in

  !--------------------------------!
  !     NON-HOMOGENEOUS VECTOR     !
  !--------------------------------!
  B_OG = 0.0d0
  B_OG(1) = 0.0d0
  B_OG(2) = 0.0d0
  B_OG(3) = -gamma_in

  !---------------------------------------------!
  !     RESONANCES (wj) AND COUPLINGS (E_j)     !
  !---------------------------------------------!
  ! Allocate array of Delta and gka values
  wl = 0.0d0
  gkl = 0.0d0
  DO j = -N_in, N_in
    IF (N_in == 0) THEN
      ! Cavity A
      wl(j, f) = w0(f)
      gkl(j, f) = DSQRT(epsilon_in * gamma_in * kappa(f))

      ! Cavity B
      wl(j, g) = w0(g)
      gkl(j, g) = DSQRT(epsilon_in * gamma_in * kappa(g))
    ELSE
      ! Blackman window coefficient
      blackman = 1.0d0
      ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
      !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))

      ! Cavity A
      wl(j, f) = w0(f) + DBLE(j) * dw(f)
      ! Mode dependent phase difference
      gkl(j, f) = DSQRT((epsilon_in / DBLE(2*N_in + 1)) * gamma_in * kappa(f)) * blackman * &
                & EXP(i * DBLE(phase_in) * DBLE(j) * pi / DBLE(N_in))

      ! Cavity B
      wl(j, g) = w0(g) + DBLE(j) * dw(g)
      ! Mode dependent phase difference
      gkl(j, g) = DSQRT((epsilon_in / DBLE(2*N_in + 1)) * gamma_in * kappa(g)) * blackman * &
                & EXP(i * DBLE(phase_in) * DBLE(j) * pi / DBLE(N_in))
    END IF
  END DO

  !------------------------------------------!
  !     INITALISE OPERATOR MOMENT ARRAYS     !
  !------------------------------------------!
  ! Steady states
  ! First-order: Cavity
  f1_out = 0.0d0
  ! Second-order: Cavity and Atom
  f1sig_out = 0.0d0
  ! Second-order: Cavity
  f2_out = 0.0d0
  ! Third-order: Cavity and Atom
  f2sig_out = 0.0d0
  ! Third-order: Cavity
  f3_out = 0.0d0
  ! Fourth-order: Cavity and atom
  f3sig_out = 0.0d0
  ! Fourth-order: Cavity
  f4_out = 0.0d0

  photon_out = 0.0d0

  !============================================================================!
  !                       CALCULATE STEADY-STATE MOMENTS                       !
  !============================================================================!
  !---------------------------!
  !     FIRST-ORDER: ATOM     !
  !---------------------------!
  ! Calculate steady states
  CALL MatrixInverseSS(N_mat, Mat_OG, B_OG, sigma_out)

  ! Add OMP clauses
  !$OMP PARALLEL DO PRIVATE(j, Mat, B_vec) COLLAPSE(1)

  ! Cycle through modes
  DO j = -N_in, N_in
    !-----------------------------!
    !     FIRST-ORDER: CAVITY     !
    !-----------------------------!
    !-----------!
    ! < f_{j} > !
    !-----------!
    f1_out(j, f) = -gkl(j, f) * sigma_out(sm)
    f1_out(j, f) = f1_out(j, f) / &
                 & (kappa(f) + i * wl(j, f))
    !---------------------!
    ! < f^{\dagger}_{j} > !
    !---------------------!
    f1_out(j, ft) = -CONJG(gkl(j, f)) * sigma_out(sp)
    f1_out(j, ft) = f1_out(j, ft) / &
                  & (kappa(f) - i * wl(j, f))
    !-----------!
    ! < g_{j} > !
    !-----------!
    f1_out(j, g) = -gkl(j, g) * sigma_out(sm)
    f1_out(j, g) = f1_out(j, g) / &
                 & (kappa(g) + i * wl(j, g))
    !---------------------!
    ! < g^{\dagger}_{j} > !
    !---------------------!
    f1_out(j, gt) = -CONJG(gkl(j, g)) * sigma_out(sp)
    f1_out(j, gt) = f1_out(j, gt) / &
                  & (kappa(g) - i * wl(j, g))
    !---------------------------------------!
    !     SECOND-ORDER: CAVITY AND ATOM     !
    !---------------------------------------!
    !------------------!
    ! < f_{j} \sigma > !
    !------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa(f) + i * wl(j, f))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = 0.0d0
    B_vec(2) = -0.5d0 * gkl(j, f) * (sigma_out(sz) + 1.0d0)
    B_vec(3) = -gamma_in * f1_out(j, f) + &
             & gkl(j, f) * sigma_out(sm)

    ! Calculate steady states
    CALL MatrixInverseSS(N_mat, Mat, B_vec, f1sig_out(j, f, :))

    !----------------------------!
    ! < f^{\dagger}_{j} \sigma > !
    !----------------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa(f) - i * wl(j, f))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -0.5d0 * CONJG(gkl(j, f)) * (sigma_out(sz) + 1.0d0)
    B_vec(2) = 0.0d0
    B_vec(3) = -gamma_in * f1_out(j, ft) + &
             & CONJG(gkl(j, f)) * sigma_out(sp)

    ! Calculate steady states
    CALL MatrixInverseSS(N_mat, Mat, B_vec, f1sig_out(j, ft, :))

    !------------------!
    ! < g_{j} \sigma > !
    !------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa(g) + i * wl(j, g))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = 0.0d0
    B_vec(2) = -0.5d0 * gkl(j, g) * (sigma_out(sz) + 1.0d0)
    B_vec(3) = -gamma_in * f1_out(j, g) + &
             & gkl(j, g) * sigma_out(sm)

    ! Calculate steady states
    CALL MatrixInverseSS(N_mat, Mat, B_vec, f1sig_out(j, g, :))

    !----------------------------!
    ! < g^{\dagger}_{j} \sigma > !
    !----------------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa(g) - i * wl(j, g))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -0.5d0 * CONJG(gkl(j, g)) * (sigma_out(sz) + 1.0d0)
    B_vec(2) = 0.0d0
    B_vec(3) = -gamma_in * f1_out(j, gt) + &
             & CONJG(gkl(j, g)) * sigma_out(sp)

    ! Calculate steady states
    CALL MatrixInverseSS(N_mat, Mat, B_vec, f1sig_out(j, gt, :))

    ! Close j loop
  END DO
  !$OMP END PARALLEL DO

  moment_out = 0.0d0
  moment_out2 = 0.0d0

  ! Add OMP clauses
  !$OMP PARALLEL DO PRIVATE(j, k, Mat, B_vec) REDUCTION(+:moment_out, moment_out2) COLLAPSE(2)

  ! Cycle through modes
  DO k = -N_in, N_in
    DO j = -N_in, N_in
      !------------------------------!
      !     SECOND-ORDER: CAVITY     !
      !------------------------------!
      !-----------------!
      ! < f_{j} g_{k} > !
      !-----------------!
      f2_out(j, k, fg) = -gkl(j, f) * f1sig_out(k, g, sm) + &
                       & (-gkl(k, g)) * f1sig_out(j, f, sm)
      f2_out(j, k, fg) = f2_out(j, k, fg) / &
                       & (kappa(f) + kappa(g) + i * (wl(j, f) + wl(k, g)))

      !-------------------------------------!
      ! < f^{\dagger}_{j} g^{\dagger}_{k} > !
      !-------------------------------------!
      f2_out(j, k, ftgt) = -CONJG(gkl(j, f)) * f1sig_out(k, gt, sp) + &
                         & (-CONJG(gkl(k, g))) * f1sig_out(j, ft, sp)
      f2_out(j, k, ftgt) = f2_out(j, k, ftgt) / &
                         & (kappa(f) + kappa(g) - i * (wl(j, f) + wl(k, g)))

      !---------------------------!
      ! < f^{\dagger}_{j} f_{k} > !
      !---------------------------!
      f2_out(j, k, ftf) = -CONJG(gkl(j, f)) * f1sig_out(k, f, sp) + &
                        & (-gkl(k, f)) * f1sig_out(j, ft, sm)
      f2_out(j, k, ftf) = f2_out(j, k, ftf) / &
                        & ((2.0d0 * kappa(f)) - i * (wl(j, f) - wl(k, f)))
      ! Update photon number
      moment_out = moment_out + f2_out(j, k, ftf)

      !---------------------------!
      ! < g^{\dagger}_{j} g_{k} > !
      !---------------------------!
      f2_out(j, k, gtg) = -CONJG(gkl(j, g)) * f1sig_out(k, g, sp) + &
                        & (-gkl(k, g)) * f1sig_out(j, gt, sm)
      f2_out(j, k, gtg) = f2_out(j, k, gtg) / &
                        & ((2.0d0 * kappa(g)) - i * (wl(j, g) - wl(k, g)))
      ! Update photon number
      moment_out2 = moment_out2 + f2_out(j, k, gtg)

      !---------------------------!
      ! < f^{\dagger}_{j} g_{k} > !
      !---------------------------!
      f2_out(j, k, ftg) = -CONJG(gkl(j, f)) * f1sig_out(k, g, sp) + &
                        & (-gkl(k, g)) * f1sig_out(j, ft, sm)
      f2_out(j, k, ftg) = f2_out(j, k, ftg) / &
                        & (kappa(f) + kappa(g) - i * (wl(j, f) - wl(k, g)))

      !---------------------------!
      ! < g^{\dagger}_{j} f_{k} > !
      !---------------------------!
      f2_out(j, k, gtf) = -CONJG(gkl(j, g)) * f1sig_out(k, f, sp) + &
                        & (-gkl(k, f)) * f1sig_out(j, gt, sm)
      f2_out(j, k, gtf) = f2_out(j, k, gtf) / &
                        & (kappa(f) + kappa(g) - i * (wl(j, g) - wl(k, f)))

      !--------------------------------------!
      !     THIRD-ORDER: CAVITY AND ATOM     !
      !--------------------------------------!
      !------------------------!
      ! < f_{j} g_{k} \sigma > !
      !------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - (kappa(f) + kappa(g) + i * (wl(j, f) + wl(k, g)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = 0.0d0
      B_vec(2) = -0.5d0 * gkl(j, f) * f1sig_out(k, g, sz) + &
               & (-0.5d0) * gkl(j, f) * f1_out(k, g) + &
               & (-0.5d0) * gkl(k, g) * f1sig_out(j, f, sz) + &
               & (-0.5d0) * gkl(k, g) * f1_out(j, f)
      B_vec(3) = -gamma_in * f2_out(j, k, fg) + &
               & gkl(j, f) * f1sig_out(k, g, sm) + &
               & gkl(k, g) * f1sig_out(j, f, sm)

      ! Calculate steady states
      CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, fg, :))

      !--------------------------------------------!
      ! < f^{\dagger}_{j} g^{\dagger}_{k} \sigma > !
      !--------------------------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - (kappa(f) + kappa(g) - i * (wl(j, f) + wl(k, g)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -0.5d0 * CONJG(gkl(j, f)) * f1sig_out(k, gt, sz) + &
               & (-0.5d0) * CONJG(gkl(j, f)) * f1_out(k, gt) + &
               & (-0.5d0) * CONJG(gkl(k, g)) * f1sig_out(j, ft, sz) + &
               & (-0.5d0) * CONJG(gkl(k, g)) * f1_out(j, ft)
      B_vec(2) = 0.0d0
      B_vec(3) = -gamma_in * f2_out(j, k, ftgt) + &
               & CONJG(gkl(j, f)) * f1sig_out(k, gt, sp) + &
               & CONJG(gkl(k, g)) * f1sig_out(j, ft, sp)

      ! Calculate steady states
      CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, ftgt, :))

      !----------------------------------!
      ! < f^{\dagger}_{j} f_{k} \sigma > !
      !----------------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa(f)) - i * (wl(j, f) - wl(k, f)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -0.5d0 * CONJG(gkl(j, f)) * f1sig_out(k, f, sz) + &
               & (-0.5d0) * CONJG(gkl(j, f)) * f1_out(k, f)
      B_vec(2) = -0.5d0 * gkl(k, f) * f1sig_out(j, ft, sz) + &
               & (-0.5d0) * gkl(k, f) * f1_out(j, ft)
      B_vec(3) = -gamma_in * f2_out(j, k, ftf) + &
               & CONJG(gkl(j, f)) * f1sig_out(k, f, sp) + &
               & gkl(k, f) * f1sig_out(j, ft, sm)

      ! Calculate steady states
      CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, ftf, :))

      !----------------------------------!
      ! < g^{\dagger}_{j} g_{k} \sigma > !
      !----------------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa(g)) - i * (wl(j, g) - wl(k, g)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -0.5d0 * CONJG(gkl(j, g)) * f1sig_out(k, g, sz) + &
               & (-0.5d0) * CONJG(gkl(j, g)) * f1_out(k, g)
      B_vec(2) = -0.5d0 * gkl(k, g) * f1sig_out(j, gt, sz) + &
               & (-0.5d0) * gkl(k, g) * f1_out(j, gt)
      B_vec(3) = -gamma_in * f2_out(j, k, gtg) + &
               & CONJG(gkl(j, g)) * f1sig_out(k, g, sp) + &
               & gkl(k, g) * f1sig_out(j, gt, sm)

      ! Calculate steady states
      CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, gtg, :))

      !----------------------------------!
      ! < f^{\dagger}_{j} g_{k} \sigma > !
      !----------------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - (kappa(f) + kappa(g) - i * (wl(j, f) - wl(k, g)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -0.5d0 * CONJG(gkl(j, f)) * f1sig_out(k, g, sz) + &
               & (-0.5d0) * CONJG(gkl(j, f)) * f1_out(k, g)
      B_vec(2) = -0.5d0 * gkl(k, g) * f1sig_out(j, ft, sz) + &
               & (-0.5d0) * gkl(k, g) * f1_out(j, ft)
      B_vec(3) = -gamma_in * f2_out(j, k, ftg) + &
               & CONJG(gkl(j, f)) * f1sig_out(k, g, sp) + &
               & gkl(k, g) * f1sig_out(j, ft, sm)

      ! Calculate steady states
      CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, ftg, :))

      !----------------------------------!
      ! < g^{\dagger}_{j} f_{k} \sigma > !
      !----------------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - (kappa(f) + kappa(g) - i * (wl(j, g) - wl(k, f)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -0.5d0 * CONJG(gkl(j, g)) * f1sig_out(k, f, sz) + &
               & (-0.5d0) * CONJG(gkl(j, g)) * f1_out(k, f)
      B_vec(2) = -0.5d0 * gkl(k, f) * f1sig_out(j, gt, sz) + &
               & (-0.5d0) * gkl(k, f) * f1_out(j, gt)
      B_vec(3) = -gamma_in * f2_out(j, k, gtf) + &
               & CONJG(gkl(j, g)) * f1sig_out(k, f, sp) + &
               & gkl(k, f) * f1sig_out(j, gt, sm)

      ! Calculate steady states
      CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, gtf, :))

      ! Close j loop
    END DO
    ! Close k loop
  END DO
  !$OMP END PARALLEL DO

  ! Update steady state photon number
  photon_out = 0.0d0
  photon_out(f) = REAL(moment_out)
  photon_out(g) = REAL(moment_out2)

  ! Add OMP clauses
  !$OMP PARALLEL DO PRIVATE(j, k, l, Mat, B_vec) COLLAPSE(3)

  ! Cycle through modes
  DO l = -N_in, N_in
    DO k = -N_in, N_in
      DO j = -N_in, N_in
        !-----------------------------!
        !     THIRD-ORDER: CAVITY     !
        !-----------------------------!
        !---------------------------------!
        ! < f^{\dagger}_{j} f_{k} g_{l} > !
        !---------------------------------!
        f3_out(j, k, l, ftfg) = -CONJG(gkl(j, f)) * f2sig_out(k, l, fg, sp) + &
                              & (-gkl(k, f)) * f2sig_out(j, l, ftg, sm) + &
                              & (-gkl(l, g)) * f2sig_out(j, k, ftf, sm)
        f3_out(j, k, l, ftfg) = f3_out(j, k, l, ftfg) / &
                              & ((2.0d0 * kappa(f)) + kappa(g) - i * (wl(j, f) - wl(k, f) - wl(l, g)))

        !---------------------------------!
        ! < g^{\dagger}_{j} g_{k} f_{l} > !
        !---------------------------------!
        f3_out(j, k, l, gtgf) = -CONJG(gkl(j, g)) * f2sig_out(l, k, fg, sp) + &
                              & (-gkl(k, g)) * f2sig_out(j, l, gtf, sm) + &
                              & (-gkl(l, f)) * f2sig_out(j, k, gtg, sm)
        f3_out(j, k, l, gtgf) = f3_out(j, k, l, gtgf) / &
                              & (kappa(f) + (2.0d0 * kappa(g)) - i * (wl(j, g) - wl(k, g) - wl(l, f)))

        !-------------------------------------------!
        ! < g^{\dagger}_{j} f^{\dagger}_{k} f_{l} > !
        !-------------------------------------------!
        f3_out(j, k, l, gtftf) = -CONJG(gkl(j, g)) * f2sig_out(k, l, ftf, sp) + &
                               & (-CONJG(gkl(k, f))) * f2sig_out(j, l, gtf, sp) + &
                               & (-gkl(l, f)) * f2sig_out(k, j, ftgt, sm)
        f3_out(j, k, l, gtftf) = f3_out(j, k, l, gtftf) / &
                               & ((2.0d0 * kappa(f)) + kappa(g) - i * (wl(j, g) + wl(k, f) - wl(l, f)))

        !-------------------------------------------!
        ! < f^{\dagger}_{j} g^{\dagger}_{k} b{l} > !
        !-------------------------------------------!
        f3_out(j, k, l, ftgtg) = -CONJG(gkl(j, f)) * f2sig_out(k, l, gtg, sp) + &
                               & (-CONJG(gkl(k, g))) * f2sig_out(j, l, ftg, sp) + &
                               & (-gkl(l, g)) * f2sig_out(j, k, ftgt, sm)
        f3_out(j, k, l, ftgtg) = f3_out(j, k, l, ftgtg) / &
                               & (kappa(f) + (2.0d0 * kappa(g)) - i * (wl(j, f) + wl(k, g) - wl(l, g)))

        !--------------------------------------!
        !     FOURTH-ORDER: CAVITY AND ATOM    !
        !--------------------------------------!
        !----------------------------------------!
        ! < f^{\dagger}_{j} f_{k} g_{l} \sigma > !
        !----------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa(f)) + kappa(g) - i * (wl(j, f) - wl(k, f) - wl(l, g)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -0.5d0 * CONJG(gkl(j, f)) * f2sig_out(k, l, fg, sz) + &
                 & (-0.5d0) * CONJG(gkl(j, f)) * f2_out(k, l, fg)
        B_vec(2) = -0.5d0 * gkl(k, f) * f2sig_out(j, l, ftg, sz) + &
                 & (-0.5d0) * gkl(k, f) * f2_out(j, l, ftg) + &
                 & (-0.5d0) * gkl(l, g) * f2sig_out(j, k, ftf, sz) + &
                 & (-0.5d0) * gkl(l, g) * f2_out(j, k, ftf)
        B_vec(3) = -gamma_in * f3_out(j, k, l, ftfg) + &
                 & CONJG(gkl(j, f)) * f2sig_out(k, l, fg, sp) + &
                 & gkl(k, f) * f2sig_out(j, l, ftg, sm) + &
                 & gkl(l, g) * f2sig_out(j, k, ftf, sm)

        ! Calculate steady states
        CALL MatrixInverseSS(N_mat, Mat, B_vec, f3sig_out(j, k, l, ftfg, :))

        !----------------------------------------!
        ! < g^{\dagger}_{j} g_{k} f_{l} \sigma > !
        !----------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - (kappa(f) + (2.0d0 * kappa(g)) - i * (wl(j, g) - wl(k, g) - wl(l, f)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -0.5d0 * CONJG(gkl(j, g)) * f2sig_out(l, k, fg, sz) + &
                 & (-0.5d0) * CONJG(gkl(j, g)) * f2_out(l, k, fg)
        B_vec(2) = -0.5d0 * gkl(k, g) * f2sig_out(j, l, gtf, sz) + &
                 & (-0.5d0) * gkl(k, g) * f2_out(j, l, gtf) + &
                 & (-0.5d0) * gkl(l, f) * f2sig_out(j, k, gtg, sz) + &
                 & (-0.5d0) * gkl(l, f) * f2_out(j, k, gtg)
        B_vec(3) = -gamma_in * f3_out(j, k, l, gtgf) + &
                 & CONJG(gkl(j, g)) * f2sig_out(l, k, fg, sp) + &
                 & gkl(k, g) * f2sig_out(j, l, gtf, sm) + &
                 & gkl(l, f) * f2sig_out(j, k, gtg, sm)

        ! Calculate steady states
        CALL MatrixInverseSS(N_mat, Mat, B_vec, f3sig_out(j, k, l, gtgf, :))

        !-------------------------------------------!
        ! < g^{\dagger}_{j} f^{\dagger}_{k} f_{l} > !
        !-------------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa(f)) + kappa(g) - i * (wl(j, g) + wl(k, f) - wl(l, f)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -0.5d0 * CONJG(gkl(j, g)) * f2sig_out(k, l, ftf, sz) + &
                 & (-0.5d0) * CONJG(gkl(j, g)) * f2_out(k, l, ftf) + &
                 & (-0.5d0) * CONJG(gkl(k, f)) * f2sig_out(j, l, gtf, sz) + &
                 & (-0.5d0) * CONJG(gkl(k, f)) * f2_out(j, l, gtf)
        B_vec(2) = -0.5d0 * gkl(l, f) * f2sig_out(k, j, ftgt, sz) + &
                 & (-0.5d0) * gkl(l, f) * f2_out(k, j, ftgt)
        B_vec(3) = -gamma_in * f3_out(j, k, l, gtftf) + &
                 & CONJG(gkl(j, g)) * f2sig_out(k, l, ftf, sp) + &
                 & CONJG(gkl(k, f)) * f2sig_out(j, l, gtf, sp) + &
                 & gkl(l, f) * f2sig_out(k, j, ftgt, sm)

        ! Calculate steady states
        CALL MatrixInverseSS(N_mat, Mat, B_vec, f3sig_out(j, k, l, gtftf, :))

        !-------------------------------------------!
        ! < f^{\dagger}_{j} g^{\dagger}_{k} g_{l} > !
        !-------------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - (kappa(f) + (2.0d0 * kappa(g)) - i * (wl(j, f) + wl(k, g) - wl(l, g)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -0.5d0 * CONJG(gkl(j, f)) * f2sig_out(k, l, gtg, sz) + &
                 & (-0.5d0) * CONJG(gkl(j, f)) * f2_out(k, l, gtg) + &
                 & (-0.5d0) * CONJG(gkl(k, g)) * f2sig_out(j, l, ftg, sz) + &
                 & (-0.5d0) * CONJG(gkl(k, g)) * f2_out(j, l, ftg)
        B_vec(2) = -0.5d0 * gkl(l, g) * f2sig_out(j, k, ftgt, sz) + &
                 & (-0.5d0) * gkl(l, g) * f2_out(j, k, ftgt)
        B_vec(3) = -gamma_in * f3_out(j, k, l, ftgtg) + &
                 & CONJG(gkl(j, f)) * f2sig_out(k, l, gtg, sp) + &
                 & CONJG(gkl(k, g)) * f2sig_out(j, l, ftg, sp) + &
                 & gkl(l, g) * f2sig_out(j, k, ftgt, sm)

        ! Calculate steady states
        CALL MatrixInverseSS(N_mat, Mat, B_vec, f3sig_out(j, k, l, ftgtg, :))

        ! Close j loop
      END DO
      ! Close k loop
    END DO
    ! Close l loop
  END DO
  !$OMP END PARALLEL DO

  ! Add OMP clauses
  !$OMP PARALLEL DO PRIVATE(j, k, l, m) COLLAPSE(4)

  ! Cycle through modes
  DO m = -N_in, N_in
    DO l = -N_in, N_in
      DO k = -N_in, N_in
        DO j = -N_in, N_in
          !------------------------------!
          !     FOURTH-ORDER: CAVITY     !
          !------------------------------!
          !-------------------------------------------------!
          ! < f^{\dagger}_{j} g^{\dagger}_{k} g_{l} f_{m} > !
          !-------------------------------------------------!
          f4_out(j, k, l, m) = -CONJG(gkl(j, f)) * f3sig_out(k, l, m, gtgf, sp) + &
                             & (-CONJG(gkl(k, g))) * f3sig_out(j, m, l, ftfg, sp) + &
                             & (-gkl(l, g)) * f3sig_out(k, j, m, gtftf, sm) + &
                             & (-gkl(m, f)) * f3sig_out(j, k, l, ftgtg, sm)
          f4_out(j, k, l, m) = f4_out(j, k, l, m) / &
                             & (((2.0d0 * kappa(f)) + (2.0d0 * kappa(g))) - i * (wl(j, f) + wl(k, g)) + i * (wl(l, g) + wl(m, f)))
          ! Close j loop
        END DO
        ! Close k loop
      END DO
      ! Close l loop
    END DO
    ! Close m loop
  END DO
  !$OMP END PARALLEL DO

END SUBROUTINE SteadyStateMoments

!==============================================================================!
!                          G2 CORRELATION SUBROUTINES                          !
!==============================================================================!
! Subroutine to calculate the initial conditions for the auto-correlations
SUBROUTINE G2_InitialConditions(gamma_in, Omega_in, &
                              & epsilon_in, N_in, phase_in, &
                              & w0a_in, kappaa_in, dwa_in, &
                              & w0b_in, kappab_in, dwb_in, &
                              & photon_ss_out, B_OG_out, &
                              & sigma_out, f1_out, f1sig_out, f2_out)

  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!

  IMPLICIT NONE

  !---------------!
  !     INPUT     !
  !---------------!
  ! Atomic decay rate
  REAL(KIND=8), INTENT(IN)                   :: gamma_in
  ! Driving amplitude
  REAL(KIND=8), INTENT(IN)                   :: Omega_in

  ! Filter parameter stuff
  ! Percentage of fluorecence aimed at cavity
  REAL(KIND=8), INTENT(IN)                   :: epsilon_in
  ! Number of mode either side of w0, 2N + 1 total mode
  INTEGER, INTENT(IN)                        :: N_in
  ! Phase modulation of mode coupling
  REAL(KIND=8), INTENT(IN)                   :: phase_in

  ! Central mode frequency of the filter cavity, with N mode frequencies either side
  REAL(KIND=8), INTENT(IN)                   :: w0a_in, w0b_in
  ! Cavity linewidth/transmission of cavity mode
  REAL(KIND=8), INTENT(IN)                   :: kappaa_in, kappab_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                   :: dwa_in, dwb_in

  !------------------------------------!
  !     MOMENT EQUATION ARRAY STUFF    !
  !------------------------------------!
  ! Dimension of M matrix
  INTEGER, PARAMETER                         :: N_mat = 3
  ! M matrix (filled as transpose)
  COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)   :: Mat, Mat_OG

  ! Integer indices for sigma operators
  INTEGER, PARAMETER                         :: sm = 1, sp = 2, sz = 3
  ! Integer indices for: f, f^{\dagger}, f^{2}, f^{\dagger}^{2} ... etc
  INTEGER, PARAMETER                         :: f = 1, ft = 2
  INTEGER, PARAMETER                         :: fg = 1, ftgt = 2, ftf = 3, gtg = 4, ftg = 5, gtf = 6
  INTEGER, PARAMETER                         :: ftfg = 1, gtgf = 2, gtftf = 3, ftgtg = 4

  ! Steady state arrays
  ! First-order moments: Atomic equations (< \sigma >)
  COMPLEX(KIND=8), DIMENSION(N_mat)                                 :: sigma_ss
  ! First-order moments: Cavity (< a >, < f^{\dagger} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 4)                         :: f1_ss
  ! Second-order moments: Cavity and atom (< a \sigma >, < f^{\dagger} \sigma >
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 4, N_mat)                  :: f1sig_ss
  ! Second-order moments: Cavity (< f^{\dagger} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 6)             :: f2_ss
  ! Third-order moments: Cavity and atom (< f^{2} \sigma >, < f^{\dagger 2} \sigma >, < f^{\dagger} a \sigma >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 6, N_mat)      :: f2sig_ss
  ! Third-order moments: Cavity (< f^{2} f^{\dagger} >, < f^{\dagger 2} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, 4) :: f3_ss
  ! Fourth-order moments: Cavity and atom ( < f^{\dagger} f^{2} \sigma >, < f^{\dagger 2} a \sigma >)
  COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE            :: f3sig_ss
  ! Fourth-order moments: Cavity (< f^{\dagger 2} f^{2} >)
  COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE               :: f4_ss

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Steady state photon number
  REAL(KIND=8), DIMENSION(2), INTENT(OUT)                           :: photon_ss_out
  ! Non-homogeneous vector
  COMPLEX(KIND=8), DIMENSION(N_mat), INTENT(OUT)                    :: B_OG_out

  ! Time integration arrays
  ! First-order moments: Atomic equations (< \sigma >)
  COMPLEX(KIND=8), DIMENSION(N_mat), INTENT(OUT)                    :: sigma_out
  ! First-order moments: Cavity (< f >, < f^{\dagger} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2), INTENT(OUT)            :: f1_out
  ! Second-order moments: Cavity and atom (< a \sigma >, < f^{\dagger} \sigma >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N_mat), INTENT(OUT)     :: f1sig_out
  ! Second-order moments: Cavity (< f^{\dagger} f >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in), INTENT(OUT)   :: f2_out

  !----------------------------!
  !     OTHER USEFUL STUFF     !
  !----------------------------!
  ! Integer counter
  INTEGER                               :: j, k, l, m
  ! Imaginary i
  COMPLEX(KIND=8), PARAMETER            :: i = CMPLX(0.0d0, 1.0d0, 8)

  !------------------------------------------!
  !     INITALISE OPERATOR MOMENT ARRAYS     !
  !------------------------------------------!
  ! Steady states
  ! First-order: Cavity
  f1_ss = 0.0d0
  ! Second-order: Cavity and Atom
  f1sig_ss = 0.0d0
  ! Second-order: Cavity
  f2_ss = 0.0d0
  ! Third-order: Cavity and Atom
  f2sig_ss = 0.0d0
  ! Third-order: Cavity
  f3_ss = 0.0d0
  ! Fourth-order: Cavity and atom
  ALLOCATE(f3sig_ss(-N_in:N_in, -N_in:N_in, -N_in:N_in, 4, N_mat)); f3sig_ss = 0.0d0
  ! Fourth-order: Cavity
  ALLOCATE(f4_ss(-N_in:N_in, -N_in:N_in, -N_in:N_in, -N_in:N_in)); f4_ss = 0.0d0

  !============================================================================!
  !                       CALCULATE STEADY-STATE MOMENTS                       !
  !============================================================================!
  CALL SteadyStateMoments(gamma_in, Omega_in, &
                        & epsilon_in, N_in, phase_in, &
                        & w0a_in, kappaa_in, dwa_in, &
                        & w0b_in, kappab_in, dwb_in, &
                        & photon_ss_out, sigma_ss, &
                        & f1_ss, f1sig_ss, &
                        & f2_ss, f2sig_ss, &
                        & f3_ss, f3sig_ss, &
                        & f4_ss)

  !============================================================================!
  !       CALCULATE SECOND-ORDER CORRELATION FUNCTION INITIAL CONDITIONS       !
  !============================================================================!
  ! Set initial conditions and non-homogeneous vector
  ! < f^{\dagger}_{j}(0) \sigma(\tau = 0) f_{m}(0) > =
  !                          < f^{\dagger}_{j} f_{m} \sigma >_{ss},
  ! < f^{\dagger}_{j}(0) f^{\dagger}_{k}(\tau = 0) f_{m}(0) > =
  !                          < f^{\dagger}_{j} f^{\dagger}_{k} f_{m} >_{ss},
  ! < f^{\dagger}_{j}(0) f_{l}(\tau = 0) f_{m}(0) > =
  !                          < f^{\dagger}_{j} f_{l} f_{m} >_{ss},
  ! < f^{\dagger}_{j}(0) f^{\dagger}_{k} \sigma(\tau = 0) f_{m}(0) =
  !                         < f^{\dagger}_{j} f^{\dagger}_{k} f_{m} \sigma >_{ss},
  ! < f^{\dagger}_{j}(0) f_{l} \sigma(\tau = 0) f_{m}(0) =
  !                         < f^{\dagger}_{j} f_{l} f_{m} \sigma >_{ss},
  ! and
  ! < f^{\dagger}_{j}(0) f^{\dagger}_{k} f_{l}(\tau = 0)  f_{m}(0) > =
  !                          < f^{\dagger}_{j} f^{\dagger}_{k} f_{l} f_{m} >_{ss}.

  sigma_out = 0.0d0
  f1_out = 0.0d0
  f1sig_out = 0.0d0
  f2_out = 0.0d0
  B_OG_out = 0.0d0
  ! Cycle through modes
  DO m = -N_in, N_in
    DO j = -N_in, N_in
      ! First-order: Atom
      sigma_out(:) = sigma_out(:) + f2sig_ss(j, m, ftf, :)

      DO k = -N_in, N_in
        ! First-order: Cavity
        f1_out(k, f) = f1_out(k, f) + f3_ss(j, m, k, ftfg)
        f1_out(k, ft) = f1_out(k, ft) + f3_ss(k, j, m, gtftf)
        ! Second-order: Cavity and atom
        f1sig_out(k, f, :) = f1sig_out(k, f, :) + f3sig_ss(j, m, k, ftfg, :)
        f1sig_out(k, ft, :) = f1sig_out(k, ft, :) + f3sig_ss(k, j, m, gtftf, :)

        DO l = -N_in, N_in
          ! Second-order: cavity
          f2_out(k, l) = f2_out(k, l) + f4_ss(j, k, l, m)

          ! Close l loop
        END DO
        ! Close k loop
      END DO
      ! Non homogeneous vector
      B_OG_out(3) = B_OG_out(3) - gamma_in * f2_ss(j, m, ftf)

      ! Close j loop
    END DO
    ! Close m loop
  END DO

END SUBROUTINE G2_InitialConditions

! Subroutine to calculate the initial value of the correlation function
SUBROUTINE G2_InitialValue(gamma_in, Omega_in, &
                         & epsilon_in, N_in, phase_in, &
                         & w0a_in, kappaa_in, dwa_in, &
                         & w0b_in, kappab_in, dwb_in, &
                         & g2_initial_out)

  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!

  IMPLICIT NONE

  !---------------!
  !     INPUT     !
  !---------------!
  ! Atomic decay rate
  REAL(KIND=8), INTENT(IN)                   :: gamma_in
  ! Driving amplitude
  REAL(KIND=8), INTENT(IN)                   :: Omega_in

  ! Filter parameter stuff
  ! Percentage of fluorecence aimed at cavity
  REAL(KIND=8), INTENT(IN)                   :: epsilon_in
  ! Number of mode either side of w0, 2N + 1 total mode
  INTEGER, INTENT(IN)                        :: N_in
  ! Phase modulation of mode coupling
  REAL(KIND=8), INTENT(IN)                   :: phase_in
  ! Central mode frequency of the filter cavity, with N mode frequencies either side
  REAL(KIND=8), INTENT(IN)                   :: w0a_in, w0b_in
  ! Cavity linewidth/transmission of cavity mode
  REAL(KIND=8), INTENT(IN)                   :: kappaa_in, kappab_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                   :: dwa_in, dwb_in

  !------------------------------------!
  !     MOMENT EQUATION ARRAY STUFF    !
  !------------------------------------!
  ! Dimension of M matrix
  INTEGER, PARAMETER                         :: N_mat = 3
  ! M matrix (filled as transpose)
  COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)   :: Mat, Mat_OG

  ! Integer indices for sigma operators
  INTEGER, PARAMETER                         :: sm = 1, sp = 2, sz = 3
  ! Integer indices for: f, f^{\dagger}, f^{2}, f^{\dagger}^{2} ... etc
  INTEGER, PARAMETER                         :: f = 1, ft = 2
  INTEGER, PARAMETER                         :: fg = 1, ftgt = 2, ftf = 3, gtg = 4, ftg = 5, gtf = 6
  INTEGER, PARAMETER                         :: ftfg = 1, gtgf = 2, gtftf = 3, ftgtg = 4

  ! Steady state arrays
  ! First-order moments: Atomic equations (< \sigma >)
  COMPLEX(KIND=8), DIMENSION(N_mat)                                 :: sigma_ss
  ! First-order moments: Cavity (< a >, < f^{\dagger} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 4)                         :: f1_ss
  ! Second-order moments: Cavity and atom (< a \sigma >, < f^{\dagger} \sigma >
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 4, N_mat)                  :: f1sig_ss
  ! Second-order moments: Cavity (< f^{\dagger} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 6)             :: f2_ss
  ! Third-order moments: Cavity and atom (< f^{2} \sigma >, < f^{\dagger 2} \sigma >, < f^{\dagger} a \sigma >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 6, N_mat)      :: f2sig_ss
  ! Third-order moments: Cavity (< f^{2} f^{\dagger} >, < f^{\dagger 2} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, 4) :: f3_ss
  ! Fourth-order moments: Cavity and atom ( < f^{\dagger} f^{2} \sigma >, < f^{\dagger 2} a \sigma >)
  COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE            :: f3sig_ss
  ! Fourth-order moments: Cavity (< f^{\dagger 2} f^{2} >)
  COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE               :: f4_ss

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Initial correlation value g^{(2)}(\tau = 0)
  REAL(KIND=8), DIMENSION(2)                                         :: photon_ss
  ! Initial correlation value g^{(2)}(\tau = 0)
  REAL(KIND=8), INTENT(OUT)                                          :: g2_initial_out

  !----------------------------!
  !     OTHER USEFUL STUFF     !
  !----------------------------!
  ! Integer counter
  INTEGER                               :: j, k, l, m
  ! Imaginary i
  COMPLEX(KIND=8), PARAMETER            :: i = CMPLX(0.0d0, 1.0d0, 8)
  ! Temporary value
  COMPLEX(KIND=8)                       :: moment_out

  !------------------------------------------!
  !     INITALISE OPERATOR MOMENT ARRAYS     !
  !------------------------------------------!
  ! Steady states
  ! First-order: Cavity
  f1_ss = 0.0d0
  ! Second-order: Cavity and Atom
  f1sig_ss = 0.0d0
  ! Second-order: Cavity
  f2_ss = 0.0d0
  ! Third-order: Cavity and Atom
  f2sig_ss = 0.0d0
  ! Third-order: Cavity
  f3_ss = 0.0d0
  ! Fourth-order: Cavity and atom
  ALLOCATE(f3sig_ss(-N_in:N_in, -N_in:N_in, -N_in:N_in, 4, N_mat)); f3sig_ss = 0.0d0
  ! Fourth-order: Cavity
  ALLOCATE(f4_ss(-N_in:N_in, -N_in:N_in, -N_in:N_in, -N_in:N_in)); f4_ss = 0.0d0

  !============================================================================!
  !                        CALCULATE STEADY-STATE MOMENTS                       !
  !===========================================================================!
  CALL SteadyStateMoments(gamma_in, Omega_in, &
                        & epsilon_in, N_in, phase_in, &
                        & w0a_in, kappaa_in, dwa_in, &
                        & w0b_in, kappab_in, dwb_in, &
                        & photon_ss, sigma_ss, &
                        & f1_ss, f1sig_ss, &
                        & f2_ss, f2sig_ss, &
                        & f3_ss, f3sig_ss, &
                        & f4_ss)

  !============================================================================!
  !          CALCULATE SECOND-ORDER CORRELATION FUNCTION INITIAL VALUE         !
  !============================================================================!
  ! Initialise values
  moment_out = 0.0d0
  g2_initial_out = 0.0d0

  ! Cycle through modes
  DO m = -N_in, N_in
    DO l = -N_in, N_in
      DO k = -N_in, N_in
        DO j = -N_in, N_in
          ! Sum over all fourth-order moments
          moment_out = moment_out + f4_ss(j, k, l, m)

          ! Close j loop
        END DO
        ! Close k loop
      END DO
      ! Close l loop
    END DO
    ! Close m loop
  END DO

  ! Normalise by steady state photon number
  g2_initial_out = REAL(moment_out) / (photon_ss(1) * photon_ss(2))

END SUBROUTINE G2_InitialValue

! Subroutine to calculate the time evolution of the g2 correlation
SUBROUTINE G2_CalculateRK4(gamma_in, Omega_in, &
                         & epsilon_in, N_in, phase_in, &
                         & w0a_in, kappaa_in, dwa_in, &
                         & w0b_in, kappab_in, dwb_in, &
                         & dt_in, tau_steps_in, &
                         & g2_array_out, WRITE_DATA_IN, filename_data_in)

  !============================================================================!
  !                   DEFINING AND DECLARING VARIABLES/ARRAYS                  !
  !============================================================================!

  IMPLICIT NONE

  !---------------!
  !     INPUT     !
  !---------------!
  ! Atomic decay rate
  REAL(KIND=8), INTENT(IN)                             :: gamma_in
  ! Driving amplitude
  REAL(KIND=8), INTENT(IN)                             :: Omega_in

  ! Filter parameter stuff
  ! Percentage of fluorecence aimed at cavity
  REAL(KIND=8), INTENT(IN)                             :: epsilon_in
  ! Number of mode either side of w0, 2N + 1 total mode
  INTEGER, INTENT(IN)                                  :: N_in
  ! Phase modulation of mode coupling
  REAL(KIND=8), INTENT(IN)                             :: phase_in
  
  ! Central mode frequency of the filter cavity, with N mode frequencies either side
  REAL(KIND=8), INTENT(IN)                             :: w0a_in, w0b_in
  ! Cavity linewidth/transmission of cavity mode
  REAL(KIND=8), INTENT(IN)                             :: kappaa_in, kappab_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                             :: dwa_in, dwb_in

  ! Time stuff
  ! Time step
  REAL(KIND=8), INTENT(IN)                             :: dt_in
  ! Maxi_inmum number of steps to integrate for
  INTEGER, INTENT(IN)                                  :: tau_steps_in

  ! Data stuff
  ! Boolean for writing data
  LOGICAL, INTENT(IN)                                  :: WRITE_DATA_IN
  ! Filename for writing data to
  CHARACTER(LEN=*), INTENT(IN)                         :: filename_data_in

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Data array
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: g2_array_out

  !------------------------------------!
  !     MOMENT EQUATION ARRAY STUFF    !
  !------------------------------------!
  ! Dimension of M matrix
  INTEGER, PARAMETER                                   :: N_mat = 3
  ! M matrix (filled as transpose)
  COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)             :: Mat, Mat_OG
  ! Non-homogeneous vector
  COMPLEX(KIND=8), DIMENSION(N_mat)                    :: B_OG, B_vec

  ! Integer indices for sigma operators
  INTEGER, PARAMETER                                   :: sm = 1, sp = 2, sz = 3
  ! Integer indices for: f, f^{\dagger}, f^{2}, f^{\dagger}^{2} ... etc
  INTEGER, PARAMETER                                   :: f = 1, ft = 2

  ! Time integration arrays
  ! First-order moments: Atomic equations (< \sigma >)
  COMPLEX(KIND=8), DIMENSION(N_mat)                    :: sigma
  COMPLEX(KIND=8), DIMENSION(N_mat)                    :: k1_sigma, k2_sigma, k3_sigma, k4_sigma
  ! First-order moments: Cavity (< a >, < f^{\dagger} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2)            :: f1
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2)            :: k1_f1, k2_f1, k3_f1, k4_f1
  ! Second-order moments: Cavity and atom (< a \sigma >, < f^{\dagger} \sigma >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N_mat)     :: f1sig
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N_mat)     :: k1_f1sig, k2_f1sig, k3_f1sig, k4_f1sig
  ! Second-order moments: Cavity (< f^{\dagger} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in)   :: f2
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in)   :: k1_f2, k2_f2, k3_f2, k4_f2

  !----------------------------!
  !     OTHER USEFUL STUFF     !
  !----------------------------!
  ! Time step integer
  INTEGER                                              :: t
  ! Integer counter
  INTEGER                                              :: j, k, l, m, x
  ! Sample rate for state populations
  INTEGER                                              :: sample_rate
  ! Imaginary i
  COMPLEX(KIND=8), PARAMETER                           :: i = CMPLX(0.0d0, 1.0d0, 8)
  ! pi
  REAL(KIND=8), PARAMETER                              :: pi = 3.1415926535897932384d0
  ! 1 / 6
  REAL(KIND=8), PARAMETER                              :: xis = 1.0d0 / 6.0d0
  ! kappa, w0, and dw values
  REAL(KIND=8)                                         :: kappa, dw, w0
  ! List of Delta values
  REAL(KIND=8), DIMENSION(-N_in:N_in)                  :: wl
  ! List of mode dependent cascade coupling values
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in)               :: gkl
  ! Blackman window coefficient
  REAL(KIND=8)                                         :: blackman
  ! Steady state photon number
  REAL(KIND=8), DIMENSION(2)                           :: photon_ss
  ! Complex data
  COMPLEX(KIND=8)                                      :: moment_out

  !============================================================================!
  !               DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS               !
  !============================================================================!
  kappa = kappab_in
  w0 = w0b_in
  dw = dwb_in

  !------------------------!
  !     BLOCH MATRIX M     !
  !------------------------!
  Mat_OG = 0.0d0
  ! Row 1: d/dt |g><e|
  Mat_OG(1, 1) = -0.5d0 * gamma_in
  Mat_OG(1, 2) = 0.0d0
  Mat_OG(1, 3) = i * 0.5d0 * Omega_in
  ! Row 2: d/dt |e><g|
  Mat_OG(2, 1) = 0.0d0
  Mat_OG(2, 2) = -0.5d0 * gamma_in
  Mat_OG(2, 3) = -i * 0.5d0 * Omega_in
  ! Row 3: d/dt |e><e| - |g><g|
  Mat_OG(3, 1) = i * Omega_in
  Mat_OG(3, 2) = -i * Omega_in
  Mat_OG(3, 3) = -gamma_in

  !--------------------------------!
  !     NON-HOMOGENEOUS VECTOR     !
  !--------------------------------!
  B_OG = 0.0d0
  B_OG(1) = 0.0d0
  B_OG(2) = 0.0d0
  B_OG(3) = -gamma_in

  !---------------------------------------------!
  !     RESONANCES (wj) AND COUPLINGS (E_j)     !
  !---------------------------------------------!
  ! Allocate array of Delta and gka values
  wl = 0.0d0
  gkl = 0.0d0
  DO j = -N_in, N_in
    IF (N_in == 0) THEN
      wl(j) = w0
      gkl(j) = DSQRT(epsilon_in * gamma_in * kappa)
    ELSE
      wl(j) = w0 + DBLE(j) * dw
      ! Blackman window coefficient
      blackman = 1.0d0
      ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
      !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))
      ! Mode dependent phase difference
      gkl(j) = DSQRT((epsilon_in / DBLE(2*N_in + 1)) * gamma_in * kappa) * blackman * &
             & EXP(i * DBLE(phase_in) * DBLE(j) * pi / DBLE(N_in))
    END IF
  END DO

  !------------------------------------------!
  !     INITALISE OPERATOR MOMENT ARRAYS     !
  !------------------------------------------!

  ! Time integration
  ! First-order: Cavity
  f1 = 0.0d0
  k1_f1 = 0.0d0; k2_f1 = 0.0d0; k3_f1 = 0.0d0; k4_f1 = 0.0d0
  ! Second-order: Cavity and Atom
  f1sig = 0.0d0
  k1_f1sig = 0.0d0; k2_f1sig = 0.0d0; k3_f1sig = 0.0d0; k4_f1sig = 0.0d0
  ! Second-order: Cavity
  f2 = 0.0d0
  k1_f2 = 0.0d0; k2_f2 = 0.0d0; k3_f2 = 0.0d0; k4_f2 = 0.0d0

  ! Data
  ALLOCATE(g2_array_out(0:tau_steps_in)); g2_array_out = 0.0d0
  photon_ss = 0.0d0

  !============================================================================!
  !                        CALCULATE INITIAL CONDITIONS                        !
  !============================================================================!
  CALL G2_InitialConditions(gamma_in, Omega_in, &
                          & epsilon_in, N_in, phase_in, &
                          & w0a_in, kappaa_in, dwa_in, &
                          & w0b_in, kappab_in, dwb_in, &
                          & photon_ss, B_OG, &
                          & sigma, f1, f1sig, f2)

  !============================================================================!
  !                 CALCULATE SECOND-ORDER CORRELATION FUNCTION                !
  !============================================================================!
  ! Calculate the sample rate for writing data to the file
  IF (tau_steps_in > 100000) THEN
    sample_rate = NINT(DBLE(tau_steps_in) / 1d5)
  ELSE
    sample_rate = 1
  END IF

  ! If WRITE_DATA_IN is TRUE, open file to write data to
  IF (WRITE_DATA_IN .EQV. .TRUE.) THEN
    ! Open file to write time and data to
    OPEN(UNIT=4, FILE=filename_data_in, STATUS='REPLACE', ACTION='WRITE', RECL=4000)
  END IF

  ! Cycle through time steps
  DO t = 0, tau_steps_in
    !============================================================================!
    !                          CALCULATE AND WRITE DATA                          !
    !============================================================================!
    !-------------------------------!
    !     CALCULATE DATA STATES     !
    !-------------------------------!
    ! Grab correlation value
    moment_out = 0.0d0
    DO k = -N_in, N_in
      DO j = -N_in, N_in
        moment_out = moment_out + f2(j, k)
      END DO
    END DO

    ! Normalise correlation by steady-state photon number
    IF (photon_ss(1) .NE. 0.0 .AND. photon_ss(2) .NE. 0.0) THEN
      moment_out = moment_out / (photon_ss(1) * photon_ss(2))
    END IF

    !-----------------------!
    !     WRITE TO FILE     !
    !-----------------------!
    g2_array_out(t) = REAL(MOMENT_OUT)
    ! Second-order correlation function
    ! If WRITE_DATA_IN is TRUE, write data to file
    IF (WRITE_DATA_IN .EQV. .TRUE.) THEN
      ! If t_max is really big, only take a sample of results to write to file
      ! so file size isn't huge-mongous.
      IF (MOD(t, sample_rate) == 0) THEN
        WRITE(4, *) dt_in * DBLE(t), REAL(moment_out)
      END IF
    END IF

    !============================================================================!
    !                  CALCULATE USING FOURTH-ORDER RUNGE-KUTTA                  !
    !============================================================================!
    !----------------------------------------!
    !     INITIALISE RUNGE-KUTTA VECTORS     !
    !----------------------------------------!
    k1_sigma = 0.0d0; k2_sigma = 0.0d0; k3_sigma = 0.0d0; k4_sigma = 0.0d0
    k1_f1 = 0.0d0; k2_f1 = 0.0d0; k3_f1 = 0.0d0; k4_f1 = 0.0d0
    ! Second-order
    k1_f1sig = 0.0d0; k2_f1sig = 0.0d0; k3_f1sig = 0.0d0; k4_f1sig = 0.0d0
    k1_f2 = 0.0d0; k2_f2 = 0.d0; k3_f2 = 0.0d0; k4_f2 = 0.0d0

    !---------------------------!
    !     FIRST-ORDER: ATOM     !
    !---------------------------!
    ! Calculate Runge-Kutta vectors
    k1_sigma = dt_in * (MATMUL(Mat_OG, sigma) + B_OG)
    k2_sigma = dt_in * (MATMUL(Mat_OG, (sigma + 0.5d0 * k1_sigma)) + B_OG)
    k3_sigma = dt_in * (MATMUL(Mat_OG, (sigma + 0.5d0 * k2_sigma)) + B_OG)
    k4_sigma = dt_in * (MATMUL(Mat_OG, (sigma + k3_sigma)) + B_OG)

    ! Cycle through modes
    DO j = -N_in, N_in
      !-----------------------------!
      !     FIRST-ORDER: CAVITY     !
      !-----------------------------!
      !-----------!
      ! < f_{j} > !
      !-----------!
      k1_f1(j, f) = -dt_in * (kappa + i * wl(j)) * f1(j, f) + &
                  & (-dt_in) * gkl(j) * sigma(sm)
      k2_f1(j, f) = -dt_in * (kappa + i * wl(j)) * (f1(j, f) + 0.5d0 * k1_f1(j, f)) + &
                  & (-dt_in) * gkl(j) * (sigma(sm) + 0.5d0 * k1_sigma(sm))
      k3_f1(j, f) = -dt_in * (kappa + i * wl(j)) * (f1(j, f) + 0.5d0 * k2_f1(j, f)) + &
                  & (-dt_in) * gkl(j) * (sigma(sm) + 0.5d0 * k2_sigma(sm))
      k4_f1(j, f) = -dt_in * (kappa + i * wl(j)) * (f1(j, f) + k3_f1(j, f)) + &
                  & (-dt_in) * gkl(j) * (sigma(sm) + k3_sigma(sm))

      !---------------------!
      ! < f^{\dagger}_{j} > !
      !---------------------!
      k1_f1(j, ft) = -dt_in * (kappa - i * wl(j)) * f1(j, ft) + &
                  & (-dt_in) * CONJG(gkl(j)) * sigma(sp)
      k2_f1(j, ft) = -dt_in * (kappa - i * wl(j)) * (f1(j, ft) + 0.5d0 * k1_f1(j, ft)) + &
                  & (-dt_in) * CONJG(gkl(j)) * (sigma(sp) + 0.5d0 * k1_sigma(sp))
      k3_f1(j, ft) = -dt_in * (kappa - i * wl(j)) * (f1(j, ft) + 0.5d0 * k2_f1(j, ft)) + &
                  & (-dt_in) * CONJG(gkl(j)) * (sigma(sp) + 0.5d0 * k2_sigma(sp))
      k4_f1(j, ft) = -dt_in * (kappa - i * wl(j)) * (f1(j, ft) + k3_f1(j, ft)) + &
                  & (-dt_in) * CONJG(gkl(j)) * (sigma(sp) + k3_sigma(sp))

      !---------------------------------------!
      !     SECOND-ORDER: CAVITY AND ATOM     !
      !---------------------------------------!
      ! Using matrix multiplication, we add to the Lindblad matrix M_OG and the
      ! non-homogeneous vector, then compute the new moments

      !------------------!
      ! < f_{j} \sigma > !
      !------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - (kappa + i * wl(j))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = 0.0d0
      B_vec(2) = -0.5d0 * gkl(j) * (sigma(sz) + photon_ss(f))
      B_vec(3) = -gamma_in * f1(j, f) + &
               & gkl(j) * sigma(sm)
      ! Calculate k1
      k1_f1sig(j, f, :) = dt_in * (MATMUL(Mat, f1sig(j, f, :)) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = 0.0d0
      B_vec(2) = -0.5d0 * gkl(j) * ((sigma(sz) + 0.5d0 * k1_sigma(sz)) + photon_ss(f))
      B_vec(3) = -gamma_in * (f1(j, f) + 0.5d0 * k1_f1(j, f)) + &
               & gkl(j) * (sigma(sm) + 0.5d0 * k1_sigma(sm))
      ! Calculate k2
      k2_f1sig(j, f, :) = dt_in * (MATMUL(Mat, (f1sig(j, f, :) + 0.5d0 * k1_f1sig(j, f, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = 0.0d0
      B_vec(2) = -0.5d0 * gkl(j) * ((sigma(sz) + 0.5d0 * k2_sigma(sz)) + photon_ss(f))
      B_vec(3) = -gamma_in * (f1(j, f) + 0.5d0 * k2_f1(j, f)) + &
               & gkl(j) * (sigma(sm) + 0.5d0 * k2_sigma(sm))
      ! Calculate k3
      k3_f1sig(j, f, :) = dt_in * (MATMUL(Mat, (f1sig(j, f, :) + 0.5d0 * k2_f1sig(j, f, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = 0.0d0
      B_vec(2) = -0.5d0 * gkl(j) * ((sigma(sz) + k3_sigma(sz)) + photon_ss(f))
      B_vec(3) = -gamma_in * (f1(j, f) + k3_f1(j, f)) + &
               & gkl(j) * (sigma(sm) + k3_sigma(sm))
      ! Calculate k4
      k4_f1sig(j, f, :) = dt_in * (MATMUL(Mat, (f1sig(j, f, :) + k3_f1sig(j, f, :))) + B_vec)

      !----------------------------!
      ! < f^{\dagger}_{j} \sigma > !
      !----------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - (kappa - i * wl(j))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -0.5d0 * CONJG(gkl(j)) * (sigma(sz) + photon_ss(f))
      B_vec(2) = 0.0d0
      B_vec(3) = -gamma_in * f1(j, ft) + &
               & CONJG(gkl(j)) * sigma(sp)
      ! Calculate k1
      k1_f1sig(j, ft, :) = dt_in * (MATMUL(Mat, f1sig(j, ft, :)) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -0.5d0 * CONJG(gkl(j)) * ((sigma(sz) + 0.5d0 * k1_sigma(sz)) + photon_ss(f))
      B_vec(2) = 0.0d0
      B_vec(3) = -gamma_in * (f1(j, ft) + 0.5d0 * k1_f1(j, ft)) + &
               & CONJG(gkl(j)) * (sigma(sp) + 0.5d0 * k1_sigma(sp))
      ! Calculate k2
      k2_f1sig(j, ft, :) = dt_in * (MATMUL(Mat, (f1sig(j, ft, :) + 0.5d0 * k1_f1sig(j, ft, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -0.5d0 * CONJG(gkl(j)) * ((sigma(sz) + 0.5d0 * k2_sigma(sz)) + photon_ss(f))
      B_vec(2) = 0.0d0
      B_vec(3) = -gamma_in * (f1(j, ft) + 0.5d0 * k2_f1(j, ft)) + &
               & CONJG(gkl(j)) * (sigma(sp) + 0.5d0 * k2_sigma(sp))
      ! Calculate k3
      k3_f1sig(j, ft, :) = dt_in * (MATMUL(Mat, (f1sig(j, ft, :) + 0.5d0 * k2_f1sig(j, ft, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -0.5d0 * CONJG(gkl(j)) * ((sigma(sz) + k3_sigma(sz)) + photon_ss(f))
      B_vec(2) = 0.0d0
      B_vec(3) = -gamma_in * (f1(j, ft) + k3_f1(j, ft)) + &
               & CONJG(gkl(j)) * (sigma(sp) + k3_sigma(sp))
      ! Calculate k4
      k4_f1sig(j, ft, :) = dt_in * (MATMUL(Mat, (f1sig(j, ft, :) + k3_f1sig(j, ft, :))) + B_vec)

      ! Close j loop
    END DO

    ! Cycle through modes
    DO j = -N_in, N_in
      DO k = -N_in, N_in
        !------------------------------!
        !     SECOND-ORDER: CAVITY     !
        !------------------------------!
        !---------------------------!
        ! < f^{\dagger}_{j} f_{k} > !
        !---------------------------!
        k1_f2(j, k) = -dt_in * (2.0d0 * kappa - i * (wl(j) - wl(k))) * f2(j, k) + &
                    & (-dt_in) * CONJG(gkl(j)) * f1sig(k, f, sp) + &
                    & (-dt_in) * gkl(k) * f1sig(j, ft, sm)
        k2_f2(j, k) = -dt_in * (2.0d0 * kappa - i * (wl(j) - wl(k))) * (f2(j, k) + 0.5d0 * k1_f2(j, k)) + &
                    & (-dt_in) * CONJG(gkl(j)) * (f1sig(k, f, sp) + 0.5d0 * k1_f1sig(k, f, sp)) + &
                    & (-dt_in) * gkl(k) * (f1sig(j, ft, sm) + 0.5d0 * k1_f1sig(j, ft, sm))
        k3_f2(j, k) = -dt_in * (2.0d0 * kappa - i * (wl(j) - wl(k))) * (f2(j, k) + 0.5d0 * k2_f2(j, k)) + &
                    & (-dt_in) * CONJG(gkl(j)) * (f1sig(k, f, sp) + 0.5d0 * k2_f1sig(k, f, sp)) + &
                    & (-dt_in) * gkl(k) * (f1sig(j, ft, sm) + 0.5d0 * k2_f1sig(j, ft, sm))
        k4_f2(j, k) = -dt_in * (2.0d0 * kappa - i * (wl(j) - wl(k))) * (f2(j, k) + k3_f2(j, k)) + &
                    & (-dt_in) * CONJG(gkl(j)) * (f1sig(k, f, sp) + k3_f1sig(k, f, sp)) + &
                    & (-dt_in) * gkl(k) * (f1sig(j, ft, sm) + k3_f1sig(j, ft, sm))

        ! Close k loop
      END DO
      ! Close j loop
    END DO


    !============================================================================!
    !                   UPDATE ARRAYS FROM RUNGE-KUTTA ARRAYS                    !
    !============================================================================!
    ! First-order
    sigma = sigma + xis * (k1_sigma + 2.0d0 * (k2_sigma + k3_sigma) + k4_sigma)
    f1 = f1 + xis * (k1_f1 + 2.0d0 * (k2_f1 + k3_f1) + k4_f1)
    ! Second-order
    f1sig = f1sig + xis * (k1_f1sig + 2.0d0 * (k2_f1sig + k3_f1sig) + k4_f1sig)
    f2 = f2 + xis * (k1_f2 + 2.0d0 * (k2_f2 + k3_f2) + k4_f2)

    ! Close t loop
  END DO

  ! If WRITE_DATA_IN is TRUE, close the file
  IF (WRITE_DATA_IN .EQV. .TRUE.) THEN
    ! Close file
    CLOSE(4)
  END IF

END SUBROUTINE G2_CalculateRK4

END MODULE TWO_FILTER_SUBROUTINES