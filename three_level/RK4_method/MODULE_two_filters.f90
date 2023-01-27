! This module contains the subroutines used in any of the two-filter
! programs [g2_cross_RK4.f90, g2_cross_twotime_RK4.f90].

! This module contains the following subroutines:
! - SquareMatrixInverse: Calculates the inverse of an input matrix.
!
! - SquareMatrixZeroEigenvalue: Calculates the eigenvector of a matrix
!                               corresponding to the zero-valued eigenvalue.
!
! - SteadyStateMoments: Calculates the steady states of the various operator
!                       moment equations for the atom-filter coupled system.
!
! - MatrixInverseSS: Uses SquareMatrixInverse to return the steady state array
!                    without having to type out the multiplication.
!
! - G2_InitialConditions: Calculates the initial conditions for the second-
!                         order cross correlation function for two filters.
!
! - G2_CalculateRK4: Calculates the time evolution of the second-order
!                    cross correlation function for two filters using
!                    Runge-Kutta 4th Order.
!
! - MatrixInverseSS : Calculates the steady state for a system of coupled
!                     differential equations using the inverse matrix method.

! This file must be added to the compilation command when compiling any of the 
! single-filter programs. Eg, with Intel oneAPI or GFORTRAN:
!     (IFORT): ifort -qmkl ./[FILENAME].f90 ./MODULE_single_filter.f90
!  (GFORTRAN): gfortran ./[FILENAME].f90 ./MODULE_single_filter.f90
!                -I/path/to/LAPACK -L/path/to/LAPACK -llapack -lblas

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

  ! Query optimal work space
  ! LWORK = -1
  ! CALL zGETRI(N_in, MatrixInv_out, N_in, IPIV, WORK, LWORK, INFO)
  ! ! Set optimal work space and run again
  ! LWORK = MIN(LWMAX, INT(WORK(1)))

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
SUBROUTINE SteadyStateMoments(Gamma_in, Omega_in, alpha_in, delta_in, xi_in, &
                            & epsilon_in, N_in, phase_in, &
                            & w0a_in, kappaa_in, dwa_in, &
                            & w0b_in, kappab_in, dwb_in, &
                            & photon_out, sigma_out, &
                            & f1_out, f1sig_out, &
                            & f2_out, f2sig_out, &
                            & f3_out, f3sig_out, &
                            & f4_out)

  !==============================================================================!
  !                    DEFINING AND DECLARING VARIABLES/ARRAYS                   !
  !==============================================================================!

  IMPLICIT NONE

  !---------------!
  !     INPUT     !
  !---------------!
  ! Atomic decay rate
  REAL(KIND=8), INTENT(IN)                   :: Gamma_in
  ! Driving amplitude
  REAL(KIND=8), INTENT(IN)                   :: Omega_in
  ! Atomic anharmonicity
  REAL(KIND=8), INTENT(IN)                   :: alpha_in
  ! Drive detuning from two-photon resonance
  REAL(KIND=8), INTENT(IN)                   :: delta_in
  ! Dipole moment ratio
  REAL(KIND=8), INTENT(IN)                   :: xi_in

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
  INTEGER, PARAMETER                        :: N_mat = 8
  ! M matrix (filled as transpose)
  COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)  :: Mat, Mat_OG, Mat_inv
  ! Non-homogeneous vector
  COMPLEX(KIND=8), DIMENSION(N_mat)         :: B_vec, B_OG

  ! Integer indices for sigma operators
  INTEGER, PARAMETER                        :: gg = 1, ge = 2, eg = 3
  INTEGER, PARAMETER                        :: ee = 4, ef = 5, fe = 6
  INTEGER, PARAMETER                        :: gf = 7, fg = 8
  ! Integer indices for: f, f^{\dagger}, f^{2}, f^{\dagger}^{2} ... etc
  INTEGER, PARAMETER                        :: f = 1, ft = 3, g = 2, gt = 4
  INTEGER, PARAMETER                        :: fg2 = 1, ftgt = 2, ftf = 3, gtg = 4, ftg = 5, gtf = 6
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
  INTEGER                                    :: j, k, l, m, x
  ! Imaginary i
  COMPLEX(KIND=8), PARAMETER                 :: i = CMPLX(0.0d0, 1.0d0, 8)
  ! pi
  REAL(KIND=8), PARAMETER                    :: pi = 3.1415926535897932384d0
  ! List of Delta values
  REAL(KIND=8), DIMENSION(-N_in:N_in, 2)     :: wl
  ! List of mode dependent cascade coupling values
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2)  :: gkl
  ! Blackman window coefficient
  REAL(KIND=8)                               :: blackman
  ! Parameters
  REAL(KIND=8), DIMENSION(2)                 :: kappa, w0, dw
  ! Temporary values
  COMPLEX(KIND=8)                            :: moment_out, moment_out2

  !==============================================================================!
  !                DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS                !
  !==============================================================================!
  kappa(f) = kappaa_in; kappa(g) = kappab_in
  w0(f) = w0a_in; w0(g) = w0b_in
  dw(f) = dwa_in; dw(g) = dwb_in

  !------------------------!
  !     BLOCH MATRIX M     !
  !------------------------!
  Mat_OG = 0.0d0
  ! Row 1: d/dt |g><g|
  Mat_OG(1, 2) = -i * 0.5d0 * Omega_in
  Mat_OG(1, 3) = i * 0.5d0 * Omega_in
  Mat_OG(1, 4) = Gamma_in
  ! Row 2: d/dt |g><e|
  Mat_OG(2, 1) = -i * 0.5d0 * Omega_in
  Mat_OG(2, 2) = -(0.5d0 * Gamma_in - i * ((0.5d0 * alpha_in) + delta_in))
  Mat_OG(2, 4) = i * 0.5d0 * Omega_in
  Mat_OG(2, 5) = Gamma_in * xi_in
  Mat_OG(2, 7) = -i * xi_in * 0.5d0 * Omega_in
  ! Row 3: d/dt |e><g|
  Mat_OG(3, 1) = i * 0.5d0 * Omega_in
  Mat_OG(3, 3) = -(0.5d0 * Gamma_in + i * ((0.5d0 * alpha_in) + delta_in))
  Mat_OG(3, 4) = -i * 0.5d0 * Omega_in
  Mat_OG(3, 6) = Gamma_in * xi_in
  Mat_OG(3, 8) = i * xi_in * 0.5d0 * Omega_in
  ! Row 4: d/dt |e><e|
  Mat_OG(4, 1) = -Gamma_in * (xi_in ** 2)
  Mat_OG(4, 2) = i * 0.5d0 * Omega_in
  Mat_OG(4, 3) = -i * 0.5d0 * Omega_in
  Mat_OG(4, 4) = -Gamma_in * (1.0d0 + (xi_in ** 2))
  Mat_OG(4, 5) = -i * xi_in * 0.5d0 * Omega_in
  Mat_OG(4, 6) = i * xi_in * 0.5d0 * Omega_in
  ! Row 5: d/dt |e><f|
  Mat_OG(5, 1) = -i * xi_in * 0.5d0 * Omega_in
  Mat_OG(5, 4) = -i * xi_in * Omega_in
  Mat_OG(5, 5) = -(0.5d0 * Gamma_in * (1.0d0 + (xi_in ** 2)) + i * ((0.5d0 * alpha_in) - delta_in))
  Mat_OG(5, 7) = i * 0.5d0 * Omega_in
  ! Row 6: d/dt |f><e|
  Mat_OG(6, 1) = i * xi_in * 0.5d0 * Omega_in
  Mat_OG(6, 4) = i * xi_in * Omega_in
  Mat_OG(6, 6) = -(0.5d0 * Gamma_in * (1.0d0 + (xi_in ** 2)) - i * ((0.5d0 * alpha_in) - delta_in))
  Mat_OG(6, 8) = -i * 0.5d0 * Omega_in
  ! Row 7: d/dt |g><f|
  Mat_OG(7, 2) = -i * xi_in * 0.5d0 * Omega_in
  Mat_OG(7, 5) = i * 0.5d0 * Omega_in
  Mat_OG(7, 7) = -(0.5d0 * Gamma_in * (xi_in ** 2) - 2.0d0 * i * delta_in)
  ! Row 8: d/dt |g><f|
  Mat_OG(8, 3) = i * xi_in * 0.5d0 * Omega_in
  Mat_OG(8, 6) = -i * 0.5d0 * Omega_in
  Mat_OG(8, 8) = -(0.5d0 * Gamma_in * (xi_in ** 2) + 2.0d0 * i * delta_in)

  !--------------------------------!
  !     NON-HOMOGENEOUS VECTOR     !
  !--------------------------------!
  B_OG = 0.0d0
  B_OG(4) = Gamma_in * (xi_in ** 2)
  B_OG(5) = i * xi_in * 0.5d0 * Omega_in
  B_OG(6) = -i * xi_in * 0.5d0 * Omega_in

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
      gkl(j, f) = DSQRT(epsilon_in * Gamma_in * kappa(f))

      ! Cavity B
      wl(j, g) = w0(g)
      gkl(j, g) = DSQRT(epsilon_in * Gamma_in * kappa(g))
    ELSE
      ! Blackman window coefficient
      blackman = 1.0d0
      ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
      !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))

      ! Cavity A
      wl(j, f) = w0(f) + DBLE(j) * dw(f)
      ! Mode dependent phase difference
      gkl(j, f) = DSQRT((epsilon_in / DBLE(2*N_in + 1)) * Gamma_in * kappa(f)) * blackman * &
                & EXP(i * DBLE(phase_in) * DBLE(j) * pi / DBLE(N_in))

      ! Cavity B
      wl(j, g) = w0(g) + DBLE(j) * dw(g)
      ! Mode dependent phase difference
      gkl(j, g) = DSQRT((epsilon_in / DBLE(2*N_in + 1)) * Gamma_in * kappa(g)) * blackman * &
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

  !==============================================================================!
  !                        CALCULATE STEADY-STATE MOMENTS                        !
  !==============================================================================!
  !---------------------------!
  !     FIRST-ORDER: ATOM     !
  !---------------------------!
  IF (xi_in .NE. 0.0d0) THEN
    ! Calculate steady states
    CALL MatrixInverseSS(N_mat, Mat_OG, B_OG, sigma_out)

  ELSE IF (xi_in .EQ. 0.0d0) THEN
    !------------------------------------------------!
    !     CALCULATE EIGENVALUES AND EIGENVECTORS     !
    !------------------------------------------------!
    ! Set Lindblad matrix
    Mat = Mat_OG

    ! Calculate steady state from eigenvectors
    CALL SquareMatrixZeroEigenvalue(N_mat, Mat, sigma_out)
    ! Normalise sigma_out so |g><g| + |e><e| = 1
    sigma_out = sigma_out / (REAL(sigma_out(1)) + REAL(sigma_out(4)))
  END IF

  ! Cycle through modes
  DO j = -N_in, N_in
    !-----------------------------!
    !     FIRST-ORDER: CAVITY     !
    !-----------------------------!
    !-----------!
    ! < f_{j} > !
    !-----------!
    f1_out(j, f) = -gkl(j, f) * sigma_out(ge) + &
                 & -gkl(j, f) * xi_in * sigma_out(ef)
    f1_out(j, f) = f1_out(j, f) / &
                 & (kappa(f) + i * wl(j, f))

    !---------------------!
    ! < f^{\dagger}_{j} > !
    !---------------------!
    f1_out(j, ft) = -CONJG(gkl(j, f)) * sigma_out(eg) + &
                  & -CONJG(gkl(j, f)) * xi_in * sigma_out(fe)
    f1_out(j, ft) = f1_out(j, ft) / &
                  & (kappa(f) - i * wl(j, f))

    !-----------!
    ! < g_{j} > !
    !-----------!
    f1_out(j, g) = -gkl(j, g) * sigma_out(ge) + &
                 & -gkl(j, g) * xi_in * sigma_out(ef)
    f1_out(j, g) = f1_out(j, g) / &
                 & (kappa(g) + i * wl(j, g))

    !---------------------!
    ! < g^{\dagger}_{j} > !
    !---------------------!
    f1_out(j, gt) = -CONJG(gkl(j, g)) * sigma_out(eg) + &
                  & -CONJG(gkl(j, g)) * xi_in * sigma_out(fe)
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
    B_vec(1) = -gkl(j, f) * sigma_out(ge)
    B_vec(2) = -gkl(j, f) * xi_in * sigma_out(gf)
    B_vec(3) = -gkl(j, f) * sigma_out(ee)
    B_vec(4) = Gamma_in * (xi_in ** 2) * f1_out(j, f) + &
             & -gkl(j, f) * xi_in * sigma_out(ef)
    B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f1_out(j, f)
    B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f1_out(j, f) + &
             & gkl(j, f) * xi_in * sigma_out(gg) + &
             & gkl(j, f) * xi_in * sigma_out(ee) + &
             & -gkl(j, f) * xi_in
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j, f) * sigma_out(fe)

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
    B_vec(1) = -CONJG(gkl(j, f)) * sigma_out(eg)
    B_vec(2) = -CONJG(gkl(j, f)) * sigma_out(ee)
    B_vec(3) = -CONJG(gkl(j, f)) * xi_in * sigma_out(fg)
    B_vec(4) = Gamma_in * (xi_in ** 2) * f1_out(j, ft) + &
             & -CONJG(gkl(j, f)) * xi_in * sigma_out(fe)
    B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f1_out(j, ft) + &
             & CONJG(gkl(j, f)) * xi_in * sigma_out(gg) + &
             & CONJG(gkl(j, f)) * xi_in * sigma_out(ee) + &
             & -CONJG(gkl(j, f)) * xi_in
    B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f1_out(j, ft)
    B_vec(7) = -CONJG(gkl(j, f)) * sigma_out(ef)
    B_vec(8) = 0.0d0

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
    B_vec(1) = -gkl(j, g) * sigma_out(ge)
    B_vec(2) = -gkl(j, g) * xi_in * sigma_out(gf)
    B_vec(3) = -gkl(j, g) * sigma_out(ee)
    B_vec(4) = Gamma_in * (xi_in ** 2) * f1_out(j, g) + &
             & -gkl(j, g) * xi_in * sigma_out(ef)
    B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f1_out(j, g)
    B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f1_out(j, g) + &
             & gkl(j, g) * xi_in * sigma_out(gg) + &
             & gkl(j, g) * xi_in * sigma_out(ee) + &
             & -gkl(j, g) * xi_in
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j, g) * sigma_out(fe)

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
    B_vec(1) = -CONJG(gkl(j, g)) * sigma_out(eg)
    B_vec(2) = -CONJG(gkl(j, g)) * sigma_out(ee)
    B_vec(3) = -CONJG(gkl(j, g)) * xi_in * sigma_out(fg)
    B_vec(4) = Gamma_in * (xi_in ** 2) * f1_out(j, gt) + &
             & -CONJG(gkl(j, g)) * xi_in * sigma_out(fe)
    B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f1_out(j, gt) + &
             & CONJG(gkl(j, g)) * xi_in * sigma_out(gg) + &
             & CONJG(gkl(j, g)) * xi_in * sigma_out(ee) + &
             & -CONJG(gkl(j, g)) * xi_in
    B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f1_out(j, gt)
    B_vec(7) = -CONJG(gkl(j, g)) * sigma_out(ef)
    B_vec(8) = 0.0d0

    ! Calculate steady states
    CALL MatrixInverseSS(N_mat, Mat, B_vec, f1sig_out(j, gt, :))

    ! Close j loop
  END DO

  moment_out = 0.0d0
  moment_out2 = 0.0d0
  ! Cycle through modes
  DO k = -N_in, N_in
    DO j = -N_in, N_in
      !------------------------------!
      !     SECOND-ORDER: CAVITY     !
      !------------------------------!
      !-----------------!
      ! < f_{j} g_{k} > !
      !-----------------!
      f2_out(j, k, fg2) = -gkl(j, f) * f1sig_out(k, g, ge) + &
                        & -gkl(j, f) * xi_in * f1sig_out(k, g, ef) + &
                        & -gkl(k, g) * f1sig_out(j, f, ge) + &
                        & -gkl(k, g) * xi_in * f1sig_out(j, f, ef)
      f2_out(j, k, fg2) = f2_out(j, k, fg2) / &
                        & ((kappa(f) + kappa(g)) + i * (wl(j, f) + wl(k, g)))

      !-------------------------------------!
      ! < f^{\dagger}_{j} g^{\dagger}_{k} > !
      !-------------------------------------!
      f2_out(j, k, ftgt) = -CONJG(gkl(j, f)) * f1sig_out(k, gt, eg) + &
                         & -CONJG(gkl(j, f)) * xi_in * f1sig_out(k, gt, fe) + &
                         & -CONJG(gkl(k, g)) * f1sig_out(j, ft, eg) + &
                         & -CONJG(gkl(k, g)) * xi_in * f1sig_out(j, ft, fe)
      f2_out(j, k, ftgt) = f2_out(j, k, ftgt) / &
                         & ((kappa(f) + kappa(g)) - i * (wl(j, f) + wl(k, g)))

      !---------------------------!
      ! < f^{\dagger}_{j} f_{k} > !
      !---------------------------!
      f2_out(j, k, ftf) = -CONJG(gkl(j, f)) * f1sig_out(k, f, eg) + &
                        & -CONJG(gkl(j, f)) * xi_in * f1sig_out(k, f, fe) + &
                        & -gkl(k, f) * f1sig_out(j, ft, ge) + &
                        & -gkl(k, f) * xi_in * f1sig_out(j, ft, ef)
      f2_out(j, k, ftf) = f2_out(j, k, ftf) / &
                        & (2.0d0 * kappa(f) - i * (wl(j, f) - wl(k, f)))

      ! Update photon number
      moment_out = moment_out + f2_out(j, k, ftf)

      !---------------------------!
      ! < g^{\dagger}_{j} g_{k} > !
      !---------------------------!
      f2_out(j, k, gtg) = -CONJG(gkl(j, g)) * f1sig_out(k, g, eg) + &
                        & -CONJG(gkl(j, g)) * xi_in * f1sig_out(k, g, fe) + &
                        & -gkl(k, g) * f1sig_out(j, gt, ge) + &
                        & -gkl(k, g) * xi_in * f1sig_out(j, gt, ef)
      f2_out(j, k, gtg) = f2_out(j, k, gtg) / &
                        & (2.0d0 * kappa(g) - i * (wl(j, g) - wl(k, g)))

      ! Update photon number
      moment_out2 = moment_out2 + f2_out(j, k, gtg)

      !---------------------------!
      ! < f^{\dagger}_{j} g_{k} > !
      !---------------------------!
      f2_out(j, k, ftg) = -CONJG(gkl(j, f)) * f1sig_out(k, g, eg) + &
                        & -CONJG(gkl(j, f)) * xi_in * f1sig_out(k, g, fe) + &
                        & -gkl(k, g) * f1sig_out(j, ft, ge) + &
                        & -gkl(k, g) * xi_in * f1sig_out(j, ft, ef)
      f2_out(j, k, ftg) = f2_out(j, k, ftg) / &
                        & ((kappa(f) + kappa(g)) - i * (wl(j, f) - wl(k, g)))

      !---------------------------!
      ! < g^{\dagger}_{j} f_{k} > !
      !---------------------------!
      f2_out(j, k, gtf) = -CONJG(gkl(j, g)) * f1sig_out(k, f, eg) + &
                        & -CONJG(gkl(j, g)) * xi_in * f1sig_out(k, f, fe) + &
                        & -gkl(k, f) * f1sig_out(j, gt, ge) + &
                        & -gkl(k, f) * xi_in * f1sig_out(j, gt, ef)
      f2_out(j, k, gtf) = f2_out(j, k, gtf) / &
                        & ((kappa(f) + kappa(g)) - i * (wl(j, g) - wl(k, f)))

      !--------------------------------------!
      !     THIRD-ORDER: CAVITY AND ATOM     !
      !--------------------------------------!
      !------------------------!
      ! < f_{j} g_{k} \sigma > !
      !------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - ((kappa(f) + kappa(g)) + i * (wl(j, f) + wl(k, g)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -gkl(j, f) * f1sig_out(k, g, ge) + &
               & -gkl(k, g) * f1sig_out(j, f, ge)
      B_vec(2) = -gkl(j, f) * xi_in * f1sig_out(k, g, gf) + &
               & -gkl(k, g) * xi_in * f1sig_out(j, f, gf)
      B_vec(3) = -gkl(j, f) * f1sig_out(k, g, ee) + &
               & -gkl(k, g) * f1sig_out(j, f, ee)
      B_vec(4) = Gamma_in * (xi_in ** 2) * f2_out(j, k, fg2) + &
               & -gkl(j, f) * xi_in * f1sig_out(k, g, ef) + &
               & -gkl(k, g) * xi_in * f1sig_out(j, f, ef)
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, fg2)
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, fg2) + &
               & gkl(j, f) * xi_in * f1sig_out(k, g, gg) + &
               & gkl(j, f) * xi_in * f1sig_out(k, g, ee) + &
               & -gkl(j, f) * xi_in * f1_out(k, g) + &
               & gkl(k, g) * xi_in * f1sig_out(j, f, gg) + &
               & gkl(k, g) * xi_in * f1sig_out(j, f, ee) + &
               & -gkl(k, g) * xi_in * f1_out(j, f)
      B_vec(7) = 0.0d0
      B_vec(8) = -gkl(j, f) * f1sig_out(k, g, fe) + &
               & -gkl(k, g) * f1sig_out(j, f, fe)

      ! Calculate steady states
      CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, fg2, :))

      !--------------------------------------------!
      ! < f^{\dagger}_{j} g^{\dagger}_{k} \sigma > !
      !--------------------------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - ((kappa(f) + kappa(g)) - i * (wl(j, f) + wl(k, g)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j, f)) * f1sig_out(k, gt, eg) + &
               & -CONJG(gkl(k, g)) * f1sig_out(j, ft, eg)
      B_vec(2) = -CONJG(gkl(j, f)) * f1sig_out(k, gt, ee) + &
               & -CONJG(gkl(k, g)) * f1sig_out(j, ft, ee)
      B_vec(3) = -CONJG(gkl(j, f)) * xi_in * f1sig_out(k, gt, fg) + &
               & -CONJG(gkl(k, g)) * xi_in * f1sig_out(j, ft, fg)
      B_vec(4) = Gamma_in * (xi_in ** 2) * f2_out(j, k, ftgt) + &
               & -CONJG(gkl(j, f)) * xi_in * f1sig_out(k, gt, fe) + &
               & -CONJG(gkl(k, g)) * xi_in * f1sig_out(j, ft, fe)
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, ftgt) + &
               & CONJG(gkl(j, f)) * xi_in * f1sig_out(k, gt, gg) + &
               & CONJG(gkl(j, f)) * xi_in * f1sig_out(k, gt, ee) + &
               & -CONJG(gkl(j, f)) * xi_in * f1_out(k, gt) + &
               & CONJG(gkl(k, g)) * xi_in * f1sig_out(j, ft, gg) + &
               & CONJG(gkl(k, g)) * xi_in * f1sig_out(j, ft, ee) + &
               & -CONJG(gkl(k, g)) * xi_in * f1_out(j, ft)
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, ftgt)
      B_vec(7) = -CONJG(gkl(j, f)) * f1sig_out(k, gt, ef) + &
               & -CONJG(gkl(k, g)) * f1sig_out(j, ft, ef)
      B_vec(8) = 0.0d0

      ! Calculate steady states
      CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, ftgt, :))

      !----------------------------------!
      ! < f^{\dagger}_{j} f_{k} \sigma > !
      !----------------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - (2.0d0 * kappa(f) - i * (wl(j, f) - wl(k, f)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j, f)) * f1sig_out(k, f, eg) + &
               & -gkl(k, f) * f1sig_out(j, ft, ge)
      B_vec(2) = -CONJG(gkl(j, f)) * f1sig_out(k, f, ee) + &
               & -gkl(k, f) * xi_in * f1sig_out(j, ft, gf)
      B_vec(3) = -CONJG(gkl(j, f)) * xi_in * f1sig_out(k, f, fg) + &
               & -gkl(k, f) * f1sig_out(j, ft, ee)
      B_vec(4) = Gamma_in * (xi_in ** 2) * f2_out(j, k, ftf) + &
               & -CONJG(gkl(j, f)) * xi_in * f1sig_out(k, f, fe) + &
               & -gkl(k, f) * xi_in * f1sig_out(j, ft, ef)
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, ftf) + &
               & CONJG(gkl(j, f)) * xi_in * f1sig_out(k, f, gg) + &
               & CONJG(gkl(j, f)) * xi_in * f1sig_out(k, f, ee) + &
               & -CONJG(gkl(j, f)) * xi_in * f1_out(k, f)
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, ftf) + &
               & gkl(k, f) * xi_in * f1sig_out(j, ft, gg) + &
               & gkl(k, f) * xi_in * f1sig_out(j, ft, ee) + &
               & -gkl(k, f) * xi_in * f1_out(j, ft)
      B_vec(7) = -CONJG(gkl(j, f)) * f1sig_out(k, f, ef)
      B_vec(8) = - gkl(k, f) * f1sig_out(j, ft, fe)

      ! Calculate steady states
      CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, ftf, :))

      !----------------------------------!
      ! < g^{\dagger}_{j} g_{k} \sigma > !
      !----------------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - (2.0d0 * kappa(g) - i * (wl(j, g) - wl(k, g)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j, g)) * f1sig_out(k, g, eg) + &
               & -gkl(k, g) * f1sig_out(j, gt, ge)
      B_vec(2) = -CONJG(gkl(j, g)) * f1sig_out(k, g, ee) + &
               & -gkl(k, g) * xi_in * f1sig_out(j, gt, gf)
      B_vec(3) = -CONJG(gkl(j, g)) * xi_in * f1sig_out(k, g, fg) + &
               & -gkl(k, g) * f1sig_out(j, gt, ee)
      B_vec(4) = Gamma_in * (xi_in ** 2) * f2_out(j, k, gtg) + &
               & -CONJG(gkl(j, g)) * xi_in * f1sig_out(k, g, fe) + &
               & -gkl(k, g) * xi_in * f1sig_out(j, gt, ef)
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, gtg) + &
               & CONJG(gkl(j, g)) * xi_in * f1sig_out(k, g, gg) + &
               & CONJG(gkl(j, g)) * xi_in * f1sig_out(k, g, ee) + &
               & -CONJG(gkl(j, g)) * xi_in * f1_out(k, g)
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, gtg) + &
               & gkl(k, g) * xi_in * f1sig_out(j, gt, gg) + &
               & gkl(k, g) * xi_in * f1sig_out(j, gt, ee) + &
               & -gkl(k, g) * xi_in * f1_out(j, gt)
      B_vec(7) = -CONJG(gkl(j, g)) * f1sig_out(k, g, ef)
      B_vec(8) = - gkl(k, g) * f1sig_out(j, gt, fe)

      ! Calculate steady states
      CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, gtg, :))

      !----------------------------------!
      ! < f^{\dagger}_{j} g_{k} \sigma > !
      !----------------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - ((kappa(f) + kappa(g)) - i * (wl(j, f) - wl(k, g)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j, f)) * f1sig_out(k, g, eg) + &
               & -gkl(k, g) * f1sig_out(j, ft, ge)
      B_vec(2) = -CONJG(gkl(j, f)) * f1sig_out(k, g, ee) + &
               & -gkl(k, g) * xi_in * f1sig_out(j, ft, gf)
      B_vec(3) = -CONJG(gkl(j, f)) * xi_in * f1sig_out(k, g, fg) + &
               & -gkl(k, g) * f1sig_out(j, ft, ee)
      B_vec(4) = Gamma_in * (xi_in ** 2) * f2_out(j, k, ftg) + &
               & -CONJG(gkl(j, f)) * xi_in * f1sig_out(k, g, fe) + &
               & -gkl(k, g) * xi_in * f1sig_out(j, ft, ef)
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, ftg) + &
               & CONJG(gkl(j, f)) * xi_in * f1sig_out(k, g, gg) + &
               & CONJG(gkl(j, f)) * xi_in * f1sig_out(k, g, ee) + &
               & -CONJG(gkl(j, f)) * xi_in * f1_out(k, g)
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, ftg) + &
               & gkl(k, g) * xi_in * f1sig_out(j, ft, gg) + &
               & gkl(k, g) * xi_in * f1sig_out(j, ft, ee) + &
               & -gkl(k, g) * xi_in * f1_out(j, ft)
      B_vec(7) = -CONJG(gkl(j, f)) * f1sig_out(k, g, ef)
      B_vec(8) = - gkl(k, g) * f1sig_out(j, ft, fe)

      ! Calculate steady states
      CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, ftg, :))

      !----------------------------------!
      ! < g^{\dagger}_{j} f_{k} \sigma > !
      !----------------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - ((kappa(f) + kappa(g)) - i * (wl(j, g) - wl(k, f)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j, g)) * f1sig_out(k, f, eg) + &
              & -gkl(k, f) * f1sig_out(j, gt, ge)
      B_vec(2) = -CONJG(gkl(j, g)) * f1sig_out(k, f, ee) + &
              & -gkl(k, f) * xi_in * f1sig_out(j, gt, gf)
      B_vec(3) = -CONJG(gkl(j, g)) * xi_in * f1sig_out(k, f, fg) + &
              & -gkl(k, f) * f1sig_out(j, gt, ee)
      B_vec(4) = Gamma_in * (xi_in ** 2) * f2_out(j, k, gtf) + &
              & -CONJG(gkl(j, g)) * xi_in * f1sig_out(k, f, fe) + &
              & -gkl(k, f) * xi_in * f1sig_out(j, gt, ef)
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, gtf) + &
              & CONJG(gkl(j, g)) * xi_in * f1sig_out(k, f, gg) + &
              & CONJG(gkl(j, g)) * xi_in * f1sig_out(k, f, ee) + &
              & -CONJG(gkl(j, g)) * xi_in * f1_out(k, f)
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, gtf) + &
              & gkl(k, f) * xi_in * f1sig_out(j, gt, gg) + &
              & gkl(k, f) * xi_in * f1sig_out(j, gt, ee) + &
              & -gkl(k, f) * xi_in * f1_out(j, gt)
      B_vec(7) = -CONJG(gkl(j, g)) * f1sig_out(k, f, ef)
      B_vec(8) = - gkl(k, f) * f1sig_out(j, gt, fe)

      ! Calculate steady states
      CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, gtf, :))

      ! Close j loop
    END DO
    ! Close k loop
  END DO

  ! Update steady state photon number
  photon_out = 0.0d0
  photon_out(f) = REAL(moment_out)
  photon_out(g) = REAL(moment_out2)

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
        f3_out(j, k, l, ftfg) = -CONJG(gkl(j, f)) * f2sig_out(k, l, fg2, eg) + &
                              & -CONJG(gkl(j, f)) * xi_in * f2sig_out(k, l, fg2, fe) + &
                              & -gkl(k, f) * f2sig_out(j, l, ftg, ge) + &
                              & -gkl(k, f) * xi_in * f2sig_out(j, l, ftg, ef) + &
                              & -gkl(l, g) * f2sig_out(j, k, ftf, ge) + &
                              & -gkl(l, g) * xi_in * f2sig_out(j, k, ftf, ef)
        f3_out(j, k, l, ftfg) = f3_out(j, k, l, ftfg) / &
                              & ((2.d0 * kappa(f) + kappa(g)) - i * (wl(j, f) - wl(k, f) - wl(l, g)))

        !---------------------------------!
        ! < g^{\dagger}_{j} g_{k} f_{l} > !
        !---------------------------------!
        f3_out(j, k, l, gtgf) = -CONJG(gkl(j, g)) * f2sig_out(l, k, fg2, eg) + &
                              & -CONJG(gkl(j, g)) * xi_in * f2sig_out(l, k, fg2, fe) + &
                              & -gkl(k, g) * f2sig_out(j, l, gtf, ge) + &
                              & -gkl(k, g) * xi_in * f2sig_out(j, l, gtf, ef) + &
                              & -gkl(l, f) * f2sig_out(j, k, gtg, ge) + &
                              & -gkl(l, f) * xi_in * f2sig_out(j, k, gtg, ef)
        f3_out(j, k, l, gtgf) = f3_out(j, k, l, gtgf) / &
                              & ((kappa(f) + 2.0d0 * kappa(g)) - i * (wl(j, g) - wl(k, g) - wl(l, f)))

        !-------------------------------------------!
        ! < g^{\dagger}_{j} f^{\dagger}_{k} f_{l} > !
        !-------------------------------------------!
        f3_out(j, k, l, gtftf) = -CONJG(gkl(j, g)) * f2sig_out(k, l, ftf, eg) + &
                               & -CONJG(gkl(j, g)) * xi_in * f2sig_out(k, l, ftf, fe) + &
                               & -CONJG(gkl(k, f)) * f2sig_out(j, l, gtf, eg) + &
                               & -CONJG(gkl(k, f)) * xi_in * f2sig_out(j, l, gtf, fe) + &
                               & -gkl(l, f) * f2sig_out(k, j, ftgt, ge) + &
                               & -gkl(l, f) * xi_in * f2sig_out(k, j, ftgt, ef)
        f3_out(j, k, l, gtftf) = f3_out(j, k, l, gtftf) / &
                               & ((2.d0 * kappa(f) + kappa(g)) - i * (wl(j, g) + wl(k, f) - wl(l, f)))

        !-------------------------------------------!
        ! < f^{\dagger}_{j} g^{\dagger}_{k} g_{l} > !
        !-------------------------------------------!
        f3_out(j, k, l, ftgtg) = -CONJG(gkl(j, f)) * f2sig_out(k, l, gtg, eg) + &
                               & -CONJG(gkl(j, f)) * xi_in * f2sig_out(k, l, gtg, fe) + &
                               & -CONJG(gkl(k, g)) * f2sig_out(j, l, ftg, eg) + &
                               & -CONJG(gkl(k, g)) * xi_in * f2sig_out(j, l, ftg, fe) + &
                               & -gkl(l, g) * f2sig_out(j, k, ftgt, ge) + &
                               & -gkl(l, g) * xi_in * f2sig_out(j, k, ftgt, ef)
        f3_out(j, k, l, ftgtg) = f3_out(j, k, l, ftgtg) / &
                               & ((kappa(f) + 2.0d0 * kappa(g)) - i * (wl(j, f) + wl(k, g) - wl(l, g)))

        !----------------------------------------!
        ! < f^{\dagger}_{j} f_{k} g_{l} \sigma > !
        !----------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((2.d0 * kappa(f) + kappa(g)) - i * (wl(j, f) - wl(k, f) - wl(l, g)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j, f)) * f2sig_out(k, l, fg2, eg) + &
                 & -gkl(k, f) * f2sig_out(j, l, ftg, ge) + &
                 & -gkl(l, g) * f2sig_out(j, k, ftf, ge)
        B_vec(2) = -CONJG(gkl(j, f)) * f2sig_out(k, l, fg2, ee) + &
                 & -gkl(k, f) * xi_in * f2sig_out(j, l, ftg, gf) + &
                 & -gkl(l, g) * xi_in * f2sig_out(j, k, ftf, gf)
        B_vec(3) = -CONJG(gkl(j, f)) * xi_in * f2sig_out(k, l, fg2, fg) + &
                 & -gkl(k, f) * f2sig_out(j, l, ftg, ee) + &
                 & -gkl(l, g) * f2sig_out(j, k, ftf, ee)
        B_vec(4) = Gamma_in * (xi_in ** 2) * f3_out(j, k, l, ftfg) + &
                 & -CONJG(gkl(j, f)) * xi_in * f2sig_out(k, l, fg2, fe) + &
                 & -gkl(k, f) * xi_in * f2sig_out(j, l, ftg, ef) + &
                 & -gkl(l, g) * xi_in * f2sig_out(j, k, ftf, ef)
        B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f3_out(j, k, l, ftfg) + &
                 & CONJG(gkl(j, f)) * xi_in * f2sig_out(k, l, fg2, gg) + &
                 & CONJG(gkl(j, f)) * xi_in * f2sig_out(k, l, fg2, ee) + &
                 & -CONJG(gkl(j, f)) * xi_in * f2_out(k, l, fg2)
        B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f3_out(j, k, l, ftfg) + &
                 & gkl(k, f) * xi_in * f2sig_out(j, l, ftg, gg) + &
                 & gkl(k, f) * xi_in * f2sig_out(j, l, ftg, ee) + &
                 & -gkl(k, f) * xi_in * f2_out(j, l, ftg) + &
                 & gkl(l, g) * xi_in * f2sig_out(j, k, ftf, gg) + &
                 & gkl(l, g) * xi_in * f2sig_out(j, k, ftf, ee) + &
                 & -gkl(l, g) * xi_in * f2_out(j, k, ftf)
        B_vec(7) = -CONJG(gkl(j, f)) * f2sig_out(k, l, fg2, ef)
        B_vec(8) = -gkl(k, f) * f2sig_out(j, l, ftg, fe) + &
                 & -gkl(l, g) * f2sig_out(j, k, ftf, fe)

        ! Calculate steady states
        CALL MatrixInverseSS(N_mat, Mat, B_vec, f3sig_out(j, k, l, ftfg, :))

        !----------------------------------------!
        ! < g^{\dagger}_{j} g_{k} f_{l} \sigma > !
        !----------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((kappa(f) + 2.0d0 * kappa(g)) - i * (wl(j, g) - wl(k, g) - wl(l, f)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j, g)) * f2sig_out(l, k, fg2, eg) + &
                 & -gkl(k, g) * f2sig_out(j, l, gtf, ge) + &
                 & -gkl(l, f) * f2sig_out(j, k, gtg, ge)
        B_vec(2) = -CONJG(gkl(j, g)) * f2sig_out(l, k, fg2, ee) + &
                 & -gkl(k, g) * xi_in * f2sig_out(j, l, gtf, gf) + &
                 & -gkl(l, f) * xi_in * f2sig_out(j, k, gtg, gf)
        B_vec(3) = -CONJG(gkl(j, g)) * xi_in * f2sig_out(l, k, fg2, fg) + &
                 & -gkl(k, g) * f2sig_out(j, l, gtf, ee) + &
                 & -gkl(l, f) * f2sig_out(j, k, gtg, ee)
        B_vec(4) = Gamma_in * (xi_in ** 2) * f3_out(j, k, l, gtgf) + &
                 & -CONJG(gkl(j, g)) * xi_in * f2sig_out(l, k, fg2, fe) + &
                 & -gkl(k, g) * xi_in * f2sig_out(j, l, gtf, ef) + &
                 & -gkl(l, f) * xi_in * f2sig_out(j, k, gtg, ef)
        B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f3_out(j, k, l, gtgf) + &
                 & CONJG(gkl(j, g)) * xi_in * f2sig_out(l, k, fg2, gg) + &
                 & CONJG(gkl(j, g)) * xi_in * f2sig_out(l, k, fg2, ee) + &
                 & -CONJG(gkl(j, g)) * xi_in * f2_out(l, k, fg2)
        B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f3_out(j, k, l, gtgf) + &
                 & gkl(k, g) * xi_in * f2sig_out(j, l, gtf, gg) + &
                 & gkl(k, g) * xi_in * f2sig_out(j, l, gtf, ee) + &
                 & -gkl(k, g) * xi_in * f2_out(j, l, gtf) + &
                 & gkl(l, f) * xi_in * f2sig_out(j, k, gtg, gg) + &
                 & gkl(l, f) * xi_in * f2sig_out(j, k, gtg, ee) + &
                 & -gkl(l, f) * xi_in * f2_out(j, k, gtg)
        B_vec(7) = -CONJG(gkl(j, g)) * f2sig_out(l, k, fg2, ef)
        B_vec(8) = -gkl(k, g) * f2sig_out(j, l, gtf, fe) + &
                 & -gkl(l, f) * f2sig_out(j, k, gtg, fe)

        ! Calculate steady states
        CALL MatrixInverseSS(N_mat, Mat, B_vec, f3sig_out(j, k, l, gtgf, :))

        !--------------------------------------------------!
        ! < g^{\dagger}_{j} f^{\dagger}_{k} f_{l} \sigma > !
        !--------------------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((2.d0 * kappa(f) + kappa(g)) - i * (wl(j, g) + wl(k, f) - wl(l, f)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j, g)) * f2sig_out(k, l, ftf, eg) + &
                 & -CONJG(gkl(k, f)) * f2sig_out(j, l, gtf, eg) + &
                 & -gkl(l, f) * f2sig_out(k, j, ftgt, ge)
        B_vec(2) = -CONJG(gkl(j, g)) * f2sig_out(k, l, ftf, ee) + &
                 & -CONJG(gkl(k, f)) * f2sig_out(j, l, gtf, ee) + &
                 & -gkl(l, f) * xi_in * f2sig_out(k, j, ftgt, gf)
        B_vec(3) = -CONJG(gkl(j, g)) * xi_in * f2sig_out(k, l, ftf, fg) + &
                 & -CONJG(gkl(k, f)) * xi_in * f2sig_out(j, l, gtf, fg) + &
                 & -gkl(l, f) * f2sig_out(k, j, ftgt, ee)
        B_vec(4) = Gamma_in * (xi_in ** 2) * f3_out(j, k, l, gtftf) + &
                 & -CONJG(gkl(j, g)) * xi_in * f2sig_out(k, l, ftf, fe) + &
                 & -CONJG(gkl(k, f)) * xi_in * f2sig_out(j, l, gtf, fe) + &
                 & -gkl(l, f) * xi_in * f2sig_out(k, j, ftgt, ef)
        B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f3_out(j, k, l, gtftf) + &
                 & CONJG(gkl(j, g)) * xi_in * f2sig_out(k, l, ftf, gg) + &
                 & CONJG(gkl(j, g)) * xi_in * f2sig_out(k, l, ftf, ee) + &
                 & -CONJG(gkl(j, g)) * xi_in * f2_out(k, l, ftf) + &
                 & CONJG(gkl(k, f)) * xi_in * f2sig_out(j, l, gtf, gg) + &
                 & CONJG(gkl(k, f)) * xi_in * f2sig_out(j, l, gtf, ee) + &
                 & -CONJG(gkl(k, f)) * xi_in * f2_out(j, l, gtf)
        B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f3_out(j, k, l, gtftf) + &
                 & gkl(l, f) * xi_in * f2sig_out(k, j, ftgt, gg) + &
                 & gkl(l, f) * xi_in * f2sig_out(k, j, ftgt, ee) + &
                 & -gkl(l, f) * xi_in * f2_out(k, j, ftgt)
        B_vec(7) = -CONJG(gkl(j, g)) * f2sig_out(k, l, ftf, ef) + &
                 & -CONJG(gkl(k, f)) * f2sig_out(j, l, gtf, ef)
        B_vec(8) = -gkl(l, f) * f2sig_out(k, j, ftgt, fe)

        ! Calculate steady states
        CALL MatrixInverseSS(N_mat, Mat, B_vec, f3sig_out(j, k, l, gtftf, :))

        !--------------------------------------------------!
        ! < f^{\dagger}_{j} g^{\dagger}_{k} g_{l} \sigma > !
        !--------------------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((kappa(f) + 2.0d0 * kappa(g)) - i * (wl(j, f) + wl(k, g) - wl(l, g)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j, f)) * f2sig_out(k, l, gtg, eg) + &
                 & -CONJG(gkl(k, g)) * f2sig_out(j, l, ftg, eg) + &
                 & -gkl(l, g) * f2sig_out(j, k, ftgt, ge)
        B_vec(2) = -CONJG(gkl(j, f)) * f2sig_out(k, l, gtg, ee) + &
                 & -CONJG(gkl(k, g)) * f2sig_out(j, l, ftg, ee) + &
                 & -gkl(l, g) * xi_in * f2sig_out(j, k, ftgt, gf)
        B_vec(3) = -CONJG(gkl(j, f)) * xi_in * f2sig_out(k, l, gtg, fg) + &
                 & -CONJG(gkl(k, g)) * xi_in * f2sig_out(j, l, ftg, fg) + &
                 & -gkl(l, g) * f2sig_out(j, k, ftgt, ee)
        B_vec(4) = Gamma_in * (xi_in ** 2) * f3_out(j, k, l, ftgtg) + &
                 & -CONJG(gkl(j, f)) * xi_in * f2sig_out(k, l, gtg, fe) + &
                 & -CONJG(gkl(k, g)) * xi_in * f2sig_out(j, l, ftg, fe) + &
                 & -gkl(l, g) * xi_in * f2sig_out(j, k, ftgt, ef)
        B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f3_out(j, k, l, ftgtg) + &
                 & CONJG(gkl(j, f)) * xi_in * f2sig_out(k, l, gtg, gg) + &
                 & CONJG(gkl(j, f)) * xi_in * f2sig_out(k, l, gtg, ee) + &
                 & -CONJG(gkl(j, f)) * xi_in * f2_out(k, l, gtg) + &
                 & CONJG(gkl(k, g)) * xi_in * f2sig_out(j, l, ftg, gg) + &
                 & CONJG(gkl(k, g)) * xi_in * f2sig_out(j, l, ftg, ee) + &
                 & -CONJG(gkl(k, g)) * xi_in * f2_out(j, l, ftg)
        B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f3_out(j, k, l, ftgtg) + &
                 & gkl(l, g) * xi_in * f2sig_out(j, k, ftgt, gg) + &
                 & gkl(l, g) * xi_in * f2sig_out(j, k, ftgt, ee) + &
                 & -gkl(l, g) * xi_in * f2_out(j, k, ftgt)
        B_vec(7) = -CONJG(gkl(j, f)) * f2sig_out(k, l, gtg, ef) + &
                 & -CONJG(gkl(k, g)) * f2sig_out(j, l, ftg, ef)
        B_vec(8) = -gkl(l, g) * f2sig_out(j, k, ftgt, fe)

        ! Calculate steady states
        CALL MatrixInverseSS(N_mat, Mat, B_vec, f3sig_out(j, k, l, ftgtg, :))

        ! Close j loop
      END DO
      ! Close k loop
    END DO
    ! Close l loop
  END DO

  ! Cycle through modes
  DO m = -N_in, N_in
    DO l = -N_in, N_in
      DO k = -N_in, N_in
        DO j = -N_in, N_in
          !------------------------------!
          !     FOURTH-ORDER: CAVITY     !
          !------------------------------!
          !-------------------------------------------------!
          ! < a^{\dagger}_{j} b^{\dagger}_{k} b_{l} a_{m} > !
          !-------------------------------------------------!
          f4_out(j, k, l, m) = -CONJG(gkl(j, f)) * f3sig_out(k, l, m, gtgf, eg) + &
                             & -CONJG(gkl(j, f)) * xi_in * f3sig_out(k, l, m, gtgf, fe) + &
                             & -CONJG(gkl(k, g)) * f3sig_out(j, m, l, ftfg, eg) + &
                             & -CONJG(gkl(k, g)) * xi_in * f3sig_out(j, m, l, ftfg, fe) + &
                             & -gkl(l, g) * f3sig_out(k, j, m, gtftf, ge) + &
                             & -gkl(l, g) * xi_in * f3sig_out(k, j, m, gtftf, ef) + &
                             & -gkl(m, f) * f3sig_out(j, k, l, ftgtg, ge) + &
                             & -gkl(m, f) * xi_in * f3sig_out(j, k, l, ftgtg, ef)
          f4_out(j, k, l, m) = f4_out(j, k, l, m) / &
                             & ((2.0d0 * kappa(f) + 2.0d0 * kappa(g)) - &
                             & i * (wl(j, f) + wl(k, g)) + i * (wl(l, g) + wl(m, f)))
          ! Close j loop
        END DO
        ! Close k loop
      END DO
      ! Close l loop
    END DO
    ! Close m loop
  END DO

END SUBROUTINE SteadyStateMoments

!==============================================================================!
!                          G2 CORRELATION SUBROUTINES                          !
!==============================================================================!

! Subroutine to calculate the initial conditions for the auto-correlations
SUBROUTINE G2_InitialConditions(Gamma_in, Omega_in, alpha_in, delta_in, xi_in, &
                              & epsilon_in, N_in, phase_in, &
                              & w0a_in, kappaa_in, dwa_in, &
                              & w0b_in, kappab_in, dwb_in, &
                              & photon_ss_out, B_OG_out, &
                              & sigma_out, f1_out, f1sig_out, f2_out)

  !==============================================================================!
  !                    DEFINING AND DECLARING VARIABLES/ARRAYS                   !
  !==============================================================================!

  IMPLICIT NONE

  !---------------!
  !     INPUT     !
  !---------------!
  ! Atomic decay rate
  REAL(KIND=8), INTENT(IN)                   :: Gamma_in
  ! Driving amplitude
  REAL(KIND=8), INTENT(IN)                   :: Omega_in
  ! Atomic anharmonicity
  REAL(KIND=8), INTENT(IN)                   :: alpha_in
  ! Drive detuning from two-photon resonance
  REAL(KIND=8), INTENT(IN)                   :: delta_in
  ! Dipole moment ratio
  REAL(KIND=8), INTENT(IN)                   :: xi_in

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
  INTEGER, PARAMETER                         :: N_mat = 8
  ! M matrix (filled as transpose)
  COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)   :: Mat, Mat_OG, Mat_inv

  ! Integer indices for sigma operators
  INTEGER, PARAMETER                         :: gg = 1, ge = 2, eg = 3
  INTEGER, PARAMETER                         :: ee = 4, ef = 5, fe = 6
  INTEGER, PARAMETER                         :: gf = 7, fg = 8
  ! Integer indices for: f, f^{\dagger}, f^{2}, f^{\dagger}^{2} ... etc
  INTEGER, PARAMETER                         :: f = 1, ft = 2
  INTEGER, PARAMETER                         :: ftf = 3
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

  !==============================================================================!
  !                        CALCULATE STEADY-STATE MOMENTS                        !
  !==============================================================================!
  CALL SteadyStateMoments(Gamma_in, Omega_in, alpha_in, delta_in, xi_in, &
                        & epsilon_in, N_in, phase_in, &
                        & w0a_in, kappaa_in, dwa_in, &
                        & w0b_in, kappab_in, dwb_in, &
                        & photon_ss_out, sigma_ss, &
                        & f1_ss, f1sig_ss, &
                        & f2_ss, f2sig_ss, &
                        & f3_ss, f3sig_ss, &
                        & f4_ss)

  !==============================================================================!
  !        CALCULATE SECOND-ORDER CORRELATION FUNCTION INITIAL CONDITIONS        !
  !==============================================================================!
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
      B_OG_out(4) = B_OG_out(4) + Gamma_in * (xi_in ** 2) * f2_ss(j, m, ftf)
      B_OG_out(5) = B_OG_out(5) + i * xi_in * 0.5d0 * Omega_in * f2_ss(j, m, ftf)
      B_OG_out(6) = B_OG_out(6) - i * xi_in * 0.5d0 * Omega_in * f2_ss(j, m, ftf)

      ! Close j loop
    END DO
    ! Close m loop
  END DO

END SUBROUTINE G2_InitialConditions

! Subroutine to calculate the time evolution of the g2 correlation
SUBROUTINE G2_CalculateRK4(Gamma_in, Omega_in, alpha_in, delta_in, xi_in, &
                         & epsilon_in, N_in, phase_in, &
                         & w0a_in, kappaa_in, dwa_in, &
                         & w0b_in, kappab_in, dwb_in, &
                         & dt_in, tau_steps_in, &
                         & g2_array_out, WRITE_DATA_IN, filename_data_in)

  !==============================================================================!
  !                    DEFINING AND DECLARING VARIABLES/ARRAYS                   !
  !==============================================================================!

  IMPLICIT NONE

  !---------------!
  !     INPUT     !
  !---------------!
  ! Atomic decay rate
  REAL(KIND=8), INTENT(IN)                           :: Gamma_in
  ! Driving amplitude
  REAL(KIND=8), INTENT(IN)                           :: Omega_in
  ! Atomic anharmonicity
  REAL(KIND=8), INTENT(IN)                           :: alpha_in
  ! Drive detuning from two-photon resonance
  REAL(KIND=8), INTENT(IN)                           :: delta_in
  ! Dipole moment ratio
  REAL(KIND=8), INTENT(IN)                           :: xi_in

  ! Filter parameter stuff
  ! Percentage of fluorecence aimed at cavity
  REAL(KIND=8), INTENT(IN)                           :: epsilon_in
  ! Number of mode either side of w0, 2N + 1 total mode
  INTEGER, INTENT(IN)                                :: N_in
  ! Phase modulation of mode coupling
  REAL(KIND=8), INTENT(IN)                           :: phase_in
  
  ! Central mode frequency of the filter cavity, with N mode frequencies either side
  REAL(KIND=8), INTENT(IN)                           :: w0a_in, w0b_in
  ! Cavity linewidth/transmission of cavity mode
  REAL(KIND=8), INTENT(IN)                           :: kappaa_in, kappab_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                           :: dwa_in, dwb_in

  ! Time stuff
  ! Time step
  REAL(KIND=8), INTENT(IN)                           :: dt_in
  ! Maxi_inmum number of steps to integrate for
  INTEGER, INTENT(IN)                                :: tau_steps_in

  ! Data stuff
  ! Boolean for writing data
  LOGICAL, INTENT(IN)                                :: WRITE_DATA_IN
  ! Filename for writing data to
  CHARACTER(LEN=*), INTENT(IN)                       :: filename_data_in

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Data array
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: g2_array_out

  !------------------------------------!
  !     MOMENT EQUATION ARRAY STUFF    !
  !------------------------------------!
  ! Dimension of M matrix
  INTEGER, PARAMETER                                 :: N_mat = 8
  ! M matrix (filled as transpose)
  COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)           :: Mat, Mat_OG
  ! Non-homogeneous vector
  COMPLEX(KIND=8), DIMENSION(N_mat)                  :: B_OG, B_vec

  ! Integer indices for sigma operators
  INTEGER, PARAMETER                                 :: gg = 1, ge = 2, eg = 3
  INTEGER, PARAMETER                                 :: ee = 4, ef = 5, fe = 6
  INTEGER, PARAMETER                                 :: gf = 7, fg = 8
  ! Integer indices for: f, f^{\dagger}, f^{2}, f^{\dagger}^{2} ... etc
  INTEGER, PARAMETER                                 :: f = 1, ft = 2

  ! Time integration arrays
  ! First-order moments: Atomic equations (< \sigma >)
  COMPLEX(KIND=8), DIMENSION(N_mat)                  :: sigma
  COMPLEX(KIND=8), DIMENSION(N_mat)                  :: k1_sigma, k2_sigma, k3_sigma, k4_sigma
  ! First-order moments: Cavity (< a >, < f^{\dagger} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2)          :: f1
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2)          :: k1_f1, k2_f1, k3_f1, k4_f1
  ! Second-order moments: Cavity and atom (< a \sigma >, < f^{\dagger} \sigma >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N_mat)   :: f1sig
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N_mat)   :: k1_f1sig, k2_f1sig, k3_f1sig, k4_f1sig
  ! Second-order moments: Cavity (< f^{\dagger} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in) :: f2
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in) :: k1_f2, k2_f2, k3_f2, k4_f2

  !----------------------------!
  !     OTHER USEFUL STUFF     !
  !----------------------------!
  ! Time step integer
  INTEGER                                            :: t
  ! Integer counter
  INTEGER                                            :: j, k, l, m, x
  ! Sample rate for state populations
  INTEGER                                            :: sample_rate
  ! Imaginary i
  COMPLEX(KIND=8), PARAMETER                         :: i = CMPLX(0.0d0, 1.0d0, 8)
  ! pi
  REAL(KIND=8), PARAMETER                            :: pi = 3.1415926535897932384d0
  ! 1 / 6
  REAL(KIND=8), PARAMETER                            :: xis = 1.0d0 / 6.0d0
  ! kappa, w0, and dw values
  REAL(KIND=8)                                       :: kappa, dw, w0
  ! List of Delta values
  REAL(KIND=8), DIMENSION(-N_in:N_in)                :: wl
  ! List of mode dependent cascade coupling values
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in)             :: gkl
  ! Blackman window coefficient
  REAL(KIND=8)                                       :: blackman
  ! Steady state photon number
  REAL(KIND=8), DIMENSION(2)                         :: photon_ss
  ! Complex data
  COMPLEX(KIND=8)                                    :: moment_out

  !==============================================================================!
  !                DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS                !
  !==============================================================================!
  kappa = kappab_in
  w0 = w0b_in
  dw = dwb_in

  !------------------------!
  !     BLOCH MATRIX M     !
  !------------------------!
  Mat_OG = 0.0d0
  ! Row 1: d/dt |g><g|
  Mat_OG(1, 2) = -i * 0.5d0 * Omega_in
  Mat_OG(1, 3) = i * 0.5d0 * Omega_in
  Mat_OG(1, 4) = Gamma_in
  ! Row 2: d/dt |g><e|
  Mat_OG(2, 1) = -i * 0.5d0 * Omega_in
  Mat_OG(2, 2) = -(0.5d0 * Gamma_in - i * ((0.5d0 * alpha_in) + delta_in))
  Mat_OG(2, 4) = i * 0.5d0 * Omega_in
  Mat_OG(2, 5) = Gamma_in * xi_in
  Mat_OG(2, 7) = -i * xi_in * 0.5d0 * Omega_in
  ! Row 3: d/dt |e><g|
  Mat_OG(3, 1) = i * 0.5d0 * Omega_in
  Mat_OG(3, 3) = -(0.5d0 * Gamma_in + i * ((0.5d0 * alpha_in) + delta_in))
  Mat_OG(3, 4) = -i * 0.5d0 * Omega_in
  Mat_OG(3, 6) = Gamma_in * xi_in
  Mat_OG(3, 8) = i * xi_in * 0.5d0 * Omega_in
  ! Row 4: d/dt |e><e|
  Mat_OG(4, 1) = -Gamma_in * (xi_in ** 2)
  Mat_OG(4, 2) = i * 0.5d0 * Omega_in
  Mat_OG(4, 3) = -i * 0.5d0 * Omega_in
  Mat_OG(4, 4) = -Gamma_in * (1.0d0 + (xi_in ** 2))
  Mat_OG(4, 5) = -i * xi_in * 0.5d0 * Omega_in
  Mat_OG(4, 6) = i * xi_in * 0.5d0 * Omega_in
  ! Row 5: d/dt |e><f|
  Mat_OG(5, 1) = -i * xi_in * 0.5d0 * Omega_in
  Mat_OG(5, 4) = -i * xi_in * Omega_in
  Mat_OG(5, 5) = -(0.5d0 * Gamma_in * (1.0d0 + (xi_in ** 2)) + i * ((0.5d0 * alpha_in) - delta_in))
  Mat_OG(5, 7) = i * 0.5d0 * Omega_in
  ! Row 6: d/dt |f><e|
  Mat_OG(6, 1) = i * xi_in * 0.5d0 * Omega_in
  Mat_OG(6, 4) = i * xi_in * Omega_in
  Mat_OG(6, 6) = -(0.5d0 * Gamma_in * (1.0d0 + (xi_in ** 2)) - i * ((0.5d0 * alpha_in) - delta_in))
  Mat_OG(6, 8) = -i * 0.5d0 * Omega_in
  ! Row 7: d/dt |g><f|
  Mat_OG(7, 2) = -i * xi_in * 0.5d0 * Omega_in
  Mat_OG(7, 5) = i * 0.5d0 * Omega_in
  Mat_OG(7, 7) = -(0.5d0 * Gamma_in * (xi_in ** 2) - 2.0d0 * i * delta_in)
  ! Row 8: d/dt |g><f|
  Mat_OG(8, 3) = i * xi_in * 0.5d0 * Omega_in
  Mat_OG(8, 6) = -i * 0.5d0 * Omega_in
  Mat_OG(8, 8) = -(0.5d0 * Gamma_in * (xi_in ** 2) + 2.0d0 * i * delta_in)

  !--------------------------------!
  !     NON-HOMOGENEOUS VECTOR     !
  !--------------------------------!
  B_OG = 0.0d0
  B_OG(4) = Gamma_in * (xi_in ** 2)
  B_OG(5) = i * xi_in * 0.5d0 * Omega_in
  B_OG(6) = -i * xi_in * 0.5d0 * Omega_in

  !---------------------------------------------!
  !     RESONANCES (wj) AND COUPLINGS (E_j)     !
  !---------------------------------------------!
  ! Allocate array of Delta and gka values
  wl = 0.0d0
  gkl = 0.0d0
  DO j = -N_in, N_in
    IF (N_in == 0) THEN
      wl(j) = w0
      gkl(j) = DSQRT(epsilon_in * Gamma_in * kappa)
    ELSE
      wl(j) = w0 + DBLE(j) * dw
      ! Blackman window coefficient
      blackman = 1.0d0
      ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
      !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))
      ! Mode dependent phase difference
      gkl(j) = DSQRT((epsilon_in / DBLE(2*N_in + 1)) * Gamma_in * kappa) * blackman * &
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

  !==============================================================================!
  !                         CALCULATE INITIAL CONDITIONS                         !
  !==============================================================================!
  CALL G2_InitialConditions(Gamma_in, Omega_in, alpha_in, delta_in, xi_in, &
                          & epsilon_in, N_in, phase_in, &
                          & w0a_in, kappaa_in, dwa_in, &
                          & w0b_in, kappab_in, dwb_in, &
                          & photon_ss, B_OG, &
                          & sigma, f1, f1sig, f2)

  !==============================================================================!
  !                  CALCULATE SECOND-ORDER CORRELATION FUNCTION                 !
  !==============================================================================!
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
                  & -dt_in * gkl(j) * sigma(ge) + &
                  & -dt_in * gkl(j) * xi_in * sigma(ef)

      k2_f1(j, f) = -dt_in * (kappa + i * wl(j)) * (f1(j, f) + 0.5d0 * k1_f1(j, f)) + &
                  & -dt_in * gkl(j) * (sigma(ge) + 0.5d0 * k1_sigma(ge)) + &
                  & -dt_in * gkl(j) * xi_in * (sigma(ef) + 0.5d0 * k1_sigma(ef))

      k3_f1(j, f) = -dt_in * (kappa + i * wl(j)) * (f1(j, f) + 0.5d0 * k2_f1(j, f)) + &
                  & -dt_in * gkl(j) * (sigma(ge) + 0.5d0 * k2_sigma(ge)) + &
                  & -dt_in * gkl(j) * xi_in * (sigma(ef) + 0.5d0 * k2_sigma(ef))

      k4_f1(j, f) = -dt_in * (kappa + i * wl(j)) * (f1(j, f) + k3_f1(j, f)) + &
                  & -dt_in * gkl(j) * (sigma(ge) + k3_sigma(ge)) + &
                  & -dt_in * gkl(j) * xi_in * (sigma(ef) + k3_sigma(ef))

      !---------------------!
      ! < f^{\dagger}_{j} > !
      !---------------------!
      k1_f1(j, ft) = -dt_in * (kappa - i * wl(j)) * f1(j, ft) + &
                   & -dt_in * CONJG(gkl(j)) * sigma(eg) + &
                   & -dt_in * CONJG(gkl(j)) * xi_in * sigma(fe)

      k2_f1(j, ft) = -dt_in * (kappa - i * wl(j)) * (f1(j, ft) + 0.5d0 * k1_f1(j, ft)) + &
                   & -dt_in * CONJG(gkl(j)) * (sigma(eg) + 0.5d0 * k1_sigma(eg)) + &
                   & -dt_in * CONJG(gkl(j)) * xi_in * (sigma(fe) + 0.5d0 * k1_sigma(fe))

      k3_f1(j, ft) = -dt_in * (kappa - i * wl(j)) * (f1(j, ft) + 0.5d0 * k2_f1(j, ft)) + &
                   & -dt_in * CONJG(gkl(j)) * (sigma(eg) + 0.5d0 * k2_sigma(eg)) + &
                   & -dt_in * CONJG(gkl(j)) * xi_in * (sigma(fe) + 0.5d0 * k2_sigma(fe))

      k4_f1(j, ft) = -dt_in * (kappa - i * wl(j)) * (f1(j, ft) + k3_f1(j, ft)) + &
                   & -dt_in * CONJG(gkl(j)) * (sigma(eg) + k3_sigma(eg)) + &
                   & -dt_in * CONJG(gkl(j)) * xi_in * (sigma(fe) + k3_sigma(fe))

      !---------------------------------------!
      !     SECOND-ORDER: CAVITY AND ATOM     !
      !---------------------------------------!
      ! Using matrix multiplication, we add to the Lindblad matrix Mat_OG and the
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
      B_vec(1) = -gkl(j) * sigma(ge)
      B_vec(2) = -gkl(j) * xi_in * sigma(gf)
      B_vec(3) = -gkl(j) * sigma(ee)
      B_vec(4) = Gamma_in * (xi_in ** 2) * f1(j, f) + &
               & -gkl(j) * xi_in * sigma(ef)
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f1(j, f)
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f1(j, f) + &
               & gkl(j) * xi_in * sigma(gg) + &
               & gkl(j) * xi_in * sigma(ee) + &
               & -gkl(j) * xi_in * photon_ss(f)
      B_vec(7) = 0.0d0
      B_vec(8) = -gkl(j) * sigma(fe)
      ! Calculate k1
      k1_f1sig(j, f, :) = dt_in * (MATMUL(Mat, f1sig(j, f, :)) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -gkl(j) * (sigma(ge) + 0.5d0 * k1_sigma(ge))
      B_vec(2) = -gkl(j) * xi_in * (sigma(gf) + 0.5d0 * k1_sigma(gf))
      B_vec(3) = -gkl(j) * (sigma(ee) + 0.5d0 * k1_sigma(ee))
      B_vec(4) = Gamma_in * (xi_in ** 2) * (f1(j, f) + 0.5d0 * k1_f1(j, f)) + &
               & -gkl(j) * xi_in * (sigma(ef) + 0.5d0 * k1_sigma(ef))
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * (f1(j, f) + 0.5d0 * k1_f1(j, f))
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * (f1(j, f) + 0.5d0 * k1_f1(j, f)) + &
               & gkl(j) * xi_in * (sigma(gg) + 0.5d0 * k1_sigma(gg)) + &
               & gkl(j) * xi_in * (sigma(ee) + 0.5d0 * k1_sigma(ee)) + &
               & -gkl(j) * xi_in * photon_ss(f)
      B_vec(7) = 0.0d0
      B_vec(8) = -gkl(j) * (sigma(fe) + 0.5d0 * k1_sigma(fe))
      ! Calculate k2
      k2_f1sig(j, f, :) = dt_in * (MATMUL(Mat, (f1sig(j, f, :) + 0.5d0 * k1_f1sig(j, f, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -gkl(j) * (sigma(ge) + 0.5d0 * k2_sigma(ge))
      B_vec(2) = -gkl(j) * xi_in * (sigma(gf) + 0.5d0 * k2_sigma(gf))
      B_vec(3) = -gkl(j) * (sigma(ee) + 0.5d0 * k2_sigma(ee))
      B_vec(4) = Gamma_in * (xi_in ** 2) * (f1(j, f) + 0.5d0 * k2_f1(j, f)) + &
               & -gkl(j) * xi_in * (sigma(ef) + 0.5d0 * k2_sigma(ef))
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * (f1(j, f) + 0.5d0 * k2_f1(j, f))
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * (f1(j, f) + 0.5d0 * k2_f1(j, f)) + &
               & gkl(j) * xi_in * (sigma(gg) + 0.5d0 * k2_sigma(gg)) + &
               & gkl(j) * xi_in * (sigma(ee) + 0.5d0 * k2_sigma(ee)) + &
               & -gkl(j) * xi_in * photon_ss(f)
      B_vec(7) = 0.0d0
      B_vec(8) = -gkl(j) * (sigma(fe) + 0.5d0 * k2_sigma(fe))
      ! Calculate k3
      k3_f1sig(j, f, :) = dt_in * (MATMUL(Mat, (f1sig(j, f, :) + 0.5d0 * k2_f1sig(j, f, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -gkl(j) * (sigma(ge) + k3_sigma(ge))
      B_vec(2) = -gkl(j) * xi_in * (sigma(gf) + k3_sigma(gf))
      B_vec(3) = -gkl(j) * (sigma(ee) + k3_sigma(ee))
      B_vec(4) = Gamma_in * (xi_in ** 2) * (f1(j, f) + k3_f1(j, f)) + &
               & -gkl(j) * xi_in * (sigma(ef) + k3_sigma(ef))
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * (f1(j, f) + k3_f1(j, f))
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * (f1(j, f) + k3_f1(j, f)) + &
               & gkl(j) * xi_in * (sigma(gg) + k3_sigma(gg)) + &
               & gkl(j) * xi_in * (sigma(ee) + k3_sigma(ee)) + &
               & -gkl(j) * xi_in * photon_ss(f)
      B_vec(7) = 0.0d0
      B_vec(8) = -gkl(j) * (sigma(fe) + k3_sigma(fe))
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
      B_vec(1) = -CONJG(gkl(j)) * sigma(eg)
      B_vec(2) = -CONJG(gkl(j)) * sigma(ee)
      B_vec(3) = -CONJG(gkl(j)) * xi_in * sigma(fg)
      B_vec(4) = Gamma_in * (xi_in ** 2) * f1(j, ft) + &
               & -CONJG(gkl(j)) * xi_in * sigma(fe)
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f1(j, ft) + &
               & CONJG(gkl(j)) * xi_in * sigma(gg) + &
               & CONJG(gkl(j)) * xi_in * sigma(ee) + &
               & -CONJG(gkl(j)) * xi_in * photon_ss(f)
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f1(j, ft)
      B_vec(7) = -CONJG(gkl(j)) * sigma(ef)
      B_vec(8) = 0.0d0
      ! Calculate k1
      k1_f1sig(j, ft, :) = dt_in * (MATMUL(Mat, f1sig(j, ft, :)) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j)) * (sigma(eg) + 0.5d0 * k1_sigma(eg))
      B_vec(2) = -CONJG(gkl(j)) * (sigma(ee) + 0.5d0 * k1_sigma(ee))
      B_vec(3) = -CONJG(gkl(j)) * xi_in * (sigma(fg) + 0.5d0 * k1_sigma(fg))
      B_vec(4) = Gamma_in * (xi_in ** 2) * (f1(j, ft) + 0.5d0 * k1_f1(j, ft)) + &
               & -CONJG(gkl(j)) * xi_in * (sigma(fe) + 0.5d0 * k1_sigma(fe))
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * (f1(j, ft) + 0.5d0 * k1_f1(j, ft)) + &
               & CONJG(gkl(j)) * xi_in * (sigma(gg) + 0.5d0 * k1_sigma(gg)) + &
               & CONJG(gkl(j)) * xi_in * (sigma(ee) + 0.5d0 * k1_sigma(ee)) + &
               & -CONJG(gkl(j)) * xi_in * photon_ss(f)
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * (f1(j, ft) + 0.5d0 * k1_f1(j, ft))
      B_vec(7) = -CONJG(gkl(j)) * (sigma(ef) + 0.5d0 * k1_sigma(ef))
      B_vec(8) = 0.0d0
      ! Calculate k2
      k2_f1sig(j, ft, :) = dt_in * (MATMUL(Mat, (f1sig(j, ft, :) + 0.5d0 * k1_f1sig(j, ft, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j)) * (sigma(eg) + 0.5d0 * k2_sigma(eg))
      B_vec(2) = -CONJG(gkl(j)) * (sigma(ee) + 0.5d0 * k2_sigma(ee))
      B_vec(3) = -CONJG(gkl(j)) * xi_in * (sigma(fg) + 0.5d0 * k2_sigma(fg))
      B_vec(4) = Gamma_in * (xi_in ** 2) * (f1(j, ft) + 0.5d0 * k2_f1(j, ft)) + &
               & -CONJG(gkl(j)) * xi_in * (sigma(fe) + 0.5d0 * k2_sigma(fe))
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * (f1(j, ft) + 0.5d0 * k2_f1(j, ft)) + &
               & CONJG(gkl(j)) * xi_in * (sigma(gg) + 0.5d0 * k2_sigma(gg)) + &
               & CONJG(gkl(j)) * xi_in * (sigma(ee) + 0.5d0 * k2_sigma(ee)) + &
               & -CONJG(gkl(j)) * xi_in * photon_ss(f)
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * (f1(j, ft) + 0.5d0 * k2_f1(j, ft))
      B_vec(7) = -CONJG(gkl(j)) * (sigma(ef) + 0.5d0 * k2_sigma(ef))
      B_vec(8) = 0.0d0
      ! Calculate k3
      k3_f1sig(j, ft, :) = dt_in * (MATMUL(Mat, (f1sig(j, ft, :) + 0.5d0 * k2_f1sig(j, ft, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j)) * (sigma(eg) + k3_sigma(eg))
      B_vec(2) = -CONJG(gkl(j)) * (sigma(ee) + k3_sigma(ee))
      B_vec(3) = -CONJG(gkl(j)) * xi_in * (sigma(fg) + k3_sigma(fg))
      B_vec(4) = Gamma_in * (xi_in ** 2) * (f1(j, ft) + k3_f1(j, ft)) + &
               & -CONJG(gkl(j)) * xi_in * (sigma(fe) + k3_sigma(fe))
      B_vec(5) = i * xi_in * 0.5d0 * Omega_in * (f1(j, ft) + k3_f1(j, ft)) + &
               & CONJG(gkl(j)) * xi_in * (sigma(gg) + k3_sigma(gg)) + &
               & CONJG(gkl(j)) * xi_in * (sigma(ee) + k3_sigma(ee)) + &
               & -CONJG(gkl(j)) * xi_in * photon_ss(f)
      B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * (f1(j, ft) + k3_f1(j, ft))
      B_vec(7) = -CONJG(gkl(j)) * (sigma(ef) + k3_sigma(ef))
      B_vec(8) = 0.0d0
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
                    & -dt_in * CONJG(gkl(j)) * f1sig(k, f, eg) + &
                    & -dt_in * CONJG(gkl(j)) * xi_in * f1sig(k, f, fe) + &
                    & -dt_in * gkl(k) * f1sig(j, ft, ge) + &
                    & -dt_in * gkl(k) * xi_in * f1sig(j, ft, ef)

        k2_f2(j, k) = -dt_in * (2.0d0 * kappa - i * (wl(j) - wl(k))) * (f2(j, k) + 0.5d0 * k1_f2(j, k)) + &
                    & -dt_in * CONJG(gkl(j)) * (f1sig(k, f, eg) + 0.5d0 * k1_f1sig(k, f, eg)) + &
                    & -dt_in * CONJG(gkl(j)) * xi_in * (f1sig(k, f, fe) + 0.5d0 * k1_f1sig(k, f, fe)) + &
                    & -dt_in * gkl(k) * (f1sig(j, ft, ge) + 0.5d0 * k1_f1sig(j, ft, ge)) + &
                    & -dt_in * gkl(k) * xi_in * (f1sig(j, ft, ef) + 0.5d0 * k1_f1sig(j, ft, ef))

        k3_f2(j, k) = -dt_in * (2.0d0 * kappa - i * (wl(j) - wl(k))) * (f2(j, k) + 0.5d0 * k2_f2(j, k)) + &
                    & -dt_in * CONJG(gkl(j)) * (f1sig(k, f, eg) + 0.5d0 * k2_f1sig(k, f, eg)) + &
                    & -dt_in * CONJG(gkl(j)) * xi_in * (f1sig(k, f, fe) + 0.5d0 * k2_f1sig(k, f, fe)) + &
                    & -dt_in * gkl(k) * (f1sig(j, ft, ge) + 0.5d0 * k2_f1sig(j, ft, ge)) + &
                    & -dt_in * gkl(k) * xi_in * (f1sig(j, ft, ef) + 0.5d0 * k2_f1sig(j, ft, ef))

        k4_f2(j, k) = -dt_in * (2.0d0 * kappa - i * (wl(j) - wl(k))) * (f2(j, k) + k3_f2(j, k)) + &
                    & -dt_in * CONJG(gkl(j)) * (f1sig(k, f, eg) + k3_f1sig(k, f, eg)) + &
                    & -dt_in * CONJG(gkl(j)) * xi_in * (f1sig(k, f, fe) + k3_f1sig(k, f, fe)) + &
                    & -dt_in * gkl(k) * (f1sig(j, ft, ge) + k3_f1sig(j, ft, ge)) + &
                    & -dt_in * gkl(k) * xi_in * (f1sig(j, ft, ef) + k3_f1sig(j, ft, ef))

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