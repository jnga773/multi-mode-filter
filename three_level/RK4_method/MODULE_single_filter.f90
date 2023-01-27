! This module contains the subroutines used in any of the single-filter
! programs [just_steady_states.f90, g1_RK4.f90, g2_RK4.f90].

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
! - G1_InitialConditions: Calculates the initial conditions for the first-
!                         order correlation function.
!
! - G1_CalculateRK4: Calculates the time evolution of the first-order
!                    correlation function using Runge-Kutta 4th Order.
!
! - G2_InitialConditions: Calculates the initial conditions for the second-
!                         order correlation function.
!
! - G2_CalculateRK4: Calculates the time evolution of the second-order
!                    correlation function using Runge-Kutta 4th Order.
!
! - MatrixInverseSS : Calculates the steady state for a system of coupled
!                     differential equations using the inverse matrix method.

! This file must be added to the compilation command when compiling any of the 
! single-filter programs. Eg, with Intel oneAPI or GFORTRAN:
!     (IFORT): ifort -qmkl ./[FILENAME].f90 ./MODULE_single_filter.f90
!  (GFORTRAN): gfortran ./[FILENAME].f90 ./MODULE_single_filter.f90
!                -I/path/to/LAPACK -L/path/to/LAPACK -llapack -lblas

MODULE SINGLE_FILTER_SUBROUTINES

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
                            & w0_in, kappa_in, dw_in, &
                            & photon_out, sigma_out, G2_SS, &
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
  REAL(KIND=8), INTENT(IN)                 :: Gamma_in
  ! Driving amplitude
  REAL(KIND=8), INTENT(IN)                 :: Omega_in
  ! Atomic anharmonicity
  REAL(KIND=8), INTENT(IN)                 :: alpha_in
  ! Drive detuning from two-photon resonance
  REAL(KIND=8), INTENT(IN)                 :: delta_in
  ! Dipole moment ratio
  REAL(KIND=8), INTENT(IN)                 :: xi_in

  ! Filter parameter stuff
  ! Percentage of fluorecence aimed at cavity
  REAL(KIND=8), INTENT(IN)                 :: epsilon_in
  ! Number of mode either side of w0, 2N + 1 total mode
  INTEGER, INTENT(IN)                      :: N_in
  ! Phase modulation of mode coupling
  REAL(KIND=8), INTENT(IN)                 :: phase_in
  ! Central mode frequency of the filter cavity, with N mode frequencies either side
  REAL(KIND=8), INTENT(IN)                 :: w0_in
  ! Cavity linewidth/transmission of cavity mode
  REAL(KIND=8), INTENT(IN)                 :: kappa_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                 :: dw_in

  ! Calculate steady states for G1 or G2
  LOGICAL, INTENT(IN)                      :: G2_SS

  !------------------------------------!
  !     MOMENT EQUATION ARRAY STUFF    !
  !------------------------------------!
  ! Dimension of M matrix
  INTEGER, PARAMETER                       :: N_mat = 8
  ! M matrix (filled as transpose)
  COMPLEX(KIND=8), DIMENSION(N_mat, N_mat) :: Mat, Mat_OG
  ! Non-homogeneous vector
  COMPLEX(KIND=8), DIMENSION(N_mat)        :: B_vec, B_OG

  ! Integer indices for sigma operators
  INTEGER, PARAMETER                       :: gg = 1, ge = 2, eg = 3
  INTEGER, PARAMETER                       :: ee = 4, ef = 5, fe = 6
  INTEGER, PARAMETER                       :: gf = 7, fg = 8
  ! Integer indices for: a, f^{\dagger}, f^{\dagger} a
  INTEGER, PARAMETER                       :: f = 1, ft = 2
  INTEGER, PARAMETER                       :: ff = 1, ftf = 2, ft2 = 3
  INTEGER, PARAMETER                       :: ftf2 = 1, ft2f = 2

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Steady state photon number
  REAL(KIND=8), INTENT(OUT)                                                               :: photon_out

  ! Steady state arrays
  ! First-order moments: Atomic equations (< \sigma >)
  COMPLEX(KIND=8), DIMENSION(N_mat), INTENT(OUT)                                          :: sigma_out
  ! First-order moments: Cavity (< a >, < f^{\dagger} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2), INTENT(OUT)                                  :: f1_out
  ! Second-order moments: Cavity and atom (< a \sigma >, < f^{\dagger} \sigma >
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N_mat), INTENT(OUT)                           :: f1sig_out
  ! Second-order moments: Cavity (< f^{\dagger} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3), INTENT(OUT)                      :: f2_out
  ! Third-order moments: Cavity and atom (< a^{2} \sigma >, < a^{\dagger 2} \sigma >, < f^{\dagger} a \sigma >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3, N_mat), INTENT(OUT)               :: f2sig_out
  ! Third-order moments: Cavity (< a^{2} f^{\dagger} >, < a^{\dagger 2} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, 2), INTENT(OUT)          :: f3_out
  ! Fourth-order moments: Cavity and atom ( < f^{\dagger} a^{2} \sigma >, < a^{\dagger 2} a \sigma >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, 2, N_mat), INTENT(OUT)   :: f3sig_out
  ! Fourth-order moments: Cavity (< a^{\dagger 2} a^{2} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, -N_in:N_in), INTENT(OUT) :: f4_out

  !----------------------------!
  !     OTHER USEFUL STUFF     !
  !----------------------------!
  ! ! Integer counter
  INTEGER                                  :: j, k, l, m, x
  ! Imaginary i
  COMPLEX(KIND=8), PARAMETER               :: i = CMPLX(0.0d0, 1.0d0, 8)
  ! pi
  REAL(KIND=8), PARAMETER                  :: pi = 3.1415926535897932384d0
  ! List of Delta values
  REAL(KIND=8), DIMENSION(-N_in:N_in)      :: wl
  ! List of mode dependent cascade coupling values
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in)   :: gkl
  ! Blackman window coefficient
  REAL(KIND=8)                             :: blackman
  ! Complex temporary values
  COMPLEX(KIND=8)                          :: moment_out

  !==============================================================================!
  !                DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS                !
  !==============================================================================!
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
      wl(j) = w0_in
      gkl(j) = DSQRT(epsilon_in * Gamma_in * kappa_in)
    ELSE
      wl(j) = w0_in + DBLE(j) * dw_in
      ! Blackman window coefficient
      blackman = 1.0d0
      ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
      !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))
      ! Mode dependent phase difference
      gkl(j) = DSQRT((epsilon_in / DBLE(2*N_in + 1)) * Gamma_in * kappa_in) * blackman * &
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
    f1_out(j, f) = -gkl(j) * sigma_out(ge) + &
                 & -gkl(j) * xi_in * sigma_out(ef)
    f1_out(j, f) = f1_out(j, f) / &
                 & (kappa_in + i * wl(j))

    !---------------------!
    ! < f^{\dagger}_{j} > !
    !---------------------!
    f1_out(j, ft) = -CONJG(gkl(j)) * sigma_out(eg) + &
                  & -CONJG(gkl(j)) * xi_in * sigma_out(fe)
    f1_out(j, ft) = f1_out(j, ft) / &
                  & (kappa_in - i * wl(j))

    !---------------------------------------!
    !     SECOND-ORDER: CAVITY AND ATOM     !
    !---------------------------------------!
    !------------------!
    ! < f_{j} \sigma > !
    !------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa_in + i * wl(j))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -gkl(j) * sigma_out(ge)
    B_vec(2) = -gkl(j) * xi_in * sigma_out(gf)
    B_vec(3) = -gkl(j) * sigma_out(ee)
    B_vec(4) = Gamma_in * (xi_in ** 2) * f1_out(j, f) + &
             & -gkl(j) * xi_in * sigma_out(ef)
    B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f1_out(j, f)
    B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f1_out(j, f) + &
             & gkl(j) * xi_in * sigma_out(gg) + &
             & gkl(j) * xi_in * sigma_out(ee) + &
             & -gkl(j) * xi_in
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j) * sigma_out(fe)

    ! Calculate steady states
    CALL MatrixInverseSS(N_mat, Mat, B_vec, f1sig_out(j, f, :))

    !----------------------------!
    ! < f^{\dagger}_{j} \sigma > !
    !----------------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa_in - i * wl(j))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -CONJG(gkl(j)) * sigma_out(eg)
    B_vec(2) = -CONJG(gkl(j)) * sigma_out(ee)
    B_vec(3) = -CONJG(gkl(j)) * xi_in * sigma_out(fg)
    B_vec(4) = Gamma_in * (xi_in ** 2) * f1_out(j, ft) + &
             & -CONJG(gkl(j)) * xi_in * sigma_out(fe)
    B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f1_out(j, ft) + &
             & CONJG(gkl(j)) * xi_in * sigma_out(gg) + &
             & CONJG(gkl(j)) * xi_in * sigma_out(ee) + &
             & -CONJG(gkl(j)) * xi_in
    B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f1_out(j, ft)
    B_vec(7) = -CONJG(gkl(j)) * sigma_out(ef)
    B_vec(8) = 0.0d0

    ! Calculate steady states
    CALL MatrixInverseSS(N_mat, Mat, B_vec, f1sig_out(j, ft, :))

    ! Close j loop
  END DO

  moment_out = 0.0d0
  ! Cycle through modes
  DO k = -N_in, N_in
    DO j = -N_in, N_in
      !------------------------------!
      !     SECOND-ORDER: CAVITY     !
      !------------------------------!
      !-----------------!
      ! < f_{j} f_{k} > !
      !-----------------!
      f2_out(j, k, ff) = -gkl(j) * f1sig_out(k, f, ge) + &
                       & -gkl(j) * xi_in * f1sig_out(k, f, ef) + &
                       & -gkl(k) * f1sig_out(j, f, ge) + &
                       & -gkl(k) * xi_in * f1sig_out(j, f, ef)
      f2_out(j, k, ff) = f2_out(j, k, ff) / &
                       & (2.0d0 * kappa_in + i * (wl(j) + wl(k)))

      !-------------------------------------!
      ! < f^{\dagger}_{j} f^{\dagger}_{k} > !
      !-------------------------------------!
      f2_out(j, k, ft2) = -CONJG(gkl(j)) * f1sig_out(k, ft, eg) + &
                        & -CONJG(gkl(j)) * xi_in * f1sig_out(k, ft, fe) + &
                        & -CONJG(gkl(k)) * f1sig_out(j, ft, eg) + &
                        & -CONJG(gkl(k)) * xi_in * f1sig_out(j, ft, fe)
      f2_out(j, k, ft2) = f2_out(j, k, ft2) / &
                        & (2.0d0 * kappa_in - i * (wl(j) + wl(k)))

      !---------------------------!
      ! < f^{\dagger}_{j} f_{k} > !
      !---------------------------!
      f2_out(j, k, ftf) = -CONJG(gkl(j)) * f1sig_out(k, f, eg) + &
                        & -CONJG(gkl(j)) * xi_in * f1sig_out(k, f, fe) + &
                        & -gkl(k) * f1sig_out(j, ft, ge) + &
                        & -gkl(k) * xi_in * f1sig_out(j, ft, ef)
      f2_out(j, k, ftf) = f2_out(j, k, ftf) / &
                        & (2.0d0 * kappa_in - i * (wl(j) - wl(k)))

      ! Update photon number
      moment_out = moment_out + f2_out(j, k, ftf)

      IF (G2_SS .EQV. .TRUE.) THEN
        !--------------------------------------!
        !     THIRD-ORDER: CAVITY AND ATOM     !
        !--------------------------------------!
        !------------------------!
        ! < f_{j} f_{k} \sigma > !
        !------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa_in) + i * (wl(j) + wl(k)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -gkl(j) * f1sig_out(k, f, ge) + &
                & -gkl(k) * f1sig_out(j, f, ge)
        B_vec(2) = -gkl(j) * xi_in * f1sig_out(k, f, gf) + &
                & -gkl(k) * xi_in * f1sig_out(j, f, gf)
        B_vec(3) = -gkl(j) * f1sig_out(k, f, ee) + &
                & -gkl(k) * f1sig_out(j, f, ee)
        B_vec(4) = Gamma_in * (xi_in ** 2) * f2_out(j, k, ff) + &
                & -gkl(j) * xi_in * f1sig_out(k, f, ef) + &
                & -gkl(k) * xi_in * f1sig_out(j, f, ef)
        B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, ff)
        B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, f) + &
                & gkl(j) * xi_in * f1sig_out(k, f, gg) + &
                & gkl(j) * xi_in * f1sig_out(k, f, ee) + &
                & -gkl(j) * xi_in * f1_out(k, f) + &
                & gkl(k) * xi_in * f1sig_out(j, f, gg) + &
                & gkl(k) * xi_in * f1sig_out(j, f, ee) + &
                & -gkl(k) * xi_in * f1_out(j, f)
        B_vec(7) = 0.0d0
        B_vec(8) = -gkl(j) * f1sig_out(k, f, fe) + &
                & -gkl(k) * f1sig_out(j, f, fe)

        ! Calculate steady states
        CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, ff, :))

        !--------------------------------------------!
        ! < f^{\dagger}_{j} f^{\dagger}_{k} \sigma > !
        !--------------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa_in) - i * (wl(j) + wl(k)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j)) * f1sig_out(k, ft, eg) + &
                & -CONJG(gkl(k)) * f1sig_out(j, ft, eg)
        B_vec(2) = -CONJG(gkl(j)) * f1sig_out(k, ft, ee) + &
                & -CONJG(gkl(k)) * f1sig_out(j, ft, ee)
        B_vec(3) = -CONJG(gkl(j)) * xi_in * f1sig_out(k, ft, fg) + &
                & -CONJG(gkl(k)) * xi_in * f1sig_out(j, ft, fg)
        B_vec(4) = Gamma_in * (xi_in ** 2) * f2_out(j, k, ft2) + &
                & -CONJG(gkl(j)) * xi_in * f1sig_out(k, ft, fe) + &
                & -CONJG(gkl(k)) * xi_in * f1sig_out(j, ft, fe)
        B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, ft2) + &
                & CONJG(gkl(j)) * xi_in * f1sig_out(k, ft, gg) + &
                & CONJG(gkl(j)) * xi_in * f1sig_out(k, ft, ee) + &
                & -CONJG(gkl(j)) * xi_in * f1_out(k, ft) + &
                & CONJG(gkl(k)) * xi_in * f1sig_out(j, ft, gg) + &
                & CONJG(gkl(k)) * xi_in * f1sig_out(j, ft, ee) + &
                & -CONJG(gkl(k)) * xi_in * f1_out(j, ft)
        B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, ft2)
        B_vec(7) = -CONJG(gkl(j)) * f1sig_out(k, ft, ef) + &
                & -CONJG(gkl(k)) * f1sig_out(j, ft, ef)
        B_vec(8) = 0.0d0

        ! Calculate steady states
        CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, ft2, :))

        !----------------------------------!
        ! < f^{\dagger}_{j} f_{k} \sigma > !
        !----------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa_in) - i * (wl(j) - wl(k)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j)) * f1sig_out(k, f, eg) + &
                & -gkl(k) * f1sig_out(j, ft, ge)
        B_vec(2) = -CONJG(gkl(j)) * f1sig_out(k, f, ee) + &
                & -gkl(k) * xi_in * f1sig_out(j, ft, gf)
        B_vec(3) = -CONJG(gkl(j)) * xi_in * f1sig_out(k, f, fg) + &
                & -gkl(k) * f1sig_out(j, ft, ee)
        B_vec(4) = Gamma_in * (xi_in ** 2) * f2_out(j, k, ftf) + &
                & -CONJG(gkl(j)) * xi_in * f1sig_out(k, f, fe) + &
                & -gkl(k) * xi_in * f1sig_out(j, ft, ef)
        B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, ftf) + &
                & CONJG(gkl(j)) * xi_in * f1sig_out(k, f, gg) + &
                & CONJG(gkl(j)) * xi_in * f1sig_out(k, f, ee) + &
                & -CONJG(gkl(j)) * xi_in * f1_out(k, f)
        B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f2_out(j, k, ftf) + &
                & gkl(k) * xi_in * f1sig_out(j, ft, gg) + &
                & gkl(k) * xi_in * f1sig_out(j, ft, ee) + &
                & -gkl(k) * xi_in * f1_out(j, ft)
        B_vec(7) = -CONJG(gkl(j)) * f1sig_out(k, f, ef)
        B_vec(8) = - gkl(k) * f1sig_out(j, ft, fe)

        ! Calculate steady states
        CALL MatrixInverseSS(N_mat, Mat, B_vec, f2sig_out(j, k, ftf, :))

      END IF

      ! Close j loop
    END DO
    ! Close k loop
  END DO

  photon_out = REAL(moment_out)

  IF (G2_SS .EQV. .TRUE.) THEN
    ! Cycle through modes
    DO l = -N_in, N_in
      DO k = -N_in, N_in
        DO j = -N_in, N_in
          !-----------------------------!
          !     THIRD-ORDER: CAVITY     !
          !-----------------------------!
          !---------------------------------!
          ! < f^{\dagger}_{j} f_{k} f_{l} > !
          !---------------------------------!
          f3_out(j, k, l, ftf2) = -CONJG(gkl(j)) * f2sig_out(k, l, ff, eg) + &
                                & -CONJG(gkl(j)) * xi_in * f2sig_out(k, l, ff, fe) + &
                                & -gkl(k) * f2sig_out(j, l, ftf, ge) + &
                                & -gkl(k) * xi_in * f2sig_out(j, l, ftf, ef) + &
                                & -gkl(l) * f2sig_out(j, k, ftf, ge) + &
                                & -gkl(l) * xi_in * f2sig_out(j, k, ftf, ef)
          f3_out(j, k, l, ftf2) = f3_out(j, k, l, ftf2) / &
                                & (3.0d0 * kappa_in - i * (wl(j) - wl(k) - wl(l)))

          !-------------------------------------------!
          ! < f^{\dagger}_{j} f^{\dagger}_{k} f_{l} > !
          !-------------------------------------------!
          f3_out(j, k, l, ft2f) = -CONJG(gkl(j)) * f2sig_out(k, l, ftf, eg) + &
                                & -CONJG(gkl(j)) * xi_in * f2sig_out(k, l, ftf, fe) + &
                                & -CONJG(gkl(k)) * f2sig_out(j, l, ftf, eg) + &
                                & -CONJG(gkl(k)) * xi_in * f2sig_out(j, l, ftf, fe) + &
                                & -gkl(l) * f2sig_out(j, k, ft2, ge) + &
                                & -gkl(l) * xi_in * f2sig_out(j, k, ft2, ef)
          f3_out(j, k, l, ft2f) = f3_out(j, k, l, ft2f) / &
                                & (3.0d0 * kappa_in - i * (wl(j) + wl(k) - wl(l)))

          !--------------------------------------!
          !     FOURTH-ORDER: CAVITY AND ATOM    !
          !--------------------------------------!
          !----------------------------------------!
          ! < f^{\dagger}_{j} f_{k} f_{l} \sigma > !
          !----------------------------------------!
          ! Set the diagonal matrix elements for M
          Mat = Mat_OG
          DO x = 1, N_mat
            Mat(x, x) = Mat(x, x) - ((3.0d0 * kappa_in) - i * (wl(j) - wl(k) - wl(l)))
          END DO

          ! Set the non-homogeneous vector
          B_vec = 0.0d0
          B_vec(1) = -CONJG(gkl(j)) * f2sig_out(k, l, ff, eg) + &
                  & -gkl(k) * f2sig_out(j, l, ftf, ge) + &
                  & -gkl(l) * f2sig_out(j, k, ftf, ge)
          B_vec(2) = -CONJG(gkl(j)) * f2sig_out(k, l, ff, ee) + &
                  & -gkl(k) * xi_in * f2sig_out(j, l, ftf, gf) + &
                  & -gkl(l) * xi_in * f2sig_out(j, k, ftf, gf)
          B_vec(3) = -CONJG(gkl(j)) * xi_in * f2sig_out(k, l, ff, fg) + &
                  & -gkl(k) * f2sig_out(j, l, ftf, ee) + &
                  & -gkl(l) * f2sig_out(j, k, ftf, ee)
          B_vec(4) = Gamma_in * (xi_in ** 2) * f3_out(j, k, l, ftf2) + &
                  & -CONJG(gkl(j)) * xi_in * f2sig_out(k, l, ff, fe) + &
                  & -gkl(k) * xi_in * f2sig_out(j, l, ftf, ef) + &
                  & -gkl(l) * xi_in * f2sig_out(j, k, ftf, ef)
          B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f3_out(j, k, l, ftf2) + &
                  & CONJG(gkl(j)) * xi_in * f2sig_out(k, l, ff, gg) + &
                  & CONJG(gkl(j)) * xi_in * f2sig_out(k, l, ff, ee) + &
                  & -CONJG(gkl(j)) * xi_in * f2_out(k, l, ff)
          B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f3_out(j, k, l, ftf2) + &
                  & gkl(k) * xi_in * f2sig_out(j, l, ftf, gg) + &
                  & gkl(k) * xi_in * f2sig_out(j, l, ftf, ee) + &
                  & -gkl(k) * xi_in * f2_out(j, l, ftf) + &
                  & gkl(l) * xi_in * f2sig_out(j, k, ftf, gg) + &
                  & gkl(l) * xi_in * f2sig_out(j, k, ftf, ee) + &
                  & -gkl(l) * xi_in * f2_out(j, k, ftf)
          B_vec(7) = -CONJG(gkl(j)) * f2sig_out(k, l, ff, ef)
          B_vec(8) = -gkl(k) * f2sig_out(j, l, ftf, fe) + &
                  & -gkl(l) * f2sig_out(j, k, ftf, fe)

          ! Calculate steady states
          CALL MatrixInverseSS(N_mat, Mat, B_vec, f3sig_out(j, k, l, ftf2, :))

          !--------------------------------------------------!
          ! < f^{\dagger}_{j} f^{\dagger}_{k} f_{l} \sigma > !
          !--------------------------------------------------!
          ! Set the diagonal matrix elements for M
          Mat = Mat_OG
          DO x = 1, N_mat
            Mat(x, x) = Mat(x, x) - ((3.0d0 * kappa_in) - i * (wl(j) + wl(k) - wl(l)))
          END DO

          ! Set the non-homogeneous vector
          B_vec = 0.0d0
          B_vec(1) = -CONJG(gkl(j)) * f2sig_out(k, l, ftf, eg) + &
                  & -CONJG(gkl(k)) * f2sig_out(j, l, ftf, eg) + &
                  & -gkl(l) * f2sig_out(j, k, ft2, ge)
          B_vec(2) = -CONJG(gkl(j)) * f2sig_out(k, l, ftf, ee) + &
                  & -CONJG(gkl(k)) * f2sig_out(j, l, ftf, ee) + &
                  & -gkl(l) * xi_in * f2sig_out(j, k, ft2, gf)
          B_vec(3) = -CONJG(gkl(j)) * xi_in * f2sig_out(k, l, ftf, fg) + &
                  & -CONJG(gkl(k)) * xi_in * f2sig_out(j, l, ftf, fg) + &
                  & -gkl(l) * f2sig_out(j, k, ft2, ee)
          B_vec(4) = Gamma_in * (xi_in ** 2) * f3_out(j, k, l, ft2f) + &
                  & -CONJG(gkl(j)) * xi_in * f2sig_out(k, l, ftf, fe) + &
                  & -CONJG(gkl(k)) * xi_in * f2sig_out(j, l, ftf, fe) + &
                  & -gkl(l) * xi_in * f2sig_out(j, k, ft2, ef)
          B_vec(5) = i * xi_in * 0.5d0 * Omega_in * f3_out(j, k, l, ft2f) + &
                  & CONJG(gkl(j)) * xi_in * f2sig_out(k, l, ftf, gg) + &
                  & CONJG(gkl(j)) * xi_in * f2sig_out(k, l, ftf, ee) + &
                  & -CONJG(gkl(j)) * xi_in * f2_out(k, l, ftf) + &
                  & CONJG(gkl(k)) * xi_in * f2sig_out(j, l, ftf, gg) + &
                  & CONJG(gkl(k)) * xi_in * f2sig_out(j, l, ftf, ee) + &
                  & -CONJG(gkl(k)) * xi_in * f2_out(j, l, ftf)
          B_vec(6) = -i * xi_in * 0.5d0 * Omega_in * f3_out(j, k, l, ft2f) + &
                  & gkl(l) * xi_in * f2sig_out(j, k, ft2, gg) + &
                  & gkl(l) * xi_in * f2sig_out(j, k, ft2, ee) + &
                  & -gkl(l) * xi_in * f2_out(j, k, ft2)
          B_vec(7) = -CONJG(gkl(j)) * f2sig_out(k, l, ftf, ef) + &
                  & -CONJG(gkl(k)) * f2sig_out(j, l, ftf, ef)
          B_vec(8) = -gkl(l) * f2sig_out(j, k, ft2, fe)

          ! Calculate steady states
          CALL MatrixInverseSS(N_mat, Mat, B_vec, f3sig_out(j, k, l, ft2f, :))

          ! Close j loop
        END DO
        ! Close k loop
      END DO
      ! Close l loop
    END DO
  END IF

  IF (G2_SS .EQV. .TRUE.) THEN
    ! Cycle through modes
    DO m = -N_in, N_in
      DO l = -N_in, N_in
        DO k = -N_in, N_in
          DO j = -N_in, N_in
            !------------------------------!
            !     FOURTH-ORDER: CAVITY     !
            !------------------------------!
            !-------------------------------------------------!
            ! < f^{\dagger}_{j} f^{\dagger}_{k} f_{l} f_{m} > !
            !-------------------------------------------------!
            f4_out(j, k, l, m) = -CONJG(gkl(j)) * f3sig_out(k, l, m, ftf2, eg) + &
                              & -CONJG(gkl(j)) * xi_in * f3sig_out(k, l, m, ftf2, fe) + &
                              & -CONJG(gkl(k)) * f3sig_out(j, l, m, ftf2, eg) + &
                              & -CONJG(gkl(k)) * xi_in * f3sig_out(j, l, m, ftf2, fe) + &
                              & -gkl(l) * f3sig_out(j, k, m, ft2f, ge) + &
                              & -gkl(l) * xi_in * f3sig_out(j, k, m, ft2f, ef) + &
                              & -gkl(m) * f3sig_out(j, k, l, ft2f, ge) + &
                              & -gkl(m) * xi_in * f3sig_out(j, k, l, ft2f, ef)
            f4_out(j, k, l, m) = f4_out(j, k, l, m) / &
                              & (4.0d0 * kappa_in - i * (wl(j) + wl(k)) + i * (wl(l) + wl(m)))

            ! Close j loop
          END DO
          ! Close k loop
        END DO
        ! Close l loop
      END DO
      ! Close m loop
    END DO
  END IF

END SUBROUTINE SteadyStateMoments

!==============================================================================!
!                          G1 CORRELATION SUBROUTINES                          !
!==============================================================================!

! Subroutine to calculate the initial conditions for the auto-correlations
SUBROUTINE G1_InitialConditions(Gamma_in, Omega_in, alpha_in, delta_in, xi_in, &
                              & epsilon_in, N_in, phase_in, &
                              & w0_in, kappa_in, dw_in, &
                              & photon_ss_out, B_OG_out, &
                              & sigma_out, f1_out)

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
  REAL(KIND=8), INTENT(IN)                   :: w0_in
  ! Cavity linewidth/transmission of cavity mode
  REAL(KIND=8), INTENT(IN)                   :: kappa_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                   :: dw_in

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
  ! Integer indices for: a, f^{\dagger}, f^{\dagger} a
  INTEGER, PARAMETER                         :: f = 1, ft = 2
  INTEGER, PARAMETER                         :: ff = 1, ftf = 2, ft2 = 3
  INTEGER, PARAMETER                         :: ftf2 = 1, ft2f = 2

  ! Steady state arrays
  ! First-order moments: Atomic equations (< \sigma >)
  COMPLEX(KIND=8), DIMENSION(N_mat)                                  :: sigma_ss
  ! First-order moments: Cavity (< a >, < f^{\dagger} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2)                          :: f1_ss
  ! Second-order moments: Cavity and atom (< a \sigma >, < f^{\dagger} \sigma >
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N_mat)                   :: f1sig_ss
  ! Second-order moments: Cavity (< f^{\dagger} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3)              :: f2_ss
  ! Third-order moments: Cavity and atom (< a^{2} \sigma >, < a^{\dagger 2} \sigma >, < f^{\dagger} a \sigma >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3, N_mat)       :: f2sig_ss
  ! Third-order moments: Cavity (< a^{2} f^{\dagger} >, < a^{\dagger 2} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, 2)  :: f3_ss
  ! Fourth-order moments: Cavity and atom ( < f^{\dagger} a^{2} \sigma >, < a^{\dagger 2} a \sigma >)
  COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE             :: f3sig_ss
  ! Fourth-order moments: Cavity (< a^{\dagger 2} a^{2} >)
  COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE                :: f4_ss

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Steady state photon number
  REAL(KIND=8), INTENT(OUT)                                          :: photon_ss_out
  ! Non-homogeneous vector
  COMPLEX(KIND=8), DIMENSION(N_mat), INTENT(OUT)                     :: B_OG_out

  ! Time integration arrays
  ! First-order moments: Atomic equations (< \sigma >)
  COMPLEX(KIND=8), DIMENSION(N_mat), INTENT(OUT)                     :: sigma_out
  ! First-order moments: Cavity (< f >, < f^{\dagger} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2), INTENT(OUT)             :: f1_out

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
  ALLOCATE(f3sig_ss(-N_in:N_in, -N_in:N_in, -N_in:N_in, 2, N_mat)); f3sig_ss = 0.0d0
  ! Fourth-order: Cavity
  ALLOCATE(f4_ss(-N_in:N_in, -N_in:N_in, -N_in:N_in, -N_in:N_in)); f4_ss = 0.0d0

  !==============================================================================!
  !                        CALCULATE STEADY-STATE MOMENTS                        !
  !==============================================================================!
  CALL SteadyStateMoments(Gamma_in, Omega_in, alpha_in, delta_in, xi_in, &
                        & epsilon_in, N_in, phase_in, &
                        & w0_in, kappa_in, dw_in, &
                        & photon_ss_out, sigma_ss, .FALSE., &
                        & f1_ss, f1sig_ss, &
                        & f2_ss, f2sig_ss, &
                        & f3_ss, f3sig_ss, &
                        & f4_ss)

  !==============================================================================!
  !         CALCULATE FIRST-ORDER CORRELATION FUNCTION INITIAL CONDITIONS        !
  !==============================================================================!
  ! Set initial conditions and non-homogeneous vector
  ! < \sigma (\tau = 0) a_{k} (0) > = < a_{k} \sigma >_{ss},
  ! and
  ! < a^{\dagger}_{j} (\tau = 0) a_{k}(0) > = < a^{\dagger}_{j} a_{k} >_{ss}.

  sigma_out = 0.0d0
  f1_out = 0.0d0
  B_OG_out = 0.0d0
  DO k = -N_in, N_in
    ! Initial conditions
    sigma_out(:) = sigma_out(:) + f1sig_ss(k, f, :)
    DO j = -N_in, N_in
      f1_out(j, ft) = f1_out(j, ft) + f2_ss(j, k, ftf)
    END DO
    ! non homogeneous vector
    B_OG_out(4) = B_OG_out(4) + Gamma_in * (xi_in ** 2) * f1_ss(k, f)
    B_OG_out(5) = B_OG_out(5) + i * xi_in * 0.5d0 * Omega_in * f1_ss(k, f)
    B_OG_out(6) = B_OG_out(6) - i * xi_in * 0.5d0 * Omega_in * f1_ss(k, f)
  END DO

END SUBROUTINE G1_InitialConditions

! Subroutine to calculate the time evolution of the g2 correlation
SUBROUTINE G1_CalculateRK4(Gamma_in, Omega_in, alpha_in, delta_in, xi_in, &
                         & epsilon_in, N_in, phase_in, &
                         & w0_in, kappa_in, dw_in, &
                         & dt_in, tau_steps_in, &
                         & g1_array_out, WRITE_DATA_IN, filename_data_in)

  !==============================================================================!
  !                    DEFINING AND DECLARING VARIABLES/ARRAYS                   !
  !==============================================================================!

  IMPLICIT NONE

  !---------------!
  !     INPUT     !
  !---------------!
  ! Atomic decay rate
  REAL(KIND=8), INTENT(IN)                  :: Gamma_in
  ! Driving amplitude
  REAL(KIND=8), INTENT(IN)                  :: Omega_in
  ! Atomic anharmonicity
  REAL(KIND=8), INTENT(IN)                  :: alpha_in
  ! Drive detuning from two-photon resonance
  REAL(KIND=8), INTENT(IN)                  :: delta_in
  ! Dipole moment ratio
  REAL(KIND=8), INTENT(IN)                  :: xi_in

  ! Filter parameter stuff
  ! Percentage of fluorecence aimed at cavity
  REAL(KIND=8), INTENT(IN)                  :: epsilon_in
  ! Number of mode either side of w0, 2N + 1 total mode
  INTEGER, INTENT(IN)                       :: N_in
  ! Phase modulation of mode coupling
  REAL(KIND=8), INTENT(IN)                  :: phase_in
  ! Central mode frequency of the filter cavity, with N mode frequencies either side
  REAL(KIND=8), INTENT(IN)                  :: w0_in
  ! Cavity linewidth/transmission of cavity mode
  REAL(KIND=8), INTENT(IN)                  :: kappa_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                  :: dw_in

  ! Time stuff
  ! Time step
  REAL(KIND=8), INTENT(IN)                  :: dt_in
  ! Maxi_inmum number of steps to integrate for
  INTEGER, INTENT(IN)                       :: tau_steps_in

  ! Data stuff
  ! Boolean for writing data
  LOGICAL, INTENT(IN)                       :: WRITE_DATA_IN
  ! Filename for writing data to
  CHARACTER(LEN=*), INTENT(IN)              :: filename_data_in

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Data array
  COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: g1_array_out

  !------------------------------------!
  !     MOMENT EQUATION ARRAY STUFF    !
  !------------------------------------!
  ! Dimension of M matrix
  INTEGER, PARAMETER                        :: N_mat = 8
  ! M matrix (filled as transpose)
  COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)  :: Mat, Mat_OG
  ! Non-homogeneous vector
  COMPLEX(KIND=8), DIMENSION(N_mat)         :: B_OG, B_vec

  ! Integer indices for sigma operators
  INTEGER, PARAMETER                        :: gg = 1, ge = 2, eg = 3
  INTEGER, PARAMETER                        :: ee = 4, ef = 5, fe = 6
  INTEGER, PARAMETER                        :: gf = 7, fg = 8
  ! Integer indices for: a, f^{\dagger}, f^{\dagger} a
  INTEGER, PARAMETER                        :: f = 1, ft = 2
  INTEGER, PARAMETER                        :: ff = 1, ftf = 2, ft2 = 3
  INTEGER, PARAMETER                        :: ftf2 = 1, ft2f = 2

  ! Time integration arrays
  ! First-order moments: Atomic equations (< \sigma >)
  COMPLEX(KIND=8), DIMENSION(N_mat)         :: sigma
  COMPLEX(KIND=8), DIMENSION(N_mat)         :: k1_sigma, k2_sigma, k3_sigma, k4_sigma
  ! First-order moments: Cavity (< a >, < f^{\dagger} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2) :: f1
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2) :: k1_f1, k2_f1, k3_f1, k4_f1

  !----------------------------!
  !     OTHER USEFUL STUFF     !
  !----------------------------!
  ! Time step integer
  INTEGER                                   :: t
  ! Integer counter
  INTEGER                                   :: j, k, l, m, x
  ! Sample rate for state populations
  INTEGER                                   :: sample_rate
  ! Imaginary i
  COMPLEX(KIND=8), PARAMETER                :: i = CMPLX(0.0d0, 1.0d0, 8)
  ! pi
  REAL(KIND=8), PARAMETER                   :: pi = 3.1415926535897932384d0
  ! 1 / 6
  REAL(KIND=8), PARAMETER                   :: xis = 1.0d0 / 6.0d0
  ! List of Delta values
  REAL(KIND=8), DIMENSION(-N_in:N_in)       :: wl
  ! List of mode dependent cascade coupling values
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in)    :: gkl
  ! Blackman window coefficient
  REAL(KIND=8)                              :: blackman
  ! Steady state photon number
  REAL(KIND=8)                              :: photon_ss
  ! Complex data
  COMPLEX(KIND=8)                           :: moment_out

  !==============================================================================!
  !                DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS                !
  !==============================================================================!
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
      wl(j) = w0_in
      gkl(j) = DSQRT(epsilon_in * Gamma_in * kappa_in)
    ELSE
      wl(j) = w0_in + DBLE(j) * dw_in
      ! Blackman window coefficient
      blackman = 1.0d0
      ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
      !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))
      ! Mode dependent phase difference
      gkl(j) = DSQRT((epsilon_in / DBLE(2*N_in + 1)) * Gamma_in * kappa_in) * blackman * &
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

  ! Data
  ALLOCATE(g1_array_out(0:tau_steps_in)); g1_array_out = 0.0d0
  photon_ss = 0.0d0

  !==============================================================================!
  !                         CALCULATE INITIAL CONDITIONS                         !
  !==============================================================================!
  CALL G1_InitialConditions(Gamma_in, Omega_in, alpha_in, delta_in, xi_in, &
                          & epsilon_in, N_in, phase_in, &
                          & w0_in, kappa_in, dw_in, &
                          & photon_ss, B_OG, &
                          & sigma, f1)

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
    DO j = -N_in, N_in
      moment_out = moment_out + f1(j, ft)
    END DO

    ! Normalise correlation by steady-state photon number
    IF (photon_ss .NE. 0.0) THEN
      moment_out = moment_out / (photon_ss)
    END IF

    !-----------------------!
    !     WRITE TO FILE     !
    !-----------------------!
    g1_array_out(t) = MOMENT_OUT
    ! Second-order correlation function
    ! If WRITE_DATA_IN is TRUE, write data to file
    IF (WRITE_DATA_IN .EQV. .TRUE.) THEN
      ! If t_max is really big, only take a sample of results to write to file
      ! so file size isn't huge-mongous.
      IF (MOD(t, sample_rate) == 0) THEN
        WRITE(4, *) dt_in * DBLE(t), REAL(moment_out), AIMAG(moment_out)
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
      !---------------------!
      ! < f^{\dagger}_{j} > !
      !---------------------!
      k1_f1(j, ft) = -dt_in * (kappa_in - i * wl(j)) * f1(j, ft) + &
                   & -dt_in * CONJG(gkl(j)) * sigma(eg) + &
                   & -dt_in * CONJG(gkl(j)) * xi_in * sigma(fe)

      k2_f1(j, ft) = -dt_in * (kappa_in - i * wl(j)) * (f1(j, ft) + 0.5d0 * k1_f1(j, ft)) + &
                   & -dt_in * CONJG(gkl(j)) * (sigma(eg) + 0.5d0 * k1_sigma(eg)) + &
                   & -dt_in * CONJG(gkl(j)) * xi_in * (sigma(fe) + 0.5d0 * k1_sigma(fe))

      k3_f1(j, ft) = -dt_in * (kappa_in - i * wl(j)) * (f1(j, ft) + 0.5d0 * k2_f1(j, ft)) + &
                   & -dt_in * CONJG(gkl(j)) * (sigma(eg) + 0.5d0 * k2_sigma(eg)) + &
                   & -dt_in * CONJG(gkl(j)) * xi_in * (sigma(fe) + 0.5d0 * k2_sigma(fe))

      k4_f1(j, ft) = -dt_in * (kappa_in - i * wl(j)) * (f1(j, ft) + k3_f1(j, ft)) + &
                   & -dt_in * CONJG(gkl(j)) * (sigma(eg) + k3_sigma(eg)) + &
                   & -dt_in * CONJG(gkl(j)) * xi_in * (sigma(fe) + k3_sigma(fe))

      ! Close j loop
    END DO


    !============================================================================!
    !                   UPDATE ARRAYS FROM RUNGE-KUTTA ARRAYS                    !
    !============================================================================!
    ! First-order
    sigma = sigma + xis * (k1_sigma + 2.0d0 * (k2_sigma + k3_sigma) + k4_sigma)
    f1 = f1 + xis * (k1_f1 + 2.0d0 * (k2_f1 + k3_f1) + k4_f1)

    ! Close t loop
  END DO

  ! If WRITE_DATA_IN is TRUE, close the file
  IF (WRITE_DATA_IN .EQV. .TRUE.) THEN
    ! Close file
    CLOSE(4)
  END IF

END SUBROUTINE G1_CalculateRK4

!==============================================================================!
!                          G2 CORRELATION SUBROUTINES                          !
!==============================================================================!

! Subroutine to calculate the initial conditions for the auto-correlations
SUBROUTINE G2_InitialConditions(Gamma_in, Omega_in, alpha_in, delta_in, xi_in, &
                              & epsilon_in, N_in, phase_in, &
                              & w0_in, kappa_in, dw_in, &
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
  REAL(KIND=8), INTENT(IN)                   :: w0_in
  ! Cavity linewidth/transmission of cavity mode
  REAL(KIND=8), INTENT(IN)                   :: kappa_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                   :: dw_in

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
  ! Integer indices for: a, f^{\dagger}, f^{\dagger} a
  INTEGER, PARAMETER                         :: f = 1, ft = 2
  INTEGER, PARAMETER                         :: ff = 1, ftf = 2, ft2 = 3
  INTEGER, PARAMETER                         :: ftf2 = 1, ft2f = 2

  ! Steady state arrays
  ! First-order moments: Atomic equations (< \sigma >)
  COMPLEX(KIND=8), DIMENSION(N_mat)                                  :: sigma_ss
  ! First-order moments: Cavity (< a >, < f^{\dagger} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2)                          :: f1_ss
  ! Second-order moments: Cavity and atom (< a \sigma >, < f^{\dagger} \sigma >
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N_mat)                   :: f1sig_ss
  ! Second-order moments: Cavity (< f^{\dagger} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3)              :: f2_ss
  ! Third-order moments: Cavity and atom (< a^{2} \sigma >, < a^{\dagger 2} \sigma >, < f^{\dagger} a \sigma >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3, N_mat)       :: f2sig_ss
  ! Third-order moments: Cavity (< a^{2} f^{\dagger} >, < a^{\dagger 2} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, 2)  :: f3_ss
  ! Fourth-order moments: Cavity and atom ( < f^{\dagger} a^{2} \sigma >, < a^{\dagger 2} a \sigma >)
  COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE             :: f3sig_ss
  ! Fourth-order moments: Cavity (< a^{\dagger 2} a^{2} >)
  COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE                :: f4_ss

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Steady state photon number
  REAL(KIND=8), INTENT(OUT)                                          :: photon_ss_out
  ! Non-homogeneous vector
  COMPLEX(KIND=8), DIMENSION(N_mat), INTENT(OUT)                     :: B_OG_out

  ! Time integration arrays
  ! First-order moments: Atomic equations (< \sigma >)
  COMPLEX(KIND=8), DIMENSION(N_mat), INTENT(OUT)                     :: sigma_out
  ! First-order moments: Cavity (< a >, < f^{\dagger} >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2), INTENT(OUT)             :: f1_out
  ! Second-order moments: Cavity and atom (< a \sigma >, < f^{\dagger} \sigma >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N_mat), INTENT(OUT)      :: f1sig_out
  ! Second-order moments: Cavity (< f^{\dagger} a >)
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in), INTENT(OUT)    :: f2_out

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
  ALLOCATE(f3sig_ss(-N_in:N_in, -N_in:N_in, -N_in:N_in, 2, N_mat)); f3sig_ss = 0.0d0
  ! Fourth-order: Cavity
  ALLOCATE(f4_ss(-N_in:N_in, -N_in:N_in, -N_in:N_in, -N_in:N_in)); f4_ss = 0.0d0

  !==============================================================================!
  !                        CALCULATE STEADY-STATE MOMENTS                        !
  !==============================================================================!
  CALL SteadyStateMoments(Gamma_in, Omega_in, alpha_in, delta_in, xi_in, &
                        & epsilon_in, N_in, phase_in, &
                        & w0_in, kappa_in, dw_in, &
                        & photon_ss_out, sigma_ss, .TRUE., &
                        & f1_ss, f1sig_ss, &
                        & f2_ss, f2sig_ss, &
                        & f3_ss, f3sig_ss, &
                        & f4_ss)

  !==============================================================================!
  !                  CALCULATE FIRST-ORDER CORRELATION FUNCTION                  !
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
        f1_out(k, f) = f1_out(k, f) + f3_ss(j, k, m, ftf2)
        f1_out(k, ft) = f1_out(k, ft) + f3_ss(j, k, m, ft2f)
        ! Second-order: Cavity and atom
        f1sig_out(k, f, :) = f1sig_out(k, f, :) + f3sig_ss(j, k, m, ftf2, :)
        f1sig_out(k, ft, :) = f1sig_out(k, ft, :) + f3sig_ss(j, k, m, ft2f, :)

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

  !--------------------------------!
  !     NON-HOMOGENEOUS VECTOR     !
  !--------------------------------!
  ! B_OG_out = 0.0d0
  ! B_OG_out(4) = Gamma_in * (xi_in ** 2) * photon_ss_out
  ! B_OG_out(5) = i * xi_in * 0.5d0 * Omega_in * photon_ss_out
  ! B_OG_out(6) = -i * xi_in * 0.5d0 * Omega_in * photon_ss_out

END SUBROUTINE G2_InitialConditions

! Subroutine to calculate the time evolution of the g2 correlation
SUBROUTINE G2_CalculateRK4(Gamma_in, Omega_in, alpha_in, delta_in, xi_in, &
                         & epsilon_in, N_in, phase_in, &
                         & w0_in, kappa_in, dw_in, &
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
  REAL(KIND=8), INTENT(IN)                 :: Gamma_in
  ! Driving amplitude
  REAL(KIND=8), INTENT(IN)                 :: Omega_in
  ! Atomic anharmonicity
  REAL(KIND=8), INTENT(IN)                 :: alpha_in
  ! Drive detuning from two-photon resonance
  REAL(KIND=8), INTENT(IN)                 :: delta_in
  ! Dipole moment ratio
  REAL(KIND=8), INTENT(IN)                 :: xi_in

  ! Filter parameter stuff
  ! Percentage of fluorecence aimed at cavity
  REAL(KIND=8), INTENT(IN)                 :: epsilon_in
  ! Number of mode either side of w0, 2N + 1 total mode
  INTEGER, INTENT(IN)                      :: N_in
  ! Phase modulation of mode coupling
  REAL(KIND=8), INTENT(IN)                 :: phase_in
  ! Central mode frequency of the filter cavity, with N mode frequencies either side
  REAL(KIND=8), INTENT(IN)                 :: w0_in
  ! Cavity linewidth/transmission of cavity mode
  REAL(KIND=8), INTENT(IN)                 :: kappa_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                 :: dw_in

  ! Time stuff
  ! Time step
  REAL(KIND=8), INTENT(IN)                 :: dt_in
  ! Maxi_inmum number of steps to integrate for
  INTEGER, INTENT(IN)                      :: tau_steps_in

  ! Data stuff
  ! Boolean for writing data
  LOGICAL, INTENT(IN)                      :: WRITE_DATA_IN
  ! Filename for writing data to
  CHARACTER(LEN=*), INTENT(IN)             :: filename_data_in

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Data array
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: g2_array_out

  !------------------------------------!
  !     MOMENT EQUATION ARRAY STUFF    !
  !------------------------------------!
  ! Dimension of M matrix
  INTEGER, PARAMETER                       :: N_mat = 8
  ! M matrix (filled as transpose)
  COMPLEX(KIND=8), DIMENSION(N_mat, N_mat) :: Mat, Mat_OG
  ! Non-homogeneous vector
  COMPLEX(KIND=8), DIMENSION(N_mat)        :: B_OG, B_vec

  ! Integer indices for sigma operators
  INTEGER, PARAMETER                       :: gg = 1, ge = 2, eg = 3
  INTEGER, PARAMETER                       :: ee = 4, ef = 5, fe = 6
  INTEGER, PARAMETER                       :: gf = 7, fg = 8
  ! Integer indices for: a, f^{\dagger}, f^{\dagger} a
  INTEGER, PARAMETER                       :: f = 1, ft = 2
  INTEGER, PARAMETER                       :: ff = 1, ftf = 2, ft2 = 3
  INTEGER, PARAMETER                       :: ftf2 = 1, ft2f = 2

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
  INTEGER                                  :: t
  ! Integer counter
  INTEGER                                  :: j, k, l, m, x
  ! Sample rate for state populations
  INTEGER                                  :: sample_rate
  ! Imaginary i
  COMPLEX(KIND=8), PARAMETER               :: i = CMPLX(0.0d0, 1.0d0, 8)
  ! pi
  REAL(KIND=8), PARAMETER                  :: pi = 3.1415926535897932384d0
  ! 1 / 6
  REAL(KIND=8), PARAMETER                  :: xis = 1.0d0 / 6.0d0
  ! List of Delta values
  REAL(KIND=8), DIMENSION(-N_in:N_in)      :: wl
  ! List of mode dependent cascade coupling values
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in)   :: gkl
  ! Blackman window coefficient
  REAL(KIND=8)                             :: blackman
  ! Steady state photon number
  REAL(KIND=8)                             :: photon_ss
  ! Complex data
  COMPLEX(KIND=8)                          :: moment_out

  !==============================================================================!
  !                DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS                !
  !==============================================================================!
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
      wl(j) = w0_in
      gkl(j) = DSQRT(epsilon_in * Gamma_in * kappa_in)
    ELSE
      wl(j) = w0_in + DBLE(j) * dw_in
      ! Blackman window coefficient
      blackman = 1.0d0
      ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
      !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))
      ! Mode dependent phase difference
      gkl(j) = DSQRT((epsilon_in / DBLE(2*N_in + 1)) * Gamma_in * kappa_in) * blackman * &
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
                          & w0_in, kappa_in, dw_in, &
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
    IF (photon_ss .NE. 0.0) THEN
      moment_out = moment_out / (photon_ss ** 2)
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
      k1_f1(j, f) = -dt_in * (kappa_in + i * wl(j)) * f1(j, f) + &
                  & -dt_in * gkl(j) * sigma(ge) + &
                  & -dt_in * gkl(j) * xi_in * sigma(ef)

      k2_f1(j, f) = -dt_in * (kappa_in + i * wl(j)) * (f1(j, f) + 0.5d0 * k1_f1(j, f)) + &
                  & -dt_in * gkl(j) * (sigma(ge) + 0.5d0 * k1_sigma(ge)) + &
                  & -dt_in * gkl(j) * xi_in * (sigma(ef) + 0.5d0 * k1_sigma(ef))

      k3_f1(j, f) = -dt_in * (kappa_in + i * wl(j)) * (f1(j, f) + 0.5d0 * k2_f1(j, f)) + &
                  & -dt_in * gkl(j) * (sigma(ge) + 0.5d0 * k2_sigma(ge)) + &
                  & -dt_in * gkl(j) * xi_in * (sigma(ef) + 0.5d0 * k2_sigma(ef))

      k4_f1(j, f) = -dt_in * (kappa_in + i * wl(j)) * (f1(j, f) + k3_f1(j, f)) + &
                  & -dt_in * gkl(j) * (sigma(ge) + k3_sigma(ge)) + &
                  & -dt_in * gkl(j) * xi_in * (sigma(ef) + k3_sigma(ef))

      !---------------------!
      ! < f^{\dagger}_{j} > !
      !---------------------!
      k1_f1(j, ft) = -dt_in * (kappa_in - i * wl(j)) * f1(j, ft) + &
                   & -dt_in * CONJG(gkl(j)) * sigma(eg) + &
                   & -dt_in * CONJG(gkl(j)) * xi_in * sigma(fe)

      k2_f1(j, ft) = -dt_in * (kappa_in - i * wl(j)) * (f1(j, ft) + 0.5d0 * k1_f1(j, ft)) + &
                   & -dt_in * CONJG(gkl(j)) * (sigma(eg) + 0.5d0 * k1_sigma(eg)) + &
                   & -dt_in * CONJG(gkl(j)) * xi_in * (sigma(fe) + 0.5d0 * k1_sigma(fe))

      k3_f1(j, ft) = -dt_in * (kappa_in - i * wl(j)) * (f1(j, ft) + 0.5d0 * k2_f1(j, ft)) + &
                   & -dt_in * CONJG(gkl(j)) * (sigma(eg) + 0.5d0 * k2_sigma(eg)) + &
                   & -dt_in * CONJG(gkl(j)) * xi_in * (sigma(fe) + 0.5d0 * k2_sigma(fe))

      k4_f1(j, ft) = -dt_in * (kappa_in - i * wl(j)) * (f1(j, ft) + k3_f1(j, ft)) + &
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
        Mat(x, x) = Mat(x, x) - (kappa_in + i * wl(j))
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
               & -gkl(j) * xi_in * photon_ss
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
               & -gkl(j) * xi_in * photon_ss
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
               & -gkl(j) * xi_in * photon_ss
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
               & -gkl(j) * xi_in * photon_ss
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
        Mat(x, x) = Mat(x, x) - (kappa_in - i * wl(j))
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
               & -CONJG(gkl(j)) * xi_in * photon_ss
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
               & -CONJG(gkl(j)) * xi_in * photon_ss
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
               & -CONJG(gkl(j)) * xi_in * photon_ss
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
               & -CONJG(gkl(j)) * xi_in * photon_ss
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
        k1_f2(j, k) = -dt_in * (2.0d0 * kappa_in - i * (wl(j) - wl(k))) * f2(j, k) + &
                    & -dt_in * CONJG(gkl(j)) * f1sig(k, f, eg) + &
                    & -dt_in * CONJG(gkl(j)) * xi_in * f1sig(k, f, fe) + &
                    & -dt_in * gkl(k) * f1sig(j, ft, ge) + &
                    & -dt_in * gkl(k) * xi_in * f1sig(j, ft, ef)

        k2_f2(j, k) = -dt_in * (2.0d0 * kappa_in - i * (wl(j) - wl(k))) * (f2(j, k) + 0.5d0 * k1_f2(j, k)) + &
                    & -dt_in * CONJG(gkl(j)) * (f1sig(k, f, eg) + 0.5d0 * k1_f1sig(k, f, eg)) + &
                    & -dt_in * CONJG(gkl(j)) * xi_in * (f1sig(k, f, fe) + 0.5d0 * k1_f1sig(k, f, fe)) + &
                    & -dt_in * gkl(k) * (f1sig(j, ft, ge) + 0.5d0 * k1_f1sig(j, ft, ge)) + &
                    & -dt_in * gkl(k) * xi_in * (f1sig(j, ft, ef) + 0.5d0 * k1_f1sig(j, ft, ef))

        k3_f2(j, k) = -dt_in * (2.0d0 * kappa_in - i * (wl(j) - wl(k))) * (f2(j, k) + 0.5d0 * k2_f2(j, k)) + &
                    & -dt_in * CONJG(gkl(j)) * (f1sig(k, f, eg) + 0.5d0 * k2_f1sig(k, f, eg)) + &
                    & -dt_in * CONJG(gkl(j)) * xi_in * (f1sig(k, f, fe) + 0.5d0 * k2_f1sig(k, f, fe)) + &
                    & -dt_in * gkl(k) * (f1sig(j, ft, ge) + 0.5d0 * k2_f1sig(j, ft, ge)) + &
                    & -dt_in * gkl(k) * xi_in * (f1sig(j, ft, ef) + 0.5d0 * k2_f1sig(j, ft, ef))

        k4_f2(j, k) = -dt_in * (2.0d0 * kappa_in - i * (wl(j) - wl(k))) * (f2(j, k) + k3_f2(j, k)) + &
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

END MODULE SINGLE_FILTER_SUBROUTINES