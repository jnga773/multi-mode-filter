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
SUBROUTINE SteadyStateMoments(kappa_p_in, Delta_in, lambda_in, &
                            & epsilon_in, N_in, phase_in, &
                            & w0_in, kappa_f_in, dw_in, &
                            & photon_f_out, &
                            & p2_out, p4_out, &
                            & f1p1_out, f1p3_out, &
                            & f2_out, f2p2_out, &
                            & f3p1_out, &
                            & f4_out)

  !==============================================================================!
  !                    DEFINING AND DECLARING VARIABLES/ARRAYS                   !
  !==============================================================================!

  IMPLICIT NONE

  !---------------!
  !     INPUT     !
  !---------------!
  ! Parametric cavity decay rate
  REAL(KIND=8)                              :: kappa_p_in
  ! Drive detuning from cavity resonance
  REAL(KIND=8)                              :: Delta_in
  ! Driving amplitude
  REAL(KIND=8)                              :: lambda_in

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
  REAL(KIND=8), INTENT(IN)                   :: kappa_f_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                   :: dw_in

  ! ! Calculate steady states for G1 or G2
  ! LOGICAL, INTENT(IN)                        :: G2_SS

  !------------------------------------!
  !     MOMENT EQUATION ARRAY STUFF    !
  !------------------------------------!
  ! Parametric Oscillator: Evolution matrices
  ! First order matrix
  INTEGER, PARAMETER                         :: N1 = 2
  COMPLEX(KIND=8), DIMENSION(N1, N1)         :: Mat1, Mat1_OG
  COMPLEX(KIND=8), DIMENSION(N1)             :: B1, B1_OG
  ! Second-order matrix
  INTEGER, PARAMETER                         :: N2 = 3
  COMPLEX(KIND=8), DIMENSION(N2, N2)         :: Mat2, Mat2_OG
  COMPLEX(KIND=8), DIMENSION(N2)             :: B2, B2_OG
  ! Third-order matrix
  INTEGER, PARAMETER                         :: N3 = 4
  COMPLEX(KIND=8), DIMENSION(N3, N3)         :: Mat3, Mat3_OG
  COMPLEX(KIND=8), DIMENSION(N3)             :: B3, B3_OG
  ! Second-order matrix
  INTEGER, PARAMETER                         :: N4 = 5
  COMPLEX(KIND=8), DIMENSION(N4, N4)         :: Mat4_OG
  COMPLEX(KIND=8), DIMENSION(N4)             :: B4, B4_OG
  

  ! Integer indices for sigma operators! Integer indices for sigma operators
  INTEGER, PARAMETER                         :: a = 1, at = 2
  INTEGER, PARAMETER                         :: a2 = 1, ata = 2, at2 = 3
  INTEGER, PARAMETER                         :: a3 = 1, ata2 = 2, at2a = 3, at3 = 4
  INTEGER, PARAMETER                         :: a4 = 1, ata3 = 2, at2a2 = 3, at3a = 4, at4 = 5
  ! Integer indices for: a, f^{\dagger}, f^{\dagger} a
  INTEGER, PARAMETER                         :: f = 1, ft = 2
  INTEGER, PARAMETER                         :: ff = 1, ftf = 2, ft2 = 3
  INTEGER, PARAMETER                         :: ftf2 = 1, ft2f = 2

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Steady state photon number
  REAL(KIND=8), INTENT(OUT)                                                               :: photon_f_out

  ! Steady state arrays
  ! Parametric oscillator moments
  ! Second-order: < a^{2} >, < a^{\dagger} a >, < a^{\dagger}^{2} >
  COMPLEX(KIND=8), DIMENSION(N2), INTENT(OUT)                                             :: p2_out
  ! Fourth-order: < a^{4} >, < a^{\dagger} a^{3} >, < a^{\dagger}^{2} a^{2} >,
  !               < a^{\dagger}^{3} a >< a^{\dagger}^{4} >
  COMPLEX(KIND=8), DIMENSION(N4), INTENT(OUT)                                             :: p4_out

  ! First-order: Filter / First-order: Parametric
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N1), INTENT(OUT)                              :: f1p1_out
  ! First-order: Filter / Third-order: Parametric
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N3), INTENT(OUT)                              :: f1p3_out

  ! Second-order: Filter operators
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3), INTENT(OUT)                      :: f2_out
  ! Second-order: Filter / Second-order: Parametric
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3, N2), INTENT(OUT)                  :: f2p2_out

  ! Third-order: Filter / First-order: Parametric
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, 2, N1), INTENT(OUT)      :: f3p1_out

  ! Fourth-order: Filter
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
  REAL(KIND=8), DIMENSION(-N_in:N_in)        :: wl
  ! List of mode dependent cascade coupling values
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in)     :: gkl
  ! Blackman window coefficient
  REAL(KIND=8)                               :: blackman
  ! Parametric: Photon number
  REAL(KIND=8)                               :: photon_p_ss

  !==============================================================================!
  !                DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS                !
  !==============================================================================!
  !---------------------------------------------!
  !     RESONANCES (wj) AND COUPLINGS (E_j)     !
  !---------------------------------------------!
  ! Allocate array of Delta and gka values
  wl = 0.0d0
  gkl = 0.0d0
  DO j = -N_in, N_in
    IF (N_in == 0) THEN
      wl(j) = w0_in
      gkl(j) = DSQRT(epsilon_in * 2.0d0 * kappa_p_in * kappa_f_in)
    ELSE
      wl(j) = w0_in + DBLE(j) * dw_in
      ! Blackman window coefficient
      blackman = 1.0d0
      ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
      !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))
      ! Mode dependent phase difference
      gkl(j) = DSQRT((epsilon_in / DBLE(2*N_in + 1)) * 2.0d0 * kappa_p_in * kappa_f_in) * blackman * &
             & EXP(i * DBLE(phase_in) * DBLE(j) * pi / DBLE(N_in))
    END IF
  END DO

  !------------------------------------------!
  !     INITALISE OPERATOR MOMENT ARRAYS     !
  !------------------------------------------!
  ! Steady states
  ! First-order: Filter / First-order: Parametric
  f1p1_out = 0.0d0
  ! First-order: Filter / Third-order: Parametric
  f1p3_out = 0.0d0
  ! Second-order: Filter
  f2_out = 0.0d0
  ! Second-order: Filter / Second-order Parametric
  f2p2_out = 0.0d0
  ! Third-order: Filter / First-order Parametric
  f3p1_out = 0.0d0
  ! Fourth-order: Filter
  f4_out = 0.0d0

  photon_f_out = 0.0d0

  !----------------------------!
  !     EVOLUTION MATRICES     !
  !----------------------------!
  !--------------------!
  !     First-order    !
  !--------------------!
  ! Evolution matrix
  Mat1_OG = 0.0d0
  ! d/dt < a >
  Mat1_OG(1, 1) = -(kappa_p_in + i * Delta_in)
  Mat1_OG(1, 2) = lambda_in
  ! d/dt < a^{\dagger} >
  Mat1_OG(2, 1) = lambda_in
  Mat1_OG(2, 2) = -(kappa_p_in - i * Delta_in)

  ! Non-homogeneous vector
  B1_OG = 0.0d0

  !----------------------!
  !     Second-Order     !
  !----------------------!
  ! Evolution matrix
  Mat2_OG = 0.0d0; Mat2 = 0.0d0
  ! d/dt < a^{2} >
  Mat2_OG(1, 1) = -2.0d0 * (kappa_p_in + i * Delta_in)
  Mat2_OG(1, 2) = 2.0d0 * lambda_in
  Mat2_OG(1, 3) = 0.0d0
  ! d/dt < a^{\dagger} a>
  Mat2_OG(2, 1) = lambda_in
  Mat2_OG(2, 2) = -2.0d0 * kappa_p_in
  Mat2_OG(2, 3) = lambda_in
  ! d/dt < a^{\dagger}^{2} >
  Mat2_OG(3, 1) = 0.0d0
  Mat2_OG(3, 2) = 2.0d0 * lambda_in
  Mat2_OG(3, 3) = -2.0d0 * (kappa_p_in - i * Delta_in)

  ! Mat2_OG = Mat2

  ! Non-homogeneous vector
  B2_OG = 0.0d0; B2 = 0.0d0
  B2_OG(1) = lambda_in
  B2_OG(2) = 0.0d0
  B2_OG(3) = lambda_in

  !---------------------!
  !     Third-Order     !
  !---------------------!
  ! Evolution matrix
  Mat3_OG = 0.0d0; Mat3 = 0.0d0
  ! d/dt < a^{3} >
  Mat3_OG(1, 1) = -3.0d0 * (kappa_p_in + i * Delta_in)
  Mat3_OG(1, 2) = 3.0d0 * lambda_in
  ! d/dt < a^{\dagger} a^{2} >
  Mat3_OG(2, 1) = lambda_in
  Mat3_OG(2, 2) = -(3.0d0 * kappa_p_in + i * Delta_in)
  Mat3_OG(2, 3) = 2.0d0 * lambda_in
  ! d/dt < a^{\dagger}^{2} a >
  Mat3_OG(3, 2) = 2.0d0 * lambda_in
  Mat3_OG(3, 3) = -(3.0d0 * kappa_p_in - i * Delta_in)
  Mat3_OG(3, 4) = lambda_in
  ! d/dt < a^{\dagger}^{3} >
  Mat3_OG(4, 3) = 3.0d0 * lambda_in
  Mat3_OG(4, 4) = -3.0d0 * (kappa_p_in - i * Delta_in)

  ! Non-homogeneous vector
  B3_OG = 0.0d0

  !----------------------!
  !     Fourth-Order     !
  !----------------------!
  ! Evolution matrix
  Mat4_OG = 0.0d0
  ! d/dt < a^{4} >
  Mat4_OG(1, 1) = -4.0d0 * (kappa_p_in + i * Delta_in)
  Mat4_OG(1, 2) = 4.0d0 * lambda_in
  ! d/dt < a^{\dagger} a^{3} >
  Mat4_OG(2, 1) = lambda_in
  Mat4_OG(2, 2) = -(4.0d0 * kappa_p_in + 2.0d0 * i * Delta_in)
  Mat4_OG(2, 3) = 3.0d0 * lambda_in
  ! d/dt < a^{\dagger}^{2} a^{2} >
  Mat4_OG(3, 2) = 2.0d0 * lambda_in
  Mat4_OG(3, 3) = -4.0d0 * kappa_p_in
  Mat4_OG(3, 4) = 2.0d0 * lambda_in
  ! d/dt < a^{\dagger}^{3} a >
  Mat4_OG(4, 3) = 3.0d0 * lambda_in
  Mat4_OG(4, 4) = -(4.0d0 * kappa_p_in - 2.0d0 * i * Delta_in)
  Mat4_OG(4, 5) = lambda_in
  ! d/dt < a^{\dagger}^{4} >
  Mat4_OG(5, 4) = 4.0d0 * lambda_in
  Mat4_OG(5, 5) = -4.0d0 * (kappa_p_in - i * Delta_in)

  ! Non-homogeneous vector
  B4_OG = 0.0d0

  !==============================================================================!
  !                        CALCULATE STEADY-STATE MOMENTS                        !
  !==============================================================================!
  !---------------------------------------------!
  !     PARAMETRIC OSCILLATOR: SECOND-ORDER     !
  !---------------------------------------------!
  CALL MatrixInverseSS(N2, Mat2_OG, B2_OG, p2_out)

  ! Save steady state photon number of the parametric oscillator
  photon_p_ss = REAL(p2_out(ata))

  !---------------------------------------------!
  !     PARAMETRIC OSCILLATOR: FOURTH-ORDER     !
  !---------------------------------------------!
  ! Set the non-homogeneous vector
  B4_OG = 0.0d0
  B4_OG(1) = 6.0d0 * lambda_in * p2_out(a2)
  B4_OG(2) = 3.0d0 * lambda_in * p2_out(ata)
  B4_OG(3) = lambda_in * p2_out(a2) + &
           & lambda_in * p2_out(at2)
  B4_OG(4) = 3.0d0 * lambda_in * p2_out(ata)
  B4_OG(5) = 6.0d0 * lambda_in * p2_out(at2)

  ! Calculate steady state
  CALL MatrixInverseSS(N4, Mat4_OG, B4_OG, p4_out)

  !=============================!
  !     FIRST-ORDER: FILTER     !
  !=============================!
  DO j = -N_in, N_in
    !===============================================!
    ! First-order: Filter / First-order: Parametric !
    !===============================================!
    !-------------------!
    ! < f_{j} a^{(1)} > !
    !-------------------!
    ! Set the diagonal matrix elements for M
    Mat1 = Mat1_OG
    DO x = 1, N1
      Mat1(x, x) = Mat1(x, x) - (kappa_f_in + i * wl(j))
    END DO

    ! Set the non-homogeneous vector
    B1 = 0.0d0
    B1(1) = -gkl(j) * p2_out(a2)
    B1(2) = -gkl(j) * p2_out(ata)

    ! Calculate steady state
    CALL MatrixInverseSS(N1, Mat1, B1, f1p1_out(j, f, :))
      
    !-----------------------------!
    ! < f^{\dagger}_{j} a^{(1)} > !
    !-----------------------------!
    ! Set the diagonal matrix elements for M
    Mat1 = Mat1_OG
    DO x = 1, N1
      Mat1(x, x) = Mat1(x, x) - (kappa_f_in - i * wl(j))
    END DO

    ! Set the non-homogeneous vector
    B1 = 0.0d0
    B1(1) = -CONJG(gkl(j)) * p2_out(ata)
    B1(2) = -CONJG(gkl(j)) * p2_out(at2)

    ! Calculate steady state
    CALL MatrixInverseSS(N1, Mat1, B1, f1p1_out(j, ft, :))

    !===============================================!
    ! First-order: Filter / Third-order: Parametric !
    !===============================================!
    !-------------------!
    ! < f_{j} a^{(3)} > !
    !-------------------!
    ! Set the diagonal matrix elements for M
    Mat3 = Mat3_OG
    DO x = 1, N3
      Mat3(x, x) = Mat3(x, x) - (kappa_f_in + i * wl(j))
    END DO

    ! Set the non-homogeneous vector
    B3 = 0.0d0
    B3(1) = 3.0d0 * lambda_in * f1p1_out(j, f, a) - &
          & gkl(j) * p4_out(a4)
    B3(2) = lambda_in * f1p1_out(j, f, at) - &
          & gkl(j) * p4_out(ata3)
    B3(3) = lambda_in * f1p1_out(j, f, a) - &
          & gkl(j) * p4_out(at2a2)
    B3(4) = 3.0d0 * lambda_in * f1p1_out(j, f, at) - &
          & gkl(j) * p4_out(at3a)

    ! Calculate steady state
    CALL MatrixInverseSS(N3, Mat3, B3, f1p3_out(j, f, :))

    !-----------------------------!
    ! < f^{\dagger}_{j} a^{(3)} > !
    !-----------------------------!
    ! Set the diagonal matrix elements for M
    Mat3 = Mat3_OG
    DO x = 1, N3
      Mat3(x, x) = Mat3(x, x) - (kappa_f_in - i * wl(j))
    END DO

    ! Set the non-homogeneous vector
    B3 = 0.0d0
    B3(1) = 3.0d0 * lambda_in * f1p1_out(j, ft, a) - &
          & CONJG(gkl(j)) * p4_out(ata3)
    B3(2) = lambda_in * f1p1_out(j, ft, at) - &
          & CONJG(gkl(j)) * p4_out(at2a2)
    B3(3) = lambda_in * f1p1_out(j, ft, a) - &
          & CONJG(gkl(j)) * p4_out(at3a)
    B3(4) = 3.0d0 * lambda_in * f1p1_out(j, ft, at) - &
          & CONJG(gkl(j)) * p4_out(at4)

    ! Calculate steady state
    CALL MatrixInverseSS(N3, Mat3, B3, f1p3_out(j, ft, :))

    ! Close j DO loop
  END DO

  !==============================!
  !     SECOND-ORDER: FILTER     !
  !==============================!
  photon_f_out = 0.0d0
  DO k = -N_in, N_in
    DO j = -N_in, N_in
      !======================!
      ! Second-order: Filter !
      !======================!
      !-----------------!
      ! < f_{j} f_{k} > !
      !-----------------!
      f2_out(j, k, ff) = -(gkl(j) * f1p1_out(k, f, a) + gkl(k) * f1p1_out(j, f, a)) / &
                       & (2.0d0 * kappa_f_in + i * (wl(j) + wl(k)))

      !---------------------------!
      ! < f^{\dagger}_{j} f_{k} > !
      !---------------------------!
      f2_out(j, k, ftf) = -(CONJG(gkl(j)) * f1p1_out(k, f, at) + gkl(k) * f1p1_out(j, ft, a)) / &
                        & (2.0d0 * kappa_f_in - i * (wl(j) - wl(k)))

      ! Mean photon number in filter
      photon_f_out = photon_f_out + REAL(f2_out(j, k, ftf))

      !-------------------------------------!
      ! < f^{\dagger}_{j} f^{\dagger}_{k} > !
      !-------------------------------------!
      f2_out(j, k, ft2) = -(CONJG(gkl(j)) * f1p1_out(k, ft, at) + CONJG(gkl(k)) * f1p1_out(j, ft, at)) / &
                        & (2.0d0 * kappa_f_in - i * (wl(j) + wl(k)))
        
      !=================================================!
      ! Second-order: Filter / Second-order: Parametric !
      !=================================================!
      !-------------------------!
      ! < f_{j} f_{k} a^{(2)} > !
      !-------------------------!
      ! Set the diagonal matrix elements for M
      Mat2 = Mat2_OG
      DO x = 1, N2
        Mat2(x, x) = Mat2(x, x) - (2.0d0 * kappa_f_in + i * (wl(j) + wl(k)))
      END DO

      ! Set the non-homogeneous vector
      B2 = 0.0d0
      B2(1) = lambda_in * f2_out(j, k, ff) - &
            & gkl(j) * f1p3_out(k, f, a3) - &
            & gkl(k) * f1p3_out(j, f, a3)
      B2(2) = -gkl(j) * f1p3_out(k, f, ata2) - &
            & gkl(k) * f1p3_out(j, f, ata2)
      B2(3) = lambda_in * f2_out(j, k, ff) - &
            & gkl(j) * f1p3_out(k, f, at2a) - &
            & gkl(k) * f1p3_out(j, f, at2a)

      ! Calculate steady state
      CALL MatrixInverseSS(N2, Mat2, B2, f2p2_out(j, k, ff, :))

      !-----------------------------------!
      ! < f^{\dagger}_{j} f_{k} a^{(2)} > !
      !-----------------------------------!
      ! Set the diagonal matrix elements for M
      Mat2 = Mat2_OG
      DO x = 1, N2
        Mat2(x, x) = Mat2(x, x) - (2.0d0 * kappa_f_in - i * (wl(j) - wl(k)))
      END DO

      ! Set the non-homogeneous vector
      B2 = 0.0d0
      B2(1) = lambda_in * f2_out(j, k, ftf) - &
            & CONJG(gkl(j)) * f1p3_out(k, f, ata2) - &
            & gkl(k) * f1p3_out(j, ft, a3)
      B2(2) = -CONJG(gkl(j)) * f1p3_out(k, f, at2a) - &
            & gkl(k) * f1p3_out(j, ft, ata2)
      B2(3) = lambda_in * f2_out(j, k, ftf) - &
            & CONJG(gkl(j)) * f1p3_out(k, f, at3) - &
            & gkl(k) * f1p3_out(j, ft, at2a)

      ! Calculate steady state
      CALL MatrixInverseSS(N2, Mat2, B2, f2p2_out(j, k, ftf, :))

      !---------------------------------------------!
      ! < f^{\dagger}_{j} f^{\dagger}_{k} a^{(2)} > !
      !---------------------------------------------!
      ! Set the diagonal matrix elements for M
      Mat2 = Mat2_OG
      DO x = 1, N2
        Mat2(x, x) = Mat2(x, x) - (2.0d0 * kappa_f_in - i * (wl(j) + wl(k)))
      END DO

      ! Set the non-homogeneous vector
      B2 = 0.0d0
      B2(1) = lambda_in * f2_out(j, k, ft2) - &
            & CONJG(gkl(j)) * f1p3_out(k, ft, ata2) - &
            & CONJG(gkl(k)) * f1p3_out(j, ft, ata2)
      B2(2) = -CONJG(gkl(j)) * f1p3_out(k, ft, at2a) - &
            & CONJG(gkl(k)) * f1p3_out(j, ft, at2a)
      B2(3) = lambda_in * f2_out(j, k, ft2) - &
            & CONJG(gkl(j)) * f1p3_out(k, ft, at3) - &
            & CONJG(gkl(k)) * f1p3_out(j, ft, at3)

      ! Calculate steady state
      CALL MatrixInverseSS(N2, Mat2, B2, f2p2_out(j, k, ft2, :))

      ! Close j DO loop
    END DO
  ! Close k DO loop
  END DO

  !=============================!
  !     THIRD-ORDER: FILTER     !
  !=============================!
  DO l = -N_in, N_in
    DO k = -N_in, N_in
      DO j = -N_in, N_in
        !===============================================!
        ! Third-order: Filter / First-order: Parametric !
        !===============================================!
        !-----------------------------------------!
        ! < f^{\dagger}_{j} f_{k} f_{l} a^{(1)} > !
        !-----------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat1 = Mat1_OG
        DO x = 1, N1
          Mat1(x, x) = Mat1(x, x) - (3.0d0 * kappa_f_in - i * (wl(j) - wl(k) - wl(l)))
        END DO

        ! Set the non-homogeneous vector
        B1 = 0.0d0
        B1(1) = -CONJG(gkl(j)) * f2p2_out(k, l, ff, ata) - &
              & gkl(k) * f2p2_out(j, l, ftf, a2) - &
              & gkl(l) * f2p2_out(j, k, ftf, a2)
        B1(2) = -CONJG(gkl(j)) * f2p2_out(k, l, ff, at2) - &
              & gkl(k) * f2p2_out(j, l, ftf, ata) - &
              & gkl(l) * f2p2_out(j, k, ftf, ata)

        ! Calculate steady state
        CALL MatrixInverseSS(N1, Mat1, B1, f3p1_out(j, k, l, ftf2, :))

        !---------------------------------------------------!
        ! < f^{\dagger}_{j} f^{\dagger}_{k} f_{l} a^{(1)} > !
        !---------------------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat1 = Mat1_OG
        DO x = 1, N1
          Mat1(x, x) = Mat1(x, x) - (3.0d0 * kappa_f_in - i * (wl(j) + wl(k) - wl(l)))
        END DO

        ! Set the non-homogeneous vector
        B1 = 0.0d0
        B1(1) = -CONJG(gkl(j)) * f2p2_out(k, l, ftf, ata) - &
              & CONJG(gkl(k)) * f2p2_out(j, l, ftf, ata) - &
              & gkl(l) * f2p2_out(j, k, ft2, a2)
        B1(2) = -CONJG(gkl(j)) * f2p2_out(k, l, ftf, at2) - &
              & CONJG(gkl(k)) * f2p2_out(j, l, ftf, at2) - &
              & gkl(l) * f2p2_out(j, k, ft2, ata)

        ! Calculate steady state
        CALL MatrixInverseSS(N1, Mat1, B1, f3p1_out(j, k, l, ft2f, :))

        ! Close j DO loop
      END DO
      ! Close k DO loop
    END DO
    ! Close l DO loop
  END DO 

  !==============================!
  !     FOURTH-ORDER: FILTER     !
  !==============================!
  DO m = -N_in, N_in
    DO l = -N_in, N_in
      DO k = -N_in, N_in
        DO j = -N_in, N_in
          !======================!
          ! Fourth-order: Cavity !
          !======================!
          !-------------------------------------------------!
          ! < f^{\dagger}_{j} f^{\dagger}_{k} f_{l} f_{m} > !
          !-------------------------------------------------!
          f4_out(j, k, l, m) = -(CONJG(gkl(j)) * f3p1_out(k, l, m, ftf2, at) + &
                             &   CONJG(gkl(k)) * f3p1_out(j, l, m, ftf2, at) + &
                             &   gkl(l) * f3p1_out(j, k, m, ft2f, a) + &
                             &   gkl(m) * f3p1_out(j, k, l, ft2f, a))
          f4_out(j, k, l, m) = f4_out(j, k, l, m) / &
                             & (4.0d0 * kappa_f_in - i * (wl(j) + wl(k)) + i * (wl(l) + wl(m)))

          ! Close j DO loop
        END DO
        ! Close k DO loop
      END DO
      ! Close l DO loop
    END DO
    ! Close m DO loop
  END DO

END SUBROUTINE SteadyStateMoments

!==============================================================================!
!                          G2 CORRELATION SUBROUTINES                          !
!==============================================================================!

! Subroutine to calculate the initial conditions for the auto-correlations
SUBROUTINE G2_InitialConditions(kappa_p_in, Delta_in, lambda_in, &
                              & epsilon_in, N_in, phase_in, &
                              & w0_in, kappa_f_in, dw_in, &
                              & photon_f_out, B2_OG_out, &
                              & p2_out, f1p1_out, f2_out)

  !==============================================================================!
  !                    DEFINING AND DECLARING VARIABLES/ARRAYS                   !
  !==============================================================================!

  IMPLICIT NONE

  !---------------!
  !     INPUT     !
  !---------------!
  ! Parametric cavity decay rate
  REAL(KIND=8)                               :: kappa_p_in
  ! Drive detuning from cavity resonance
  REAL(KIND=8)                               :: Delta_in
  ! Driving amplitude
  REAL(KIND=8)                               :: lambda_in

  ! Filter parameter stuff
  ! Central mode frequency of the filter cavity, with N mode frequencies either side
  REAL(KIND=8), INTENT(IN)                   :: w0_in
  ! Cavity linewidth/transmission of cavity mode
  REAL(KIND=8), INTENT(IN)                   :: kappa_f_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                   :: dw_in
  ! Percentage of fluorecence aimed at cavity
  REAL(KIND=8), INTENT(IN)                   :: epsilon_in
  ! Number of mode either side of w0, 2N + 1 total mode
  INTEGER, INTENT(IN)                        :: N_in
  ! Phase modulation of mode coupling
  REAL(KIND=8), INTENT(IN)                   :: phase_in

  !------------------------------------!
  !     MOMENT EQUATION ARRAY STUFF    !
  !------------------------------------!
  ! Parametric Oscillator: Evolution matrices
  ! First order matrix
  INTEGER, PARAMETER                         :: N1 = 2
  ! Second-order matrix
  INTEGER, PARAMETER                         :: N2 = 3
  ! Third-order matrix
  INTEGER, PARAMETER                         :: N3 = 4
  ! Second-order matrix
  INTEGER, PARAMETER                         :: N4 = 5
  

  ! Integer indices for sigma operators! Integer indices for sigma operators
  INTEGER, PARAMETER                         :: a = 1, at = 2
  INTEGER, PARAMETER                         :: a2 = 1, ata = 2, at2 = 3
  INTEGER, PARAMETER                         :: a3 = 1, ata2 = 2, at2a = 3, at3 = 4
  INTEGER, PARAMETER                         :: a4 = 1, ata3 = 2, at2a2 = 3, at3a = 4, at4 = 5
  ! Integer indices for: a, f^{\dagger}, f^{\dagger} a
  INTEGER, PARAMETER                         :: f = 1, ft = 2
  INTEGER, PARAMETER                         :: ff = 1, ftf = 2, ft2 = 3
  INTEGER, PARAMETER                         :: ftf2 = 1, ft2f = 2

  ! Steady state arrays
  ! Parametric oscillator moments
  ! Second-order: < a^{2} >, < a^{\dagger} a >, < a^{\dagger}^{2} >
  COMPLEX(KIND=8), DIMENSION(N2)                                             :: p2_ss
  ! Fourth-order: < a^{4} >, < a^{\dagger} a^{3} >, < a^{\dagger}^{2} a^{2} >,
  !               < a^{\dagger}^{3} a >< a^{\dagger}^{4} >
  COMPLEX(KIND=8), DIMENSION(N4)                                             :: p4_ss

  ! First-order: Filter / First-order: Parametric
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N1)                              :: f1p1_ss
  ! First-order: Filter / Third-order: Parametric
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N3)                              :: f1p3_ss

  ! Second-order: Filter operators
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3)                      :: f2_ss
  ! Second-order: Filter / Second-order: Parametric
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3, N2)                  :: f2p2_ss

  ! Third-order: Filter / First-order: Parametric
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, 2, N1)      :: f3p1_ss

  ! Fourth-order: Filter
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, -N_in:N_in, -N_in:N_in) :: f4_ss

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Steady state photon number
  REAL(KIND=8), INTENT(OUT)                                                  :: photon_f_out
  ! Non-homogeneous vector
  COMPLEX(KIND=8), DIMENSION(N2), INTENT(OUT)                                :: B2_OG_out

  ! Time integration arrays
  ! Parametric oscillator moments
  ! Second-order: < a^{2} >, < a^{\dagger} a >, < a^{\dagger}^{2} >
  COMPLEX(KIND=8), DIMENSION(N2), INTENT(OUT)                                :: p2_out

  ! First-order: Filter / First-order: Parametric
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N1), INTENT(OUT)                 :: f1p1_out

  ! Second-order: Filter operators
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3), INTENT(OUT)         :: f2_out

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
  ! First-order: Filter / First-order: Parametric
  f1p1_ss = 0.0d0
  ! First-order: Filter / Third-order: Parametric
  f1p3_ss = 0.0d0
  ! Second-order: Filter
  f2_ss = 0.0d0
  ! Second-order: Filter / Second-order Parametric
  f2p2_ss = 0.0d0
  ! Third-order: Filter / First-order Parametric
  f3p1_ss = 0.0d0
  ! Fourth-order: Filter
  f4_ss = 0.0d0

  ! Steady state photon number
  photon_f_out = 0.0d0

  !==============================================================================!
  !                        CALCULATE STEADY-STATE MOMENTS                        !
  !==============================================================================!
  CALL SteadyStateMoments(kappa_p_in, Delta_in, lambda_in, &
                        & epsilon_in, N_in, phase_in, &
                        & w0_in, kappa_f_in, dw_in, &
                        & photon_f_out, &
                        & p2_ss, p4_ss, &
                        & f1p1_ss, f1p3_ss, &
                        & f2_ss, f2p2_ss, &
                        & f3p1_ss, &
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

  p2_out = 0.0d0
  f1p1_out = 0.0d0
  f2_out = 0.0d0
  B2_OG_out = 0.0d0

  ! Cycle through modes
  DO m = -N_in, N_in
    DO j = -N_in, N_in
      ! Second-order: Parametric
      p2_out(:) = p2_out(:) + f2p2_ss(j, m, ftf, :)

      DO k = -N_in, N_in
        ! First-order: Filter / First-order: Parametric
        f1p1_out(k, f, :) = f1p1_out(k, f, :) + f3p1_ss(j, k, m, ftf2, :)
        f1p1_out(k, ft, :) = f1p1_out(k, ft, :) + f3p1_ss(j, k, m, ft2f, :)

        DO l = -N_in, N_in
          ! Second-order: cavity
          f2_out(k, l, ftf) = f2_out(k, l, ftf) + f4_ss(j, k, l, m)

          ! Close l loop
        END DO
        ! Close k loop
      END DO
      ! Non homogeneous vector
      B2_OG_out(1) = B2_OG_out(1) + lambda_in * f2_ss(j, m, ftf)
      B2_OG_out(3) = B2_OG_out(3) + lambda_in * f2_ss(j, m, ftf)

      ! Close j loop
    END DO
    ! Close m loop
  END DO

END SUBROUTINE G2_InitialConditions

! Subroutine to calculate the time evolution of the g2 correlation
SUBROUTINE G2_CalculateRK4(kappa_p_in, Delta_in, lambda_in, &
                         & epsilon_in, N_in, phase_in, &
                         & w0_in, kappa_f_in, dw_in, &
                         & dt_in, tau_steps_in, &
                         & g2_array_out, WRITE_DATA_IN, filename_data_in)

  !==============================================================================!
  !                    DEFINING AND DECLARING VARIABLES/ARRAYS                   !
  !==============================================================================!

  IMPLICIT NONE

  !---------------!
  !     INPUT     !
  !---------------!
  ! Parametric cavity decay rate
  REAL(KIND=8)                                          :: kappa_p_in
  ! Drive detuning from cavity resonance
  REAL(KIND=8)                                          :: Delta_in
  ! Driving amplitude
  REAL(KIND=8)                                          :: lambda_in

  ! Filter parameter stuff
  ! Central mode frequency of the filter cavity, with N mode frequencies either side
  REAL(KIND=8), INTENT(IN)                              :: w0_in
  ! Cavity linewidth/transmission of cavity mode
  REAL(KIND=8), INTENT(IN)                              :: kappa_f_in
  ! Frequency spacing of modes
  REAL(KIND=8), INTENT(IN)                              :: dw_in
  ! Percentage of fluorecence aimed at cavity
  REAL(KIND=8), INTENT(IN)                              :: epsilon_in
  ! Number of mode either side of w0, 2N + 1 total mode
  INTEGER, INTENT(IN)                                   :: N_in
  ! Phase modulation of mode coupling
  REAL(KIND=8), INTENT(IN)                              :: phase_in

  ! Time stuff
  ! Time step
  REAL(KIND=8), INTENT(IN)                              :: dt_in
  ! Maxi_inmum number of steps to integrate for
  INTEGER, INTENT(IN)                                   :: tau_steps_in

  ! Data stuff
  ! Boolean for writing data
  LOGICAL, INTENT(IN)                                   :: WRITE_DATA_IN
  ! Filename for writing data to
  CHARACTER(LEN=*), INTENT(IN)                          :: filename_data_in

  !----------------!
  !     OUTPUT     !
  !----------------!
  ! Data array
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(OUT)  :: g2_array_out

  !------------------------------------!
  !     MOMENT EQUATION ARRAY STUFF    !
  !------------------------------------!
  ! Parametric Oscillator: Evolution matrices
  ! First order matrix
  INTEGER, PARAMETER                                    :: N1 = 2
  COMPLEX(KIND=8), DIMENSION(N1, N1)                    :: Mat1, Mat1_OG
  COMPLEX(KIND=8), DIMENSION(N1)                        :: B1, B1_OG
  ! Second-order matrix
  INTEGER, PARAMETER                                    :: N2 = 3
  COMPLEX(KIND=8), DIMENSION(N2, N2)                    :: Mat2, Mat2_OG
  COMPLEX(KIND=8), DIMENSION(N2)                        :: B2, B2_OG
  ! Third-order matrix
  INTEGER, PARAMETER                                    :: N3 = 4
  ! Second-order matrix
  INTEGER, PARAMETER                                    :: N4 = 5
  

  ! Integer indices for sigma operators! Integer indices for sigma operators
  INTEGER, PARAMETER                                    :: a = 1, at = 2
  INTEGER, PARAMETER                                    :: a2 = 1, ata = 2, at2 = 3
  INTEGER, PARAMETER                                    :: a3 = 1, ata2 = 2, at2a = 3, at3 = 4
  INTEGER, PARAMETER                                    :: a4 = 1, ata3 = 2, at2a2 = 3, at3a = 4, at4 = 5
  ! Integer indices for: a, f^{\dagger}, f^{\dagger} a
  INTEGER, PARAMETER                                    :: f = 1, ft = 2
  INTEGER, PARAMETER                                    :: ff = 1, ftf = 2, ft2 = 3
  INTEGER, PARAMETER                                    :: ftf2 = 1, ft2f = 2

  ! Time integration arrays
  ! Parametric oscillator moments
  ! Second-order: < a^{2} >, < a^{\dagger} a >, < a^{\dagger}^{2} >
  COMPLEX(KIND=8), DIMENSION(N2)                        :: p2
  COMPLEX(KIND=8), DIMENSION(N2)                        :: k1_p2, k2_p2, k3_p2, k4_p2

  ! First-order: Filter / First-order: Parametric
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N1)         :: f1p1
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, 2, N1)         :: k1_f1p1, k2_f1p1, k3_f1p1, k4_f1p1

  ! Second-order: Filter operators
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3) :: f2
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in, -N_in:N_in, 3) :: k1_f2, k2_f2, k3_f2, k4_f2

  !----------------------------!
  !     OTHER USEFUL STUFF     !
  !----------------------------!
  ! Time step integer
  INTEGER                                               :: t
  ! Integer counter
  INTEGER                                               :: j, k, l, m, x
  ! Sample rate for state populations
  INTEGER                                               :: sample_rate
  ! Imaginary i
  COMPLEX(KIND=8), PARAMETER                            :: i = CMPLX(0.0d0, 1.0d0, 8)
  ! pi
  REAL(KIND=8), PARAMETER                               :: pi = 3.1415926535897932384d0
  ! 1 / 6
  REAL(KIND=8), PARAMETER                               :: xis = 1.0d0 / 6.0d0
  ! List of Delta values
  REAL(KIND=8), DIMENSION(-N_in:N_in)                   :: wl
  ! List of mode dependent cascade coupling values
  COMPLEX(KIND=8), DIMENSION(-N_in:N_in)                :: gkl
  ! Blackman window coefficient
  REAL(KIND=8)                                          :: blackman
  ! Steady state photon number of the filter
  REAL(KIND=8)                                          :: photon_f_ss
  ! Complex data
  COMPLEX(KIND=8)                                       :: moment_out

  !==============================================================================!
  !                DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS                !
  !==============================================================================!
  !----------------------------!
  !     EVOLUTION MATRICES     !
  !----------------------------!
  !--------------------!
  !     First-order    !
  !--------------------!
  ! Evolution matrix
  Mat1_OG = 0.0d0
  ! d/dt < a >
  Mat1_OG(1, 1) = -(kappa_p_in + i * Delta_in)
  Mat1_OG(1, 2) = lambda_in
  ! d/dt < a^{\dagger} >
  Mat1_OG(2, 1) = lambda_in
  Mat1_OG(2, 2) = -(kappa_p_in - i * Delta_in)

  ! Non-homogeneous vector
  B1_OG = 0.0d0

  !----------------------!
  !     Second-Order     !
  !----------------------!
  ! Evolution matrix
  Mat2_OG = 0.0d0; Mat2 = 0.0d0
  ! d/dt < a^{2} >
  Mat2_OG(1, 1) = -2.0d0 * (kappa_p_in + i * Delta_in)
  Mat2_OG(1, 2) = 2.0d0 * lambda_in
  Mat2_OG(1, 3) = 0.0d0
  ! d/dt < a^{\dagger} a>
  Mat2_OG(2, 1) = lambda_in
  Mat2_OG(2, 2) = -2.0d0 * kappa_p_in
  Mat2_OG(2, 3) = lambda_in
  ! d/dt < a^{\dagger}^{2} >
  Mat2_OG(3, 1) = 0.0d0
  Mat2_OG(3, 2) = 2.0d0 * lambda_in
  Mat2_OG(3, 3) = -2.0d0 * (kappa_p_in - i * Delta_in)

  !---------------------------------------------!
  !     RESONANCES (wj) AND COUPLINGS (E_j)     !
  !---------------------------------------------!
  ! Allocate array of Delta and gka values
  wl = 0.0d0
  gkl = 0.0d0
  DO j = -N_in, N_in
    IF (N_in == 0) THEN
      wl(j) = w0_in
      gkl(j) = DSQRT(epsilon_in * 2.0d0 * kappa_p_in * kappa_f_in)
    ELSE
      wl(j) = w0_in + DBLE(j) * dw_in
      ! Blackman window coefficient
      blackman = 1.0d0
      ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
      !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))
      ! Mode dependent phase difference
      gkl(j) = DSQRT((epsilon_in / DBLE(2*N_in + 1)) * 2.0d0 * kappa_p_in * kappa_f_in) * blackman * &
             & EXP(i * DBLE(phase_in) * DBLE(j) * pi / DBLE(N_in))
    END IF
  END DO

  !------------------------------------------!
  !     INITALISE OPERATOR MOMENT ARRAYS     !
  !------------------------------------------!

  ! Time integration
  ! Second-order: Cavity and Atom
  f1p1 = 0.0d0
  k1_f1p1 = 0.0d0; k2_f1p1 = 0.0d0; k3_f1p1 = 0.0d0; k4_f1p1 = 0.0d0
  ! Second-order: Cavity
  f2 = 0.0d0
  k1_f2 = 0.0d0; k2_f2 = 0.0d0; k3_f2 = 0.0d0; k4_f2 = 0.0d0

  ! Data
  photon_f_ss = 0.0d0

  !==============================================================================!
  !                         CALCULATE INITIAL CONDITIONS                         !
  !==============================================================================!
  CALL G2_InitialConditions(kappa_p_in, Delta_in, lambda_in, &
                          & epsilon_in, N_in, phase_in, &
                          & w0_in, kappa_f_in, dw_in, &
                          & photon_f_ss, B2_OG, &
                          & p2, f1p1, f2)

  !==============================================================================!
  !                  CALCULATE SECOND-ORDER CORRELATION FUNCTION                 !
  !==============================================================================!
  ! Calculate the sample rate for writing data to the file
  IF (tau_steps_in > 100000) THEN
    ! Set sample rate
    sample_rate = NINT(DBLE(tau_steps_in) / 1d5)
  ELSE
    sample_rate = 1
  END IF

  ! Allocate data array
  ALLOCATE(g2_array_out(0:tau_steps_in)); g2_array_out = 0.0d0

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
        moment_out = moment_out + f2(j, k, ftf)
      END DO
    END DO

    ! Normalise correlation by steady-state photon number
    IF (photon_f_ss .NE. 0.0) THEN
      moment_out = moment_out / (photon_f_ss ** 2)
    END IF

    !-----------------------!
    !     WRITE TO FILE     !
    !-----------------------!
    ! Second-order correlation function
    ! If t_max is really big, only take a sample of results to write to file
    ! so file size isn't huge-mongous.
    IF (MOD(t, sample_rate) == 0) THEN
      ! Save to array
      g2_array_out(NINT(DBLE(t) / DBLE(sample_rate))) = REAL(MOMENT_OUT)

      ! If WRITE_DATA_IN is TRUE, write data to file
      IF (WRITE_DATA_IN .EQV. .TRUE.) THEN
        ! Write to file
        WRITE(4, *) dt_in * DBLE(t), REAL(moment_out)
      END IF
    END IF

    !============================================================================!
    !                  CALCULATE USING FOURTH-ORDER RUNGE-KUTTA                  !
    !============================================================================!
    !----------------------------------------!
    !     INITIALISE RUNGE-KUTTA VECTORS     !
    !----------------------------------------!
    k1_p2 = 0.0d0; k2_p2 = 0.0d0; k3_p2 = 0.0d0; k4_p2 = 0.0d0
    ! Second-order
    k1_f1p1 = 0.0d0; k2_f1p1 = 0.0d0; k3_f1p1 = 0.0d0; k4_f1p1 = 0.0d0
    k1_f2 = 0.0d0; k2_f2 = 0.d0; k3_f2 = 0.0d0; k4_f2 = 0.0d0

    !----------------------------------!
    !     SECOND-ORDER: PARAMETRIC     !
    !----------------------------------!
    k1_p2 = dt_in * (MATMUL(Mat2_OG, p2) + B2_OG)
    k2_p2 = dt_in * (MATMUL(Mat2_OG, (p2 + 0.5d0 * k1_p2)) + B2_OG)
    k3_p2 = dt_in * (MATMUL(Mat2_OG, (p2 + 0.5d0 * k2_p2)) + B2_OG)
    k4_p2 = dt_in * (MATMUL(Mat2_OG, (p2 + k3_p2)) + B2_OG)

    !=============================!
    !     FIRST-ORDER: FILTER     !
    !=============================!
    DO j = -N_in, N_in
      !===============================================!
      ! First-order: Filter / First-order: Parametric !
      !===============================================!
      !-------------------!
      ! < f_{j} a^{(1)} > !
      !-------------------!
      ! Set the diagonal matrix elements for M
      Mat1 = Mat1_OG
      DO x = 1, N1
        Mat1(x, x) = Mat1(x, x) - (kappa_f_in + i * wl(j))
      END DO

      ! k1
      ! Set the non-homogeneous vector
      B1 = 0.0d0
      B1(1) = -gkl(j) * p2(a2)
      B1(2) = -gkl(j) * p2(ata)
      ! Calculate k1
      k1_f1p1(j, f, :) = dt_in * (MATMUL(Mat1, f1p1(j, f, :)) + B1)

      ! k2
      ! Set the non-homogeneous vector
      B1 = 0.0d0
      B1(1) = -gkl(j) * (p2(a2) + 0.5d0 * k1_p2(a2))
      B1(2) = -gkl(j) * (p2(ata) + 0.5d0 * k1_p2(ata))
      ! Calculate k2
      k2_f1p1(j, f, :) = dt_in * (MATMUL(Mat1, (f1p1(j, f, :) + 0.5d0 * k1_f1p1(j, f, :))) + B1)

      ! k3
      ! Set the non-homogeneous vector
      B1 = 0.0d0
      B1(1) = -gkl(j) * (p2(a2) + 0.5d0 * k2_p2(a2))
      B1(2) = -gkl(j) * (p2(ata) + 0.5d0 * k2_p2(ata))
      ! Calculate k3
      k3_f1p1(j, f, :) = dt_in * (MATMUL(Mat1, (f1p1(j, f, :) + 0.5d0 * k2_f1p1(j, f, :))) + B1)

      ! k4
      ! Set the non-homogeneous vector
      B1 = 0.0d0
      B1 = 0.0d0
      B1(1) = -gkl(j) * (p2(a2) + k3_p2(a2))
      B1(2) = -gkl(j) * (p2(ata) + k3_p2(ata))
      ! Calculate k4
      k4_f1p1(j, f, :) = dt_in * (MATMUL(Mat1, (f1p1(j, f, :) + k3_f1p1(j, f, :))) + B1)
      
      !-----------------------------!
      ! < f^{\dagger}_{j} a^{(1)} > !
      !-----------------------------!
      ! Set the diagonal matrix elements for M
      Mat1 = Mat1_OG
      DO x = 1, N1
        Mat1(x, x) = Mat1(x, x) - (kappa_f_in - i * wl(j))
      END DO

      ! k1
      ! Set the non-homogeneous vector
      B1 = 0.0d0
      B1(1) = -CONJG(gkl(j)) * p2(ata)
      B1(2) = -CONJG(gkl(j)) * p2(at2)
      ! Calculate k1
      k1_f1p1(j, ft, :) = dt_in * (MATMUL(Mat1, f1p1(j, ft, :)) + B1)

      ! k2
      ! Set the non-homogeneous vector
      B1 = 0.0d0
      B1(1) = -CONJG(gkl(j)) * (p2(ata) + 0.5d0 * k1_p2(ata))
      B1(2) = -CONJG(gkl(j)) * (p2(at2) + 0.5d0 * k1_p2(at2))
      ! Calculate k1
      k2_f1p1(j, ft, :) = dt_in * (MATMUL(Mat1, (f1p1(j, ft, :) + 0.5d0 * k1_f1p1(j, ft, :))) + B1)

      ! k3
      ! Set the non-homogeneous vector
      B1 = 0.0d0
      B1(1) = -CONJG(gkl(j)) * (p2(ata) + 0.5d0 * k2_p2(ata))
      B1(2) = -CONJG(gkl(j)) * (p2(at2) + 0.5d0 * k2_p2(at2))
      ! Calculate k1
      k3_f1p1(j, ft, :) = dt_in * (MATMUL(Mat1, (f1p1(j, ft, :) + 0.5d0 * k2_f1p1(j, ft, :))) + B1)

      ! k4
      ! Set the non-homogeneous vector
      B1 = 0.0d0
      B1(1) = -CONJG(gkl(j)) * (p2(ata) + k3_p2(ata))
      B1(2) = -CONJG(gkl(j)) * (p2(at2) + k3_p2(at2))
      ! Calculate k1
      k4_f1p1(j, ft, :) = dt_in * (MATMUL(Mat1, (f1p1(j, ft, :) + k3_f1p1(j, ft, :))) + B1)

      ! Close j DO loop
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
        k1_f2(j, k, ftf) = -dt_in * (2.0d0 * kappa_f_in - i * (wl(j) - wl(k))) * f2(j, k, ftf) - &
                         & dt_in * CONJG(gkl(j)) * f1p1(k, f, at) - &
                         & dt_in * gkl(k) * f1p1(j, ft, a)
        k2_f2(j, k, ftf) = -dt_in * (2.0d0 * kappa_f_in - i * (wl(j) - wl(k))) * (f2(j, k, ftf) + 0.5d0 * k1_f2(j, k, ftf)) - &
                         & dt_in * CONJG(gkl(j)) * (f1p1(k, f, at) + 0.5d0 * k1_f1p1(k, f, at)) - &
                         & dt_in * gkl(k) * (f1p1(j, ft, a) + 0.5d0 * k1_f1p1(j, ft, a))
        k3_f2(j, k, ftf) = -dt_in * (2.0d0 * kappa_f_in - i * (wl(j) - wl(k))) * (f2(j, k, ftf) + 0.5d0 * k2_f2(j, k, ftf)) - &
                         & dt_in * CONJG(gkl(j)) * (f1p1(k, f, at) + 0.5d0 * k2_f1p1(k, f, at)) - &
                         & dt_in * gkl(k) * (f1p1(j, ft, a) + 0.5d0 * k2_f1p1(j, ft, a))
        k4_f2(j, k, ftf) = -dt_in * (2.0d0 * kappa_f_in - i * (wl(j) - wl(k))) * (f2(j, k, ftf) + k3_f2(j, k, ftf)) - &
                         & dt_in * CONJG(gkl(j)) * (f1p1(k, f, at) + k3_f1p1(k, f, at)) - &
                         & dt_in * gkl(k) * (f1p1(j, ft, a) + k3_f1p1(j, ft, a))

        ! Close k loop
      END DO
      ! Close j loop
    END DO


    !============================================================================!
    !                   UPDATE ARRAYS FROM RUNGE-KUTTA ARRAYS                    !
    !============================================================================!
    !-------------------------------!
    !     Parametric Oscillator     !
    !-------------------------------!
    p2 = p2 + xis * (k1_p2 + 2.0d0 * (k2_p2 + k3_p2) + k4_p2)

    !-----------------------------!
    !     First-order: Filter     !
    !-----------------------------!
    f1p1 = f1p1 + xis * (k1_f1p1 + 2.0d0 * (k2_f1p1 + k3_f1p1) + k4_f1p1)

    !------------------------------!
    !     Second-order: Filter     !
    !------------------------------!
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