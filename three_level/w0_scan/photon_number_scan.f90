! The system in this program is a resonantly driven three-level ladder-type atom
! that is coupled into a multi-mode array of filter cavities, as described in
! some notes somewhere (or with this program if this goes to anyone).

! After verifying the moment equations with "all_order_moments_RK4.f90" with
! QuTiP, this program calculates the steady state values for each moment by
! inverting the Lindblad matrix for the atomic and cavity-atom coupled moments.
! The steady state values are printed to the console and can be compared with
! similar outouts from the RK4 program.

! The input parameters are taken from a NameList file [filename_ParamList] which,
! by default, points to "./ParamList.nml". The code can thus be compiled once,
! and parameters can be changed in the NameList file for subsequent runs.

! For the default filenames, the NameList file "./ParamList.nml" MUST EXIST IN
! THE WORKING DIRECTORY.

! To compiled the code, I use the Intel Parallel Studio compiler IFORT with the
! command:
!       (LINUX): ifort -O3 -o ss -mkl just_steady_states.f90
!     (WINDOWS): ifort /O3 /o ss /Qmkl just_steady_states.f90
! where the -O3 (/O3) flag gives maximum optimisation, the -o (/o) ss flag
! names the executable as "ss" ("ss.exe"), and the -mkl (/Qmkl) flag links
! the program to Intel's Math Kernel Library, to make use of the LAPACK routines.

PROGRAM THREE_LEVEL_ATOM_MULTI_MODE_FILTER_PHOTON_NUMBER_SCAN

! Import the MKL Library LAPACK95, from which the eigenvalue/eigenvector and
! matrix inversion subroutines come from.
! MUST INCLUDE THE -mkl OR /Qmkl FLAG AS A COMPILER OPTION IF USING INTEL.
! Otherwise you'll have to link it yourself and I don't know how to do that :)

USE LAPACK95

! The subroutines used from LAPACK are:
! - zGETRF - Calculates LU-factorisation of a complexmatrix so it can be
!            inverted by...,
! - zGETRI - Calculates the inverse of a complex matrix,
! - zGEEV  - Calculates the eigenvalues and eigevectors of a complex matrix.

!==============================================================================!
!                    DEFINING AND DECLARING VARIABLES/ARRAYS                   !
!==============================================================================!

IMPLICIT NONE

!---------------------------------!
!     SYSTEM PARAMETERS STUFF     !
!---------------------------------!
! Atomic decay rate
REAL(KIND=8)                                           :: Gamma
! Driving amplitude
REAL(KIND=8)                                           :: Omega
! Atomic anharmonicity
REAL(KIND=8)                                           :: alpha
! Drive detuning from two-photon resonance
REAL(KIND=8)                                           :: delta
! Dipole moment ratio
REAL(KIND=8)                                           :: xi

! Filter parameter stuff
! Central mode frequency of the filter cavity, with N mode frequencies either side
REAL(KIND=8)                                           :: w0a
! REAL(KIND=8)                                           :: w0b
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                           :: kappaa
! REAL(KIND=8)                                           :: kappab
! Percentage of fluorecence aimed at cavity
REAL(KIND=8)                                           :: epsilon
! Number of mode either side of w0, 2N + 1 total mode
INTEGER                                                :: N
! Frequency spacing of modes
REAL(KIND=8)                                           :: dwa
! REAL(KIND=8)                                           :: dwb
! Phase modulation of mode coupling
INTEGER                                                :: phase
! List of Delta values
REAL(KIND=8), DIMENSION(:), ALLOCATABLE                :: wl
! List of mode dependent cascade coupling values
COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE             :: gkl
! Blackman window coefficient
REAL(KIND=8)                                           :: blackman
! kappa, w0, and dw values
REAL(KIND=8)                                           :: kappa, w0, dw

! Scanning parameters
! Number of scans either side of centre
INTEGER                                                :: no_scans
! List of central frequencies
REAL(KIND=8), DIMENSION(:), ALLOCATABLE                :: w0l
! Scanning integer
INTEGER                                                :: scan
! Max/min value of D0a [Default = 30.0]
REAL(KIND=8)                                           :: scan_max
! Scan step size [Default = 0.5]
REAL(KIND=8)                                           :: scan_step

! ! Time stuff
! Time step
! REAL(KIND=8)                                           :: dt
! Maximum time to integrate for
! REAL(KIND=8)                                           :: t_max, tau1_max, tau2_max
! Maximum number of steps to integrate for
! INTEGER                                                :: t_steps
! Time step integer
! INTEGER                                                :: t
! Runtime variables
REAL(KIND=8)                                           :: start_time, end_time

!------------------------------------!
!     MOMENT EQUATION ARRAY STUFF    !
!------------------------------------!
! Dimension of M matrix
INTEGER, PARAMETER                                     :: N_mat = 8
! M matrix (filled as transpose)
COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)               :: Mat, Mat_OG, Mat_inv
! Non-homogeneous vector
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: B_vec, B_OG

! Steady state arrays
! First-order moments: Atomic equations (< \sigma >)
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: sigma_ss
! First-order moments: Cavity (< a >, < a^{\dagger} >)
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: cav1_ss
! Second-order moments: Cavity and atom (< a \sigma >, < a^{\dagger} \sigma >
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cavsig2_ss
! Second-order moments: Cavity (< a^{\dagger} a >)
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cav2_ss

! Integer indices for sigma operators
INTEGER, PARAMETER                                     :: gg = 1, ge = 2, eg = 3
INTEGER, PARAMETER                                     :: ee = 4, ef = 5, fe = 6
INTEGER, PARAMETER                                     :: gf = 7, fg = 8
! Integer indices for: a, a^{\dagger}, a^{\dagger} a
INTEGER, PARAMETER                                     :: a = 1, at = 2, ata = 3

!----------------------------!
!     OTHER USEFUL STUFF     !
!----------------------------!
! Integer counter
INTEGER                                                :: j, k, l, m, x, y
! Imaginary i
COMPLEX(KIND=8), PARAMETER                             :: i = CMPLX(0.0d0, 1.0d0, 8)
! pi
REAL(KIND=8), PARAMETER                                :: pi = 3.1415926535897932384d0
! 1 / 6
REAL(KIND=8), PARAMETER                                :: xis = 1.0d0 / 6.0d0
! Atomic population and photon number
! REAL(KIND=8)                                           :: popg, pope, popf, photon
! Steady state photon number
REAL(KIND=8)                                           :: photon_ss
! Complex data
COMPLEX(KIND=8)                                        :: moment_out

! Run-time percentage printing
! Print completion time check
LOGICAL, PARAMETER                                     :: progress_bar = .TRUE.
! Runtime variables
REAL(KIND=8)                                           :: loop_start_time, loop_check_time
REAL(KIND=8)                                           :: loop_run_time, loop_remaining_time
! Ten percent of time steps for loop
INTEGER                                                :: ten_percent
! Percentage of time steps completed
INTEGER                                                :: percentage
CHARACTER(LEN=28), PARAMETER                           :: FMT_ss = "(T2,I3,A13,F9.2,A19,F9.2,A1)"

!--------------------------!
!     SUBROUTINE STUFF     !
!--------------------------!
! Work space dimension
INTEGER, PARAMETER                                     :: LWMAX = 300
INTEGER                                                :: LWORK
! Work space array
COMPLEX(KIND=8), DIMENSION(LWMAX)                      :: WORK
REAL(KIND=8), DIMENSION(2*N_mat)                       :: RWORK
! Eigenvalues of matrix M
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: eigval
! S, and S^{-1} matrix for diagonalising eigenvectors
COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)               :: S, Sinv
! LU-factorisation array
INTEGER, DIMENSION(N_mat)                              :: IPIV
! Info IO
INTEGER                                                :: INFO

!------------------------!
!     FILENAME STUFF     !
!------------------------!
! Paramert Name List
CHARACTER(LEN=15), PARAMETER :: filename_ParamList = "./ParamList.nml"
! Filename of parameters
CHARACTER(LEN=27), PARAMETER :: filename_parameters = "./data_files/parameters.txt"
! Filename for photon number
CHARACTER(LEN=23), PARAMETER :: filename_photon = "./data_files/photon.txt"

!==============================================================================!
!                 NAMELIST AND PARAMETERS TO BE READ FROM FILE                 !
!==============================================================================!
! NameList things
! Status and unit integers
INTEGER            :: ISTAT, IUNIT
! Line to be read from file
CHARACTER(LEN=512) :: LINE
! Namelist parameters
NAMELIST /ATOM/ Gamma, Omega, alpha, delta, xi
NAMELIST /FILTER/ epsilon, N, phase
NAMELIST /CAVITYA/ kappaa, w0a, dwa
NAMELIST /SCANLIST/ scan_max, scan_step

! Call start time from CPU_TIME
CALL CPU_TIME(start_time)

! Read the parameters from the NAMELIST file
IUNIT = 420
OPEN(IUNIT, FILE=filename_ParamList, STATUS="OLD", DELIM="QUOTE")

READ(IUNIT, NML=ATOM, IOSTAT=ISTAT)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') LINE
  CLOSE(IUNIT)
  PRINT *, "Invalid line in ATOM namelist: " // TRIM(line)
  CALL EXIT(1)
END IF

READ(IUNIT, NML=FILTER, IOSTAT=ISTAT)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') LINE
  CLOSE(IUNIT)
  PRINT *, "Invalid line in FILTER namelist: " // TRIM(line)
  CALL EXIT(1)
END IF

READ(IUNIT, NML=CAVITYA, IOSTAT=ISTAT)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') LINE
  CLOSE(IUNIT)
  PRINT *, "Invalid line in CAVITYA namelist: " // TRIM(line)
  CALL EXIT(1)
END IF

READ(IUNIT, NML=SCANLIST, IOSTAT=ISTAT)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') LINE
  CLOSE(IUNIT)
  PRINT *, "Invalid line in SCANLIST namelist: " // TRIM(line)
  CALL EXIT(1)
END IF

! Close unit
CLOSE(IUNIT)

! Number of time-steps
! t_steps = NINT(t_max / dt)

! Max number of scans
no_scans = NINT(scan_max / scan_step)

! Set system parameters
kappa = kappaa
w0 = w0a
dw = dwa

!==============================================================================!
!                DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS                !
!==============================================================================!
!------------------------!
!     BLOCH MATRIX M     !
!------------------------!
Mat_OG = 0.0d0
! Row 1: d/dt |g><g|
Mat_OG(1, 2) = -i * 0.5d0 * Omega
Mat_OG(1, 3) = i * 0.5d0 * Omega
Mat_OG(1, 4) = Gamma
! Row 2: d/dt |g><e|
Mat_OG(2, 1) = -i * 0.5d0 * Omega
Mat_OG(2, 2) = -(0.5d0 * Gamma - i * ((0.5d0 * alpha) + delta))
Mat_OG(2, 4) = i * 0.5d0 * Omega
Mat_OG(2, 5) = Gamma * xi
Mat_OG(2, 7) = -i * xi * 0.5d0 * Omega
! Row 3: d/dt |e><g|
Mat_OG(3, 1) = i * 0.5d0 * Omega
Mat_OG(3, 3) = -(0.5d0 * Gamma + i * ((0.5d0 * alpha) + delta))
Mat_OG(3, 4) = -i * 0.5d0 * Omega
Mat_OG(3, 6) = Gamma * xi
Mat_OG(3, 8) = i * xi * 0.5d0 * Omega
! Row 4: d/dt |e><e|
Mat_OG(4, 1) = -Gamma * (xi ** 2)
Mat_OG(4, 2) = i * 0.5d0 * Omega
Mat_OG(4, 3) = -i * 0.5d0 * Omega
Mat_OG(4, 4) = -Gamma * (1.0d0 + (xi ** 2))
Mat_OG(4, 5) = -i * xi * 0.5d0 * Omega
Mat_OG(4, 6) = i * xi * 0.5d0 * Omega
! Row 5: d/dt |e><f|
Mat_OG(5, 1) = -i * xi * 0.5d0 * Omega
Mat_OG(5, 4) = -i * xi * Omega
Mat_OG(5, 5) = -(0.5d0 * Gamma * (1.0d0 + (xi ** 2)) + i * ((0.5d0 * alpha) - delta))
Mat_OG(5, 7) = i * 0.5d0 * Omega
! Row 6: d/dt |f><e|
Mat_OG(6, 1) = i * xi * 0.5d0 * Omega
Mat_OG(6, 4) = i * xi * Omega
Mat_OG(6, 6) = -(0.5d0 * Gamma * (1.0d0 + (xi ** 2)) - i * ((0.5d0 * alpha) - delta))
Mat_OG(6, 8) = -i * 0.5d0 * Omega
! Row 7: d/dt |g><f|
Mat_OG(7, 2) = -i * xi * 0.5d0 * Omega
Mat_OG(7, 5) = i * 0.5d0 * Omega
Mat_OG(7, 7) = -(0.5d0 * Gamma * (xi ** 2) - 2.0d0 * i * delta)
! Row 8: d/dt |g><f|
Mat_OG(8, 3) = i * xi * 0.5d0 * Omega
Mat_OG(8, 6) = -i * 0.5d0 * Omega
Mat_OG(8, 8) = -(0.5d0 * Gamma * (xi ** 2) + 2.0d0 * i * delta)

!--------------------------------!
!     NON-HOMOGENEOUS VECTOR     !
!--------------------------------!
B_OG = 0.0d0
B_OG(4) = Gamma * (xi ** 2)
B_OG(5) = i * xi * 0.5d0 * Omega
B_OG(6) = -i * xi * 0.5d0 * Omega

!---------------------------------------------!
!     RESONANCES (wj) AND COUPLINGS (E_j)     !
!---------------------------------------------!
! Allocate array of Delta and gka values
ALLOCATE(wl(-N:N))
wl = 0.0d0
ALLOCATE(gkl(-N:N))
gkl = 0.0d0
DO j = -N, N
  IF (N == 0) THEN
    wl(j) = w0
    gkl(j) = DSQRT(epsilon * gamma * kappa)
  ELSE
    wl(j) = w0 + DBLE(j) * dw
    ! Blackman window coefficient
    blackman = 1.0d0
    ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
    !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))
    ! Mode dependent phase difference
    gkl(j) = DSQRT((epsilon / DBLE(2*N + 1)) * gamma * kappa) * blackman * EXP(i * DBLE(phase) * DBLE(j) * pi / DBLE(N))
  END IF
END DO

! Allocate central frequency list
ALLOCATE(w0l(-no_scans:no_scans))
w0l = 0.0d0
DO scan = -no_scans, no_scans
  w0l(scan) = 0.0d0 + DBLE(scan) * scan_step
END DO

!------------------------------------------!
!     INITALISE OPERATOR MOMENT ARRAYS     !
!------------------------------------------!
! Steady states
! First-order: Cavity
ALLOCATE(cav1_ss(-N:N, 2)); cav1_ss = 0.0d0
! Second-order: Cavity and Atom
ALLOCATE(cavsig2_ss(-N:N, 2, N_mat)); cavsig2_ss = 0.0d0
! Second-order: Cavity
ALLOCATE(cav2_ss(-N:N, -N:N, 3)); cav2_ss = 0.0d0

!==============================================================================!
!                           WRITE PARAMETERS TO FILE                           !
!==============================================================================!

! Open file to write time to
OPEN(UNIT=1, FILE=filename_parameters, STATUS='REPLACE', ACTION='WRITE')
! Write parameter
WRITE(1,*) "Parameters are in the following order:"
WRITE(1,"(A10,F25.15)") "Gamma =", Gamma
WRITE(1,"(A10,F25.15)") "Omega =", Omega
WRITE(1,"(A10,F25.15)") "alpha =", alpha
WRITE(1,"(A10,F25.15)") "delta =", delta
WRITE(1,"(A10,F25.15)") "xi =", xi
WRITE(1,"(A10,F25.15)") "w0 =", w0
WRITE(1,"(A10,F25.15)") "kappa =", kappa
WRITE(1,"(A11,F25.15)") "dw = ", dw
WRITE(1,"(A10,F25.15)") "epsilon =", epsilon
WRITE(1,"(A11,I9)") "N = ", N
WRITE(1,"(A11,I9)") "phase = ", phase
! WRITE(1,"(A10,F25.15)") "dt =", dt
! WRITE(1,"(A10,F25.15)") "Max t =", t_max
! WRITE(1,"(A10,F25.15)") "Max tau1 =", tau1_max
! WRITE(1,"(A10,F25.15)") "Max tau2 =", tau2_max
! Close file
CLOSE(1)

!==============================================================================!
!                        CALCULATE STEADY-STATE MOMENTS                        !
!==============================================================================!
! Calculate atomic steady states
!---------------------------!
!     FIRST-ORDER: ATOM     !
!---------------------------!
IF (xi .NE. 0.0d0) THEN
  ! Set matrix to be inverted
  Mat_inv = Mat_OG

  ! Perform LU-factorization of matrix
  CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
  ! Invert Matrix (Optimal LWORK = 8)
  LWORK = 8
  CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

  ! Calculate steady states
  sigma_ss = 0.0d0
  sigma_ss = -MATMUL(Mat_inv, B_OG)
ELSE IF (xi .EQ. 0.0d0) THEN
  !------------------------------------------------!
  !     CALCULATE EIGENVALUES AND EIGENVECTORS     !
  !------------------------------------------------!
  ! Set Lindblad matrix
  Mat = Mat_OG

  ! ! Query optimal work space
  ! LWORK = -1
  ! CALL zGEEV('N', 'V', N_mat, M, N_mat, eigval, S, N_mat, S, N_mat, WORK, LWORK, RWORK, INFO)
  ! ! Set optimal work space
  ! LWORK = MIN(LWMAX, INT(WORK(1)))

  ! Calculate eigenvalues and eigenvectors (Optimal LWORK = 264)
  LWORK = 264
  CALL zGEEV('N', 'V', N_mat, Mat, N_mat, eigval, S, N_mat, S, N_mat, WORK, LWORK, RWORK, INFO)

  ! Cycle through eigenvalues and, for the eigenvalue that = 0, use that
  ! eigenvector as the steady state
  DO x = 1, N_mat
    IF (ABS(REAL(eigval(x))) .LT. 1D-10 .AND. ABS(REAL(eigval(x))) .LT. 1D-10) THEN
      ! Save steady state eigenvector
      DO j = 1, N_mat
        sigma_ss(j) = S(j, x)
      END DO
    END IF
  END DO

  ! Normalise sigma_ss so |g><g| + |e><e| = 1
  sigma_ss = sigma_ss / (REAL(sigma_ss(1)) + REAL(sigma_ss(4)))
END IF

! Cycle through each different \omega_{0} value and calculate steady states
! Ten percent of time steps
ten_percent = NINT((DBLE(SIZE(w0l)) / 10.0d0))

! Open file to write time and data to
OPEN(UNIT=2, FILE=filename_photon, STATUS='REPLACE', ACTION='WRITE', RECL=4000)

! Start scan loop
DO scan = -no_scans, no_scans
  ! Set central resonance frequency
  w0 = w0l(scan)
  ! Initialise mode resonance frequency list
  wl = 0.0d0
  DO j = -N, N
    IF (N == 0) THEN
      wl(j) = w0
    ELSE
      wl(j) = w0 + DBLE(j) * dw
    END IF
  END DO

  ! Cycle through modes
  DO j = -N, N
    !-----------------------------!
    !     FIRST-ORDER: CAVITY     !
    !-----------------------------!
    !-----------!
    ! < a_{j} > !
    !-----------!
    cav1_ss(j, a) = -gkl(j) * sigma_ss(ge) + &
                  & -gkl(j) * xi * sigma_ss(ef)
    cav1_ss(j, a) = cav1_ss(j, a) / &
                  & (kappa + i * wl(j))

    !---------------------!
    ! < a^{\dagger}_{j} > !
    !---------------------!
    cav1_ss(j, at) = -CONJG(gkl(j)) * sigma_ss(eg) + &
                   & -CONJG(gkl(j)) * xi * sigma_ss(fe)
    cav1_ss(j, at) = cav1_ss(j, at) / &
                   & (kappa - i * wl(j))

    !---------------------------------------!
    !     SECOND-ORDER: CAVITY AND ATOM     !
    !---------------------------------------!
    !------------------!
    ! < a_{j} \sigma > !
    !------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa + i * wl(j))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -gkl(j) * sigma_ss(ge)
    B_vec(2) = -gkl(j) * xi * sigma_ss(gf)
    B_vec(3) = -gkl(j) * sigma_ss(ee)
    B_vec(4) = Gamma * (xi ** 2) * cav1_ss(j, a) + &
             & -gkl(j) * xi * sigma_ss(ef)
    B_vec(5) = i * xi * 0.5d0 * Omega * cav1_ss(j, a)
    B_vec(6) = -i * xi * 0.5d0 * Omega * cav1_ss(j, a) + &
             & gkl(j) * xi * sigma_ss(gg) + &
             & gkl(j) * xi * sigma_ss(ee) + &
             & -gkl(j) * xi
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j) * sigma_ss(fe)

    ! Set inverse matrix
    Mat_inv = Mat
    ! Perform LU-factorization of matrix
    CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
    IF (INFO .NE. 0) THEN
      WRITE(*, '(A26,I1,A15,I3)') "zGETRF M failed on cavsig(", j , ", a) :( INFO = ", INFO
      STOP
    END IF

    ! ! Query optimal work space
    ! LWORK = -1
    ! CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)
    ! ! Set optimal work space and run again
    ! LWORK = MIN(LWMAX, INT(WORK(1)))

    ! Invert Matrix (Optimal LWORK = 8)
    LWORK = 8
    CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

    ! Calculate steady state
    cavsig2_ss(j, a, :) = -MATMUL(Mat_inv, B_vec)

    !----------------------------!
    ! < a^{\dagger}_{j} \sigma > !
    !----------------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa - i * wl(j))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -CONJG(gkl(j)) * sigma_ss(eg)
    B_vec(2) = -CONJG(gkl(j)) * sigma_ss(ee)
    B_vec(3) = -CONJG(gkl(j)) * xi * sigma_ss(fg)
    B_vec(4) = Gamma * (xi ** 2) * cav1_ss(j, at) + &
             & -CONJG(gkl(j)) * xi * sigma_ss(fe)
    B_vec(5) = i * xi * 0.5d0 * Omega * cav1_ss(j, at) + &
             & CONJG(gkl(j)) * xi * sigma_ss(gg) + &
             & CONJG(gkl(j)) * xi * sigma_ss(ee) + &
             & -CONJG(gkl(j)) * xi
    B_vec(6) = -i * xi * 0.5d0 * Omega * cav1_ss(j, at)
    B_vec(7) = -CONJG(gkl(j)) * sigma_ss(ef)
    B_vec(8) = 0.0d0

    ! Set inverse matrix
    Mat_inv = Mat
    ! Perform LU-factorization of matrix
    CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
    IF (INFO .NE. 0) THEN
      WRITE(*, '(A26,I1,A16,I3)') "zGETRF M failed on cavsig(", j , ", at) :( INFO = ", INFO
      STOP
    END IF

    ! ! Query optimal work space
    ! LWORK = -1
    ! CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)
    ! ! Set optimal work space and run again
    ! LWORK = MIN(LWMAX, INT(WORK(1)))

    ! Invert Matrix (Optimal LWORK = 8)
    LWORK = 8
    CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

    ! Calculate steady state
    cavsig2_ss(j, at, :) = -MATMUL(Mat_inv, B_vec)

    ! Close j loop
  END DO

  moment_out = 0.0d0
  ! Cycle through modes
  DO k = -N, N
    DO j = -N, N
      !------------------------------!
      !     SECOND-ORDER: CAVITY     !
      !------------------------------!
      !---------------------------!
      ! < a^{\dagger}_{j} a_{k} > !
      !---------------------------!
      cav2_ss(j, k, ata) = -CONJG(gkl(j)) * cavsig2_ss(k, a, eg) + &
                         & -CONJG(gkl(j)) * xi * cavsig2_ss(k, a, fe) + &
                         & -gkl(k) * cavsig2_ss(j, at, ge) + &
                         & -gkl(k) * xi * cavsig2_ss(j, at, ef)
      cav2_ss(j, k, ata) = cav2_ss(j, k, ata) / &
                         & (2.0d0 * kappa - i * (wl(j) - wl(k)))

      ! Update photon number
      moment_out = moment_out + cav2_ss(j, k, ata)

      ! Close j loop
    END DO
    ! Close k loop
  END DO

  photon_ss = REAL(moment_out)

  !-----------------------!
  !     WRITE TO FILE     !
  !-----------------------!
  ! Write to file
  WRITE(2,*) w0, photon_ss

  ! Check percentage
  IF (progress_bar .EQV. .TRUE.) THEN
    IF (MOD(scan, ten_percent) == 0 .AND. scan /= -no_scans) THEN
      CALL CPU_TIME(loop_check_time)
      percentage = NINT((100.0d0 * DBLE(scan + no_scans)) / (1.0d0 * DBLE(SIZE(w0l))))
      loop_run_time = loop_check_time - loop_start_time
      loop_remaining_time = ((100.0d0 * (loop_check_time - loop_start_time)) / (1.0d0 * percentage)) - loop_run_time
      WRITE(*, FMT_ss) percentage, "%. Run time: ", loop_run_time, "s. Est. time left: ", &
                     & loop_remaining_time, "s"
    END IF
  END IF

  ! Close scan loop
END DO

! Close file
CLOSE(2)

!==============================================================================!
!                                END OF PROGRAM                                !
!==============================================================================!

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

END PROGRAM THREE_LEVEL_ATOM_MULTI_MODE_FILTER_PHOTON_NUMBER_SCAN
