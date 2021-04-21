! The system in this program is a resonantly driven two-level atom that is
! coupled into a multi-mode array of filter cavities, as described in some notes
! somewhere (or with this program if this goes to anyone).

! The operator moment steady states are calculated using the inverse matrix/
! analytical method. Using quantum regression theorem, we then use Runge-Kutta
! fourth-order to numerically integrate the moment equations, suing the steady
! states and initial conditions, to calculate the normalised second-order
! cross-correlation function, correlation two different-frequency photons.

! The input parameters are taken from a NameList file [filename_ParamList] which,
! by default, points to "./ParamList.nml". The code can thus be compiled once,
! and parameters can be changed in the NameList file for subsequent runs.

! The parameters for each run are written to [filename_parameters] which, by
! default, points to "./data_files/parameters.txt".

! The normalised first-order correlation function is written to [filename_g1]
! which, by default, points to "./data_files/cross_g1_corr.txt". The file has
! two columns:
!                            t     REAL(g2)

! For the default filenames, the folder "./data_files/" and NameList file
! "./ParamList.nml" MUST EXIST IN THE WORKING DIRECTORY.

! To compiled the code, I use the Intel Parallel Studio compiler IFORT with the
! command
!       (LINUX): ifort -O3 -o g2 -mkl g2_RK4.f90
!     (WINDOWS): ifort /O3 /o g2 /Qmkl g2_RK4.f90
! where the -O3 (/O3) flag gives maximum optimisation, the -o (/o) g1 flag
! names the executable as "g2" ("g2.exe"), and the -mkl (/Qmkl) flag links
! the program to Intel's Math Kernel Library, to make use of the LAPACK routines.

PROGRAM TWO_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_CROSS_G2

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
REAL(KIND=8)                                           :: gamma
! Driving amplitude
REAL(KIND=8)                                           :: Omega

! Filter parameter stuff
! Central mode frequency of the filter cavity, with N mode frequencies either side
REAL(KIND=8)                                           :: w0a
REAL(KIND=8)                                           :: w0b
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                           :: kappaa
REAL(KIND=8)                                           :: kappab
! Percentage of fluorecence aimed at cavity
REAL(KIND=8)                                           :: epsilon
! Number of mode either side of w0, 2N + 1 total mode
INTEGER                                                :: N
! Frequency spacing of modes
REAL(KIND=8)                                           :: dwa
REAL(KIND=8)                                           :: dwb
! Phase modulation of mode coupling
INTEGER                                                :: phase
! Blackman window coefficient
REAL(KIND=8)                                           :: blackman
! Kappa values for both cavities
REAL(KIND=8), DIMENSION(2)                             :: kappa
! List of omega values
REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE             :: wl
! List of mode dependent cascade coupling values
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: gkl

! Time stuff
! Time step
REAL(KIND=8)                                           :: dt
! Maximum time to integrate for
REAL(KIND=8)                                           :: t_max, tau1_max, tau2_max
! Maximum number of steps to integrate for
INTEGER                                                :: t_steps, tau_steps
! Time step integer
INTEGER                                                :: t
! Runtime variables
REAL(KIND=8)                                           :: start_time, end_time

!------------------------------------!
!     MOMENT EQUATION ARRAY STUFF    !
!------------------------------------!
! Dimension of M matrix
INTEGER, PARAMETER                                              :: N_mat = 3
! M matrix (filled as transpose)
COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)                        :: Mat, Mat_OG, Mat_inv
! Non-homogeneous vector
COMPLEX(KIND=8), DIMENSION(N_mat)                               :: B, B_OG

! First-order moments: Atomic equations (< \sigma >)
COMPLEX(KIND=8), DIMENSION(N_mat)                               :: sigma
COMPLEX(KIND=8), DIMENSION(N_mat)                               :: k1_sigma, k2_sigma, k3_sigma, k4_sigma
COMPLEX(KIND=8), DIMENSION(N_mat)                               :: sigma_ss

! First-order moments: Cavity (< a >, < a^{\dagger} >)
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE                   :: cav1
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE                   :: k1_cav1, k2_cav1, k3_cav1, k4_cav1
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE                :: cav1_ss

! Second-order moments: Cavity and atom (< a \sigma >, < a^{\dagger} \sigma >
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE                :: cavsig2
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE                :: k1_cavsig2, k2_cavsig2, k3_cavsig2, k4_cavsig2
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE             :: cavsig2_ss

! Second-order moments: Cavity (< a^{\dagger} a >)
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE                :: cav2
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE                :: k1_cav2, k2_cav2, k3_cav2, k4_cav2
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE          :: cav2_ss

! Third-order moments: Cavity and atom (< a^{2} \sigma >, < a^{\dagger 2} \sigma >, < a^{\dagger} a \sigma >)
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cavsig3
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: k1_cavsig3, k2_cavsig3, k3_cavsig3, k4_cavsig3
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :), ALLOCATABLE       :: cavsig3_ss

! Third-order moments: Cavity (< a^{2} a^{\dagger} >, < a^{\dagger 2} a >)
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cav3
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: k1_cav3, k2_cav3, k3_cav3, k4_cav3
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :, :), ALLOCATABLE    :: cav3_ss

! Fourth-order moments: Cavity and atom ( < a^{\dagger} a^{2} \sigma >, < a^{\dagger 2} a \sigma >)
! COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: cavsig4
! COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: k1_cavsig4, k2_cavsig4, k3_cavsig4, k4_cavsig4
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :, :, :), ALLOCATABLE :: cavsig4_ss

! Fourth-order moments: Cavity (< a^{\dagger 2} a^{2} >)
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cav4
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: k1_cav4, k2_cav4, k3_cav4, k4_cav4
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE :: cav4_ss

! Integer indices for sigma operators
INTEGER, PARAMETER                                        :: sm = 1, sp = 2, sz = 3
! Integer indices for: a, a^{\dagger}, a^{\dagger} a
INTEGER, PARAMETER                                        :: a = 1, at = 2, ata = 3

!----------------------------!
!     OTHER USEFUL STUFF     !
!----------------------------!
! Integer counter
INTEGER                                                :: j, k, l, m, x
! Integer cavity counters
INTEGER                                                :: cj, ck, cl, cm
! Imaginary i
COMPLEX(KIND=8), PARAMETER                             :: i = CMPLX(0.0d0, 1.0d0, 8)
! pi
REAL(KIND=8), PARAMETER                                :: pi = 3.1415926535897932384d0
! 1 / 6
REAL(KIND=8), PARAMETER                                :: xis = 1.0d0 / 6.0d0
! Atomic population
! REAL(KIND=8)                                           :: popg, pope, popf
! Photon number
! REAL(KIND=8)                                           :: photon
REAL(KIND=8)                                           :: photona_ss, photonb_ss
! Complex data
COMPLEX(KIND=8)                                        :: moment_out

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
CHARACTER(LEN=33), PARAMETER :: filename_parameters = "./data_files/cross_parameters.txt"
! Filename for state population
! CHARACTER(LEN=23), PARAMETER :: filename_states = "./data_files/states.txt"
! Filename for operators
! CHARACTER(LEN=26), PARAMETER :: filename_cavity = "./data_files/operators.txt"
! Filename for first-order correlation
! CHARACTER(LEN=24), PARAMETER :: filename_g1 = "./data_files/g1_corr.txt"
! Filename for second-order correlation
CHARACTER(LEN=30), PARAMETER :: filename_g2 = "./data_files/cross_g2_corr.txt"

!==============================================================================!
!                 NAMELIST AND PARAMETERS TO BE READ FROM FILE                 !
!==============================================================================!
! NameList things
! Status and unit integers
INTEGER :: ISTAT, IUNIT
! Line to be read from file
CHARACTER(LEN=512) :: LINE
! Namelist parameters
NAMELIST /ATOM/ Gamma, Omega
NAMELIST /FILTER/ epsilon, N, phase
NAMELIST /CAVITYA/ kappaa, w0a, dwa
NAMELIST /CAVITYB/ kappab, w0b, dwb
NAMELIST /TIME/ dt, t_max, tau1_max, tau2_max

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

READ(IUNIT, NML=CAVITYB, IOSTAT=ISTAT)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') LINE
  CLOSE(IUNIT)
  PRINT *, "Invalid line in CAVITYB namelist: " // TRIM(line)
  CALL EXIT(1)
END IF

READ(IUNIT, NML=TIME, IOSTAT=ISTAT)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') LINE
  CLOSE(IUNIT)
  PRINT *, "Invalid line in TIME namelist: " // TRIM(line)
  CALL EXIT(1)
END IF
CLOSE(IUNIT)

! Number of time-steps
t_steps = NINT(t_max / dt)

!==============================================================================!
!                DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS                !
!==============================================================================!
!------------------------!
!     BLOCH MATRIX M     !
!------------------------!
Mat_OG = 0.0d0
! Row 1: d/dt |g><e|
Mat_OG(1, 1) = -0.5d0 * gamma
Mat_OG(1, 2) = 0.0d0
Mat_OG(1, 3) = i * 0.5d0 * Omega
! Row 2: d/dt |e><g|
Mat_OG(2, 1) = 0.0d0
Mat_OG(2, 2) = -0.5d0 * gamma
Mat_OG(2, 3) = -i * 0.5d0 * Omega
! Row 3: d/dt |e><e| - |g><g|
Mat_OG(3, 1) = i * Omega
Mat_OG(3, 2) = -i * Omega
Mat_OG(3, 3) = -gamma

!--------------------------------!
!     NON-HOMOGENEOUS VECTOR     !
!--------------------------------!
B_OG = 0.0d0
B_OG(1) = 0.0d0
B_OG(2) = 0.0d0
B_OG(3) = -gamma

!---------------------------------------------!
!     RESONANCES (wj) AND COUPLINGS (E_j)     !
!---------------------------------------------!
! Kappa values
kappa = 0.0d0
kappa(1) = kappaa; kappa(2) = kappab
! Allocate array of omega and gka values
ALLOCATE(wl(-N:N, 2)); wl = 0.0d0
ALLOCATE(gkl(-N:N, 2)); gkl = 0.0d0
DO j = -N, N
  IF (N == 0) THEN
    ! Cavity a
    wl(j, 1) = w0a
    gkl(j, 1) = DSQRT(0.5d0 * epsilon * gamma * kappa(1))
    ! Cavity b
    wl(j, 2) = w0b
    gkl(j, 2) = DSQRT(0.5d0 * epsilon * gamma * kappa(2))
  ELSE
    ! Blackman window coefficient
    blackman = 1.0d0
    ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + ja) / (2.0d0 * DBLE(N))) + &
    !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + ja) / (2.0d0 * DBLE(N)))
    ! Cavity a
    wl(j, 1) = w0a + DBLE(j) * dwa
    gkl(j, 1) = DSQRT(0.5d0 * (epsilon / DBLE(2*N + 1)) * gamma * kappa(1)) * blackman * EXP(i * DBLE(phase) * DBLE(j) * pi / DBLE(N))
    ! Cavity a
    wl(j, 2) = w0b + DBLE(j) * dwb
    gkl(j, 2) = DSQRT(0.5d0 * (epsilon / DBLE(2*N + 1)) * gamma * kappa(2)) * blackman * EXP(i * DBLE(phase) * DBLE(j) * pi / DBLE(N))
  END IF
END DO

!------------------------------------------!
!     INITALISE OPERATOR MOMENT ARRAYS     !
!------------------------------------------!
! First-order: Cavity
ALLOCATE(cav1(-N:N, 2)); cav1 = 0.0d0
ALLOCATE(k1_cav1(-N:N, 2)); k1_cav1 = 0.0d0
ALLOCATE(k2_cav1(-N:N, 2)); k2_cav1 = 0.0d0
ALLOCATE(k3_cav1(-N:N, 2)); k3_cav1 = 0.0d0
ALLOCATE(k4_cav1(-N:N, 2)); k4_cav1 = 0.0d0
ALLOCATE(cav1_ss(-N:N, 2, 2)); cav1_ss = 0.0d0

! Second-order: Cavity and Atom
ALLOCATE(cavsig2(-N:N, 2, N_mat)); cavsig2 = 0.0d0
ALLOCATE(k1_cavsig2(-N:N, 2, N_mat)); k1_cavsig2 = 0.0d0
ALLOCATE(k2_cavsig2(-N:N, 2, N_mat)); k2_cavsig2 = 0.0d0
ALLOCATE(k3_cavsig2(-N:N, 2, N_mat)); k3_cavsig2 = 0.0d0
ALLOCATE(k4_cavsig2(-N:N, 2, N_mat)); k4_cavsig2 = 0.0d0
ALLOCATE(cavsig2_ss(-N:N, 2, N_mat, 2)); cavsig2_ss = 0.0d0

! Second-order: Cavity
ALLOCATE(cav2(-N:N, -N:N, 3)); cav2 = 0.0d0
ALLOCATE(k1_cav2(-N:N, -N:N, 3)); k1_cav2 = 0.0d0
ALLOCATE(k2_cav2(-N:N, -N:N, 3)); k2_cav2 = 0.0d0
ALLOCATE(k3_cav2(-N:N, -N:N, 3)); k3_cav2 = 0.0d0
ALLOCATE(k4_cav2(-N:N, -N:N, 3)); k4_cav2 = 0.0d0
ALLOCATE(cav2_ss(-N:N, -N:N, 3, 2, 2)); cav2_ss = 0.0d0

! Third-order: Cavity and Atom
! ALLOCATE(cavsig3(-N:N, -N:N, 3, N_mat)); cavsig3 = 0.0d0
! ALLOCATE(k1_cavsig3(-N:N, -N:N, 3, N_mat)); k1_cavsig3 = 0.0d0
! ALLOCATE(k2_cavsig3(-N:N, -N:N, 3, N_mat)); k2_cavsig3 = 0.0d0
! ALLOCATE(k3_cavsig3(-N:N, -N:N, 3, N_mat)); k3_cavsig3 = 0.0d0
! ALLOCATE(k4_cavsig3(-N:N, -N:N, 3, N_mat)); k4_cavsig3 = 0.0d0
ALLOCATE(cavsig3_ss(-N:N, -N:N, 3, N_mat, 2, 2)); cavsig3_ss = 0.0d0

! Third-order: Cavity
! ALLOCATE(cav3(-N:N, -N:N, -N:N, 2)); cav3 = 0.0d0
! ALLOCATE(k1_cav3(-N:N, -N:N, -N:N, 2)); k1_cav3 = 0.0d0
! ALLOCATE(k2_cav3(-N:N, -N:N, -N:N, 2)); k2_cav3 = 0.0d0
! ALLOCATE(k3_cav3(-N:N, -N:N, -N:N, 2)); k3_cav3 = 0.0d0
! ALLOCATE(k4_cav3(-N:N, -N:N, -N:N, 2)); k4_cav3 = 0.0d0
ALLOCATE(cav3_ss(-N:N, -N:N, -N:N, 2, 2, 2, 2)); cav3_ss = 0.0d0

! Fourth-order: Cavity and atom
! ALLOCATE(cavsig4(-N:N, -N:N, -N:N, 2, N_mat)); cavsig4 = 0.0d0
! ALLOCATE(k1_cavsig4(-N:N, -N:N, -N:N, 2, N_mat)); k1_cavsig4 = 0.0d0
! ALLOCATE(k2_cavsig4(-N:N, -N:N, -N:N, 2, N_mat)); k2_cavsig4 = 0.0d0
! ALLOCATE(k3_cavsig4(-N:N, -N:N, -N:N, 2, N_mat)); k3_cavsig4 = 0.0d0
! ALLOCATE(k4_cavsig4(-N:N, -N:N, -N:N, 2, N_mat)); k4_cavsig4 = 0.0d0
ALLOCATE(cavsig4_ss(-N:N, -N:N, -N:N, 2, N_mat, 2, 2, 2)); cavsig4_ss = 0.0d0

! Fourth-order: Cavity
! ALLOCATE(cav4(-N:N, -N:N, -N:N, -N:N)); cav4 = 0.0d0
! ALLOCATE(k1_cav4(-N:N, -N:N, -N:N, -N:N)); k1_cav4 = 0.0d0
! ALLOCATE(k2_cav4(-N:N, -N:N, -N:N, -N:N)); k2_cav4 = 0.0d0
! ALLOCATE(k3_cav4(-N:N, -N:N, -N:N, -N:N)); k3_cav4 = 0.0d0
! ALLOCATE(k4_cav4(-N:N, -N:N, -N:N, -N:N)); k4_cav4 = 0.0d0
ALLOCATE(cav4_ss(-N:N, -N:N, -N:N, -N:N)); cav4_ss = 0.0d0

!----------------------------!
!     INITIAL CONDITIONS     !
!----------------------------!
sigma = 0.0d0
! Atom in ground state, cjity in vacuum.
sigma(sz) = -1.0d0

!==============================================================================!
!                           WRITE PARAMETERS TO FILE                           !
!==============================================================================!

! Open file to write time to
OPEN(UNIT=1, FILE=filename_parameters, STATUS='REPLACE', ACTION='WRITE')
! Write parameter
WRITE(1,*) "Parameters are in the following order:"
WRITE(1,"(A10,F25.15)") "Gamma =", Gamma
WRITE(1,"(A10,F25.15)") "Omega =", Omega
WRITE(1,"(A10,F25.15)") "w0a =", w0a
WRITE(1,"(A10,F25.15)") "kappaa =", kappaa
WRITE(1,"(A11,F25.15)") "dwa = ", dwa
WRITE(1,"(A10,F25.15)") "w0b =", w0b
WRITE(1,"(A10,F25.15)") "kappab =", kappab
WRITE(1,"(A11,F25.15)") "dwb = ", dwb
WRITE(1,"(A10,F25.15)") "epsilon =", epsilon
WRITE(1,"(A11,I9)") "N = ", N
WRITE(1,"(A11,I9)") "phase = ", phase
! WRITE(1,"(A10,F25.15)") "dt =", dt
! WRITE(1,"(A10,F25.15)") "Max tau2 =", tau2_max
! Close file
CLOSE(1)

!==============================================================================!
!                        CALCULATE STEADY-STATE MOMENTS                        !
!==============================================================================!
!---------------------------!
!     FIRST-ORDER: ATOM     !
!---------------------------!
! Set matrix to be inverted
Mat_inv = Mat_OG

! Perform LU-factorization of matrix
CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
IF (INFO .NE. 0) THEN
  PRINT*, "zGETRF M failed :( INFO = ", INFO
  STOP
END IF

! Query optimal work space
! LWORK = -1
! CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)
! ! Set optimal work space and run again
! LWORK = MIN(LWMAX, INT(WORK(1)))

! Invert Matrix (Optimal LWORK = 3)
LWORK = 3
CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

! Calculate steady states
sigma_ss = 0.0d0
sigma_ss = -MATMUL(Mat_inv, B_OG)

! Cycle through cavity mode combinations
DO cj = 1, 2
  ! Cycle through modes
  DO j = -N, N
    !-----------------------------!
    !     FIRST-ORDER: CAVITY     !
    !-----------------------------!
    !-----------!
    ! < a_{j} > !
    !-----------!
    cav1_ss(j, a, cj) = -gkl(j, cj) * sigma_ss(sm) / &
                      & (kappa(cj) + i * wl(j, cj))
    !---------------------!
    ! < a^{\dagger}_{j} > !
    !---------------------!
    cav1_ss(j, at, cj) = -CONJG(gkl(j, cj)) * sigma_ss(sp) / &
                       & (kappa(cj) - i * wl(j, cj))
    !---------------------------------------!
    !     SECOND-ORDER: CAVITY AND ATOM     !
    !---------------------------------------!
    !------------------!
    ! < a_{j} \sigma > !
    !------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa(cj) + i * wl(j, cj))
    END DO

    ! Set the non-homogeneous vector
    B = 0.0d0
    B(1) = 0.0d0
    B(2) = -0.5d0 * gkl(j, cj) * (sigma_ss(sz) + 1.0d0)
    B(3) = -gamma * cav1_ss(j, a, cj) + &
         & gkl(j, cj) * sigma_ss(sm)

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

    ! Invert Matrix (Optimal LWORK = 3)
    LWORK = 3
    CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

    ! Calculate steady state
    cavsig2_ss(j, a, :, cj) = -MATMUL(Mat_inv, B)

    !----------------------------!
    ! < a^{\dagger}_{j} \sigma > !
    !----------------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa(cj) - i * wl(j, cj))
    END DO

    ! Set the non-homogeneous vector
    B = 0.0d0
    B(1) = -0.5d0 * CONJG(gkl(j, cj)) * (sigma_ss(sz) + 1.0d0)
    B(2) = 0.0d0
    B(3) = -gamma * cav1_ss(j, at, cj) + &
         & CONJG(gkl(j, cj)) * sigma_ss(sp)

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

    ! Invert Matrix (Optimal LWORK = 3)
    LWORK = 3
    CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

    ! Calculate steady state
    cavsig2_ss(j, at, :, cj) = -MATMUL(Mat_inv, B)

    ! Close j loop
  END DO
  ! Close cj loop
END DO

! Cycle through cavity mode combinations
DO ck = 1, 2
  DO cj = 1, 2
    ! Cycle through modes
    DO k = -N, N
      DO j = -N, N
        !------------------------------!
        !     SECOND-ORDER: CAVITY     !
        !------------------------------!
        !-----------------!
        ! < a_{j} a_{k} > !
        !-----------------!
        cav2_ss(j, k, a, cj, ck) = -gkl(j, cj) * cavsig2_ss(k, a, sm, ck) + &
                                 & -gkl(k, ck) * cavsig2_ss(j, a, sm, cj)
        cav2_ss(j, k, a, cj, ck) = cav2_ss(j, k, a, cj, ck) / &
                                 & ((kappa(cj) + kappa(ck)) + i * (wl(j, cj) + wl(k, ck)))

        !-------------------------------------!
        ! < a^{\dagger}_{j} a^{\dagger}_{k} > !
        !-------------------------------------!
        cav2_ss(j, k, at, cj, ck) = -CONJG(gkl(j, cj)) * cavsig2_ss(k, at, sp, ck) + &
                                  & -CONJG(gkl(k, ck)) * cavsig2_ss(j, at, sp, cj)
        cav2_ss(j, k, at, cj, ck) = cav2_ss(j, k, at, cj, ck) / &
                                  & ((kappa(cj) + kappa(ck)) - i * (wl(j, cj) + wl(k, ck)))

        !---------------------------!
        ! < a^{\dagger}_{j} a_{k} > !
        !---------------------------!
        cav2_ss(j, k, ata, cj, ck) = -CONJG(gkl(j, cj)) * cavsig2_ss(k, a, sp, ck) + &
                                   & -gkl(k, ck) * cavsig2_ss(j, at, sm, cj)
        cav2_ss(j, k, ata, cj, ck) = cav2_ss(j, k, ata, cj, ck) / &
                                   & ((kappa(cj) + kappa(ck)) - i * (wl(j, cj) - wl(k, ck)))

        !--------------------------------------!
        !     THIRD-ORDER: CAVITY AND ATOM     !
        !--------------------------------------!
        !------------------------!
        ! < a_{j} a_{k} \sigma > !
        !------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((kappa(cj) + kappa(ck)) + i * (wl(j, cj) + wl(k, ck)))
        END DO

        ! Set the non-homogeneous vector
        B = 0.0d0
        B(1) = 0.0d0
        B(2) = -0.5d0 * gkl(j, cj) * cavsig2_ss(k, a, sz, ck) + &
             & -0.5d0 * gkl(j, cj) * cav1_ss(k, a, ck) + &
             & -0.5d0 * gkl(k, ck) * cavsig2_ss(j, a, sz, cj) + &
             & -0.5d0 * gkl(k, ck) * cav1_ss(j, a, cj)
        B(3) = -gamma * cav2_ss(j, k, a, cj, ck) + &
             & gkl(j, cj) * cavsig2_ss(k, a, sm, ck) + &
             & gkl(k, ck) * cavsig2_ss(j, a, sm, cj)

        ! Set inverse matrix
        Mat_inv = Mat
        ! Perform LU-factorization of matrix
        CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
        IF (INFO .NE. 0) THEN
          WRITE(*, '(A27,I3,A3,I3,A15,I4)') "zGETRF M failed on cavsig3(", j , ", ", k, ", a) :( INFO = ", INFO
          STOP
        END IF

        ! ! Query optimal work space
        ! LWORK = -1
        ! CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)
        ! ! Set optimal work space and run again
        ! LWORK = MIN(LWMAX, INT(WORK(1)))

        ! Invert Matrix (Optimal LWORK = 3)
        LWORK = 3
        CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

        ! Calculate steady state
        cavsig3_ss(j, k, a, :, cj, ck) = -MATMUL(Mat_inv, B)

        !--------------------------------------------!
        ! < a^{\dagger}_{j} a^{\dagger}_{k} \sigma > !
        !--------------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((kappa(cj) + kappa(ck)) - i * (wl(j, cj) + wl(k, ck)))
        END DO

        ! Set the non-homogeneous vector
        B = 0.0d0
        B(1) = -0.5d0 * CONJG(gkl(j, cj)) * cavsig2_ss(k, at, sz, ck) + &
             & -0.5d0 * CONJG(gkl(j, cj)) * cav1_ss(k, at, ck) + &
             & -0.5d0 * CONJG(gkl(k, ck)) * cavsig2_ss(j, at, sz, cj) + &
             & -0.5d0 * CONJG(gkl(k, ck)) * cav1_ss(j, at, cj)
        B(2) = 0.0d0
        B(3) = -gamma * cav2_ss(j, k, at, cj, ck) + &
             & CONJG(gkl(j, cj)) * cavsig2_ss(k, at, sp, ck) + &
             & CONJG(gkl(k, ck)) * cavsig2_ss(j, at, sp, cj)

        ! Set inverse matrix
        Mat_inv = Mat
        ! Perform LU-factorization of matrix
        CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
        IF (INFO .NE. 0) THEN
          WRITE(*, '(A27,I3,A3,I3,A16,I4)') "zGETRF M failed on cavsig3(", j , ", ", k, ", at) :( INFO = ", INFO
          STOP
        END IF

        ! ! Query optimal work space
        ! LWORK = -1
        ! CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)
        ! ! Set optimal work space and run again
        ! LWORK = MIN(LWMAX, INT(WORK(1)))

        ! Invert Matrix (Optimal LWORK = 3)
        LWORK = 3
        CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

        ! Calculate steady state
        cavsig3_ss(j, k, at, :, cj, ck) = -MATMUL(Mat_inv, B)

        !----------------------------------!
        ! < a^{\dagger}_{j} a_{k} \sigma > !
        !----------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((kappa(cj) + kappa(ck)) - i * (wl(j, cj) - wl(k, ck)))
        END DO

        ! Set the non-homogeneous vector
        B = 0.0d0
        B(1) = -0.5d0 * CONJG(gkl(j, cj)) * cavsig2_ss(k, a, sz, ck) + &
             & -0.5d0 * CONJG(gkl(j, cj)) * cav1_ss(k, a, ck)
        B(2) = -0.5d0 * gkl(k, ck) * cavsig2_ss(j, at, sz, cj) + &
             & -0.5d0 * gkl(k, ck) * cav1_ss(j, at, cj)
        B(3) = -gamma * cav2_ss(j, k, ata, cj, ck) + &
             & CONJG(gkl(j, cj)) * cavsig2_ss(k, a, sp, ck) + &
             & gkl(k, ck) * cavsig2_ss(j, at, sm, cj)

        ! Set inverse matrix
        Mat_inv = Mat
        ! Perform LU-factorization of matrix
        CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
        IF (INFO .NE. 0) THEN
          WRITE(*, '(A27,I3,A3,I3,A17,I4)') "zGETRF M failed on cavsig3(", j , ", ", k, ", ata) :( INFO = ", INFO
          STOP
        END IF

        ! ! Query optimal work space
        ! LWORK = -1
        ! CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)
        ! ! Set optimal work space and run again
        ! LWORK = MIN(LWMAX, INT(WORK(1)))

        ! Invert Matrix (Optimal LWORK = 3)
        LWORK = 3
        CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

        ! Calculate steady state
        cavsig3_ss(j, k, ata, :, cj, ck) = -MATMUL(Mat_inv, B)

        ! Close j loop
      END DO
      ! Close k loop
    END DO
    ! Close cj loop
  END DO
  ! Close ck loop
END DO

! Cycle through cavity mode combinations
DO cl = 1, 2
  DO ck = 1, 2
    DO cj = 1, 2
      ! Cycle through modes
      DO l = -N, N
        DO k = -N, N
          DO j = -N, N
            !-----------------------------!
            !     THIRD-ORDER: CAVITY     !
            !-----------------------------!
            !---------------------------------!
            ! < a^{\dagger}_{j} a_{k} a_{l} > !
            !---------------------------------!
            cav3_ss(j, k, l, a, cj, ck, cl) = -CONJG(gkl(j, cj)) * cavsig3_ss(k, l, a, sp, ck, cl) + &
                                            & -gkl(k, ck) * cavsig3_ss(j, l, ata, sm, cj, cl) + &
                                            & -gkl(l, cl) * cavsig3_ss(j, k, ata, sm, cj, ck)
            cav3_ss(j, k, l, a, cj, ck, cl) = cav3_ss(j, k, l, a, cj, ck, cl) / &
                                            & ((kappa(cj) + kappa(ck) + kappa(cl)) - i * (wl(j, cj) - wl(k, ck) - wl(l, cl)))

            !-------------------------------------------!
            ! < a^{\dagger}_{j} a^{\dagger}_{k} a_{l} > !
            !-------------------------------------------!
            cav3_ss(j, k, l, at, cj, ck, cl) = -CONJG(gkl(j, cj)) * cavsig3_ss(k, l, ata, sp, ck, cl) + &
                                             & -CONJG(gkl(k, ck)) * cavsig3_ss(j, l, ata, sp, cj, cl) + &
                                             & -gkl(l, cl) * cavsig3_ss(j, k, at, sm, cj, ck)
            cav3_ss(j, k, l, at, cj, ck, cl) = cav3_ss(j, k, l, at, cj, ck, cl) / &
                                             & ((kappa(cj) + kappa(ck) + kappa(cl)) - i * (wl(j, cj) + wl(k, ck) - wl(l, cl)))

            !--------------------------------------!
            !     FOURTH-ORDER: CAVITY AND ATOM    !
            !--------------------------------------!
            !----------------------------------------!
            ! < a^{\dagger}_{j} a_{k} a_{l} \sigma > !
            !----------------------------------------!
            ! Set the diagonal matrix elements for M
            Mat = Mat_OG
            DO x = 1, N_mat
              Mat(x, x) = Mat(x, x) - ((kappa(cj) + kappa(ck) + kappa(cl)) - i * (wl(j, cj) - wl(k, ck) - wl(l, cl)))
            END DO

            ! Set the non-homogeneous vector
            B = 0.0d0
            B(1) = -0.5d0 * CONJG(gkl(j, cj)) * cavsig3_ss(k, l, a, sz, ck, cl) + &
                 & -0.5d0 * CONJG(gkl(j, cj)) * cav2_ss(k, l, a, ck, cl)
            B(2) = -0.5d0 * gkl(k, ck) * cavsig3_ss(j, l, ata, sz, cj, cl) + &
                 & -0.5d0 * gkl(k, ck) * cav2_ss(j, l, ata, cj, cl) + &
                 & -0.5d0 * gkl(l, cl) * cavsig3_ss(j, k, ata, sz, cj, ck) + &
                 & -0.5d0 * gkl(l, cl) * cav2_ss(j, k, ata, cj, ck)
            B(3) = -gamma * cav3_ss(j, k, l, a, cj, ck, cl) + &
                 & CONJG(gkl(j, cj)) * cavsig3_ss(k, l, a, sp, ck, cl) + &
                 & gkl(k, ck) * cavsig3_ss(j, l, ata, sm, cj, cl) + &
                 & gkl(l, cl) * cavsig3_ss(j, k, ata, sm, cj, ck)

            ! Set inverse matrix
            Mat_inv = Mat
            ! Perform LU-factorization of matrix
            CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
            IF (INFO .NE. 0) THEN
              WRITE(*, '(A27,I3,A3,I3,A3,I3,A15,I4)') "zGETRF M failed on cavsig4(", j , ", ", k, ", ", l, ", a) :( INFO = ", INFO
              STOP
            END IF

            ! ! Query optimal work space
            ! LWORK = -1
            ! CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)
            ! ! Set optimal work space and run again
            ! LWORK = MIN(LWMAX, INT(WORK(1)))

            ! Invert Matrix (Optimal LWORK = 3)
            LWORK = 3
            CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

            ! Calculate steady state
            cavsig4_ss(j, k, l, a, :, cj, ck, cl) = -MATMUL(Mat_inv, B)

            !-------------------------------------------!
            ! < a^{\dagger}_{j} a^{\dagger}_{k} a_{l} > !
            !-------------------------------------------!
            ! Set the diagonal matrix elements for M
            Mat = Mat_OG
            DO x = 1, N_mat
              Mat(x, x) = Mat(x, x) - ((kappa(cj) + kappa(ck) + kappa(cl)) - i * (wl(j, cj) + wl(k, ck) - wl(l, cl)))
            END DO

            ! Set the non-homogeneous vector
            B = 0.0d0
            B(1) = -0.5d0 * CONJG(gkl(j, cj)) * cavsig3_ss(k, l, ata, sz, ck, cl) + &
                 & -0.5d0 * CONJG(gkl(j, cj)) * cav2_ss(k, l, ata, ck, cl) + &
                 & -0.5d0 * CONJG(gkl(k, ck)) * cavsig3_ss(j, l, ata, sz, cj, cl) + &
                 & -0.5d0 * CONJG(gkl(k, ck)) * cav2_ss(j, l, ata, cj, cl)
            B(2) = -0.5d0 * gkl(l, cl) * cavsig3_ss(j, k, at, sz, cj, ck) + &
                 & -0.5d0 * gkl(l, cl) * cav2_ss(j, k, at, cj, ck)
            B(3) = -gamma * cav3_ss(j, k, l, at, cj, ck, cl) + &
                 & CONJG(gkl(j, cj)) * cavsig3_ss(k, l, ata, sp, ck, cl) + &
                 & CONJG(gkl(k, ck)) * cavsig3_ss(j, l, ata, sp, cj, cl) + &
                 & gkl(l, cl) * cavsig3_ss(j, k, at, sm, cj, ck)

            ! Set inverse matrix
            Mat_inv = Mat
            ! Perform LU-factorization of matrix
            CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
            IF (INFO .NE. 0) THEN
              WRITE(*, '(A27,I3,A3,I3,A3,I3,A16,I4)') "zGETRF M failed on cavsig4(", j , ", ", k, ", ", l, ", at) :( INFO = ", INFO
              STOP
            END IF

            ! ! Query optimal work space
            ! LWORK = -1
            ! CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)
            ! ! Set optimal work space and run again
            ! LWORK = MIN(LWMAX, INT(WORK(1)))

            ! Invert Matrix (Optimal LWORK = 3)
            LWORK = 3
            CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

            ! Calculate steady state
            cavsig4_ss(j, k, l, at, :, cj, ck, cl) = -MATMUL(Mat_inv, B)

            ! Close j loop
          END DO
          ! Close k loop
        END DO
        ! Close l loop
      END DO
      ! Close cj loop
    END DO
    ! Close ck loop
  END DO
  ! Close cl loop
END DO

! Cycle through cavity mode combinations
cj = 1; cm = 1
ck = 2; cl = 2
! Cycle through modes
DO m = -N, N
  DO l = -N, N
    DO k = -N, N
      DO j = -N, N
        !------------------------------!
        !     FOURTH-ORDER: CAVITY     !
        !------------------------------!
        !-------------------------------------------------!
        ! < a^{\dagger}_{j} a^{\dagger}_{k} a_{l} a_{m} > !
        !-------------------------------------------------!
        cav4_ss(j, k, l, m) = -CONJG(gkl(j, 1)) * cavsig4_ss(k, l, m, a, sp, 2, 2, 1) + &
                            & -CONJG(gkl(k, 2)) * cavsig4_ss(j, m, l, a, sp, 1, 1, 2) + &
                            & -gkl(l, 2) * cavsig4_ss(k, j, m, at, sm, 2, 1, 1) + &
                            & -gkl(m, 1) * cavsig4_ss(j, k, l, at, sm, 1, 2, 2)
        cav4_ss(j, k, l, m) = cav4_ss(j, k, l, m) / &
                            & (((2.0d0 * kappa(1)) + (2.0d0 * kappa(2))) - i * (wl(j, 1) + wl(k, 2)) + i * (wl(l, 2) + wl(m, 1)))
        ! Close j loop
      END DO
      ! Close k loop
    END DO
    ! Close l loop
  END DO
  ! Close m loop
END DO

! moment_out = cavsig4(1, 2)
! WRITE(*, '(A26,ES18.11E2,A3,ES18.11E2,A1)') " < at_0 a_0 a_0 sm >_ss = ", REAL(moment_out), " + ", AIMAG(moment_out), "i"
! moment_out =
! WRITE(*, '(A26,ES18.11E2,A3,ES18.11E2,A1)') "< at_0 at_0 a_0 sm >_ss = ", REAL(moment_out), " + ", AIMAG(moment_out), "i"

! Calculate photon number
photona_ss = 0.0d0
photonb_ss = 0.0d0
DO k = -N, N
  DO j = -N, N
    photona_ss = photona_ss + cav2_ss(j, k, ata, 1, 1)
    photonb_ss = photonb_ss + cav2_ss(j, k, ata, 2, 2)
  END DO
END DO

!==============================================================================!
!                  CALCULATE FIRST-ORDER CORRELATION FUNCTION                  !
!==============================================================================!
! Set initial conditions and non-homogeneous vector
! < a^{\dagger}_{j}(0) \sigma(\tau = 0) a_{m}(0) > =
!                          < a^{\dagger}_{j} a_{m} \sigma >_{ss},
! < a^{\dagger}_{j}(0) b^{\dagger}_{k}(\tau = 0) a_{m}(0) > =
!                          < a^{\dagger}_{j} b^{\dagger}_{k} a_{m} >_{ss},
! < a^{\dagger}_{j}(0) b_{l}(\tau = 0) a_{m}(0) > =
!                          < a^{\dagger}_{j} b_{l} a_{m} >_{ss},
! < a^{\dagger}_{j}(0) b^{\dagger}_{k} \sigma(\tau = 0) a_{m}(0) =
!                         < a^{\dagger}_{j} b^{\dagger}_{k} a_{m} \sigma >_{ss},
! < a^{\dagger}_{j}(0) b_{l} \sigma(\tau = 0) a_{m}(0) =
!                         < a^{\dagger}_{j} b_{l} a_{m} \sigma >_{ss},
! and
! < a^{\dagger}_{j}(0) b^{\dagger}_{k} b_{l}(\tau = 0)  a_{m}(0) > =
!                          < a^{\dagger}_{j} b^{\dagger}_{k} b_{l} a_{m} >_{ss}.

sigma = 0.0d0
cav1 = 0.0d0
cavsig2 = 0.0d0
cav2 = 0.0d0
B_OG = 0.0d0
! Set modes to calculate initial conditions for
cj = 1
ck = 2
cl = 2
cm = 1
! Cycle through modes
DO m = -N, N
  DO j = -N, N
    ! First-order: Atom
    sigma(:) = sigma(:) + cavsig3_ss(j, m, ata, :, cj, cm)

    DO k = -N, N
      ! First-order: Cavity
      cav1(k, a) = cav1(k, a) + cav3_ss(j, k, m, a, cj, ck, cm)
      cav1(k, at) = cav1(k, at) + cav3_ss(j, k, m, at, cj, ck, cm)
      ! Second-order: Cavity and atom
      cavsig2(k, a, :) = cavsig2(k, a, :) + cavsig4_ss(j, k, m, a, :, cj, ck, cm)
      cavsig2(k, at, :) = cavsig2(k, at, :) + cavsig4_ss(j, k, m, at, :, cj, ck, cm)

      DO l = -N, N
        ! Second-order: cavity
        cav2(k, l, ata) = cav2(k, l, ata) + cav4_ss(j, k, l, m)

        ! Close l loop
      END DO
      ! Close k loop
    END DO
    ! Non homogeneous vector
    B_OG(3) = B_OG(3) - gamma * cav2_ss(j, m, ata, cj, cm)

    ! Close j loop
  END DO
  ! Close m loop
END DO

! Set tau_steps
tau_steps = NINT(tau2_max / DBLE(dt))

! Open file to write time and data to
OPEN(UNIT=4, FILE=filename_g2, STATUS='REPLACE', ACTION='WRITE', RECL=4000)

! Set modes to calculate correlation for
cj = 2

! Cycle through time steps
DO t = 0, tau_steps
  !============================================================================!
  !                          CALCULATE AND WRITE DATA                          !
  !============================================================================!
  !-------------------------------!
  !     CALCULATE DATA STATES     !
  !-------------------------------!
  ! Grab correlation value
  moment_out = 0.0d0
  DO k = -N, N
    DO j = -N, N
      moment_out = moment_out + cav2(j, k, ata)
    END DO
  END DO

  ! Normalise correlation by steady-state photon number
  IF (photona_ss .NE. 0.0d0 .AND. photonb_ss .NE. 0.0d0) THEN
    moment_out = moment_out / (photona_ss * photonb_ss)
  END IF

  !-----------------------!
  !     WRITE TO FILE     !
  !-----------------------!
  ! First-order correlation function
  WRITE(4, *) DBLE(t) * dt, REAL(moment_out)

  !============================================================================!
  !                  CALCULATE USING FOURTH-ORDER RUNGE-KUTTA                  !
  !============================================================================!
  !----------------------------------------!
  !     INITIALISE RUNGE-KUTTA VECTORS     !
  !----------------------------------------!
  k1_sigma = 0.0d0; k2_sigma = 0.0d0; k3_sigma = 0.0d0; k4_sigma = 0.0d0
  k1_cav1 = 0.0d0; k2_cav1 = 0.0d0; k3_cav1 = 0.0d0; k4_cav1 = 0.0d0
  ! Second-order
  k1_cavsig2 = 0.0d0; k2_cavsig2 = 0.0d0; k3_cavsig2 = 0.0d0; k4_cavsig2 = 0.0d0
  k1_cav2 = 0.0d0; k2_cav2 = 0.d0; k3_cav2 = 0.0d0; k4_cav2 = 0.0d0

  !---------------------------!
  !     FIRST-ORDER: ATOM     !
  !---------------------------!
  ! Calculate Runge-Kutta vectors
  k1_sigma = dt * (MATMUL(Mat_OG, sigma) + B_OG)
  k2_sigma = dt * (MATMUL(Mat_OG, (sigma + 0.5d0 * k1_sigma)) + B_OG)
  k3_sigma = dt * (MATMUL(Mat_OG, (sigma + 0.5d0 * k2_sigma)) + B_OG)
  k4_sigma = dt * (MATMUL(Mat_OG, (sigma + k3_sigma)) + B_OG)

  ! Cycle through modes
  DO j = -N, N
    !-----------------------------!
    !     FIRST-ORDER: CAVITY     !
    !-----------------------------!
    !-----------!
    ! < a_{j} > !
    !-----------!
    k1_cav1(j, a) = -dt * (kappa(cj) + i * wl(j, cj)) * cav1(j, a) + &
                  & -dt * gkl(j, cj) * sigma(sm)
    k2_cav1(j, a) = -dt * (kappa(cj) + i * wl(j, cj)) * (cav1(j, a) + 0.5d0 * k1_cav1(j, a)) + &
                  & -dt * gkl(j, cj) * (sigma(sm) + 0.5d0 * k1_sigma(sm))
    k3_cav1(j, a) = -dt * (kappa(cj) + i * wl(j, cj)) * (cav1(j, a) + 0.5d0 * k2_cav1(j, a)) + &
                  & -dt * gkl(j, cj) * (sigma(sm) + 0.5d0 * k2_sigma(sm))
    k4_cav1(j, a) = -dt * (kappa(cj) + i * wl(j, cj)) * (cav1(j, a) + k3_cav1(j, a)) + &
                  & -dt * gkl(j, cj) * (sigma(sm) + k3_sigma(sm))

    !---------------------!
    ! < a^{\dagger}_{j} > !
    !---------------------!
    k1_cav1(j, at) = -dt * (kappa(cj) - i * wl(j, cj)) * cav1(j, at) + &
                   & -dt * CONJG(gkl(j, cj)) * sigma(sp)
    k2_cav1(j, at) = -dt * (kappa(cj) - i * wl(j, cj)) * (cav1(j, at) + 0.5d0 * k1_cav1(j, at)) + &
                   & -dt * CONJG(gkl(j, cj)) * (sigma(sp) + 0.5d0 * k1_sigma(sp))
    k3_cav1(j, at) = -dt * (kappa(cj) - i * wl(j, cj)) * (cav1(j, at) + 0.5d0 * k2_cav1(j, at)) + &
                   & -dt * CONJG(gkl(j, cj)) * (sigma(sp) + 0.5d0 * k2_sigma(sp))
    k4_cav1(j, at) = -dt * (kappa(cj) - i * wl(j, cj)) * (cav1(j, at) + k3_cav1(j, at)) + &
                   & -dt * CONJG(gkl(j, cj)) * (sigma(sp) + k3_sigma(sp))

    !---------------------------------------!
    !     SECOND-ORDER: CAVITY AND ATOM     !
    !---------------------------------------!
    ! Using matrix multiplication, we add to the Lindblad matrix M_OG and the
    ! non-homogeneous vector, then compute the new moments

    !------------------!
    ! < a_{j} \sigma > !
    !------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa(cj) + i * wl(j, cj))
    END DO

    ! Set the non-homogeneous vector
    B = 0.0d0
    B(1) = 0.0d0
    B(2) = -0.5d0 * gkl(j, cj) * (sigma(sz) + photonb_ss)
    B(3) = -gamma * cav1(j, a) + &
         & gkl(j, cj) * sigma(sm)
    ! Calculate k1
    k1_cavsig2(j, a, :) = dt * (MATMUL(Mat, cavsig2(j, a, :)) + B)

    ! Set the non-homogeneous vector
    B = 0.0d0
    B(1) = 0.0d0
    B(2) = -0.5d0 * gkl(j, cj) * ((sigma(sz) + 0.5d0 * k1_sigma(sz)) + photonb_ss)
    B(3) = -gamma * (cav1(j, a) + 0.5d0 * k1_cav1(j, a)) + &
         & gkl(j, cj) * (sigma(sm) + 0.5d0 * k1_sigma(sm))
    ! Calculate k2
    k2_cavsig2(j, a, :) = dt * (MATMUL(Mat, (cavsig2(j, a, :) + 0.5d0 * k1_cavsig2(j, a, :))) + B)

    ! Set the non-homogeneous vector
    B = 0.0d0
    B(1) = 0.0d0
    B(2) = -0.5d0 * gkl(j, cj) * ((sigma(sz) + 0.5d0 * k2_sigma(sz)) + photonb_ss)
    B(3) = -gamma * (cav1(j, a) + 0.5d0 * k2_cav1(j, a)) + &
         & gkl(j, cj) * (sigma(sm) + 0.5d0 * k2_sigma(sm))
    ! Calculate k3
    k3_cavsig2(j, a, :) = dt * (MATMUL(Mat, (cavsig2(j, a, :) + 0.5d0 * k2_cavsig2(j, a, :))) + B)

    ! Set the non-homogeneous vector
    B = 0.0d0
    B(1) = 0.0d0
    B(2) = -0.5d0 * gkl(j, cj) * ((sigma(sz) + k3_sigma(sz)) + photonb_ss)
    B(3) = -gamma * (cav1(j, a) + k3_cav1(j, a)) + &
         & gkl(j, cj) * (sigma(sm) + k3_sigma(sm))
    ! Calculate k4
    k4_cavsig2(j, a, :) = dt * (MATMUL(Mat, (cavsig2(j, a, :) + k3_cavsig2(j, a, :))) + B)

    !----------------------------!
    ! < a^{\dagger}_{j} \sigma > !
    !----------------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa(cj) - i * wl(j, cj))
    END DO

    ! Set the non-homogeneous vector
    B = 0.0d0
    B(1) = -0.5d0 * CONJG(gkl(j, cj)) * (sigma(sz) + photonb_ss)
    B(2) = 0.0d0
    B(3) = -gamma * cav1(j, at) + &
         & CONJG(gkl(j, cj)) * sigma(sp)
    ! Calculate k1
    k1_cavsig2(j, at, :) = dt * (MATMUL(Mat, cavsig2(j, at, :)) + B)

    ! Set the non-homogeneous vector
    B = 0.0d0
    B(1) = -0.5d0 * CONJG(gkl(j, cj)) * ((sigma(sz) + 0.5d0 * k1_sigma(sz)) + photonb_ss)
    B(2) = 0.0d0
    B(3) = -gamma * (cav1(j, at) + 0.5d0 * k1_cav1(j, at)) + &
         & CONJG(gkl(j, cj)) * (sigma(sp) + 0.5d0 * k1_sigma(sp))
    ! Calculate k2
    k2_cavsig2(j, at, :) = dt * (MATMUL(Mat, (cavsig2(j, at, :) + 0.5d0 * k1_cavsig2(j, at, :))) + B)

    ! Set the non-homogeneous vector
    B = 0.0d0
    B(1) = -0.5d0 * CONJG(gkl(j, cj)) * ((sigma(sz) + 0.5d0 * k2_sigma(sz)) + photonb_ss)
    B(2) = 0.0d0
    B(3) = -gamma * (cav1(j, at) + 0.5d0 * k2_cav1(j, at)) + &
         & CONJG(gkl(j, cj)) * (sigma(sp) + 0.5d0 * k2_sigma(sp))
    ! Calculate k3
    k3_cavsig2(j, at, :) = dt * (MATMUL(Mat, (cavsig2(j, at, :) + 0.5d0 * k2_cavsig2(j, at, :))) + B)

    ! Set the non-homogeneous vector
    B = 0.0d0
    B(1) = -0.5d0 * CONJG(gkl(j, cj)) * ((sigma(sz) + k3_sigma(sz)) + photonb_ss)
    B(2) = 0.0d0
    B(3) = -gamma * (cav1(j, at) + k3_cav1(j, at)) + &
         & CONJG(gkl(j, cj)) * (sigma(sp) + k3_sigma(sp))
    ! Calculate k4
    k4_cavsig2(j, at, :) = dt * (MATMUL(Mat, (cavsig2(j, at, :) + k3_cavsig2(j, at, :))) + B)

    ! Close j loop
  END DO

  ! Cycle through modes
  DO j = -N, N
    DO k = -N, N
      !------------------------------!
      !     SECOND-ORDER: CAVITY     !
      !------------------------------!
      !---------------------------!
      ! < a^{\dagger}_{j} a_{k} > !
      !---------------------------!
      k1_cav2(j, k, ata) = -dt * (2.0d0 * kappa(cj) - i * (wl(j, cj) - wl(k, cj))) * cav2(j, k, ata) + &
                         & -dt * CONJG(gkl(j, cj)) * cavsig2(k, a, sp) + &
                         & -dt * gkl(k, cj) * cavsig2(j, at, sm)
      k2_cav2(j, k, ata) = -dt * (2.0d0 * kappa(cj) - i * (wl(j, cj) - wl(k, cj))) * (cav2(j, k, ata) + 0.5d0 * k1_cav2(j, k, ata)) + &
                         & -dt * CONJG(gkl(j, cj)) * (cavsig2(k, a, sp) + 0.5d0 * k1_cavsig2(k, a, sp)) + &
                         & -dt * gkl(k, cj) * (cavsig2(j, at, sm) + 0.5d0 * k1_cavsig2(j, at, sm))
      k3_cav2(j, k, ata) = -dt * (2.0d0 * kappa(cj) - i * (wl(j, cj) - wl(k, cj))) * (cav2(j, k, ata) + 0.5d0 * k2_cav2(j, k, ata)) + &
                         & -dt * CONJG(gkl(j, cj)) * (cavsig2(k, a, sp) + 0.5d0 * k2_cavsig2(k, a, sp)) + &
                         & -dt * gkl(k, cj) * (cavsig2(j, at, sm) + 0.5d0 * k2_cavsig2(j, at, sm))
      k4_cav2(j, k, ata) = -dt * (2.0d0 * kappa(cj) - i * (wl(j, cj) - wl(k, cj))) * (cav2(j, k, ata) + k3_cav2(j, k, ata)) + &
                         & -dt * CONJG(gkl(j, cj)) * (cavsig2(k, a, sp) + k3_cavsig2(k, a, sp)) + &
                         & -dt * gkl(k, cj) * (cavsig2(j, at, sm) + k3_cavsig2(j, at, sm))

      ! Close k loop
    END DO
    ! Close j loop
  END DO


  !============================================================================!
  !                   UPDATE ARRAYS FROM RUNGE-KUTTA ARRAYS                    !
  !============================================================================!
  ! First-order
  sigma = sigma + xis * (k1_sigma + 2.0d0 * (k2_sigma + k3_sigma) + k4_sigma)
  cav1 = cav1 + xis * (k1_cav1 + 2.0d0 * (k2_cav1 + k3_cav1) + k4_cav1)
  ! Second-order
  cavsig2 = cavsig2 + xis * (k1_cavsig2 + 2.0d0 * (k2_cavsig2 + k3_cavsig2) + k4_cavsig2)
  cav2 = cav2 + xis * (k1_cav2 + 2.0d0 * (k2_cav2 + k3_cav2) + k4_cav2)

  ! Close t loop
END DO

! Close file
CLOSE(4)

!==============================================================================!
!                                END OF PROGRAM                                !
!==============================================================================!

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

END PROGRAM TWO_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_CROSS_G2
