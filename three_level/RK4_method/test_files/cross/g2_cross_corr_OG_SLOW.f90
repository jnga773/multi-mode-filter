! The system in this program is a resonantly driven three-level ladder-type atom
! that is coupled into a multi-mode array of filter cavities, as described in
! some notes somewhere (or with this program if this goes to anyone).

! The operator moment steady states are calculated using the inverse matrix/
! analytical method. Using quantum regression theorem, we then use Runge-Kutta
! fourth-order to numerically integrate the moment equations, using the steady
! states and initial conditions, to calculate the normalised second-order
! cross correlation function between filters A and B.

! The input parameters are taken from a NameList file [filename_ParamList] which,
! by default, points to "./ParamList.nml". The code can thus be compiled once,
! and parameters can be changed in the NameList file for subsequent runs.

! The parameters for each run are written to [filename_parameters] which, by
! default, points to "./data_files/parameters.txt".

! The normalised first-order correlation function is written to [filename_g1]
! which, by default, points to "./data_files/g1_corr.txt". The file has two
! columns:
!                            t     REAL(g2_cross)

! For the default filenames, the folder "./data_files/" and NameList file
! "./ParamList.nml" MUST EXIST IN THE WORKING DIRECTORY.

! To compiled the code, I use the Intel Parallel Studio compiler IFORT with the
! command
!       (LINUX): ifort -O3 -o g2 -mkl g2_RK4.f90
!     (WINDOWS): ifort /O3 /o g2 /Qmkl g2_RK4.f90
! where the -O3 (/O3) flag gives maximum optimisation, the -o (/o) g1 flag
! names the executable as "g2" ("g2.exe"), and the -mkl (/Qmkl) flag links
! the program to Intel's Math Kernel Library, to make use of the LAPACK routines.

PROGRAM THREE_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_G2_CROSS

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
! Percentage of fluorecence aimed at cavity
REAL(KIND=8)                                           :: epsilon
! Number of mode either side of w0, 2N + 1 total mode
INTEGER                                                :: N
! Phase modulation of mode coupling
INTEGER                                                :: phase
! Central mode frequency of the filter cavity, with N mode frequencies either side
REAL(KIND=8)                                           :: w0a
REAL(KIND=8)                                           :: w0b
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                           :: kappaa
REAL(KIND=8)                                           :: kappab
! Frequency spacing of modes
REAL(KIND=8)                                           :: dwa
REAL(KIND=8)                                           :: dwb
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
INTEGER, PARAMETER                                     :: N_mat = 8
! M matrix (filled as transpose)
COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)               :: Mat, Mat_OG, Mat_inv
! Non-homogeneous vector
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: B_vec, B_OG

! First-order moments: Atomic equations (< \sigma >)
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: sigma
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: k1_sigma, k2_sigma, k3_sigma, k4_sigma
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: sigma_ss

! First-order moments: Cavity (< a >, < a^{\dagger} >)
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: cav1
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: k1_cav1, k2_cav1, k3_cav1, k4_cav1
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cav1_ss

! Second-order moments: Cavity and atom (< a \sigma >, < a^{\dagger} \sigma >
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cavsig2
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: k1_cavsig2, k2_cavsig2, k3_cavsig2, k4_cavsig2
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cavsig2_ss

! Second-order moments: Cavity (< a^{\dagger} a >)
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cav2
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: k1_cav2, k2_cav2, k3_cav2, k4_cav2
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: cav2_ss

! Third-order moments: Cavity and atom (< a^{2} \sigma >, < a^{\dagger 2} \sigma >, < a^{\dagger} a \sigma >)
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cavsig3
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: k1_cavsig3, k2_cavsig3, k3_cavsig3, k4_cavsig3
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: cavsig3_ss

! Third-order moments: Cavity (< a^{2} a^{\dagger} >, < a^{\dagger 2} a >)
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cav3
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: k1_cav3, k2_cav3, k3_cav3, k4_cav3
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :, :), ALLOCATABLE :: cav3_ss

! Fourth-order moments: Cavity and atom ( < a^{\dagger} a^{2} \sigma >, < a^{\dagger 2} a \sigma >)
! COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: cavsig4
! COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: k1_cavsig4, k2_cavsig4, k3_cavsig4, k4_cavsig4
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :, :, :), ALLOCATABLE :: cavsig4_ss

! Fourth-order moments: Cavity (< a^{\dagger 2} a^{2} >)
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cav4
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: k1_cav4, k2_cav4, k3_cav4, k4_cav4
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :, :, :), ALLOCATABLE :: cav4_ss


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
NAMELIST /ATOM/ Gamma, Omega, alpha, delta, xi
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
ALLOCATE(cav4_ss(-N:N, -N:N, -N:N, -N:N, 2, 2, 2, 2)); cav4_ss = 0.0d0

!----------------------------!
!     INITIAL CONDITIONS     !
!----------------------------!
! sigma = 0.0d0
! ! Atom in ground state, cavity in vacuum.
! sigma(sz) = -1.0d0

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
! WRITE(1,"(A10,F25.15)") "Max t =", t_max
! WRITE(1,"(A10,F25.15)") "Max tau1 =", tau1_max
! WRITE(1,"(A10,F25.15)") "Max tau2 =", tau2_max
! Close file
CLOSE(1)

!==============================================================================!
!                        CALCULATE STEADY-STATE MOMENTS                        !
!==============================================================================!
!---------------------------!
!     FIRST-ORDER: ATOM     !
!---------------------------!
IF (xi .NE. 0.0d0) THEN
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
  ! Check convergence
  IF (INFO .GT. 0) THEN
     PRINT*, "zGEEV failed ON eigenvalues/vectors of Mat_OG"
     STOP
  END IF

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

! Cycle through cavity operators
DO cj = 1, 2
  ! Cycle through modes
  DO j = -N, N
    !-----------------------------!
    !     FIRST-ORDER: CAVITY     !
    !-----------------------------!
    !-----------!
    ! < a_{j} > !
    !-----------!
    cav1_ss(j, a, cj) = -gkl(j, cj) * sigma_ss(ge) + &
                      & -gkl(j, cj) * xi * sigma_ss(ef)
    cav1_ss(j, a, cj) = cav1_ss(j, a, cj) / &
                      & (kappa(cj) + i * wl(j, cj))

    !---------------------!
    ! < a^{\dagger}_{j} > !
    !---------------------!
    cav1_ss(j, at, cj) = -CONJG(gkl(j, cj)) * sigma_ss(eg) + &
                       & -CONJG(gkl(j, cj)) * xi * sigma_ss(fe)
    cav1_ss(j, at, cj) = cav1_ss(j, at, cj) / &
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
    B_vec = 0.0d0
    B_vec(1) = -gkl(j, cj) * sigma_ss(ge)
    B_vec(2) = -gkl(j, cj) * xi * sigma_ss(gf)
    B_vec(3) = -gkl(j, cj) * sigma_ss(ee)
    B_vec(4) = Gamma * (xi ** 2) * cav1_ss(j, a, cj) + &
             & -gkl(j, cj) * xi * sigma_ss(ef)
    B_vec(5) = i * xi * 0.5d0 * Omega * cav1_ss(j, a, cj)
    B_vec(6) = -i * xi * 0.5d0 * Omega * cav1_ss(j, a, cj) + &
             & gkl(j, cj) * xi * sigma_ss(gg) + &
             & gkl(j, cj) * xi * sigma_ss(ee) + &
             & -gkl(j, cj) * xi
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j, cj) * sigma_ss(fe)

    ! Set inverse matrix
    Mat_inv = Mat
    ! Perform LU-factorization of matrix
    CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
    ! Invert Matrix (Optimal LWORK = 8)
    LWORK = 8
    CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

    ! Calculate steady state
    cavsig2_ss(j, a, :, cj) = -MATMUL(Mat_inv, B_vec)

    !----------------------------!
    ! < a^{\dagger}_{j} \sigma > !
    !----------------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa(cj) - i * wl(j, cj))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -CONJG(gkl(j, cj)) * sigma_ss(eg)
    B_vec(2) = -CONJG(gkl(j, cj)) * sigma_ss(ee)
    B_vec(3) = -CONJG(gkl(j, cj)) * xi * sigma_ss(fg)
    B_vec(4) = Gamma * (xi ** 2) * cav1_ss(j, at, cj) + &
             & -CONJG(gkl(j, cj)) * xi * sigma_ss(fe)
    B_vec(5) = i * xi * 0.5d0 * Omega * cav1_ss(j, at, cj) + &
             & CONJG(gkl(j, cj)) * xi * sigma_ss(gg) + &
             & CONJG(gkl(j, cj)) * xi * sigma_ss(ee) + &
             & -CONJG(gkl(j, cj)) * xi
    B_vec(6) = -i * xi * 0.5d0 * Omega * cav1_ss(j, at, cj)
    B_vec(7) = -CONJG(gkl(j, cj)) * sigma_ss(ef)
    B_vec(8) = 0.0d0

    ! Set inverse matrix
    Mat_inv = Mat
    ! Perform LU-factorization of matrix
    CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
    ! Invert Matrix (Optimal LWORK = 8)
    LWORK = 8
    CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

    ! Calculate steady state
    cavsig2_ss(j, at, :, cj) = -MATMUL(Mat_inv, B_vec)

    ! Close j loop
  END DO
  ! Close cj loop
END DO

! Cycle through cavity operators
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
        cav2_ss(j, k, a, cj, ck) = -gkl(j, cj) * cavsig2_ss(k, a, ge, ck) + &
                                 & -gkl(j, cj) * xi * cavsig2_ss(k, a, ef, ck) + &
                                 & -gkl(k, ck) * cavsig2_ss(j, a, ge, cj) + &
                                 & -gkl(k, ck) * xi * cavsig2_ss(j, a, ef, cj)
        cav2_ss(j, k, a, cj, ck) = cav2_ss(j, k, a, cj, ck) / &
                                 & ((kappa(cj) + kappa(ck)) + i * (wl(j, cj) + wl(k, ck)))

        !-------------------------------------!
        ! < a^{\dagger}_{j} a^{\dagger}_{k} > !
        !-------------------------------------!
        cav2_ss(j, k, at, cj, ck) = -CONJG(gkl(j, cj)) * cavsig2_ss(k, at, eg, ck) + &
                                  & -CONJG(gkl(j, cj)) * xi * cavsig2_ss(k, at, fe, ck) + &
                                  & -CONJG(gkl(k, ck)) * cavsig2_ss(j, at, eg, cj) + &
                                  & -CONJG(gkl(k, ck)) * xi * cavsig2_ss(j, at, fe, cj)
        cav2_ss(j, k, at, cj, ck) = cav2_ss(j, k, at, cj, ck) / &
                                  & ((kappa(cj) + kappa(ck)) - i * (wl(j, cj) + wl(k, ck)))

        !---------------------------!
        ! < a^{\dagger}_{j} a_{k} > !
        !---------------------------!
        cav2_ss(j, k, ata, cj, ck) = -CONJG(gkl(j, cj)) * cavsig2_ss(k, a, eg, ck) + &
                                   & -CONJG(gkl(j, cj)) * xi * cavsig2_ss(k, a, fe, ck) + &
                                   & -gkl(k, ck) * cavsig2_ss(j, at, ge, cj) + &
                                   & -gkl(k, ck) * xi * cavsig2_ss(j, at, ef, cj)
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
          Mat(x, x) = Mat(x, x) - (((kappa(cj) + kappa(ck))) + i * (wl(j, cj) + wl(k, ck)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -gkl(j, cj) * cavsig2_ss(k, a, ge, ck) + &
                 & -gkl(k, ck) * cavsig2_ss(j, a, ge, cj)
        B_vec(2) = -gkl(j, cj) * xi * cavsig2_ss(k, a, gf, ck) + &
                 & -gkl(k, ck) * xi * cavsig2_ss(j, a, gf, cj)
        B_vec(3) = -gkl(j, cj) * cavsig2_ss(k, a, ee, ck) + &
                 & -gkl(k, ck) * cavsig2_ss(j, a, ee, cj)
        B_vec(4) = Gamma * (xi ** 2) * cav2_ss(j, k, a, cj, ck) + &
                 & -gkl(j, cj) * xi * cavsig2_ss(k, a, ef, ck) + &
                 & -gkl(k, ck) * xi * cavsig2_ss(j, a, ef, cj)
        B_vec(5) = i * xi * 0.5d0 * Omega * cav2_ss(j, k, a, cj, ck)
        B_vec(6) = -i * xi * 0.5d0 * Omega * cav2_ss(j, k, a, cj, ck) + &
                 & gkl(j, cj) * xi * cavsig2_ss(k, a, gg, ck) + &
                 & gkl(j, cj) * xi * cavsig2_ss(k, a, ee, ck) + &
                 & -gkl(j, cj) * xi * cav1_ss(k, a, ck) + &
                 & gkl(k, ck) * xi * cavsig2_ss(j, a, gg, cj) + &
                 & gkl(k, ck) * xi * cavsig2_ss(j, a, ee, cj) + &
                 & -gkl(k, ck) * xi * cav1_ss(j, a, cj)
        B_vec(7) = 0.0d0
        B_vec(8) = -gkl(j, cj) * cavsig2_ss(k, a, fe, ck) + &
                 & -gkl(k, ck) * cavsig2_ss(j, a, fe, cj)

        ! Set inverse matrix
        Mat_inv = Mat
        ! Perform LU-factorization of matrix
        CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
        ! Invert Matrix (Optimal LWORK = 8)
        LWORK = 8
        CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

        ! Calculate steady state
        cavsig3_ss(j, k, a, :, cj, ck) = -MATMUL(Mat_inv, B_vec)

        !--------------------------------------------!
        ! < a^{\dagger}_{j} a^{\dagger}_{k} \sigma > !
        !--------------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - (((kappa(cj) + kappa(ck))) - i * (wl(j, cj) + wl(k, ck)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j, cj)) * cavsig2_ss(k, at, eg, ck) + &
                 & -CONJG(gkl(k, ck)) * cavsig2_ss(j, at, eg, cj)
        B_vec(2) = -CONJG(gkl(j, cj)) * cavsig2_ss(k, at, ee, ck) + &
                 & -CONJG(gkl(k, ck)) * cavsig2_ss(j, at, ee, cj)
        B_vec(3) = -CONJG(gkl(j, cj)) * xi * cavsig2_ss(k, at, fg, ck) + &
                 & -CONJG(gkl(k, ck)) * xi * cavsig2_ss(j, at, fg, cj)
        B_vec(4) = Gamma * (xi ** 2) * cav2_ss(j, k, at, cj, ck) + &
                 & -CONJG(gkl(j, cj)) * xi * cavsig2_ss(k, at, fe, ck) + &
                 & -CONJG(gkl(k, ck)) * xi * cavsig2_ss(j, at, fe, cj)
        B_vec(5) = i * xi * 0.5d0 * Omega * cav2_ss(j, k, at, cj, ck) + &
                 & CONJG(gkl(j, cj)) * xi * cavsig2_ss(k, at, gg, ck) + &
                 & CONJG(gkl(j, cj)) * xi * cavsig2_ss(k, at, ee, cj) + &
                 & -CONJG(gkl(j, cj)) * xi * cav1_ss(k, at, ck) + &
                 & CONJG(gkl(k, ck)) * xi * cavsig2_ss(j, at, gg, cj) + &
                 & CONJG(gkl(k, ck)) * xi * cavsig2_ss(j, at, ee, cj) + &
                 & -CONJG(gkl(k, ck)) * xi * cav1_ss(j, at, cj)
        B_vec(6) = -i * xi * 0.5d0 * Omega * cav2_ss(j, k, at, cj, ck)
        B_vec(7) = -CONJG(gkl(j, cj)) * cavsig2_ss(k, at, ef, ck) + &
                 & -CONJG(gkl(k, ck)) * cavsig2_ss(j, at, ef, cj)
        B_vec(8) = 0.0d0

        ! Set inverse matrix
        Mat_inv = Mat
        ! Perform LU-factorization of matrix
        CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
        ! Invert Matrix (Optimal LWORK = 8)
        LWORK = 8
        CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

        ! Calculate steady state
        cavsig3_ss(j, k, at, :, cj, ck) = -MATMUL(Mat_inv, B_vec)

        !----------------------------------!
        ! < a^{\dagger}_{j} a_{k} \sigma > !
        !----------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - (((kappa(cj) + kappa(ck))) - i * (wl(j, cj) - wl(k, ck)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j, cj)) * cavsig2_ss(k, a, eg, ck) + &
                 & -gkl(k, ck) * cavsig2_ss(j, at, ge, cj)
        B_vec(2) = -CONJG(gkl(j, cj)) * cavsig2_ss(k, a, ee, ck) + &
                 & -gkl(k, ck) * xi * cavsig2_ss(j, at, gf, cj)
        B_vec(3) = -CONJG(gkl(j, cj)) * xi * cavsig2_ss(k, a, fg, ck) + &
                 & -gkl(k, ck) * cavsig2_ss(j, at, ee, cj)
        B_vec(4) = Gamma * (xi ** 2) * cav2_ss(j, k, ata, cj, ck) + &
                 & -CONJG(gkl(j, cj)) * xi * cavsig2_ss(k, a, fe, ck) + &
                 & -gkl(k, ck) * xi * cavsig2_ss(j, at, ef, cj)
        B_vec(5) = i * xi * 0.5d0 * Omega * cav2_ss(j, k, ata, cj, ck) + &
                 & CONJG(gkl(j, cj)) * xi * cavsig2_ss(k, a, gg, ck) + &
                 & CONJG(gkl(j, cj)) * xi * cavsig2_ss(k, a, ee, ck) + &
                 & -CONJG(gkl(j, cj)) * xi * cav1_ss(k, a, ck)
        B_vec(6) = -i * xi * 0.5d0 * Omega * cav2_ss(j, k, ata, cj, ck) + &
                 & gkl(k, ck) * xi * cavsig2_ss(j, at, gg, cj) + &
                 & gkl(k, ck) * xi * cavsig2_ss(j, at, ee, cj) + &
                 & -gkl(k, ck) * xi * cav1_ss(j, at, cj)
        B_vec(7) = -CONJG(gkl(j, cj)) * cavsig2_ss(k, a, ef, ck)
        B_vec(8) = - gkl(k, ck) * cavsig2_ss(j, at, fe, cj)

        ! Set inverse matrix
        Mat_inv = Mat
        ! Perform LU-factorization of matrix
        CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
        ! Invert Matrix (Optimal LWORK = 8)
        LWORK = 8
        CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

        ! Calculate steady state
        cavsig3_ss(j, k, ata, :, cj, ck) = -MATMUL(Mat_inv, B_vec)

        ! Close j loop
      END DO
      ! Close k loop
    END DO
    ! Close cj loop
  END DO
  ! Close ck loop
END DO

! Cycle through cavity operators
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
            cav3_ss(j, k, l, a, cj, ck, cl) = -CONJG(gkl(j, cj)) * cavsig3_ss(k, l, a, eg, ck, cl) + &
                                            & -CONJG(gkl(j, cj)) * xi * cavsig3_ss(k, l, a, fe, ck, cl) + &
                                            & -gkl(k, ck) * cavsig3_ss(j, l, ata, ge, cj, cl) + &
                                            & -gkl(k, ck) * xi * cavsig3_ss(j, l, ata, ef, cj, cl) + &
                                            & -gkl(l, cl) * cavsig3_ss(j, k, ata, ge, cj, ck) + &
                                            & -gkl(l, cl) * xi * cavsig3_ss(j, k, ata, ef, cj, ck)
            cav3_ss(j, k, l, a, cj, ck, cl) = cav3_ss(j, k, l, a, cj, ck, cl) / &
                                            & ((kappa(cj) + kappa(ck) + kappa(cl)) - i * (wl(j, cj) - wl(k, ck) - wl(l, cl)))

            !-------------------------------------------!
            ! < a^{\dagger}_{j} a^{\dagger}_{k} a_{l} > !
            !-------------------------------------------!
            cav3_ss(j, k, l, at, cj, ck, cl) = -CONJG(gkl(j, cj)) * cavsig3_ss(k, l, ata, eg, ck, cl) + &
                                             & -CONJG(gkl(j, cj)) * xi * cavsig3_ss(k, l, ata, fe, ck, cl) + &
                                             & -CONJG(gkl(k, ck)) * cavsig3_ss(j, l, ata, eg, cj, cl) + &
                                             & -CONJG(gkl(k, ck)) * xi * cavsig3_ss(j, l, ata, fe, cj, cl) + &
                                             & -gkl(l, cl) * cavsig3_ss(j, k, at, ge, cj, ck) + &
                                             & -gkl(l, cl) * xi * cavsig3_ss(j, k, at, ef, cj, ck)
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
            B_vec = 0.0d0
            B_vec(1) = -CONJG(gkl(j, cj)) * cavsig3_ss(k, l, a, eg, ck, cl) + &
                     & -gkl(k, ck) * cavsig3_ss(j, l, ata, ge, cj, cl) + &
                       & -gkl(l, cl) * cavsig3_ss(j, k, ata, ge, cj, ck)
            B_vec(2) = -CONJG(gkl(j, cj)) * cavsig3_ss(k, l, a, ee, ck, cl) + &
                     & -gkl(k, ck) * xi * cavsig3_ss(j, l, ata, gf, cj, cl) + &
                     & -gkl(l, cl) * xi * cavsig3_ss(j, k, ata, gf, cj, ck)
            B_vec(3) = -CONJG(gkl(j, cj)) * xi * cavsig3_ss(k, l, a, fg, ck, cl) + &
                     & -gkl(k, ck) * cavsig3_ss(j, l, ata, ee, cj, cl) + &
                     & -gkl(l, cl) * cavsig3_ss(j, k, ata, ee, cj, ck)
            B_vec(4) = Gamma * (xi ** 2) * cav3_ss(j, k, l, a, cj, ck, cl) + &
                     & -CONJG(gkl(j, cj)) * xi * cavsig3_ss(k, l, a, fe, ck, cl) + &
                     & -gkl(k, ck) * xi * cavsig3_ss(j, l, ata, ef, cj, cl) + &
                     & -gkl(l, cl) * xi * cavsig3_ss(j, k, ata, ef, cj, ck)
            B_vec(5) = i * xi * 0.5d0 * Omega * cav3_ss(j, k, l, a, cj, ck, cl) + &
                     & CONJG(gkl(j, cj)) * xi * cavsig3_ss(k, l, a, gg, ck, cl) + &
                     & CONJG(gkl(j, cj)) * xi * cavsig3_ss(k, l, a, ee, ck, cl) + &
                     & -CONJG(gkl(j, cj)) * xi * cav2_ss(k, l, a, ck, cl)
            B_vec(6) = -i * xi * 0.5d0 * Omega * cav3_ss(j, k, l, a, cj, ck, cl) + &
                     & gkl(k, ck) * xi * cavsig3_ss(j, l, ata, gg, cj, cl) + &
                     & gkl(k, ck) * xi * cavsig3_ss(j, l, ata, ee, cj, cl) + &
                     & -gkl(k, ck) * xi * cav2_ss(j, l, ata, cj, cl) + &
                     & gkl(l, cl) * xi * cavsig3_ss(j, k, ata, gg, cj, ck) + &
                     & gkl(l, cl) * xi * cavsig3_ss(j, k, ata, ee, cj, ck) + &
                     & -gkl(l, cl) * xi * cav2_ss(j, k, ata, cj, ck)
            B_vec(7) = -CONJG(gkl(j, cj)) * cavsig3_ss(k, l, a, ef, ck, cl)
            B_vec(8) = -gkl(k, ck) * cavsig3_ss(j, l, ata, fe, cj, cl) + &
                     & -gkl(l, cl) * cavsig3_ss(j, k, ata, fe, cj, ck)

            ! Set inverse matrix
            Mat_inv = Mat
            ! Perform LU-factorization of matrix
            CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
            ! Invert Matrix (Optimal LWORK = 8)
            LWORK = 8
            CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

            ! Calculate steady state
            cavsig4_ss(j, k, l, a, :, cj, ck, cl) = -MATMUL(Mat_inv, B_vec)

            !-------------------------------------------!
            ! < a^{\dagger}_{j} a^{\dagger}_{k} a_{l} > !
            !-------------------------------------------!
            ! Set the diagonal matrix elements for M
            Mat = Mat_OG
            DO x = 1, N_mat
              Mat(x, x) = Mat(x, x) - ((kappa(cj) + kappa(ck) + kappa(cl)) - i * (wl(j, cj) + wl(k, ck) - wl(l, cl)))
            END DO

            ! Set the non-homogeneous vector
            B_vec = 0.0d0
            B_vec(1) = -CONJG(gkl(j, cj)) * cavsig3_ss(k, l, ata, eg, ck, cl) + &
                     & -CONJG(gkl(k, ck)) * cavsig3_ss(j, l, ata, eg, cj, cl) + &
                     & -gkl(l, cl) * cavsig3_ss(j, k, at, ge, cj, ck)
            B_vec(2) = -CONJG(gkl(j, cj)) * cavsig3_ss(k, l, ata, ee, ck, cl) + &
                     & -CONJG(gkl(k, ck)) * cavsig3_ss(j, l, ata, ee, cj, cl) + &
                     & -gkl(l, cl) * xi * cavsig3_ss(j, k, at, gf, cj, ck)
            B_vec(3) = -CONJG(gkl(j, cj)) * xi * cavsig3_ss(k, l, ata, fg, ck, cl) + &
                     & -CONJG(gkl(k, ck)) * xi * cavsig3_ss(j, l, ata, fg, cj, cl) + &
                     & -gkl(l, cl) * cavsig3_ss(j, k, at, ee, cj, ck)
            B_vec(4) = Gamma * (xi ** 2) * cav3_ss(j, k, l, at, cj, ck, cl) + &
                     & -CONJG(gkl(j, cj)) * xi * cavsig3_ss(k, l, ata, fe, ck, cl) + &
                     & -CONJG(gkl(k, ck)) * xi * cavsig3_ss(j, l, ata, fe, cj, cl) + &
                     & -gkl(l, cl) * xi * cavsig3_ss(j, k, at, ef, cj, ck)
            B_vec(5) = i * xi * 0.5d0 * Omega * cav3_ss(j, k, l, at, cj, ck, cl) + &
                     & CONJG(gkl(j, cj)) * xi * cavsig3_ss(k, l, ata, gg, ck, cl) + &
                     & CONJG(gkl(j, cj)) * xi * cavsig3_ss(k, l, ata, ee, ck, cl) + &
                     & -CONJG(gkl(j, cj)) * xi * cav2_ss(k, l, ata, ck, cl) + &
                     & CONJG(gkl(k, ck)) * xi * cavsig3_ss(j, l, ata, gg, cj, cl) + &
                     & CONJG(gkl(k, ck)) * xi * cavsig3_ss(j, l, ata, ee, cj, cl) + &
                     & -CONJG(gkl(k, ck)) * xi * cav2_ss(j, l, ata, cj, cl)
            B_vec(6) = -i * xi * 0.5d0 * Omega * cav3_ss(j, k, l, at, cj, ck, cl) + &
                     & gkl(l, cl) * xi * cavsig3_ss(j, k, at, gg, cj, ck) + &
                     & gkl(l, cl) * xi * cavsig3_ss(j, k, at, ee, cj, ck) + &
                     & -gkl(l, cl) * xi * cav2_ss(j, k, at, cj, ck)
            B_vec(7) = -CONJG(gkl(j, cj)) * cavsig3_ss(k, l, ata, ef, ck, cl) + &
                     & -CONJG(gkl(k, ck)) * cavsig3_ss(j, l, ata, ef, cj, cl)
            B_vec(8) = -gkl(l, cl) * cavsig3_ss(j, k, at, fe, cj, ck)

            ! Set inverse matrix
            Mat_inv = Mat
            ! Perform LU-factorization of matrix
            CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
            ! Invert Matrix (Optimal LWORK = 8)
            LWORK = 8
            CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

            ! Calculate steady state
            cavsig4_ss(j, k, l, at, :, cj, ck, cl) = -MATMUL(Mat_inv, B_vec)

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

! Cycle through cavity operators
DO cm = 1, 2
  DO cl = 1, 2
    DO ck = 1, 2
      DO cj = 1, 2
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
                cav4_ss(j, k, l, m, cj, ck, cl, cm) = -CONJG(gkl(j, cj)) * cavsig4_ss(k, l, m, a, eg, ck, cl, cm) + &
                                                    & -CONJG(gkl(j, cj)) * xi * cavsig4_ss(k, l, m, a, fe, ck, cl, cm) + &
                                                    & -CONJG(gkl(k, ck)) * cavsig4_ss(j, l, m, a, eg, cj, cl, cm) + &
                                                    & -CONJG(gkl(k, ck)) * xi * cavsig4_ss(j, l, m, a, fe, cj, cl, cm) + &
                                                    & -gkl(l, cl) * cavsig4_ss(j, k, m, at, ge, cj, ck, cm) + &
                                                    & -gkl(l, cl) * xi * cavsig4_ss(j, k, m, at, ef, cj, ck, cm) + &
                                                    & -gkl(m, cm) * cavsig4_ss(j, k, l, at, ge, cj, ck, cl) + &
                                                    & -gkl(m, cm) * xi * cavsig4_ss(j, k, l, at, ef, cj, ck, cl)
                cav4_ss(j, k, l, m, cj, ck, cl, cm) = cav4_ss(j, k, l, m, cj, ck, cl, cm) / &
                                                    & ((kappa(cj) + kappa(ck) + kappa(cl) + kappa(cm)) - i * (wl(j, cj) + wl(k, ck)) + i * (wl(l, cl) + wl(m, cm)))
                ! Close j loop
              END DO
              ! Close k loop
            END DO
            ! Close l loop
          END DO
          ! Close m loop
        END DO
        ! Close cj loop
      END DO
      ! Close ck loop
    END DO
    ! Close cl loop
  END DO
  ! Close cm loop
END DO

! Calculate photon number
photona_ss = 0.0d0
photonb_ss = 0.0d0
DO k = -N, N
  DO j = -N, N
    photona_ss = photona_ss + cav2_ss(j, k, ata, 1, 1)
    photonb_ss = photonb_ss + cav2_ss(j, k, ata, 2, 2)
  END DO
END DO

WRITE(*, *) "=============================================="
WRITE(*, *) "THIRD-ORDER: CAVITY"
moment_out = cav3_ss(0, 0, 0, a, 1, 1, 2)
WRITE(*, '(A23,ES18.11E2,A3,ES18.11E2,A1)') " < at_0 a_0 b_0 >_ss = ", REAL(moment_out), " + ", AIMAG(moment_out), "i"
moment_out = cav3_ss(0, 0, 0, a, 2, 2, 1)
WRITE(*, '(A23,ES18.11E2,A3,ES18.11E2,A1)') " < bt_0 b_0 a_0 >_ss = ", REAL(moment_out), " + ", AIMAG(moment_out), "i"
moment_out = cav3_ss(0, 0, 0, at, 2, 1, 1)
WRITE(*, '(A23,ES18.11E2,A3,ES18.11E2,A1)') "< bt_0 at_0 a_0 >_ss = ", REAL(moment_out), " + ", AIMAG(moment_out), "i"
moment_out = cav3_ss(0, 0, 0, at, 1, 2, 2)
WRITE(*, '(A23,ES18.11E2,A3,ES18.11E2,A1)') "< at_0 bt_0 b_0 >_ss = ", REAL(moment_out), " + ", AIMAG(moment_out), "i"
WRITE(*, *) "=============================================="


! WRITE(*, *) "=============================================="
! WRITE(*, *) "FOURTH-ORDER: CAVITY AND ATOM"
! moment_out = cavsig4_ss(0, 1, 2, a, fe, 1, 1, 2)
! WRITE(*, '(A26,ES18.11E2,A3,ES18.11E2,A1)') " < at_0 a_1 b_2 fe >_ss = ", REAL(moment_out), " + ", AIMAG(moment_out), "i"
! moment_out = cavsig4_ss(0, 1, 2, a, fe, 2, 2, 1)
! WRITE(*, '(A26,ES18.11E2,A3,ES18.11E2,A1)') " < bt_0 b_1 a_2 fe >_ss = ", REAL(moment_out), " + ", AIMAG(moment_out), "i"
! moment_out = cavsig4_ss(0, 1, 2, at, fe, 2, 1, 1)
! WRITE(*, '(A26,ES18.11E2,A3,ES18.11E2,A1)') "< bt_0 at_1 a_2 fe >_ss = ", REAL(moment_out), " + ", AIMAG(moment_out), "i"
! moment_out = cavsig4_ss(0, 1, 2, at, fe, 1, 2, 2)
! WRITE(*, '(A26,ES18.11E2,A3,ES18.11E2,A1)') "< at_0 bt_1 b_2 fe >_ss = ", REAL(moment_out), " + ", AIMAG(moment_out), "i"
! WRITE(*, *) "=============================================="

!==============================================================================!
!                  CALCULATE FIRST-ORDER CORRELATION FUNCTION                  !
!==============================================================================!
! Set initial conditions and non-homogeneous vector
! < a^{\dagger}_{j}(0) \sigma(\tau = 0) a_{m}(0) > =
!                         < a^{\dagger}_{j} a_{m} \sigma >_{ss},
! < a^{\dagger}_{j}(0) a^{\dagger}_{k}(\tau = 0) a_{m}(0) > =
!                         < a^{\dagger}_{j} a^{\dagger}_{k} a_{m} >_{ss},
! < a^{\dagger}_{j}(0) a_{l}(\tau = 0) a_{m}(0) > =
!                         < a^{\dagger}_{j} a_{l} a_{m} >_{ss},
! < a^{\dagger}_{j}(0) a^{\dagger}_{k} \sigma(\tau = 0) a_{m}(0) =
!                         < a^{\dagger}_{j} a^{\dagger}_{k} a_{m} \sigma >_{ss},
! < a^{\dagger}_{j}(0) a_{l} \sigma(\tau = 0) a_{m}(0) =
!                         < a^{\dagger}_{j} a_{l} a_{m} \sigma >_{ss},
! and
! < a^{\dagger}_{j}(0) a^{\dagger}_{k} a_{l}(\tau = 0)  a_{m}(0) > =
!                         < a^{\dagger}_{j} a^{\dagger}_{k} a_{l} a_{m} >_{ss}.

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
        cav2(k, l, ata) = cav2(k, l, ata) + cav4_ss(j, k, l, m, cj, ck, cl, cm)

        ! Close l loop
      END DO
      ! Close k loop
    END DO
    ! Non homogeneous vector
    B_OG(4) = B_OG(4) + Gamma * (xi ** 2) * cav2_ss(j, m, ata, cj, cm)
    B_OG(5) = B_OG(5) + i * xi * 0.5d0 * Omega * cav2_ss(j, m, ata, cj, cm)
    B_OG(6) = B_OG(6) - i * xi * 0.5d0 * Omega * cav2_ss(j, m, ata, cj, cm)

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
ck = 2

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
  IF (cj .EQ. 1 .AND. ck .EQ. 1) THEN
    IF (photona_ss .NE. 0.0d0) THEN
      moment_out = moment_out / (photona_ss ** 2)
    END IF
  ELSE IF (cj .EQ. 2 .AND. ck .EQ. 2) THEN
    IF (photona_ss .NE. 0.0d0 .AND. photonb_ss .NE. 0.0d0) THEN
      moment_out = moment_out / (photona_ss * photonb_ss)
    END IF
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
  ! First-order
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
                  & -dt * gkl(j, cj) * sigma(ge) + &
                  & -dt * gkl(j, cj) * xi * sigma(ef)

    k2_cav1(j, a) = -dt * (kappa(cj) + i * wl(j, cj)) * (cav1(j, a) + 0.5d0 * k1_cav1(j, a)) + &
                  & -dt * gkl(j, cj) * (sigma(ge) + 0.5d0 * k1_sigma(ge)) + &
                  & -dt * gkl(j, cj) * xi * (sigma(ef) + 0.5d0 * k1_sigma(ef))

    k3_cav1(j, a) = -dt * (kappa(cj) + i * wl(j, cj)) * (cav1(j, a) + 0.5d0 * k2_cav1(j, a)) + &
                  & -dt * gkl(j, cj) * (sigma(ge) + 0.5d0 * k2_sigma(ge)) + &
                  & -dt * gkl(j, cj) * xi * (sigma(ef) + 0.5d0 * k2_sigma(ef))

    k4_cav1(j, a) = -dt * (kappa(cj) + i * wl(j, cj)) * (cav1(j, a) + k3_cav1(j, a)) + &
                  & -dt * gkl(j, cj) * (sigma(ge) + k3_sigma(ge)) + &
                  & -dt * gkl(j, cj) * xi * (sigma(ef) + k3_sigma(ef))

    !---------------------!
    ! < a^{\dagger}_{j} > !
    !---------------------!
    k1_cav1(j, at) = -dt * (kappa(cj) - i * wl(j, cj)) * cav1(j, at) + &
                   & -dt * CONJG(gkl(j, cj)) * sigma(eg) + &
                   & -dt * CONJG(gkl(j, cj)) * xi * sigma(fe)

    k2_cav1(j, at) = -dt * (kappa(cj) - i * wl(j, cj)) * (cav1(j, at) + 0.5d0 * k1_cav1(j, at)) + &
                   & -dt * CONJG(gkl(j, cj)) * (sigma(eg) + 0.5d0 * k1_sigma(eg)) + &
                   & -dt * CONJG(gkl(j, cj)) * xi * (sigma(fe) + 0.5d0 * k1_sigma(fe))

    k3_cav1(j, at) = -dt * (kappa(cj) - i * wl(j, cj)) * (cav1(j, at) + 0.5d0 * k2_cav1(j, at)) + &
                   & -dt * CONJG(gkl(j, cj)) * (sigma(eg) + 0.5d0 * k2_sigma(eg)) + &
                   & -dt * CONJG(gkl(j, cj)) * xi * (sigma(fe) + 0.5d0 * k2_sigma(fe))

    k4_cav1(j, at) = -dt * (kappa(cj) - i * wl(j, cj)) * (cav1(j, at) + k3_cav1(j, at)) + &
                   & -dt * CONJG(gkl(j, cj)) * (sigma(eg) + k3_sigma(eg)) + &
                   & -dt * CONJG(gkl(j, cj)) * xi * (sigma(fe) + k3_sigma(fe))

    !---------------------------------------!
    !     SECOND-ORDER: CAVITY AND ATOM     !
    !---------------------------------------!
    ! Using matrix multiplication, we add to the Lindblad matrix Mat_OG and the
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
    B_vec = 0.0d0
    B_vec(1) = -gkl(j, cj) * sigma(ge)
    B_vec(2) = -gkl(j, cj) * xi * sigma(gf)
    B_vec(3) = -gkl(j, cj) * sigma(ee)
    B_vec(4) = Gamma * (xi ** 2) * cav1(j, a) + &
         & -gkl(j, cj) * xi * sigma(ef)
    B_vec(5) = i * xi * 0.5d0 * Omega * cav1(j, a)
    B_vec(6) = -i * xi * 0.5d0 * Omega * cav1(j, a) + &
         & gkl(j, cj) * xi * sigma(gg) + &
         & gkl(j, cj) * xi * sigma(ee) + &
         & -gkl(j, cj) * xi * photona_ss
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j, cj) * sigma(fe)
    ! Calculate k1
    k1_cavsig2(j, a, :) = dt * (MATMUL(Mat, cavsig2(j, a, :)) + B_vec)

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -gkl(j, cj) * (sigma(ge) + 0.5d0 * k1_sigma(ge))
    B_vec(2) = -gkl(j, cj) * xi * (sigma(gf) + 0.5d0 * k1_sigma(gf))
    B_vec(3) = -gkl(j, cj) * (sigma(ee) + 0.5d0 * k1_sigma(ee))
    B_vec(4) = Gamma * (xi ** 2) * (cav1(j, a) + 0.5d0 * k1_cav1(j, a)) + &
         & -gkl(j, cj) * xi * (sigma(ef) + 0.5d0 * k1_sigma(ef))
    B_vec(5) = i * xi * 0.5d0 * Omega * (cav1(j, a) + 0.5d0 * k1_cav1(j, a))
    B_vec(6) = -i * xi * 0.5d0 * Omega * (cav1(j, a) + 0.5d0 * k1_cav1(j, a)) + &
         & gkl(j, cj) * xi * (sigma(gg) + 0.5d0 * k1_sigma(gg)) + &
         & gkl(j, cj) * xi * (sigma(ee) + 0.5d0 * k1_sigma(ee)) + &
         & -gkl(j, cj) * xi * photona_ss
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j, cj) * (sigma(fe) + 0.5d0 * k1_sigma(fe))
    ! Calculate k2
    k2_cavsig2(j, a, :) = dt * (MATMUL(Mat, (cavsig2(j, a, :) + 0.5d0 * k1_cavsig2(j, a, :))) + B_vec)

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -gkl(j, cj) * (sigma(ge) + 0.5d0 * k2_sigma(ge))
    B_vec(2) = -gkl(j, cj) * xi * (sigma(gf) + 0.5d0 * k2_sigma(gf))
    B_vec(3) = -gkl(j, cj) * (sigma(ee) + 0.5d0 * k2_sigma(ee))
    B_vec(4) = Gamma * (xi ** 2) * (cav1(j, a) + 0.5d0 * k2_cav1(j, a)) + &
         & -gkl(j, cj) * xi * (sigma(ef) + 0.5d0 * k2_sigma(ef))
    B_vec(5) = i * xi * 0.5d0 * Omega * (cav1(j, a) + 0.5d0 * k2_cav1(j, a))
    B_vec(6) = -i * xi * 0.5d0 * Omega * (cav1(j, a) + 0.5d0 * k2_cav1(j, a)) + &
         & gkl(j, cj) * xi * (sigma(gg) + 0.5d0 * k2_sigma(gg)) + &
         & gkl(j, cj) * xi * (sigma(ee) + 0.5d0 * k2_sigma(ee)) + &
         & -gkl(j, cj) * xi * photona_ss
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j, cj) * (sigma(fe) + 0.5d0 * k2_sigma(fe))
    ! Calculate k3
    k3_cavsig2(j, a, :) = dt * (MATMUL(Mat, (cavsig2(j, a, :) + 0.5d0 * k2_cavsig2(j, a, :))) + B_vec)

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -gkl(j, cj) * (sigma(ge) + k3_sigma(ge))
    B_vec(2) = -gkl(j, cj) * xi * (sigma(gf) + k3_sigma(gf))
    B_vec(3) = -gkl(j, cj) * (sigma(ee) + k3_sigma(ee))
    B_vec(4) = Gamma * (xi ** 2) * (cav1(j, a) + k3_cav1(j, a)) + &
         & -gkl(j, cj) * xi * (sigma(ef) + k3_sigma(ef))
    B_vec(5) = i * xi * 0.5d0 * Omega * (cav1(j, a) + k3_cav1(j, a))
    B_vec(6) = -i * xi * 0.5d0 * Omega * (cav1(j, a) + k3_cav1(j, a)) + &
         & gkl(j, cj) * xi * (sigma(gg) + k3_sigma(gg)) + &
         & gkl(j, cj) * xi * (sigma(ee) + k3_sigma(ee)) + &
         & -gkl(j, cj) * xi * photona_ss
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j, cj) * (sigma(fe) + k3_sigma(fe))
    ! Calculate k4
    k4_cavsig2(j, a, :) = dt * (MATMUL(Mat, (cavsig2(j, a, :) + k3_cavsig2(j, a, :))) + B_vec)

    !----------------------------!
    ! < a^{\dagger}_{j} \sigma > !
    !----------------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - (kappa(cj) - i * wl(j, cj))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -CONJG(gkl(j, cj)) * sigma(eg)
    B_vec(2) = -CONJG(gkl(j, cj)) * sigma(ee)
    B_vec(3) = -CONJG(gkl(j, cj)) * xi * sigma(fg)
    B_vec(4) = Gamma * (xi ** 2) * cav1(j, at) + &
         & -CONJG(gkl(j, cj)) * xi * sigma(fe)
    B_vec(5) = i * xi * 0.5d0 * Omega * cav1(j, at) + &
         & CONJG(gkl(j, cj)) * xi * sigma(gg) + &
         & CONJG(gkl(j, cj)) * xi * sigma(ee) + &
         & -CONJG(gkl(j, cj)) * xi * photona_ss
    B_vec(6) = -i * xi * 0.5d0 * Omega * cav1(j, at)
    B_vec(7) = -CONJG(gkl(j, cj)) * sigma(ef)
    B_vec(8) = 0.0d0
    ! Calculate k1
    k1_cavsig2(j, at, :) = dt * (MATMUL(Mat, cavsig2(j, at, :)) + B_vec)

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -CONJG(gkl(j, cj)) * (sigma(eg) + 0.5d0 * k1_sigma(eg))
    B_vec(2) = -CONJG(gkl(j, cj)) * (sigma(ee) + 0.5d0 * k1_sigma(ee))
    B_vec(3) = -CONJG(gkl(j, cj)) * xi * (sigma(fg) + 0.5d0 * k1_sigma(fg))
    B_vec(4) = Gamma * (xi ** 2) * (cav1(j, at) + 0.5d0 * k1_cav1(j, at)) + &
         & -CONJG(gkl(j, cj)) * xi * (sigma(fe) + 0.5d0 * k1_sigma(fe))
    B_vec(5) = i * xi * 0.5d0 * Omega * (cav1(j, at) + 0.5d0 * k1_cav1(j, at)) + &
         & CONJG(gkl(j, cj)) * xi * (sigma(gg) + 0.5d0 * k1_sigma(gg)) + &
         & CONJG(gkl(j, cj)) * xi * (sigma(ee) + 0.5d0 * k1_sigma(ee)) + &
         & -CONJG(gkl(j, cj)) * xi * photona_ss
    B_vec(6) = -i * xi * 0.5d0 * Omega * (cav1(j, at) + 0.5d0 * k1_cav1(j, at))
    B_vec(7) = -CONJG(gkl(j, cj)) * (sigma(ef) + 0.5d0 * k1_sigma(ef))
    B_vec(8) = 0.0d0
    ! Calculate k2
    k2_cavsig2(j, at, :) = dt * (MATMUL(Mat, (cavsig2(j, at, :) + 0.5d0 * k1_cavsig2(j, at, :))) + B_vec)

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -CONJG(gkl(j, cj)) * (sigma(eg) + 0.5d0 * k2_sigma(eg))
    B_vec(2) = -CONJG(gkl(j, cj)) * (sigma(ee) + 0.5d0 * k2_sigma(ee))
    B_vec(3) = -CONJG(gkl(j, cj)) * xi * (sigma(fg) + 0.5d0 * k2_sigma(fg))
    B_vec(4) = Gamma * (xi ** 2) * (cav1(j, at) + 0.5d0 * k2_cav1(j, at)) + &
         & -CONJG(gkl(j, cj)) * xi * (sigma(fe) + 0.5d0 * k2_sigma(fe))
    B_vec(5) = i * xi * 0.5d0 * Omega * (cav1(j, at) + 0.5d0 * k2_cav1(j, at)) + &
         & CONJG(gkl(j, cj)) * xi * (sigma(gg) + 0.5d0 * k2_sigma(gg)) + &
         & CONJG(gkl(j, cj)) * xi * (sigma(ee) + 0.5d0 * k2_sigma(ee)) + &
         & -CONJG(gkl(j, cj)) * xi * photona_ss
    B_vec(6) = -i * xi * 0.5d0 * Omega * (cav1(j, at) + 0.5d0 * k2_cav1(j, at))
    B_vec(7) = -CONJG(gkl(j, cj)) * (sigma(ef) + 0.5d0 * k2_sigma(ef))
    B_vec(8) = 0.0d0
    ! Calculate k3
    k3_cavsig2(j, at, :) = dt * (MATMUL(Mat, (cavsig2(j, at, :) + 0.5d0 * k2_cavsig2(j, at, :))) + B_vec)

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -CONJG(gkl(j, cj)) * (sigma(eg) + k3_sigma(eg))
    B_vec(2) = -CONJG(gkl(j, cj)) * (sigma(ee) + k3_sigma(ee))
    B_vec(3) = -CONJG(gkl(j, cj)) * xi * (sigma(fg) + k3_sigma(fg))
    B_vec(4) = Gamma * (xi ** 2) * (cav1(j, at) + k3_cav1(j, at)) + &
         & -CONJG(gkl(j, cj)) * xi * (sigma(fe) + k3_sigma(fe))
    B_vec(5) = i * xi * 0.5d0 * Omega * (cav1(j, at) + k3_cav1(j, at)) + &
         & CONJG(gkl(j, cj)) * xi * (sigma(gg) + k3_sigma(gg)) + &
         & CONJG(gkl(j, cj)) * xi * (sigma(ee) + k3_sigma(ee)) + &
         & -CONJG(gkl(j, cj)) * xi * photona_ss
    B_vec(6) = -i * xi * 0.5d0 * Omega * (cav1(j, at) + k3_cav1(j, at))
    B_vec(7) = -CONJG(gkl(j, cj)) * (sigma(ef) + k3_sigma(ef))
    B_vec(8) = 0.0d0
    ! Calculate k4
    k4_cavsig2(j, at, :) = dt * (MATMUL(Mat, (cavsig2(j, at, :) + k3_cavsig2(j, at, :))) + B_vec)

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
      k1_cav2(j, k, ata) = -dt * ((kappa(cj) + kappa(ck)) - i * (wl(j, cj) - wl(k, ck))) * cav2(j, k, ata) + &
                         & -dt * CONJG(gkl(j, cj)) * cavsig2(k, a, eg) + &
                         & -dt * CONJG(gkl(j, cj)) * xi * cavsig2(k, a, fe) + &
                         & -dt * gkl(k, ck) * cavsig2(j, at, ge) + &
                         & -dt * gkl(k, ck) * xi * cavsig2(j, at, ef)

      k2_cav2(j, k, ata) = -dt * ((kappa(cj) + kappa(ck)) - i * (wl(j, cj) - wl(k, ck))) * (cav2(j, k, ata) + 0.5d0 * k1_cav2(j, k, ata)) + &
                         & -dt * CONJG(gkl(j, cj)) * (cavsig2(k, a, eg) + 0.5d0 * k1_cavsig2(k, a, eg)) + &
                         & -dt * CONJG(gkl(j, cj)) * xi * (cavsig2(k, a, fe) + 0.5d0 * k1_cavsig2(k, a, fe)) + &
                         & -dt * gkl(k, ck) * (cavsig2(j, at, ge) + 0.5d0 * k1_cavsig2(j, at, ge)) + &
                         & -dt * gkl(k, ck) * xi * (cavsig2(j, at, ef) + 0.5d0 * k1_cavsig2(j, at, ef))

      k3_cav2(j, k, ata) = -dt * ((kappa(cj) + kappa(ck)) - i * (wl(j, cj) - wl(k, ck))) * (cav2(j, k, ata) + 0.5d0 * k2_cav2(j, k, ata)) + &
                         & -dt * CONJG(gkl(j, cj)) * (cavsig2(k, a, eg) + 0.5d0 * k2_cavsig2(k, a, eg)) + &
                         & -dt * CONJG(gkl(j, cj)) * xi * (cavsig2(k, a, fe) + 0.5d0 * k2_cavsig2(k, a, fe)) + &
                         & -dt * gkl(k, ck) * (cavsig2(j, at, ge) + 0.5d0 * k2_cavsig2(j, at, ge)) + &
                         & -dt * gkl(k, ck) * xi * (cavsig2(j, at, ef) + 0.5d0 * k2_cavsig2(j, at, ef))

      k4_cav2(j, k, ata) = -dt * ((kappa(cj) + kappa(ck)) - i * (wl(j, cj) - wl(k, ck))) * (cav2(j, k, ata) + k3_cav2(j, k, ata)) + &
                         & -dt * CONJG(gkl(j, cj)) * (cavsig2(k, a, eg) + k3_cavsig2(k, a, eg)) + &
                         & -dt * CONJG(gkl(j, cj)) * xi * (cavsig2(k, a, fe) + k3_cavsig2(k, a, fe)) + &
                         & -dt * gkl(k, ck) * (cavsig2(j, at, ge) + k3_cavsig2(j, at, ge)) + &
                         & -dt * gkl(k, ck) * xi * (cavsig2(j, at, ef) + k3_cavsig2(j, at, ef))

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

END PROGRAM THREE_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_G2_CROSS
