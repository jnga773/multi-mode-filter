! The system in this program is a resonantly driven three-level ladder-type atom
! that is coupled into a multi-mode array of filter cavities, as described in
! some notes somewhere (or with this program if this goes to anyone).

! The operator moment steady states are calculated using the inverse matrix/
! analytical method. Using quantum regression theorem, we then use Runge-Kutta
! fourth-order to numerically integrate the moment equations, suing the steady
! states and initial conditions, to calculate the normalised second-order
! correlation function.

! The input parameters are taken from a NameList file [filename_ParamList] which,
! by default, points to "./ParamList.nml". The code can thus be compiled once,
! and parameters can be changed in the NameList file for subsequent runs.

! The parameters for each run are written to [filename_parameters] which, by
! default, points to "./data_files/parameters.txt".

! The normalised first-order correlation function is written to [filename_g1]
! which, by default, points to "./data_files/g1_corr.txt". The file has two
! columns:
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

PROGRAM THREE_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_G2

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
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                           :: kappaa
! Frequency spacing of modes
REAL(KIND=8)                                           :: dwa
! List of Delta values
REAL(KIND=8), DIMENSION(:), ALLOCATABLE                :: wl
! List of mode dependent cascade coupling values
COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE             :: gkl
! Blackman window coefficient
REAL(KIND=8)                                           :: blackman
! kappa, w0, and dw values
REAL(KIND=8)                                           :: kappa, w0, dw

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

! Time integration arrays
! First-order moments: Atomic equations (< \sigma >)
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: sigma
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: k1_sigma, k2_sigma, k3_sigma, k4_sigma
! First-order moments: Cavity (< a >, < a^{\dagger} >)
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: cav1
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: k1_cav1, k2_cav1, k3_cav1, k4_cav1
! Second-order moments: Cavity and atom (< a \sigma >, < a^{\dagger} \sigma >)
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cavsig2
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: k1_cavsig2, k2_cavsig2, k3_cavsig2, k4_cavsig2
! Second-order moments: Cavity (< a^{\dagger} a >)
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cav2
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: k1_cav2, k2_cav2, k3_cav2, k4_cav2

! Steady state arrays
! First-order moments: Atomic equations (< \sigma >)
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: sigma_ss
! First-order moments: Cavity (< a >, < a^{\dagger} >)
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: cav1_ss
! Second-order moments: Cavity and atom (< a \sigma >, < a^{\dagger} \sigma >
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cavsig2_ss
! Second-order moments: Cavity (< a^{\dagger} a >)
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cav2_ss
! Third-order moments: Cavity and atom (< a^{2} \sigma >, < a^{\dagger 2} \sigma >, < a^{\dagger} a \sigma >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cavsig3_ss
! Third-order moments: Cavity (< a^{2} a^{\dagger} >, < a^{\dagger 2} a >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cav3_ss
! Fourth-order moments: Cavity and atom ( < a^{\dagger} a^{2} \sigma >, < a^{\dagger 2} a \sigma >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: cavsig4_ss
! Fourth-order moments: Cavity (< a^{\dagger 2} a^{2} >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cav4_ss


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
REAL(KIND=8)                                           :: photon_ss
! Complex data
COMPLEX(KIND=8)                                        :: moment_out

!------------------------!
!     FILENAME STUFF     !
!------------------------!
! Paramert Name List
CHARACTER(LEN=15), PARAMETER :: filename_ParamList = "./ParamList.nml"
! Filename of parameters
CHARACTER(LEN=27), PARAMETER :: filename_parameters = "./data_files/parameters.txt"
! Filename for second-order correlation
CHARACTER(LEN=24), PARAMETER :: filename_g2 = "./data_files/g2_corr.txt"

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
    gkl(j) = DSQRT(epsilon * Gamma * kappa)
  ELSE
    wl(j) = w0 + DBLE(j) * dw
    ! Blackman window coefficient
    blackman = 1.0d0
    ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
    !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))
    ! Mode dependent phase difference
    gkl(j) = DSQRT((epsilon / DBLE(2*N + 1)) * Gamma * kappa) * blackman * EXP(i * DBLE(phase) * DBLE(j) * pi / DBLE(N))
  END IF
END DO

!------------------------------------------!
!     INITALISE OPERATOR MOMENT ARRAYS     !
!------------------------------------------!
! Time Integration
! First-order: Cavity
ALLOCATE(cav1(-N:N, 2)); cav1 = 0.0d0
ALLOCATE(k1_cav1(-N:N, 2)); k1_cav1 = 0.0d0
ALLOCATE(k2_cav1(-N:N, 2)); k2_cav1 = 0.0d0
ALLOCATE(k3_cav1(-N:N, 2)); k3_cav1 = 0.0d0
ALLOCATE(k4_cav1(-N:N, 2)); k4_cav1 = 0.0d0
! Second-order: Cavity and Atom
ALLOCATE(cavsig2(-N:N, 2, N_mat)); cavsig2 = 0.0d0
ALLOCATE(k1_cavsig2(-N:N, 2, N_mat)); k1_cavsig2 = 0.0d0
ALLOCATE(k2_cavsig2(-N:N, 2, N_mat)); k2_cavsig2 = 0.0d0
ALLOCATE(k3_cavsig2(-N:N, 2, N_mat)); k3_cavsig2 = 0.0d0
ALLOCATE(k4_cavsig2(-N:N, 2, N_mat)); k4_cavsig2 = 0.0d0
! Second-order: Cavity
ALLOCATE(cav2(-N:N, -N:N, 3)); cav2 = 0.0d0
ALLOCATE(k1_cav2(-N:N, -N:N, 3)); k1_cav2 = 0.0d0
ALLOCATE(k2_cav2(-N:N, -N:N, 3)); k2_cav2 = 0.0d0
ALLOCATE(k3_cav2(-N:N, -N:N, 3)); k3_cav2 = 0.0d0
ALLOCATE(k4_cav2(-N:N, -N:N, 3)); k4_cav2 = 0.0d0

! Steady states
! First-order: Cavity
ALLOCATE(cav1_ss(-N:N, 2)); cav1_ss = 0.0d0
! Second-order: Cavity and Atom
ALLOCATE(cavsig2_ss(-N:N, 2, N_mat)); cavsig2_ss = 0.0d0
! Second-order: Cavity
ALLOCATE(cav2_ss(-N:N, -N:N, 3)); cav2_ss = 0.0d0
! Third-order: Cavity and Atom
ALLOCATE(cavsig3_ss(-N:N, -N:N, 3, N_mat)); cavsig3_ss = 0.0d0
! Third-order: Cavity
ALLOCATE(cav3_ss(-N:N, -N:N, -N:N, 2)); cav3_ss = 0.0d0
! Fourth-order: Cavity and atom
ALLOCATE(cavsig4_ss(-N:N, -N:N, -N:N, 2, N_mat)); cavsig4_ss = 0.0d0
! Fourth-order: Cavity
ALLOCATE(cav4_ss(-N:N, -N:N, -N:N, -N:N)); cav4_ss = 0.0d0

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
!---------------------------!
!     FIRST-ORDER: ATOM     !
!---------------------------!
IF (xi .NE. 0.0d0) THEN
  ! Set matrix to be inverted
  Mat_inv = Mat_OG
  ! Invert matrix
  CALL SquareMatrixInverse(N_mat, Mat_inv)

  ! Calculate steady states
  sigma_ss = 0.0d0
  sigma_ss = -MATMUL(Mat_inv, B_OG)
ELSE IF (xi .EQ. 0.0d0) THEN
  !------------------------------------------------!
  !     CALCULATE EIGENVALUES AND EIGENVECTORS     !
  !------------------------------------------------!
  ! Set Lindblad matrix
  Mat = Mat_OG

  ! Calculate steady state from eigenvectors
  CALL SquareMatrixZeroEigenvalue(N_mat, Mat, sigma_ss)
  ! Normalise sigma_ss so |g><g| + |e><e| = 1
  sigma_ss = sigma_ss / (REAL(sigma_ss(1)) + REAL(sigma_ss(4)))
END IF

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
  ! Invert matrix
  CALL SquareMatrixInverse(N_mat, Mat_inv)

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
  ! Invert matrix
  CALL SquareMatrixInverse(N_mat, Mat_inv)

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
    !-----------------!
    ! < a_{j} a_{k} > !
    !-----------------!
    cav2_ss(j, k, a) = -gkl(j) * cavsig2_ss(k, a, ge) + &
                     & -gkl(j) * xi * cavsig2_ss(k, a, ef) + &
                     & -gkl(k) * cavsig2_ss(j, a, ge) + &
                     & -gkl(k) * xi * cavsig2_ss(j, a, ef)
    cav2_ss(j, k, a) = cav2_ss(j, k, a) / &
                     & (2.0d0 * kappa + i * (wl(j) + wl(k)))

    !-------------------------------------!
    ! < a^{\dagger}_{j} a^{\dagger}_{k} > !
    !-------------------------------------!
    cav2_ss(j, k, at) = -CONJG(gkl(j)) * cavsig2_ss(k, at, eg) + &
                      & -CONJG(gkl(j)) * xi * cavsig2_ss(k, at, fe) + &
                      & -CONJG(gkl(k)) * cavsig2_ss(j, at, eg) + &
                      & -CONJG(gkl(k)) * xi * cavsig2_ss(j, at, fe)
    cav2_ss(j, k, at) = cav2_ss(j, k, at) / &
                      & (2.0d0 * kappa - i * (wl(j) + wl(k)))

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

    !--------------------------------------!
    !     THIRD-ORDER: CAVITY AND ATOM     !
    !--------------------------------------!
    !------------------------!
    ! < a_{j} a_{k} \sigma > !
    !------------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa) + i * (wl(j) + wl(k)))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -gkl(j) * cavsig2_ss(k, a, ge) + &
             & -gkl(k) * cavsig2_ss(j, a, ge)
    B_vec(2) = -gkl(j) * xi * cavsig2_ss(k, a, gf) + &
             & -gkl(k) * xi * cavsig2_ss(j, a, gf)
    B_vec(3) = -gkl(j) * cavsig2_ss(k, a, ee) + &
             & -gkl(k) * cavsig2_ss(j, a, ee)
    B_vec(4) = Gamma * (xi ** 2) * cav2_ss(j, k, a) + &
             & -gkl(j) * xi * cavsig2_ss(k, a, ef) + &
             & -gkl(k) * xi * cavsig2_ss(j, a, ef)
    B_vec(5) = i * xi * 0.5d0 * Omega * cav2_ss(j, k, a)
    B_vec(6) = -i * xi * 0.5d0 * Omega * cav2_ss(j, k, a) + &
             & gkl(j) * xi * cavsig2_ss(k, a, gg) + &
             & gkl(j) * xi * cavsig2_ss(k, a, ee) + &
             & -gkl(j) * xi * cav1_ss(k, a) + &
             & gkl(k) * xi * cavsig2_ss(j, a, gg) + &
             & gkl(k) * xi * cavsig2_ss(j, a, ee) + &
             & -gkl(k) * xi * cav1_ss(j, a)
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j) * cavsig2_ss(k, a, fe) + &
             & -gkl(k) * cavsig2_ss(j, a, fe)

    ! Set inverse matrix
    Mat_inv = Mat
    ! Invert matrix
    CALL SquareMatrixInverse(N_mat, Mat_inv)

    ! Calculate steady state
    cavsig3_ss(j, k, a, :) = -MATMUL(Mat_inv, B_vec)

    !--------------------------------------------!
    ! < a^{\dagger}_{j} a^{\dagger}_{k} \sigma > !
    !--------------------------------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa) - i * (wl(j) + wl(k)))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -CONJG(gkl(j)) * cavsig2_ss(k, at, eg) + &
             & -CONJG(gkl(k)) * cavsig2_ss(j, at, eg)
    B_vec(2) = -CONJG(gkl(j)) * cavsig2_ss(k, at, ee) + &
             & -CONJG(gkl(k)) * cavsig2_ss(j, at, ee)
    B_vec(3) = -CONJG(gkl(j)) * xi * cavsig2_ss(k, at, fg) + &
             & -CONJG(gkl(k)) * xi * cavsig2_ss(j, at, fg)
    B_vec(4) = Gamma * (xi ** 2) * cav2_ss(j, k, at) + &
             & -CONJG(gkl(j)) * xi * cavsig2_ss(k, at, fe) + &
             & -CONJG(gkl(k)) * xi * cavsig2_ss(j, at, fe)
    B_vec(5) = i * xi * 0.5d0 * Omega * cav2_ss(j, k, at) + &
             & CONJG(gkl(j)) * xi * cavsig2_ss(k, at, gg) + &
             & CONJG(gkl(j)) * xi * cavsig2_ss(k, at, ee) + &
             & -CONJG(gkl(j)) * xi * cav1_ss(k, at) + &
             & CONJG(gkl(k)) * xi * cavsig2_ss(j, at, gg) + &
             & CONJG(gkl(k)) * xi * cavsig2_ss(j, at, ee) + &
             & -CONJG(gkl(k)) * xi * cav1_ss(j, at)
    B_vec(6) = -i * xi * 0.5d0 * Omega * cav2_ss(j, k, at)
    B_vec(7) = -CONJG(gkl(j)) * cavsig2_ss(k, at, ef) + &
             & -CONJG(gkl(k)) * cavsig2_ss(j, at, ef)
    B_vec(8) = 0.0d0

    ! Set inverse matrix
    Mat_inv = Mat
    ! Invert matrix
    CALL SquareMatrixInverse(N_mat, Mat_inv)

    ! Calculate steady state
    cavsig3_ss(j, k, at, :) = -MATMUL(Mat_inv, B_vec)

    !----------------------------------!
    ! < a^{\dagger}_{j} a_{k} \sigma > !
    !----------------------------------!
    ! Set the diagonal matrix elements for M
    Mat = Mat_OG
    DO x = 1, N_mat
      Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa) - i * (wl(j) - wl(k)))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -CONJG(gkl(j)) * cavsig2_ss(k, a, eg) + &
             & -gkl(k) * cavsig2_ss(j, at, ge)
    B_vec(2) = -CONJG(gkl(j)) * cavsig2_ss(k, a, ee) + &
             & -gkl(k) * xi * cavsig2_ss(j, at, gf)
    B_vec(3) = -CONJG(gkl(j)) * xi * cavsig2_ss(k, a, fg) + &
             & -gkl(k) * cavsig2_ss(j, at, ee)
    B_vec(4) = Gamma * (xi ** 2) * cav2_ss(j, k, ata) + &
             & -CONJG(gkl(j)) * xi * cavsig2_ss(k, a, fe) + &
             & -gkl(k) * xi * cavsig2_ss(j, at, ef)
    B_vec(5) = i * xi * 0.5d0 * Omega * cav2_ss(j, k, ata) + &
             & CONJG(gkl(j)) * xi * cavsig2_ss(k, a, gg) + &
             & CONJG(gkl(j)) * xi * cavsig2_ss(k, a, ee) + &
             & -CONJG(gkl(j)) * xi * cav1_ss(k, a)
    B_vec(6) = -i * xi * 0.5d0 * Omega * cav2_ss(j, k, ata) + &
             & gkl(k) * xi * cavsig2_ss(j, at, gg) + &
             & gkl(k) * xi * cavsig2_ss(j, at, ee) + &
             & -gkl(k) * xi * cav1_ss(j, at)
    B_vec(7) = -CONJG(gkl(j)) * cavsig2_ss(k, a, ef)
    B_vec(8) = - gkl(k) * cavsig2_ss(j, at, fe)

    ! Set inverse matrix
    Mat_inv = Mat
    ! Invert matrix
    CALL SquareMatrixInverse(N_mat, Mat_inv)

    ! Calculate steady state
    cavsig3_ss(j, k, ata, :) = -MATMUL(Mat_inv, B_vec)

    ! Close j loop
  END DO
  ! Close k loop
END DO

photon_ss = REAL(moment_out)

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
      cav3_ss(j, k, l, a) = -CONJG(gkl(j)) * cavsig3_ss(k, l, a, eg) + &
                          & -CONJG(gkl(j)) * xi * cavsig3_ss(k, l, a, fe) + &
                          & -gkl(k) * cavsig3_ss(j, l, ata, ge) + &
                          & -gkl(k) * xi * cavsig3_ss(j, l, ata, ef) + &
                          & -gkl(l) * cavsig3_ss(j, k, ata, ge) + &
                          & -gkl(l) * xi * cavsig3_ss(j, k, ata, ef)
      cav3_ss(j, k, l, a) = cav3_ss(j, k, l, a) / &
                          & (3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l)))

      !-------------------------------------------!
      ! < a^{\dagger}_{j} a^{\dagger}_{k} a_{l} > !
      !-------------------------------------------!
      cav3_ss(j, k, l, at) = -CONJG(gkl(j)) * cavsig3_ss(k, l, ata, eg) + &
                           & -CONJG(gkl(j)) * xi * cavsig3_ss(k, l, ata, fe) + &
                           & -CONJG(gkl(k)) * cavsig3_ss(j, l, ata, eg) + &
                           & -CONJG(gkl(k)) * xi * cavsig3_ss(j, l, ata, fe) + &
                           & -gkl(l) * cavsig3_ss(j, k, at, ge) + &
                           & -gkl(l) * xi * cavsig3_ss(j, k, at, ef)
      cav3_ss(j, k, l, at) = cav3_ss(j, k, l, at) / &
                           & (3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l)))

      !--------------------------------------!
      !     FOURTH-ORDER: CAVITY AND ATOM    !
      !--------------------------------------!
      !----------------------------------------!
      ! < a^{\dagger}_{j} a_{k} a_{l} \sigma > !
      !----------------------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - ((3.0d0 * kappa) - i * (wl(j) - wl(k) - wl(l)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j)) * cavsig3_ss(k, l, a, eg) + &
               & -gkl(k) * cavsig3_ss(j, l, ata, ge) + &
               & -gkl(l) * cavsig3_ss(j, k, ata, ge)
      B_vec(2) = -CONJG(gkl(j)) * cavsig3_ss(k, l, a, ee) + &
               & -gkl(k) * xi * cavsig3_ss(j, l, ata, gf) + &
               & -gkl(l) * xi * cavsig3_ss(j, k, ata, gf)
      B_vec(3) = -CONJG(gkl(j)) * xi * cavsig3_ss(k, l, a, fg) + &
               & -gkl(k) * cavsig3_ss(j, l, ata, ee) + &
               & -gkl(l) * cavsig3_ss(j, k, ata, ee)
      B_vec(4) = Gamma * (xi ** 2) * cav3_ss(j, k, l, a) + &
               & -CONJG(gkl(j)) * xi * cavsig3_ss(k, l, a, fe) + &
               & -gkl(k) * xi * cavsig3_ss(j, l, ata, ef) + &
               & -gkl(l) * xi * cavsig3_ss(j, k, ata, ef)
      B_vec(5) = i * xi * 0.5d0 * Omega * cav3_ss(j, k, l, a) + &
               & CONJG(gkl(j)) * xi * cavsig3_ss(k, l, a, gg) + &
               & CONJG(gkl(j)) * xi * cavsig3_ss(k, l, a, ee) + &
               & -CONJG(gkl(j)) * xi * cav2_ss(k, l, a)
      B_vec(6) = -i * xi * 0.5d0 * Omega * cav3_ss(j, k, l, a) + &
               & gkl(k) * xi * cavsig3_ss(j, l, ata, gg) + &
               & gkl(k) * xi * cavsig3_ss(j, l, ata, ee) + &
               & -gkl(k) * xi * cav2_ss(j, l, ata) + &
               & gkl(l) * xi * cavsig3_ss(j, k, ata, gg) + &
               & gkl(l) * xi * cavsig3_ss(j, k, ata, ee) + &
               & -gkl(l) * xi * cav2_ss(j, k, ata)
      B_vec(7) = -CONJG(gkl(j)) * cavsig3_ss(k, l, a, ef)
      B_vec(8) = -gkl(k) * cavsig3_ss(j, l, ata, fe) + &
               & -gkl(l) * cavsig3_ss(j, k, ata, fe)

      ! Set inverse matrix
      Mat_inv = Mat
      ! Invert matrix
      CALL SquareMatrixInverse(N_mat, Mat_inv)

      ! Calculate steady state
      cavsig4_ss(j, k, l, a, :) = -MATMUL(Mat_inv, B_vec)

      !--------------------------------------------------!
      ! < a^{\dagger}_{j} a^{\dagger}_{k} a_{l} \sigma > !
      !--------------------------------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - ((3.0d0 * kappa) - i * (wl(j) + wl(k) - wl(l)))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j)) * cavsig3_ss(k, l, ata, eg) + &
               & -CONJG(gkl(k)) * cavsig3_ss(j, l, ata, eg) + &
               & -gkl(l) * cavsig3_ss(j, k, at, ge)
      B_vec(2) = -CONJG(gkl(j)) * cavsig3_ss(k, l, ata, ee) + &
               & -CONJG(gkl(k)) * cavsig3_ss(j, l, ata, ee) + &
               & -gkl(l) * xi * cavsig3_ss(j, k, at, gf)
      B_vec(3) = -CONJG(gkl(j)) * xi * cavsig3_ss(k, l, ata, fg) + &
               & -CONJG(gkl(k)) * xi * cavsig3_ss(j, l, ata, fg) + &
               & -gkl(l) * cavsig3_ss(j, k, at, ee)
      B_vec(4) = Gamma * (xi ** 2) * cav3_ss(j, k, l, at) + &
               & -CONJG(gkl(j)) * xi * cavsig3_ss(k, l, ata, fe) + &
               & -CONJG(gkl(k)) * xi * cavsig3_ss(j, l, ata, fe) + &
               & -gkl(l) * xi * cavsig3_ss(j, k, at, ef)
      B_vec(5) = i * xi * 0.5d0 * Omega * cav3_ss(j, k, l, at) + &
               & CONJG(gkl(j)) * xi * cavsig3_ss(k, l, ata, gg) + &
               & CONJG(gkl(j)) * xi * cavsig3_ss(k, l, ata, ee) + &
               & -CONJG(gkl(j)) * xi * cav2_ss(k, l, ata) + &
               & CONJG(gkl(k)) * xi * cavsig3_ss(j, l, ata, gg) + &
               & CONJG(gkl(k)) * xi * cavsig3_ss(j, l, ata, ee) + &
               & -CONJG(gkl(k)) * xi * cav2_ss(j, l, ata)
      B_vec(6) = -i * xi * 0.5d0 * Omega * cav3_ss(j, k, l, at) + &
               & gkl(l) * xi * cavsig3_ss(j, k, at, gg) + &
               & gkl(l) * xi * cavsig3_ss(j, k, at, ee) + &
               & -gkl(l) * xi * cav2_ss(j, k, at)
      B_vec(7) = -CONJG(gkl(j)) * cavsig3_ss(k, l, ata, ef) + &
               & -CONJG(gkl(k)) * cavsig3_ss(j, l, ata, ef)
      B_vec(8) = -gkl(l) * cavsig3_ss(j, k, at, fe)

      ! Set inverse matrix
      Mat_inv = Mat
      ! Invert matrix
      CALL SquareMatrixInverse(N_mat, Mat_inv)

      ! Calculate steady state
      cavsig4_ss(j, k, l, at, :) = -MATMUL(Mat_inv, B_vec)

      ! Close j loop
    END DO
    ! Close k loop
  END DO
  ! Close l loop
END DO

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
        cav4_ss(j, k, l, m) = -CONJG(gkl(j)) * cavsig4_ss(k, l, m, a, eg) + &
                            & -CONJG(gkl(j)) * xi * cavsig4_ss(k, l, m, a, fe) + &
                            & -CONJG(gkl(k)) * cavsig4_ss(j, l, m, a, eg) + &
                            & -CONJG(gkl(k)) * xi * cavsig4_ss(j, l, m, a, fe) + &
                            & -gkl(l) * cavsig4_ss(j, k, m, at, ge) + &
                            & -gkl(l) * xi * cavsig4_ss(j, k, m, at, ef) + &
                            & -gkl(m) * cavsig4_ss(j, k, l, at, ge) + &
                            & -gkl(m) * xi * cavsig4_ss(j, k, l, at, ef)
        cav4_ss(j, k, l, m) = cav4_ss(j, k, l, m) / &
                            & (4.0d0 * kappa - i * (wl(j) + wl(k)) + i * (wl(l) + wl(m)))
        ! Close j loop
      END DO
      ! Close k loop
    END DO
    ! Close l loop
  END DO
  ! Close m loop
END DO

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
! Cycle through modes
DO m = -N, N
  DO j = -N, N
    ! First-order: Atom
    sigma(:) = sigma(:) + cavsig3_ss(j, m, ata, :)

    DO k = -N, N
      ! First-order: Cavity
      cav1(k, a) = cav1(k, a) + cav3_ss(j, k, m, a)
      cav1(k, at) = cav1(k, at) + cav3_ss(j, k, m, at)
      ! Second-order: Cavity and atom
      cavsig2(k, a, :) = cavsig2(k, a, :) + cavsig4_ss(j, k, m, a, :)
      cavsig2(k, at, :) = cavsig2(k, at, :) + cavsig4_ss(j, k, m, at, :)

      DO l = -N, N
        ! Second-order: cavity
        cav2(k, l, ata) = cav2(k, l, ata) + cav4_ss(j, k, l, m)

        ! Close l loop
      END DO
      ! Close k loop
    END DO
    ! Non homogeneous vector
    B_OG(4) = B_OG(4) + Gamma * (xi ** 2) * cav2_ss(j, m, ata)
    B_OG(5) = B_OG(5) + i * xi * 0.5d0 * Omega * cav2_ss(j, m, ata)
    B_OG(6) = B_OG(6) - i * xi * 0.5d0 * Omega * cav2_ss(j, m, ata)

    ! Close j loop
  END DO
  ! Close m loop
END DO

! Set tau_steps
tau_steps = NINT(tau2_max / DBLE(dt))

! Open file to write time and data to
OPEN(UNIT=4, FILE=filename_g2, STATUS='REPLACE', ACTION='WRITE', RECL=4000)

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
  IF (photon_ss .NE. 0.0) THEN
    moment_out = moment_out / (photon_ss ** 2)
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
    k1_cav1(j, a) = -dt * (kappa + i * wl(j)) * cav1(j, a) + &
                  & -dt * gkl(j) * sigma(ge) + &
                  & -dt * gkl(j) * xi * sigma(ef)

    k2_cav1(j, a) = -dt * (kappa + i * wl(j)) * (cav1(j, a) + 0.5d0 * k1_cav1(j, a)) + &
                  & -dt * gkl(j) * (sigma(ge) + 0.5d0 * k1_sigma(ge)) + &
                  & -dt * gkl(j) * xi * (sigma(ef) + 0.5d0 * k1_sigma(ef))

    k3_cav1(j, a) = -dt * (kappa + i * wl(j)) * (cav1(j, a) + 0.5d0 * k2_cav1(j, a)) + &
                  & -dt * gkl(j) * (sigma(ge) + 0.5d0 * k2_sigma(ge)) + &
                  & -dt * gkl(j) * xi * (sigma(ef) + 0.5d0 * k2_sigma(ef))

    k4_cav1(j, a) = -dt * (kappa + i * wl(j)) * (cav1(j, a) + k3_cav1(j, a)) + &
                  & -dt * gkl(j) * (sigma(ge) + k3_sigma(ge)) + &
                  & -dt * gkl(j) * xi * (sigma(ef) + k3_sigma(ef))

    !---------------------!
    ! < a^{\dagger}_{j} > !
    !---------------------!
    k1_cav1(j, at) = -dt * (kappa - i * wl(j)) * cav1(j, at) + &
                   & -dt * CONJG(gkl(j)) * sigma(eg) + &
                   & -dt * CONJG(gkl(j)) * xi * sigma(fe)

    k2_cav1(j, at) = -dt * (kappa - i * wl(j)) * (cav1(j, at) + 0.5d0 * k1_cav1(j, at)) + &
                   & -dt * CONJG(gkl(j)) * (sigma(eg) + 0.5d0 * k1_sigma(eg)) + &
                   & -dt * CONJG(gkl(j)) * xi * (sigma(fe) + 0.5d0 * k1_sigma(fe))

    k3_cav1(j, at) = -dt * (kappa - i * wl(j)) * (cav1(j, at) + 0.5d0 * k2_cav1(j, at)) + &
                   & -dt * CONJG(gkl(j)) * (sigma(eg) + 0.5d0 * k2_sigma(eg)) + &
                   & -dt * CONJG(gkl(j)) * xi * (sigma(fe) + 0.5d0 * k2_sigma(fe))

    k4_cav1(j, at) = -dt * (kappa - i * wl(j)) * (cav1(j, at) + k3_cav1(j, at)) + &
                   & -dt * CONJG(gkl(j)) * (sigma(eg) + k3_sigma(eg)) + &
                   & -dt * CONJG(gkl(j)) * xi * (sigma(fe) + k3_sigma(fe))

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
      Mat(x, x) = Mat(x, x) - (kappa + i * wl(j))
    END DO

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -gkl(j) * sigma(ge)
    B_vec(2) = -gkl(j) * xi * sigma(gf)
    B_vec(3) = -gkl(j) * sigma(ee)
    B_vec(4) = Gamma * (xi ** 2) * cav1(j, a) + &
             & -gkl(j) * xi * sigma(ef)
    B_vec(5) = i * xi * 0.5d0 * Omega * cav1(j, a)
    B_vec(6) = -i * xi * 0.5d0 * Omega * cav1(j, a) + &
             & gkl(j) * xi * sigma(gg) + &
             & gkl(j) * xi * sigma(ee) + &
             & -gkl(j) * xi * photon_ss
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j) * sigma(fe)
    ! Calculate k1
    k1_cavsig2(j, a, :) = dt * (MATMUL(Mat, cavsig2(j, a, :)) + B_vec)

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -gkl(j) * (sigma(ge) + 0.5d0 * k1_sigma(ge))
    B_vec(2) = -gkl(j) * xi * (sigma(gf) + 0.5d0 * k1_sigma(gf))
    B_vec(3) = -gkl(j) * (sigma(ee) + 0.5d0 * k1_sigma(ee))
    B_vec(4) = Gamma * (xi ** 2) * (cav1(j, a) + 0.5d0 * k1_cav1(j, a)) + &
             & -gkl(j) * xi * (sigma(ef) + 0.5d0 * k1_sigma(ef))
    B_vec(5) = i * xi * 0.5d0 * Omega * (cav1(j, a) + 0.5d0 * k1_cav1(j, a))
    B_vec(6) = -i * xi * 0.5d0 * Omega * (cav1(j, a) + 0.5d0 * k1_cav1(j, a)) + &
             & gkl(j) * xi * (sigma(gg) + 0.5d0 * k1_sigma(gg)) + &
             & gkl(j) * xi * (sigma(ee) + 0.5d0 * k1_sigma(ee)) + &
             & -gkl(j) * xi * photon_ss
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j) * (sigma(fe) + 0.5d0 * k1_sigma(fe))
    ! Calculate k2
    k2_cavsig2(j, a, :) = dt * (MATMUL(Mat, (cavsig2(j, a, :) + 0.5d0 * k1_cavsig2(j, a, :))) + B_vec)

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -gkl(j) * (sigma(ge) + 0.5d0 * k2_sigma(ge))
    B_vec(2) = -gkl(j) * xi * (sigma(gf) + 0.5d0 * k2_sigma(gf))
    B_vec(3) = -gkl(j) * (sigma(ee) + 0.5d0 * k2_sigma(ee))
    B_vec(4) = Gamma * (xi ** 2) * (cav1(j, a) + 0.5d0 * k2_cav1(j, a)) + &
             & -gkl(j) * xi * (sigma(ef) + 0.5d0 * k2_sigma(ef))
    B_vec(5) = i * xi * 0.5d0 * Omega * (cav1(j, a) + 0.5d0 * k2_cav1(j, a))
    B_vec(6) = -i * xi * 0.5d0 * Omega * (cav1(j, a) + 0.5d0 * k2_cav1(j, a)) + &
             & gkl(j) * xi * (sigma(gg) + 0.5d0 * k2_sigma(gg)) + &
             & gkl(j) * xi * (sigma(ee) + 0.5d0 * k2_sigma(ee)) + &
             & -gkl(j) * xi * photon_ss
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j) * (sigma(fe) + 0.5d0 * k2_sigma(fe))
    ! Calculate k3
    k3_cavsig2(j, a, :) = dt * (MATMUL(Mat, (cavsig2(j, a, :) + 0.5d0 * k2_cavsig2(j, a, :))) + B_vec)

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -gkl(j) * (sigma(ge) + k3_sigma(ge))
    B_vec(2) = -gkl(j) * xi * (sigma(gf) + k3_sigma(gf))
    B_vec(3) = -gkl(j) * (sigma(ee) + k3_sigma(ee))
    B_vec(4) = Gamma * (xi ** 2) * (cav1(j, a) + k3_cav1(j, a)) + &
             & -gkl(j) * xi * (sigma(ef) + k3_sigma(ef))
    B_vec(5) = i * xi * 0.5d0 * Omega * (cav1(j, a) + k3_cav1(j, a))
    B_vec(6) = -i * xi * 0.5d0 * Omega * (cav1(j, a) + k3_cav1(j, a)) + &
             & gkl(j) * xi * (sigma(gg) + k3_sigma(gg)) + &
             & gkl(j) * xi * (sigma(ee) + k3_sigma(ee)) + &
             & -gkl(j) * xi * photon_ss
    B_vec(7) = 0.0d0
    B_vec(8) = -gkl(j) * (sigma(fe) + k3_sigma(fe))
    ! Calculate k4
    k4_cavsig2(j, a, :) = dt * (MATMUL(Mat, (cavsig2(j, a, :) + k3_cavsig2(j, a, :))) + B_vec)

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
    B_vec(1) = -CONJG(gkl(j)) * sigma(eg)
    B_vec(2) = -CONJG(gkl(j)) * sigma(ee)
    B_vec(3) = -CONJG(gkl(j)) * xi * sigma(fg)
    B_vec(4) = Gamma * (xi ** 2) * cav1(j, at) + &
             & -CONJG(gkl(j)) * xi * sigma(fe)
    B_vec(5) = i * xi * 0.5d0 * Omega * cav1(j, at) + &
             & CONJG(gkl(j)) * xi * sigma(gg) + &
             & CONJG(gkl(j)) * xi * sigma(ee) + &
             & -CONJG(gkl(j)) * xi * photon_ss
    B_vec(6) = -i * xi * 0.5d0 * Omega * cav1(j, at)
    B_vec(7) = -CONJG(gkl(j)) * sigma(ef)
    B_vec(8) = 0.0d0
    ! Calculate k1
    k1_cavsig2(j, at, :) = dt * (MATMUL(Mat, cavsig2(j, at, :)) + B_vec)

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -CONJG(gkl(j)) * (sigma(eg) + 0.5d0 * k1_sigma(eg))
    B_vec(2) = -CONJG(gkl(j)) * (sigma(ee) + 0.5d0 * k1_sigma(ee))
    B_vec(3) = -CONJG(gkl(j)) * xi * (sigma(fg) + 0.5d0 * k1_sigma(fg))
    B_vec(4) = Gamma * (xi ** 2) * (cav1(j, at) + 0.5d0 * k1_cav1(j, at)) + &
             & -CONJG(gkl(j)) * xi * (sigma(fe) + 0.5d0 * k1_sigma(fe))
    B_vec(5) = i * xi * 0.5d0 * Omega * (cav1(j, at) + 0.5d0 * k1_cav1(j, at)) + &
             & CONJG(gkl(j)) * xi * (sigma(gg) + 0.5d0 * k1_sigma(gg)) + &
             & CONJG(gkl(j)) * xi * (sigma(ee) + 0.5d0 * k1_sigma(ee)) + &
             & -CONJG(gkl(j)) * xi * photon_ss
    B_vec(6) = -i * xi * 0.5d0 * Omega * (cav1(j, at) + 0.5d0 * k1_cav1(j, at))
    B_vec(7) = -CONJG(gkl(j)) * (sigma(ef) + 0.5d0 * k1_sigma(ef))
    B_vec(8) = 0.0d0
    ! Calculate k2
    k2_cavsig2(j, at, :) = dt * (MATMUL(Mat, (cavsig2(j, at, :) + 0.5d0 * k1_cavsig2(j, at, :))) + B_vec)

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -CONJG(gkl(j)) * (sigma(eg) + 0.5d0 * k2_sigma(eg))
    B_vec(2) = -CONJG(gkl(j)) * (sigma(ee) + 0.5d0 * k2_sigma(ee))
    B_vec(3) = -CONJG(gkl(j)) * xi * (sigma(fg) + 0.5d0 * k2_sigma(fg))
    B_vec(4) = Gamma * (xi ** 2) * (cav1(j, at) + 0.5d0 * k2_cav1(j, at)) + &
             & -CONJG(gkl(j)) * xi * (sigma(fe) + 0.5d0 * k2_sigma(fe))
    B_vec(5) = i * xi * 0.5d0 * Omega * (cav1(j, at) + 0.5d0 * k2_cav1(j, at)) + &
             & CONJG(gkl(j)) * xi * (sigma(gg) + 0.5d0 * k2_sigma(gg)) + &
             & CONJG(gkl(j)) * xi * (sigma(ee) + 0.5d0 * k2_sigma(ee)) + &
             & -CONJG(gkl(j)) * xi * photon_ss
    B_vec(6) = -i * xi * 0.5d0 * Omega * (cav1(j, at) + 0.5d0 * k2_cav1(j, at))
    B_vec(7) = -CONJG(gkl(j)) * (sigma(ef) + 0.5d0 * k2_sigma(ef))
    B_vec(8) = 0.0d0
    ! Calculate k3
    k3_cavsig2(j, at, :) = dt * (MATMUL(Mat, (cavsig2(j, at, :) + 0.5d0 * k2_cavsig2(j, at, :))) + B_vec)

    ! Set the non-homogeneous vector
    B_vec = 0.0d0
    B_vec(1) = -CONJG(gkl(j)) * (sigma(eg) + k3_sigma(eg))
    B_vec(2) = -CONJG(gkl(j)) * (sigma(ee) + k3_sigma(ee))
    B_vec(3) = -CONJG(gkl(j)) * xi * (sigma(fg) + k3_sigma(fg))
    B_vec(4) = Gamma * (xi ** 2) * (cav1(j, at) + k3_cav1(j, at)) + &
             & -CONJG(gkl(j)) * xi * (sigma(fe) + k3_sigma(fe))
    B_vec(5) = i * xi * 0.5d0 * Omega * (cav1(j, at) + k3_cav1(j, at)) + &
             & CONJG(gkl(j)) * xi * (sigma(gg) + k3_sigma(gg)) + &
             & CONJG(gkl(j)) * xi * (sigma(ee) + k3_sigma(ee)) + &
             & -CONJG(gkl(j)) * xi * photon_ss
    B_vec(6) = -i * xi * 0.5d0 * Omega * (cav1(j, at) + k3_cav1(j, at))
    B_vec(7) = -CONJG(gkl(j)) * (sigma(ef) + k3_sigma(ef))
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
      k1_cav2(j, k, ata) = -dt * (2.0d0 * kappa - i * (wl(j) - wl(k))) * cav2(j, k, ata) + &
                         & -dt * CONJG(gkl(j)) * cavsig2(k, a, eg) + &
                         & -dt * CONJG(gkl(j)) * xi * cavsig2(k, a, fe) + &
                         & -dt * gkl(k) * cavsig2(j, at, ge) + &
                         & -dt * gkl(k) * xi * cavsig2(j, at, ef)

      k2_cav2(j, k, ata) = -dt * (2.0d0 * kappa - i * (wl(j) - wl(k))) * (cav2(j, k, ata) + 0.5d0 * k1_cav2(j, k, ata)) + &
                         & -dt * CONJG(gkl(j)) * (cavsig2(k, a, eg) + 0.5d0 * k1_cavsig2(k, a, eg)) + &
                         & -dt * CONJG(gkl(j)) * xi * (cavsig2(k, a, fe) + 0.5d0 * k1_cavsig2(k, a, fe)) + &
                         & -dt * gkl(k) * (cavsig2(j, at, ge) + 0.5d0 * k1_cavsig2(j, at, ge)) + &
                         & -dt * gkl(k) * xi * (cavsig2(j, at, ef) + 0.5d0 * k1_cavsig2(j, at, ef))

      k3_cav2(j, k, ata) = -dt * (2.0d0 * kappa - i * (wl(j) - wl(k))) * (cav2(j, k, ata) + 0.5d0 * k2_cav2(j, k, ata)) + &
                         & -dt * CONJG(gkl(j)) * (cavsig2(k, a, eg) + 0.5d0 * k2_cavsig2(k, a, eg)) + &
                         & -dt * CONJG(gkl(j)) * xi * (cavsig2(k, a, fe) + 0.5d0 * k2_cavsig2(k, a, fe)) + &
                         & -dt * gkl(k) * (cavsig2(j, at, ge) + 0.5d0 * k2_cavsig2(j, at, ge)) + &
                         & -dt * gkl(k) * xi * (cavsig2(j, at, ef) + 0.5d0 * k2_cavsig2(j, at, ef))

      k4_cav2(j, k, ata) = -dt * (2.0d0 * kappa - i * (wl(j) - wl(k))) * (cav2(j, k, ata) + k3_cav2(j, k, ata)) + &
                         & -dt * CONJG(gkl(j)) * (cavsig2(k, a, eg) + k3_cavsig2(k, a, eg)) + &
                         & -dt * CONJG(gkl(j)) * xi * (cavsig2(k, a, fe) + k3_cavsig2(k, a, fe)) + &
                         & -dt * gkl(k) * (cavsig2(j, at, ge) + k3_cavsig2(j, at, ge)) + &
                         & -dt * gkl(k) * xi * (cavsig2(j, at, ef) + k3_cavsig2(j, at, ef))

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

END PROGRAM THREE_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_G2

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

  ! Set optimal LWORK for N_in = 3 or N_in = 8
  IF (N_in .EQ. 3) THEN
    LWORK = 3
  ELSE IF (N_in .EQ. 8) THEN
    LWORK = 8
  END IF

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
