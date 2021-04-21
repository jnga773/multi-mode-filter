! The system in this program is a resonantly driven three-level ladder-type atom
! that is coupled into a multi-mode array of filter cavities, as described in
! some notes somewhere (or with this program if this goes to anyone).

! Using Runge-Kutta Fourth-Order method of numerically integrating differential
! equations we solve for the first four orders of cavity-atom operator moment
! equations.

! The input parameters are taken from a NameList file [filename_ParamList] which,
! by default, points to "./ParamList.nml". The code can thus be compiled once,
! and parameters can be changed in the NameList file for subsequent runs.

! The parameters for each run are written to [filename_parameters] which, by
! default, points to "./data_files/parameters.txt".

! The atomic populations of the ground |g>, intermediate |e>, and excited |f>
! states, and the mean photon number inside the cavity <A^{\dagger} A(t)> are
! calculated from the operator moments and written to [filename_states] which,
! by default, points to "./data_files/states.txt". The file is in five columns:
!      t     |g><g|(t)     |e><e|(t)     |f><f|(t)     <A^{\dagger} A(t)>

! moment_out is a complex variable that can be assigned to any of the calculated
! moments. This can be compared with results from other software, such as QuTiP
! to verify the equations and program are correct. These results are written to
! [filename_moments] which, by default, points to "./data_files/operators.txt".
! The file is in three columns:
!               t     REAL(moment_out)     IMAG(moment_out)

! For the default filenames, the folder "./data_files/" and NameList file
! "./ParamList.nml" MUST EXIST IN THE WORKING DIRECTORY.

! To compiled the code, I use the Intel Parallel Studio compiler IFORT with the
! command:
!       (LINUX): ifort -O3 -o RK4 all_order_moments_RK4.f90
!     (WINDOWS): ifort /O3 /o RK4 all_order_moments_RK4.f90
! where the -O3 (/O3) flag gives maximum optimisation, and the -o (/o) RK4 flag
! names the executable as "RK4" ("RK4.exe")

PROGRAM THREE_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_RK4

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
! REAL(KIND=8)                                           :: w0b
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                           :: kappaa
! REAL(KIND=8)                                           :: kappab
! Frequency spacing of modes
REAL(KIND=8)                                           :: dwa
! REAL(KIND=8)                                           :: dwb
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
INTEGER                                                :: t_steps
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
! Third-order moments: Cavity and atom (< a^{2} \sigma >, < a^{\dagger 2} \sigma >, < a^{\dagger} a \sigma >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cavsig3
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: k1_cavsig3, k2_cavsig3, k3_cavsig3, k4_cavsig3
! Third-order moments: Cavity (< a^{2} a^{\dagger} >, < a^{\dagger 2} a >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cav3
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: k1_cav3, k2_cav3, k3_cav3, k4_cav3
! Fourth-order moments: Cavity and atom ( < a^{\dagger} a^{2} \sigma >, < a^{\dagger 2} a \sigma >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: cavsig4
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: k1_cavsig4, k2_cavsig4, k3_cavsig4, k4_cavsig4
! Fourth-order moments: Cavity (< a^{\dagger 2} a^{2} >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cav4
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: k1_cav4, k2_cav4, k3_cav4, k4_cav4

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
INTEGER                                                :: j, k, l, m, x
! Imaginary i
COMPLEX(KIND=8), PARAMETER                             :: i = CMPLX(0.0d0, 1.0d0, 8)
! pi
REAL(KIND=8), PARAMETER                                :: pi = 3.1415926535897932384d0
! 1 / 6
REAL(KIND=8), PARAMETER                                :: xis = 1.0d0 / 6.0d0
! Atomic population
REAL(KIND=8)                                           :: popg, pope, popf
! Photon number
REAL(KIND=8)                                           :: photon, photon_ss
! Complex data
COMPLEX(KIND=8)                                        :: moment_out

!------------------------!
!     FILENAME STUFF     !
!------------------------!
! Paramert Name List
CHARACTER(LEN=15), PARAMETER :: filename_ParamList = "./ParamList.nml"
! Filename of parameters
CHARACTER(LEN=27), PARAMETER :: filename_parameters = "./data_files/parameters.txt"
! Filename for state population
CHARACTER(LEN=23), PARAMETER :: filename_states = "./data_files/states.txt"
! Filename for operators
CHARACTER(LEN=26), PARAMETER :: filename_cavity = "./data_files/operators.txt"

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
! NAMELIST /CAVITYB/ kappab, w0b, dwb
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

! READ(IUNIT, NML=CAVITYB, IOSTAT=ISTAT)
! IF (ISTAT .NE. 0) THEN
!   BACKSPACE(IUNIT)
!   READ(IUNIT, FMT='(A)') LINE
!   CLOSE(IUNIT)
!   PRINT *, "Invalid line in CAVITYB namelist: " // TRIM(line)
!   CALL EXIT(1)
! END IF

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
! Third-order: Cavity and Atom
ALLOCATE(cavsig3(-N:N, -N:N, 3, N_mat)); cavsig3 = 0.0d0
ALLOCATE(k1_cavsig3(-N:N, -N:N, 3, N_mat)); k1_cavsig3 = 0.0d0
ALLOCATE(k2_cavsig3(-N:N, -N:N, 3, N_mat)); k2_cavsig3 = 0.0d0
ALLOCATE(k3_cavsig3(-N:N, -N:N, 3, N_mat)); k3_cavsig3 = 0.0d0
ALLOCATE(k4_cavsig3(-N:N, -N:N, 3, N_mat)); k4_cavsig3 = 0.0d0
! Third-order: Cavity
ALLOCATE(cav3(-N:N, -N:N, -N:N, 2)); cav3 = 0.0d0
ALLOCATE(k1_cav3(-N:N, -N:N, -N:N, 2)); k1_cav3 = 0.0d0
ALLOCATE(k2_cav3(-N:N, -N:N, -N:N, 2)); k2_cav3 = 0.0d0
ALLOCATE(k3_cav3(-N:N, -N:N, -N:N, 2)); k3_cav3 = 0.0d0
ALLOCATE(k4_cav3(-N:N, -N:N, -N:N, 2)); k4_cav3 = 0.0d0
! Fourth-order: Cavity and atom
ALLOCATE(cavsig4(-N:N, -N:N, -N:N, 2, N_mat)); cavsig4 = 0.0d0
ALLOCATE(k1_cavsig4(-N:N, -N:N, -N:N, 2, N_mat)); k1_cavsig4 = 0.0d0
ALLOCATE(k2_cavsig4(-N:N, -N:N, -N:N, 2, N_mat)); k2_cavsig4 = 0.0d0
ALLOCATE(k3_cavsig4(-N:N, -N:N, -N:N, 2, N_mat)); k3_cavsig4 = 0.0d0
ALLOCATE(k4_cavsig4(-N:N, -N:N, -N:N, 2, N_mat)); k4_cavsig4 = 0.0d0
! Fourth-order: Cavity
ALLOCATE(cav4(-N:N, -N:N, -N:N, -N:N)); cav4 = 0.0d0
ALLOCATE(k1_cav4(-N:N, -N:N, -N:N, -N:N)); k1_cav4 = 0.0d0
ALLOCATE(k2_cav4(-N:N, -N:N, -N:N, -N:N)); k2_cav4 = 0.0d0
ALLOCATE(k3_cav4(-N:N, -N:N, -N:N, -N:N)); k3_cav4 = 0.0d0
ALLOCATE(k4_cav4(-N:N, -N:N, -N:N, -N:N)); k4_cav4 = 0.0d0

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

!----------------------------!
!     INITIAL CONDITIONS     !
!----------------------------!
sigma = 0.0d0
! Atom in ground state, cavity in vacuum.
sigma(gg) = 1.0d0

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
!                        CALCULATE TIME-DEPENDENT VALUES                       !
!==============================================================================!

! Open file to write time and data to
OPEN(UNIT=2, FILE=filename_states, STATUS='REPLACE', ACTION='WRITE', RECL=4000)
OPEN(UNIT=10, FILE=filename_cavity, STATUS='REPLACE', ACTION='WRITE', RECL=4000)

! Cycle through time steps
DO t = 0, t_steps
  !============================================================================!
  !                          CALCULATE AND WRITE DATA                          !
  !============================================================================!
  !-------------------------------!
  !     CALCULATE DATA STATES     !
  !-------------------------------!
  ! < |g><g|(t) >
  popg = REAL(sigma(gg))
  ! < |e><e|(t) >
  pope = REAL(sigma(ee))
  ! < |f><f|(t) >
  popf = 1.0d0 - (popg + pope)

  ! < A^{\dagger} A >
  moment_out = 0.0d0
  ! Cycle through modes
  DO k = -N, N
    DO j = -N, N
      moment_out = moment_out + cav2(j, k, ata)
    END DO
  END DO
  photon = REAL(moment_out)

  ! Random moment to be written to filename_cavity
  ! Set mode j
  j = 0
  ! moment_out = cav1(0, a)
  ! moment_out = cavsig2(0, at, ee)
  ! moment_out = cav2(0, 0, at)
  ! moment_out = cavsig3(0, 0, ata, fe)
  ! moment_out = cav3(0, 0, 0, at)
  ! moment_out = cavsig4(0, 0, 0, at, ee)
  moment_out = cav4(0, 0, 0, 0)

  !-----------------------!
  !     WRITE TO FILE     !
  !-----------------------!
  ! Population data
  WRITE(2, *) DBLE(t) * dt, popg, pope, popf, photon

  ! Random moment data
  WRITE(10, *) DBLE(t) * dt, REAL(moment_out), AIMAG(moment_out)

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
  ! Third-order
  k1_cavsig3 = 0.0d0; k2_cavsig3 = 0.0d0; k3_cavsig3 = 0.0d0; k4_cavsig3 = 0.0d0
  k1_cav3 = 0.0d0; k2_cav3 = 0.d0; k3_cav3 = 0.0d0; k4_cav3 = 0.0d0
  ! Fourth-order
  k1_cavsig4 = 0.0d0; k2_cavsig4 = 0.0d0; k3_cavsig4 = 0.0d0; k4_cavsig4 = 0.0d0
  k1_cav4 = 0.0d0; k2_cav4 = 0.d0; k3_cav4 = 0.0d0; k4_cav4 = 0.0d0

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
             & -gkl(j) * xi
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
             & -gkl(j) * xi
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
             & -gkl(j) * xi
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
             & -gkl(j) * xi
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
             & -CONJG(gkl(j)) * xi
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
             & -CONJG(gkl(j)) * xi
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
             & -CONJG(gkl(j)) * xi
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
             & -CONJG(gkl(j)) * xi
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
      !-----------------!
      ! < a_{j} a_{k} > !
      !-----------------!
      k1_cav2(j, k, a) = -dt * (2.0d0 * kappa + i * (wl(j) + wl(k))) * cav2(j, k, a) + &
                       & -dt * gkl(j) * cavsig2(k, a, ge) + &
                       & -dt * gkl(j) * xi * cavsig2(k, a, ef) + &
                       & -dt * gkl(k) * cavsig2(j, a, ge) + &
                       & -dt * gkl(k) * xi * cavsig2(j, a, ef)

      k2_cav2(j, k, a) = -dt * (2.0d0 * kappa + i * (wl(j) + wl(k))) * (cav2(j, k, a) + 0.5d0 * k1_cav2(j, k, a)) + &
                       & -dt * gkl(j) * (cavsig2(k, a, ge) + 0.5d0 * k1_cavsig2(k, a, ge)) + &
                       & -dt * gkl(j) * xi * (cavsig2(k, a, ef) + 0.5d0 * k1_cavsig2(k, a, ef)) + &
                       & -dt * gkl(k) * (cavsig2(j, a, ge) + 0.5d0 * k1_cavsig2(j, a, ge)) + &
                       & -dt * gkl(k) * xi * (cavsig2(j, a, ef) + 0.5d0 * k1_cavsig2(j, a, ef))

      k3_cav2(j, k, a) = -dt * (2.0d0 * kappa + i * (wl(j) + wl(k))) * (cav2(j, k, a) + 0.5d0 * k2_cav2(j, k, a)) + &
                       & -dt * gkl(j) * (cavsig2(k, a, ge) + 0.5d0 * k2_cavsig2(k, a, ge)) + &
                       & -dt * gkl(j) * xi * (cavsig2(k, a, ef) + 0.5d0 * k2_cavsig2(k, a, ef)) + &
                       & -dt * gkl(k) * (cavsig2(j, a, ge) + 0.5d0 * k2_cavsig2(j, a, ge)) + &
                       & -dt * gkl(k) * xi * (cavsig2(j, a, ef) + 0.5d0 * k2_cavsig2(j, a, ef))

      k4_cav2(j, k, a) = -dt * (2.0d0 * kappa + i * (wl(j) + wl(k))) * (cav2(j, k, a) + k3_cav2(j, k, a)) + &
                       & -dt * gkl(j) * (cavsig2(k, a, ge) + k3_cavsig2(k, a, ge)) + &
                       & -dt * gkl(j) * xi * (cavsig2(k, a, ef) + k3_cavsig2(k, a, ef)) + &
                       & -dt * gkl(k) * (cavsig2(j, a, ge) + k3_cavsig2(j, a, ge)) + &
                       & -dt * gkl(k) * xi * (cavsig2(j, a, ef) + k3_cavsig2(j, a, ef))

      !-------------------------------------!
      ! < a^{\dagger}_{j} a^{\dagger}_{k} > !
      !-------------------------------------!
      k1_cav2(j, k, at) = -dt * (2.0d0 * kappa - i * (wl(j) + wl(k))) * cav2(j, k, at) + &
                        & -dt * CONJG(gkl(j)) * cavsig2(k, at, eg) + &
                        & -dt * CONJG(gkl(j)) * xi * cavsig2(k, at, fe) + &
                        & -dt * CONJG(gkl(k)) * cavsig2(j, at, eg) + &
                        & -dt * CONJG(gkl(k)) * xi * cavsig2(j, at, fe)

      k2_cav2(j, k, at) = -dt * (2.0d0 * kappa - i * (wl(j) + wl(k))) * (cav2(j, k, at) + 0.5d0 * k1_cav2(j, k, at)) + &
                        & -dt * CONJG(gkl(j)) * (cavsig2(k, at, eg) + 0.5d0 * k1_cavsig2(k, at, eg)) + &
                        & -dt * CONJG(gkl(j)) * xi * (cavsig2(k, at, fe) + 0.5d0 * k1_cavsig2(k, at, fe)) + &
                        & -dt * CONJG(gkl(k)) * (cavsig2(j, at, eg) + 0.5d0 * k1_cavsig2(j, at, eg)) + &
                        & -dt * CONJG(gkl(k)) * xi * (cavsig2(j, at, fe) + 0.5d0 * k1_cavsig2(j, at, fe))

      k3_cav2(j, k, at) = -dt * (2.0d0 * kappa - i * (wl(j) + wl(k))) * (cav2(j, k, at) + 0.5d0 * k2_cav2(j, k, at)) + &
                        & -dt * CONJG(gkl(j)) * (cavsig2(k, at, eg) + 0.5d0 * k1_cavsig2(k, at, eg)) + &
                        & -dt * CONJG(gkl(j)) * xi * (cavsig2(k, at, fe) + 0.5d0 * k1_cavsig2(k, at, fe)) + &
                        & -dt * CONJG(gkl(k)) * (cavsig2(j, at, eg) + 0.5d0 * k2_cavsig2(j, at, eg)) + &
                        & -dt * CONJG(gkl(k)) * xi * (cavsig2(j, at, fe) + 0.5d0 * k2_cavsig2(j, at, fe))

      k4_cav2(j, k, at) = -dt * (2.0d0 * kappa - i * (wl(j) + wl(k))) * (cav2(j, k, at) + k3_cav2(j, k, at)) + &
                        & -dt * CONJG(gkl(j)) * (cavsig2(k, at, eg) + k3_cavsig2(k, at, eg)) + &
                        & -dt * CONJG(gkl(j)) * xi * (cavsig2(k, at, fe) + k3_cavsig2(k, at, fe)) + &
                        & -dt * CONJG(gkl(k)) * (cavsig2(j, at, eg) + k3_cavsig2(j, at, eg)) + &
                        & -dt * CONJG(gkl(k)) * xi * (cavsig2(j, at, fe) + k3_cavsig2(j, at, fe))

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
      B_vec(1) = -gkl(j) * cavsig2(k, a, ge) + &
               & -gkl(k) * cavsig2(j, a, ge)
      B_vec(2) = -gkl(j) * xi * cavsig2(k, a, gf) + &
               & -gkl(k) * xi * cavsig2(j, a, gf)
      B_vec(3) = -gkl(j) * cavsig2(k, a, ee) + &
               & -gkl(k) * cavsig2(j, a, ee)
      B_vec(4) = Gamma * (xi ** 2) * cav2(j, k, a) + &
               & -gkl(j) * xi * cavsig2(k, a, ef) + &
               & -gkl(k) * xi * cavsig2(j, a, ef)
      B_vec(5) = i * xi * 0.5d0 * Omega * cav2(j, k, a)
      B_vec(6) = -i * xi * 0.5d0 * Omega * cav2(j, k, a) + &
               & gkl(j) * xi * cavsig2(k, a, gg) + &
               & gkl(j) * xi * cavsig2(k, a, ee) + &
               & -gkl(j) * xi * cav1(k, a) + &
               & gkl(k) * xi * cavsig2(j, a, gg) + &
               & gkl(k) * xi * cavsig2(j, a, ee) + &
               & -gkl(k) * xi * cav1(j, a)
      B_vec(7) = 0.0d0
      B_vec(8) = -gkl(j) * cavsig2(k, a, fe) + &
               & -gkl(k) * cavsig2(j, a, fe)
      ! Calculate k1
      k1_cavsig3(j, k, a, :) = dt * (MATMUL(Mat, cavsig3(j, k, a, :)) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -gkl(j) * (cavsig2(k, a, ge) + 0.5d0 * k1_cavsig2(k, a, ge)) + &
               & -gkl(k) * (cavsig2(j, a, ge) + 0.5d0 * k1_cavsig2(j, a, ge))
      B_vec(2) = -gkl(j) * xi * (cavsig2(k, a, gf) + 0.5d0 * k1_cavsig2(k, a, gf)) + &
               & -gkl(k) * xi * (cavsig2(j, a, gf) + 0.5d0 * k1_cavsig2(j, a, gf))
      B_vec(3) = -gkl(j) * (cavsig2(k, a, ee) + 0.5d0 * k1_cavsig2(k, a, ee)) + &
               & -gkl(k) * (cavsig2(j, a, ee) + 0.5d0 * k1_cavsig2(j, a, ee))
      B_vec(4) = Gamma * (xi ** 2) * (cav2(j, k, a) + 0.5d0 * k1_cav2(j, k, a)) + &
               & -gkl(j) * xi * (cavsig2(k, a, ef) + 0.5d0 * k1_cavsig2(k, a, ef)) + &
               & -gkl(k) * xi * (cavsig2(j, a, ef) + 0.5d0 * k1_cavsig2(j, a, ef))
      B_vec(5) = i * xi * 0.5d0 * Omega * (cav2(j, k, a) + 0.5d0 * k1_cav2(j, k, a))
      B_vec(6) = -i * xi * 0.5d0 * Omega * (cav2(j, k, a) + 0.5d0 * k1_cav2(j, k, a)) + &
               & gkl(j) * xi * (cavsig2(k, a, gg) + 0.5d0 * k1_cavsig2(k, a, gg)) + &
               & gkl(j) * xi * (cavsig2(k, a, ee) + 0.5d0 * k1_cavsig2(k, a, ee)) + &
               & -gkl(j) * xi * (cav1(k, a) + 0.5d0 * k1_cav1(k, a)) + &
               & gkl(k) * xi * (cavsig2(j, a, gg) + 0.5d0 * k1_cavsig2(j, a, gg)) + &
               & gkl(k) * xi * (cavsig2(j, a, ee) + 0.5d0 * k1_cavsig2(j, a, ee)) + &
               & -gkl(k) * xi * (cav1(j, a) + 0.5d0 * k1_cav1(j, a))
      B_vec(7) = 0.0d0
      B_vec(8) = -gkl(j) * (cavsig2(k, a, fe) + 0.5d0 * k1_cavsig2(k, a, fe)) + &
               & -gkl(k) * (cavsig2(j, a, fe) + 0.5d0 * k1_cavsig2(j, a, fe))
      ! Calculate k2
      k2_cavsig3(j, k, a, :) = dt * (MATMUL(Mat, (cavsig3(j, k, a, :) + 0.5d0 * k1_cavsig3(j, k, a, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -gkl(j) * (cavsig2(k, a, ge) + 0.5d0 * k2_cavsig2(k, a, ge)) + &
               & -gkl(k) * (cavsig2(j, a, ge) + 0.5d0 * k2_cavsig2(j, a, ge))
      B_vec(2) = -gkl(j) * xi * (cavsig2(k, a, gf) + 0.5d0 * k2_cavsig2(k, a, gf)) + &
               & -gkl(k) * xi * (cavsig2(j, a, gf) + 0.5d0 * k2_cavsig2(j, a, gf))
      B_vec(3) = -gkl(j) * (cavsig2(k, a, ee) + 0.5d0 * k2_cavsig2(k, a, ee)) + &
               & -gkl(k) * (cavsig2(j, a, ee) + 0.5d0 * k2_cavsig2(j, a, ee))
      B_vec(4) = Gamma * (xi ** 2) * (cav2(j, k, a) + 0.5d0 * k2_cav2(j, k, a)) + &
               & -gkl(j) * xi * (cavsig2(k, a, ef) + 0.5d0 * k2_cavsig2(k, a, ef)) + &
               & -gkl(k) * xi * (cavsig2(j, a, ef) + 0.5d0 * k2_cavsig2(j, a, ef))
      B_vec(5) = i * xi * 0.5d0 * Omega * (cav2(j, k, a) + 0.5d0 * k2_cav2(j, k, a))
      B_vec(6) = -i * xi * 0.5d0 * Omega * (cav2(j, k, a) + 0.5d0 * k2_cav2(j, k, a)) + &
               & gkl(j) * xi * (cavsig2(k, a, gg) + 0.5d0 * k2_cavsig2(k, a, gg)) + &
               & gkl(j) * xi * (cavsig2(k, a, ee) + 0.5d0 * k2_cavsig2(k, a, ee)) + &
               & -gkl(j) * xi * (cav1(k, a) + 0.5d0 * k2_cav1(k, a)) + &
               & gkl(k) * xi * (cavsig2(j, a, gg) + 0.5d0 * k2_cavsig2(j, a, gg)) + &
               & gkl(k) * xi * (cavsig2(j, a, ee) + 0.5d0 * k2_cavsig2(j, a, ee)) + &
               & -gkl(k) * xi * (cav1(j, a) + 0.5d0 * k2_cav1(j, a))
      B_vec(7) = 0.0d0
      B_vec(8) = -gkl(j) * (cavsig2(k, a, fe) + 0.5d0 * k2_cavsig2(k, a, fe)) + &
               & -gkl(k) * (cavsig2(j, a, fe) + 0.5d0 * k2_cavsig2(j, a, fe))
      ! Calculate k3
      k3_cavsig3(j, k, a, :) = dt * (MATMUL(Mat, (cavsig3(j, k, a, :) + 0.5d0 * k2_cavsig3(j, k, a, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -gkl(j) * (cavsig2(k, a, ge) + k3_cavsig2(k, a, ge)) + &
               & -gkl(k) * (cavsig2(j, a, ge) + k3_cavsig2(j, a, ge))
      B_vec(2) = -gkl(j) * xi * (cavsig2(k, a, gf) + k3_cavsig2(k, a, gf)) + &
               & -gkl(k) * xi * (cavsig2(j, a, gf) + k3_cavsig2(j, a, gf))
      B_vec(3) = -gkl(j) * (cavsig2(k, a, ee) + k3_cavsig2(k, a, ee)) + &
               & -gkl(k) * (cavsig2(j, a, ee) + k3_cavsig2(j, a, ee))
      B_vec(4) = Gamma * (xi ** 2) * (cav2(j, k, a) + k3_cav2(j, k, a)) + &
               & -gkl(j) * xi * (cavsig2(k, a, ef) + k3_cavsig2(k, a, ef)) + &
               & -gkl(k) * xi * (cavsig2(j, a, ef) + k3_cavsig2(j, a, ef))
      B_vec(5) = i * xi * 0.5d0 * Omega * (cav2(j, k, a) + k3_cav2(j, k, a))
      B_vec(6) = -i * xi * 0.5d0 * Omega * (cav2(j, k, a) + k3_cav2(j, k, a)) + &
               & gkl(j) * xi * (cavsig2(k, a, gg) + k3_cavsig2(k, a, gg)) + &
               & gkl(j) * xi * (cavsig2(k, a, ee) + k3_cavsig2(k, a, ee)) + &
               & -gkl(j) * xi * (cav1(k, a) + k3_cav1(k, a)) + &
               & gkl(k) * xi * (cavsig2(j, a, gg) + k3_cavsig2(j, a, gg)) + &
               & gkl(k) * xi * (cavsig2(j, a, ee) + k3_cavsig2(j, a, ee)) + &
               & -gkl(k) * xi * (cav1(j, a) + k3_cav1(j, a))
      B_vec(7) = 0.0d0
      B_vec(8) = -gkl(j) * (cavsig2(k, a, fe) + k3_cavsig2(k, a, fe)) + &
               & -gkl(k) * (cavsig2(j, a, fe) + k3_cavsig2(j, a, fe))
      ! Calculate k4
      k4_cavsig3(j, k, a, :) = dt * (MATMUL(Mat, (cavsig3(j, k, a, :) + k3_cavsig3(j, k, a, :))) + B_vec)

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
      B_vec(1) = -CONJG(gkl(j)) * cavsig2(k, at, eg) + &
               & -CONJG(gkl(k)) * cavsig2(j, at, eg)
      B_vec(2) = -CONJG(gkl(j)) * cavsig2(k, at, ee) + &
               & -CONJG(gkl(k)) * cavsig2(j, at, ee)
      B_vec(3) = -CONJG(gkl(j)) * xi * cavsig2(k, at, fg) + &
               & -CONJG(gkl(k)) * xi * cavsig2(j, at, fg)
      B_vec(4) = Gamma * (xi ** 2) * cav2(j, k, at) + &
               & -CONJG(gkl(j)) * xi * cavsig2(k, at, fe) + &
               & -CONJG(gkl(k)) * xi * cavsig2(j, at, fe)
      B_vec(5) = i * xi * 0.5d0 * Omega * cav2(j, k, at) + &
               & CONJG(gkl(j)) * xi * cavsig2(k, at, gg) + &
               & CONJG(gkl(j)) * xi * cavsig2(k, at, ee) + &
               & -CONJG(gkl(j)) * xi * cav1(k, at) + &
               & CONJG(gkl(k)) * xi * cavsig2(j, at, gg) + &
               & CONJG(gkl(k)) * xi * cavsig2(j, at, ee) + &
               & -CONJG(gkl(k)) * xi * cav1(j, at)
      B_vec(6) = -i * xi * 0.5d0 * Omega * cav2(j, k, at)
      B_vec(7) = -CONJG(gkl(j)) * cavsig2(k, at, ef) + &
               & -CONJG(gkl(k)) * cavsig2(j, at, ef)
      B_vec(8) = 0.0d0
      ! Calculate k1
      k1_cavsig3(j, k, at, :) = dt * (MATMUL(Mat, cavsig3(j, k, at, :)) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j)) * (cavsig2(k, at, eg) + 0.5d0 * k1_cavsig2(k, at, eg)) + &
               & -CONJG(gkl(k)) * (cavsig2(j, at, eg) + 0.5d0 * k1_cavsig2(j, at, eg))
      B_vec(2) = -CONJG(gkl(j)) * (cavsig2(k, at, ee) + 0.5d0 * k1_cavsig2(k, at, ee)) + &
               & -CONJG(gkl(k)) * (cavsig2(j, at, ee) + 0.5d0 * k1_cavsig2(j, at, ee))
      B_vec(3) = -CONJG(gkl(j)) * xi * (cavsig2(k, at, fg) + 0.5d0 * k1_cavsig2(k, at, fg)) + &
               & -CONJG(gkl(k)) * xi * (cavsig2(j, at, fg) + 0.5d0 * k1_cavsig2(j, at, fg))
      B_vec(4) = Gamma * (xi ** 2) * (cav2(j, k, at) + 0.5d0 * k1_cav2(j, k, at)) + &
               & -CONJG(gkl(j)) * xi * (cavsig2(k, at, fe) + 0.5d0 * k1_cavsig2(k, at, fe)) + &
               & -CONJG(gkl(k)) * xi * (cavsig2(j, at, fe) + 0.5d0 * k1_cavsig2(j, at, fe))
      B_vec(5) = i * xi * 0.5d0 * Omega * (cav2(j, k, at) + 0.5d0 * k1_cav2(j, k, at)) + &
               & CONJG(gkl(j)) * xi * (cavsig2(k, at, gg) + 0.5d0 * k1_cavsig2(k, at, gg)) + &
               & CONJG(gkl(j)) * xi * (cavsig2(k, at, ee) + 0.5d0 * k1_cavsig2(k, at, ee)) + &
               & -CONJG(gkl(j)) * xi * (cav1(k, at) + 0.5d0 * k1_cav1(k, at)) + &
               & CONJG(gkl(k)) * xi * (cavsig2(j, at, gg) + 0.5d0 * k1_cavsig2(j, at, gg)) + &
               & CONJG(gkl(k)) * xi * (cavsig2(j, at, ee) + 0.5d0 * k1_cavsig2(j, at, ee)) + &
               & -CONJG(gkl(k)) * xi * (cav1(j, at) + 0.5d0 * k1_cav1(j, at))
      B_vec(6) = -i * xi * 0.5d0 * Omega * (cav2(j, k, at) + 0.5d0 * k1_cav2(j, k, at))
      B_vec(7) = -CONJG(gkl(j)) * (cavsig2(k, at, ef) + 0.5d0 * k1_cavsig2(k, at, ef)) + &
           & -CONJG(gkl(k)) * (cavsig2(j, at, ef) + 0.5d0 * k1_cavsig2(j, at, ef))
      B_vec(8) = 0.0d0
      ! Calculate k2
      k2_cavsig3(j, k, at, :) = dt * (MATMUL(Mat, (cavsig3(j, k, at, :) + 0.5d0 * k1_cavsig3(j, k, at, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j)) * (cavsig2(k, at, eg) + 0.5d0 * k2_cavsig2(k, at, eg)) + &
               & -CONJG(gkl(k)) * (cavsig2(j, at, eg) + 0.5d0 * k2_cavsig2(j, at, eg))
      B_vec(2) = -CONJG(gkl(j)) * (cavsig2(k, at, ee) + 0.5d0 * k2_cavsig2(k, at, ee)) + &
               & -CONJG(gkl(k)) * (cavsig2(j, at, ee) + 0.5d0 * k2_cavsig2(j, at, ee))
      B_vec(3) = -CONJG(gkl(j)) * xi * (cavsig2(k, at, fg) + 0.5d0 * k2_cavsig2(k, at, fg)) + &
               & -CONJG(gkl(k)) * xi * (cavsig2(j, at, fg) + 0.5d0 * k2_cavsig2(j, at, fg))
      B_vec(4) = Gamma * (xi ** 2) * (cav2(j, k, at) + 0.5d0 * k2_cav2(j, k, at)) + &
               & -CONJG(gkl(j)) * xi * (cavsig2(k, at, fe) + 0.5d0 * k2_cavsig2(k, at, fe)) + &
               & -CONJG(gkl(k)) * xi * (cavsig2(j, at, fe) + 0.5d0 * k2_cavsig2(j, at, fe))
      B_vec(5) = i * xi * 0.5d0 * Omega * (cav2(j, k, at) + 0.5d0 * k2_cav2(j, k, at)) + &
               & CONJG(gkl(j)) * xi * (cavsig2(k, at, gg) + 0.5d0 * k2_cavsig2(k, at, gg)) + &
               & CONJG(gkl(j)) * xi * (cavsig2(k, at, ee) + 0.5d0 * k2_cavsig2(k, at, ee)) + &
               & -CONJG(gkl(j)) * xi * (cav1(k, at) + 0.5d0 * k2_cav1(k, at)) + &
               & CONJG(gkl(k)) * xi * (cavsig2(j, at, gg) + 0.5d0 * k2_cavsig2(j, at, gg)) + &
               & CONJG(gkl(k)) * xi * (cavsig2(j, at, ee) + 0.5d0 * k2_cavsig2(j, at, ee)) + &
               & -CONJG(gkl(k)) * xi * (cav1(j, at) + 0.5d0 * k2_cav1(j, at))
      B_vec(6) = -i * xi * 0.5d0 * Omega * (cav2(j, k, at) + 0.5d0 * k2_cav2(j, k, at))
      B_vec(7) = -CONJG(gkl(j)) * (cavsig2(k, at, ef) + 0.5d0 * k2_cavsig2(k, at, ef)) + &
               & -CONJG(gkl(k)) * (cavsig2(j, at, ef) + 0.5d0 * k2_cavsig2(j, at, ef))
      B_vec(8) = 0.0d0
      ! Calculate k3
      k3_cavsig3(j, k, at, :) = dt * (MATMUL(Mat, (cavsig3(j, k, at, :) + 0.5d0 * k2_cavsig3(j, k, at, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j)) * (cavsig2(k, at, eg) + k3_cavsig2(k, at, eg)) + &
               & -CONJG(gkl(k)) * (cavsig2(j, at, eg) + k3_cavsig2(j, at, eg))
      B_vec(2) = -CONJG(gkl(j)) * (cavsig2(k, at, ee) + k3_cavsig2(k, at, ee)) + &
               & -CONJG(gkl(k)) * (cavsig2(j, at, ee) + k3_cavsig2(j, at, ee))
      B_vec(3) = -CONJG(gkl(j)) * xi * (cavsig2(k, at, fg) + k3_cavsig2(k, at, fg)) + &
               & -CONJG(gkl(k)) * xi * (cavsig2(j, at, fg) + k3_cavsig2(j, at, fg))
      B_vec(4) = Gamma * (xi ** 2) * (cav2(j, k, at) + k3_cav2(j, k, at)) + &
               & -CONJG(gkl(j)) * xi * (cavsig2(k, at, fe) + k3_cavsig2(k, at, fe)) + &
               & -CONJG(gkl(k)) * xi * (cavsig2(j, at, fe) + k3_cavsig2(j, at, fe))
      B_vec(5) = i * xi * 0.5d0 * Omega * (cav2(j, k, at) + k3_cav2(j, k, at)) + &
               & CONJG(gkl(j)) * xi * (cavsig2(k, at, gg) + k3_cavsig2(k, at, gg)) + &
               & CONJG(gkl(j)) * xi * (cavsig2(k, at, ee) + k3_cavsig2(k, at, ee)) + &
               & -CONJG(gkl(j)) * xi * (cav1(k, at) + k3_cav1(k, at)) + &
               & CONJG(gkl(k)) * xi * (cavsig2(j, at, gg) + k3_cavsig2(j, at, gg)) + &
               & CONJG(gkl(k)) * xi * (cavsig2(j, at, ee) + k3_cavsig2(j, at, ee)) + &
               & -CONJG(gkl(k)) * xi * (cav1(j, at) + k3_cav1(j, at))
      B_vec(6) = -i * xi * 0.5d0 * Omega * (cav2(j, k, at) + k3_cav2(j, k, at))
      B_vec(7) = -CONJG(gkl(j)) * (cavsig2(k, at, ef) + k3_cavsig2(k, at, ef)) + &
               & -CONJG(gkl(k)) * (cavsig2(j, at, ef) + k3_cavsig2(j, at, ef))
      B_vec(8) = 0.0d0
      ! Calculate k4
      k4_cavsig3(j, k, at, :) = dt * (MATMUL(Mat, (cavsig3(j, k, at, :) + k3_cavsig3(j, k, at, :))) + B_vec)

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
      B_vec(1) = -CONJG(gkl(j)) * cavsig2(k, a, eg) + &
               & -gkl(k) * cavsig2(j, at, ge)
      B_vec(2) = -CONJG(gkl(j)) * cavsig2(k, a, ee) + &
               & -gkl(k) * xi * cavsig2(j, at, gf)
      B_vec(3) = -CONJG(gkl(j)) * xi * cavsig2(k, a, fg) + &
               & -gkl(k) * cavsig2(j, at, ee)
      B_vec(4) = Gamma * (xi ** 2) * cav2(j, k, ata) + &
               & -CONJG(gkl(j)) * xi * cavsig2(k, a, fe) + &
               & -gkl(k) * xi * cavsig2(j, at, ef)
      B_vec(5) = i * xi * 0.5d0 * Omega * cav2(j, k, ata) + &
               & CONJG(gkl(j)) * xi * cavsig2(k, a, gg) + &
               & CONJG(gkl(j)) * xi * cavsig2(k, a, ee) + &
               & -CONJG(gkl(j)) * xi * cav1(k, a)
      B_vec(6) = -i * xi * 0.5d0 * Omega * cav2(j, k, ata) + &
               & gkl(k) * xi * cavsig2(j, at, gg) + &
               & gkl(k) * xi * cavsig2(j, at, ee) + &
               & -gkl(k) * xi * cav1(j, at)
      B_vec(7) = -CONJG(gkl(j)) * cavsig2(k, a, ef)
      B_vec(8) = - gkl(k) * cavsig2(j, at, fe)
      ! Calculate k1
      k1_cavsig3(j, k, ata, :) = dt * (MATMUL(Mat, cavsig3(j, k, ata, :)) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j)) * (cavsig2(k, a, eg) + 0.5d0 * k1_cavsig2(k, a, eg)) + &
               & -gkl(k) * (cavsig2(j, at, ge) + 0.5d0 * k1_cavsig2(j, at, ge))
      B_vec(2) = -CONJG(gkl(j)) * (cavsig2(k, a, ee) + 0.5d0 * k1_cavsig2(k, a, ee)) + &
               & -gkl(k) * xi * (cavsig2(j, at, gf) + 0.5d0 * k1_cavsig2(j, at, gf))
      B_vec(3) = -CONJG(gkl(j)) * xi * (cavsig2(k, a, fg) + 0.5d0 * k1_cavsig2(k, a, fg)) + &
               & -gkl(k) * (cavsig2(j, at, ee) + 0.5d0 * k1_cavsig2(j, at, ee))
      B_vec(4) = Gamma * (xi ** 2) * (cav2(j, k, ata) + 0.5d0 * k1_cav2(j, k, ata)) + &
               & -CONJG(gkl(j)) * xi * (cavsig2(k, a, fe) + 0.5d0 * k1_cavsig2(k, a, fe)) + &
               & -gkl(k) * xi * (cavsig2(j, at, ef) + 0.5d0 * k1_cavsig2(j, at, ef))
      B_vec(5) = i * xi * 0.5d0 * Omega * (cav2(j, k, ata) + 0.5d0 * k1_cav2(j, k, ata)) + &
               & CONJG(gkl(j)) * xi * (cavsig2(k, a, gg) + 0.5d0 * k1_cavsig2(k, a, gg)) + &
               & CONJG(gkl(j)) * xi * (cavsig2(k, a, ee) + 0.5d0 * k1_cavsig2(k, a, ee)) + &
               & -CONJG(gkl(j)) * xi * (cav1(k, a) + 0.5d0 * k1_cav1(k, a))
      B_vec(6) = -i * xi * 0.5d0 * Omega * (cav2(j, k, ata) + 0.5d0 * k1_cav2(j, k, ata)) + &
               & gkl(k) * xi * (cavsig2(j, at, gg) + 0.5d0 * k1_cavsig2(j, at, gg)) + &
               & gkl(k) * xi * (cavsig2(j, at, ee) + 0.5d0 * k1_cavsig2(j, at, ee)) + &
               & -gkl(k) * xi * (cav1(j, at) + 0.5d0 * k1_cav1(j, at))
      B_vec(7) = -CONJG(gkl(j)) * (cavsig2(k, a, ef) + 0.5d0 * k1_cavsig2(k, a, ef))
      B_vec(8) = - gkl(k) * (cavsig2(j, at, fe) + 0.5d0 * k1_cavsig2(j, at, fe))
      ! Calculate k2
      k2_cavsig3(j, k, ata, :) = dt * (MATMUL(Mat, (cavsig3(j, k, ata, :) + 0.5d0 * k1_cavsig3(j, k, ata, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j)) * (cavsig2(k, a, eg) + 0.5d0 * k2_cavsig2(k, a, eg)) + &
               & -gkl(k) * (cavsig2(j, at, ge) + 0.5d0 * k2_cavsig2(j, at, ge))
      B_vec(2) = -CONJG(gkl(j)) * (cavsig2(k, a, ee) + 0.5d0 * k2_cavsig2(k, a, ee)) + &
               & -gkl(k) * xi * (cavsig2(j, at, gf) + 0.5d0 * k2_cavsig2(j, at, gf))
      B_vec(3) = -CONJG(gkl(j)) * xi * (cavsig2(k, a, fg) + 0.5d0 * k2_cavsig2(k, a, fg)) + &
               & -gkl(k) * (cavsig2(j, at, ee) + 0.5d0 * k2_cavsig2(j, at, ee))
      B_vec(4) = Gamma * (xi ** 2) * (cav2(j, k, ata) + 0.5d0 * k2_cav2(j, k, ata)) + &
               & -CONJG(gkl(j)) * xi * (cavsig2(k, a, fe) + 0.5d0 * k2_cavsig2(k, a, fe)) + &
               & -gkl(k) * xi * (cavsig2(j, at, ef) + 0.5d0 * k2_cavsig2(j, at, ef))
      B_vec(5) = i * xi * 0.5d0 * Omega * (cav2(j, k, ata) + 0.5d0 * k2_cav2(j, k, ata)) + &
               & CONJG(gkl(j)) * xi * (cavsig2(k, a, gg) + 0.5d0 * k2_cavsig2(k, a, gg)) + &
               & CONJG(gkl(j)) * xi * (cavsig2(k, a, ee) + 0.5d0 * k2_cavsig2(k, a, ee)) + &
               & -CONJG(gkl(j)) * xi * (cav1(k, a) + 0.5d0 * k2_cav1(k, a))
      B_vec(6) = -i * xi * 0.5d0 * Omega * (cav2(j, k, ata) + 0.5d0 * k2_cav2(j, k, ata)) + &
               & gkl(k) * xi * (cavsig2(j, at, gg) + 0.5d0 * k2_cavsig2(j, at, gg)) + &
               & gkl(k) * xi * (cavsig2(j, at, ee) + 0.5d0 * k2_cavsig2(j, at, ee)) + &
               & -gkl(k) * xi * (cav1(j, at) + 0.5d0 * k2_cav1(j, at))
      B_vec(7) = -CONJG(gkl(j)) * (cavsig2(k, a, ef) + 0.5d0 * k2_cavsig2(k, a, ef))
      B_vec(8) = - gkl(k) * (cavsig2(j, at, fe) + 0.5d0 * k2_cavsig2(j, at, fe))
      ! Calculate k3
      k3_cavsig3(j, k, ata, :) = dt * (MATMUL(Mat, (cavsig3(j, k, ata, :) + 0.5d0 * k2_cavsig3(j, k, ata, :))) + B_vec)

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -CONJG(gkl(j)) * (cavsig2(k, a, eg) + k3_cavsig2(k, a, eg)) + &
               & -gkl(k) * (cavsig2(j, at, ge) + k3_cavsig2(j, at, ge))
      B_vec(2) = -CONJG(gkl(j)) * (cavsig2(k, a, ee) + k3_cavsig2(k, a, ee)) + &
               & -gkl(k) * xi * (cavsig2(j, at, gf) + k3_cavsig2(j, at, gf))
      B_vec(3) = -CONJG(gkl(j)) * xi * (cavsig2(k, a, fg) + k3_cavsig2(k, a, fg)) + &
               & -gkl(k) * (cavsig2(j, at, ee) + k3_cavsig2(j, at, ee))
      B_vec(4) = Gamma * (xi ** 2) * (cav2(j, k, ata) + k3_cav2(j, k, ata)) + &
               & -CONJG(gkl(j)) * xi * (cavsig2(k, a, fe) + k3_cavsig2(k, a, fe)) + &
               & -gkl(k) * xi * (cavsig2(j, at, ef) + k3_cavsig2(j, at, ef))
      B_vec(5) = i * xi * 0.5d0 * Omega * (cav2(j, k, ata) + k3_cav2(j, k, ata)) + &
               & CONJG(gkl(j)) * xi * (cavsig2(k, a, gg) + k3_cavsig2(k, a, gg)) + &
               & CONJG(gkl(j)) * xi * (cavsig2(k, a, ee) + k3_cavsig2(k, a, ee)) + &
               & -CONJG(gkl(j)) * xi * (cav1(k, a) + k3_cav1(k, a))
      B_vec(6) = -i * xi * 0.5d0 * Omega * (cav2(j, k, ata) + k3_cav2(j, k, ata)) + &
               & gkl(k) * xi * (cavsig2(j, at, gg) + k3_cavsig2(j, at, gg)) + &
               & gkl(k) * xi * (cavsig2(j, at, ee) + k3_cavsig2(j, at, ee)) + &
               & -gkl(k) * xi * (cav1(j, at) + k3_cav1(j, at))
      B_vec(7) = -CONJG(gkl(j)) * (cavsig2(k, a, ef) + k3_cavsig2(k, a, ef))
      B_vec(8) = - gkl(k) * (cavsig2(j, at, fe) + k3_cavsig2(j, at, fe))
      ! Calculate k4
      k4_cavsig3(j, k, ata, :) = dt * (MATMUL(Mat, (cavsig3(j, k, ata, :) + k3_cavsig3(j, k, ata, :))) + B_vec)

      ! Close k loop
    END DO
    ! Close j loop
  END DO

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
        k1_cav3(j, k, l, a) = -dt * (3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * cav3(j, k, l, a) + &
                            & -dt * CONJG(gkl(j)) * cavsig3(k, l, a, eg) + &
                            & -dt * CONJG(gkl(j)) * xi * cavsig3(k, l, a, fe) + &
                            & -dt * gkl(k) * cavsig3(j, l, ata, ge) + &
                            & -dt * gkl(k) * xi * cavsig3(j, l, ata, ef) + &
                            & -dt * gkl(l) * cavsig3(j, k, ata, ge) + &
                            & -dt * gkl(l) * xi * cavsig3(j, k, ata, ef)

        k2_cav3(j, k, l, a) = -dt * (3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cav3(j, k, l, a) + 0.5d0 * k1_cav3(j, k, l, a)) + &
                            & -dt * CONJG(gkl(j)) * (cavsig3(k, l, a, eg) + 0.5d0 * k1_cavsig3(k, l, a, eg)) + &
                            & -dt * CONJG(gkl(j)) * xi * (cavsig3(k, l, a, fe) + 0.5d0 * k1_cavsig3(k, l, a, fe)) + &
                            & -dt * gkl(k) * (cavsig3(j, l, ata, ge) + 0.5d0 * k1_cavsig3(j, l, ata, ge)) + &
                            & -dt * gkl(k) * xi * (cavsig3(j, l, ata, ef) + 0.5d0 * k1_cavsig3(j, l, ata, ef)) + &
                            & -dt * gkl(l) * (cavsig3(j, k, ata, ge) + 0.5d0 * k1_cavsig3(j, k, ata, ge)) + &
                            & -dt * gkl(l) * xi * (cavsig3(j, k, ata, ef) + 0.5d0 * k1_cavsig3(j, k, ata, ef))

        k3_cav3(j, k, l, a) = -dt * (3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cav3(j, k, l, a) + 0.5d0 * k2_cav3(j, k, l, a)) + &
                            & -dt * CONJG(gkl(j)) * (cavsig3(k, l, a, eg) + 0.5d0 * k2_cavsig3(k, l, a, eg)) + &
                            & -dt * CONJG(gkl(j)) * xi * (cavsig3(k, l, a, fe) + 0.5d0 * k2_cavsig3(k, l, a, fe)) + &
                            & -dt * gkl(k) * (cavsig3(j, l, ata, ge) + 0.5d0 * k2_cavsig3(j, l, ata, ge)) + &
                            & -dt * gkl(k) * xi * (cavsig3(j, l, ata, ef) + 0.5d0 * k2_cavsig3(j, l, ata, ef)) + &
                            & -dt * gkl(l) * (cavsig3(j, k, ata, ge) + 0.5d0 * k2_cavsig3(j, k, ata, ge)) + &
                            & -dt * gkl(l) * xi * (cavsig3(j, k, ata, ef) + 0.5d0 * k2_cavsig3(j, k, ata, ef))

        k4_cav3(j, k, l, a) = -dt * (3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cav3(j, k, l, a) + k3_cav3(j, k, l, a)) + &
                            & -dt * CONJG(gkl(j)) * (cavsig3(k, l, a, eg) + k3_cavsig3(k, l, a, eg)) + &
                            & -dt * CONJG(gkl(j)) * xi * (cavsig3(k, l, a, fe) + k3_cavsig3(k, l, a, fe)) + &
                            & -dt * gkl(k) * (cavsig3(j, l, ata, ge) + k3_cavsig3(j, l, ata, ge)) + &
                            & -dt * gkl(k) * xi * (cavsig3(j, l, ata, ef) + k3_cavsig3(j, l, ata, ef)) + &
                            & -dt * gkl(l) * (cavsig3(j, k, ata, ge) + k3_cavsig3(j, k, ata, ge)) + &
                            & -dt * gkl(l) * xi * (cavsig3(j, k, ata, ef) + k3_cavsig3(j, k, ata, ef))

        !-------------------------------------------!
        ! < a^{\dagger}_{j} a^{\dagger}_{k} a_{l} > !
        !-------------------------------------------!
        k1_cav3(j, k, l, at) = -dt * (3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * cav3(j, k, l, at) + &
                             & -dt * CONJG(gkl(j)) * cavsig3(k, l, ata, eg) + &
                             & -dt * xi * CONJG(gkl(j)) * cavsig3(k, l, ata, fe) + &
                             & -dt * CONJG(gkl(k)) * cavsig3(j, l, ata, eg) + &
                             & -dt * xi * CONJG(gkl(k)) * cavsig3(j, l, ata, fe) + &
                             & -dt * gkl(l) * cavsig3(j, k, at, ge) + &
                             & -dt * xi * gkl(l) * cavsig3(j, k, at, ef)

        k2_cav3(j, k, l, at) = -dt * (3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cav3(j, k, l, at) + 0.5d0 * k1_cav3(j, k, l, at)) + &
                             & -dt * CONJG(gkl(j)) * (cavsig3(k, l, ata, eg) + 0.5d0 * k1_cavsig3(k, l, ata, eg)) + &
                             & -dt * xi * CONJG(gkl(j)) * (cavsig3(k, l, ata, fe) + 0.5d0 * k1_cavsig3(k, l, ata, fe)) + &
                             & -dt * CONJG(gkl(k)) * (cavsig3(j, l, ata, eg) + 0.5d0 * k1_cavsig3(j, l, ata, eg)) + &
                             & -dt * xi * CONJG(gkl(k)) * (cavsig3(j, l, ata, fe) + 0.5d0 * k1_cavsig3(j, l, ata, fe)) + &
                             & -dt * gkl(l) * (cavsig3(j, k, at, ge) + 0.5d0 * k1_cavsig3(j, k, at, ge)) + &
                             & -dt * xi * gkl(l) * (cavsig3(j, k, at, ef) + 0.5d0 * k1_cavsig3(j, k, at, ef))

        k3_cav3(j, k, l, at) = -dt * (3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cav3(j, k, l, at) + 0.5d0 * k2_cav3(j, k, l, at)) + &
                             & -dt * CONJG(gkl(j)) * (cavsig3(k, l, ata, eg) + 0.5d0 * k2_cavsig3(k, l, ata, eg)) + &
                             & -dt * xi * CONJG(gkl(j)) * (cavsig3(k, l, ata, fe) + 0.5d0 * k2_cavsig3(k, l, ata, fe)) + &
                             & -dt * CONJG(gkl(k)) * (cavsig3(j, l, ata, eg) + 0.5d0 * k2_cavsig3(j, l, ata, eg)) + &
                             & -dt * xi * CONJG(gkl(k)) * (cavsig3(j, l, ata, fe) + 0.5d0 * k2_cavsig3(j, l, ata, fe)) + &
                             & -dt * gkl(l) * (cavsig3(j, k, at, ge) + 0.5d0 * k2_cavsig3(j, k, at, ge)) + &
                             & -dt * xi * gkl(l) * (cavsig3(j, k, at, ef) + 0.5d0 * k2_cavsig3(j, k, at, ef))

        k4_cav3(j, k, l, at) = -dt * (3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cav3(j, k, l, at) + k3_cav3(j, k, l, at)) + &
                             & -dt * CONJG(gkl(j)) * (cavsig3(k, l, ata, eg) + k3_cavsig3(k, l, ata, eg)) + &
                             & -dt * xi * CONJG(gkl(j)) * (cavsig3(k, l, ata, fe) + k3_cavsig3(k, l, ata, fe)) + &
                             & -dt * CONJG(gkl(k)) * (cavsig3(j, l, ata, eg) + k3_cavsig3(j, l, ata, eg)) + &
                             & -dt * xi * CONJG(gkl(k)) * (cavsig3(j, l, ata, fe) + k3_cavsig3(j, l, ata, fe)) + &
                             & -dt * gkl(l) * (cavsig3(j, k, at, ge) + k3_cavsig3(j, k, at, ge)) + &
                             & -dt * xi * gkl(l) * (cavsig3(j, k, at, ef) + k3_cavsig3(j, k, at, ef))

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
        B_vec(1) = -CONJG(gkl(j)) * cavsig3(k, l, a, eg) + &
                 & -gkl(k) * cavsig3(j, l, ata, ge) + &
                 & -gkl(l) * cavsig3(j, k, ata, ge)
        B_vec(2) = -CONJG(gkl(j)) * cavsig3(k, l, a, ee) + &
                 & -gkl(k) * xi * cavsig3(j, l, ata, gf) + &
                 & -gkl(l) * xi * cavsig3(j, k, ata, gf)
        B_vec(3) = -CONJG(gkl(j)) * xi * cavsig3(k, l, a, fg) + &
                 & -gkl(k) * cavsig3(j, l, ata, ee) + &
                 & -gkl(l) * cavsig3(j, k, ata, ee)
        B_vec(4) = Gamma * (xi ** 2) * cav3(j, k, l, a) + &
                 & -CONJG(gkl(j)) * xi * cavsig3(k, l, a, fe) + &
                 & -gkl(k) * xi * cavsig3(j, l, ata, ef) + &
                 & -gkl(l) * xi * cavsig3(j, k, ata, ef)
        B_vec(5) = i * xi * 0.5d0 * Omega * cav3(j, k, l, a) + &
                 & CONJG(gkl(j)) * xi * cavsig3(k, l, a, gg) + &
                 & CONJG(gkl(j)) * xi * cavsig3(k, l, a, ee) + &
                 & -CONJG(gkl(j)) * xi * cav2(k, l, a)
        B_vec(6) = -i * xi * 0.5d0 * Omega * cav3(j, k, l, a) + &
                 & gkl(k) * xi * cavsig3(j, l, ata, gg) + &
                 & gkl(k) * xi * cavsig3(j, l, ata, ee) + &
                 & -gkl(k) * xi * cav2(j, l, ata) + &
                 & gkl(l) * xi * cavsig3(j, k, ata, gg) + &
                 & gkl(l) * xi * cavsig3(j, k, ata, ee) + &
                 & -gkl(l) * xi * cav2(j, k, ata)
        B_vec(7) = -CONJG(gkl(j)) * cavsig3(k, l, a, ef)
        B_vec(8) = -gkl(k) * cavsig3(j, l, ata, fe) + &
                 & -gkl(l) * cavsig3(j, k, ata, fe)
        ! Calculate k1
        k1_cavsig4(j, k, l, a, :) = dt * (MATMUL(Mat, cavsig4(j, k, l, a, :)) + B_vec)

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j)) * (cavsig3(k, l, a, eg) + 0.5d0 * k1_cavsig3(k, l, a, eg)) + &
                 & -gkl(k) * (cavsig3(j, l, ata, ge) + 0.5d0 * k1_cavsig3(j, l, ata, ge)) + &
                 & -gkl(l) * (cavsig3(j, k, ata, ge) + 0.5d0 * k1_cavsig3(j, k, ata, ge))
        B_vec(2) = -CONJG(gkl(j)) * (cavsig3(k, l, a, ee) + 0.5d0 * k1_cavsig3(k, l, a, ee)) + &
                 & -gkl(k) * xi * (cavsig3(j, l, ata, gf) + 0.5d0 * k1_cavsig3(j, l, ata, gf)) + &
                 & -gkl(l) * xi * (cavsig3(j, k, ata, gf) + 0.5d0 * k1_cavsig3(j, k, ata, gf))
        B_vec(3) = -CONJG(gkl(j)) * xi * (cavsig3(k, l, a, fg) + 0.5d0 * k1_cavsig3(k, l, a, fg)) + &
                 & -gkl(k) * (cavsig3(j, l, ata, ee) + 0.5d0 * k1_cavsig3(j, l, ata, ee)) + &
                 & -gkl(l) * (cavsig3(j, k, ata, ee) + 0.5d0 * k1_cavsig3(j, k, ata, ee))
        B_vec(4) = Gamma * (xi ** 2) * (cav3(j, k, l, a) + 0.5d0 * k1_cav3(j, k, l, a)) + &
                 & -CONJG(gkl(j)) * xi * (cavsig3(k, l, a, fe) + 0.5d0 * k1_cavsig3(k, l, a, fe)) + &
                 & -gkl(k) * xi * (cavsig3(j, l, ata, ef) + 0.5d0 * k1_cavsig3(j, l, ata, ef)) + &
                 & -gkl(l) * xi * (cavsig3(j, k, ata, ef) + 0.5d0 * k1_cavsig3(j, k, ata, ef))
        B_vec(5) = i * xi * 0.5d0 * Omega * (cav3(j, k, l, a) + 0.5d0 * k1_cav3(j, k, l, a)) + &
                 & CONJG(gkl(j)) * xi * (cavsig3(k, l, a, gg) + 0.5d0 * k1_cavsig3(k, l, a, gg)) + &
                 & CONJG(gkl(j)) * xi * (cavsig3(k, l, a, ee) + 0.5d0 * k1_cavsig3(k, l, a, ee)) + &
                 & -CONJG(gkl(j)) * xi * (cav2(k, l, a) + 0.5d0 * k1_cav2(k, l, a))
        B_vec(6) = -i * xi * 0.5d0 * Omega * (cav3(j, k, l, a) + 0.5d0 * k1_cav3(j, k, l, a)) + &
                 & gkl(k) * xi * (cavsig3(j, l, ata, gg) + 0.5d0 * k1_cavsig3(j, l, ata, gg)) + &
                 & gkl(k) * xi * (cavsig3(j, l, ata, ee) + 0.5d0 * k1_cavsig3(j, l, ata, ee)) + &
                 & -gkl(k) * xi * (cav2(j, l, ata) + 0.5d0 * k1_cav2(j, l, ata)) + &
                 & gkl(l) * xi * (cavsig3(j, k, ata, gg) + 0.5d0 * k1_cavsig3(j, k, ata, gg)) + &
                 & gkl(l) * xi * (cavsig3(j, k, ata, ee) + 0.5d0 * k1_cavsig3(j, k, ata, ee)) + &
                 & -gkl(l) * xi * (cav2(j, k, ata) + 0.5d0 * k1_cav2(j, k, ata))
        B_vec(7) = -CONJG(gkl(j)) * (cavsig3(k, l, a, ef) + 0.5d0 * k1_cavsig3(k, l, a, ef))
        B_vec(8) = -gkl(k) * (cavsig3(j, l, ata, fe) + 0.5d0 * k1_cavsig3(j, l, ata, fe)) + &
                 & -gkl(l) * (cavsig3(j, k, ata, fe) + 0.5d0 * k1_cavsig3(j, k, ata, fe))
        ! Calculate k2
        k2_cavsig4(j, k, l, a, :) = dt * (MATMUL(Mat, (cavsig4(j, k, l, a, :) + 0.5d0 * k1_cavsig4(j, k, l, a, :))) + B_vec)

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j)) * (cavsig3(k, l, a, eg) + 0.5d0 * k2_cavsig3(k, l, a, eg)) + &
                 & -gkl(k) * (cavsig3(j, l, ata, ge) + 0.5d0 * k2_cavsig3(j, l, ata, ge)) + &
                 & -gkl(l) * (cavsig3(j, k, ata, ge) + 0.5d0 * k2_cavsig3(j, k, ata, ge))
        B_vec(2) = -CONJG(gkl(j)) * (cavsig3(k, l, a, ee) + 0.5d0 * k2_cavsig3(k, l, a, ee)) + &
                 & -gkl(k) * xi * (cavsig3(j, l, ata, gf) + 0.5d0 * k2_cavsig3(j, l, ata, gf)) + &
                 & -gkl(l) * xi * (cavsig3(j, k, ata, gf) + 0.5d0 * k2_cavsig3(j, k, ata, gf))
        B_vec(3) = -CONJG(gkl(j)) * xi * (cavsig3(k, l, a, fg) + 0.5d0 * k2_cavsig3(k, l, a, fg)) + &
                 & -gkl(k) * (cavsig3(j, l, ata, ee) + 0.5d0 * k2_cavsig3(j, l, ata, ee)) + &
                 & -gkl(l) * (cavsig3(j, k, ata, ee) + 0.5d0 * k2_cavsig3(j, k, ata, ee))
        B_vec(4) = Gamma * (xi ** 2) * (cav3(j, k, l, a) + 0.5d0 * k2_cav3(j, k, l, a)) + &
                 & -CONJG(gkl(j)) * xi * (cavsig3(k, l, a, fe) + 0.5d0 * k2_cavsig3(k, l, a, fe)) + &
                 & -gkl(k) * xi * (cavsig3(j, l, ata, ef) + 0.5d0 * k2_cavsig3(j, l, ata, ef)) + &
                 & -gkl(l) * xi * (cavsig3(j, k, ata, ef) + 0.5d0 * k2_cavsig3(j, k, ata, ef))
        B_vec(5) = i * xi * 0.5d0 * Omega * (cav3(j, k, l, a) + 0.5d0 * k2_cav3(j, k, l, a)) + &
                 & CONJG(gkl(j)) * xi * (cavsig3(k, l, a, gg) + 0.5d0 * k2_cavsig3(k, l, a, gg)) + &
                 & CONJG(gkl(j)) * xi * (cavsig3(k, l, a, ee) + 0.5d0 * k2_cavsig3(k, l, a, ee)) + &
                 & -CONJG(gkl(j)) * xi * (cav2(k, l, a) + 0.5d0 * k2_cav2(k, l, a))
        B_vec(6) = -i * xi * 0.5d0 * Omega * (cav3(j, k, l, a) + 0.5d0 * k2_cav3(j, k, l, a)) + &
                 & gkl(k) * xi * (cavsig3(j, l, ata, gg) + 0.5d0 * k2_cavsig3(j, l, ata, gg)) + &
                 & gkl(k) * xi * (cavsig3(j, l, ata, ee) + 0.5d0 * k2_cavsig3(j, l, ata, ee)) + &
                 & -gkl(k) * xi * (cav2(j, l, ata) + 0.5d0 * k2_cav2(j, l, ata)) + &
                 & gkl(l) * xi * (cavsig3(j, k, ata, gg) + 0.5d0 * k2_cavsig3(j, k, ata, gg)) + &
                 & gkl(l) * xi * (cavsig3(j, k, ata, ee) + 0.5d0 * k2_cavsig3(j, k, ata, ee)) + &
                 & -gkl(l) * xi * (cav2(j, k, ata) + 0.5d0 * k2_cav2(j, k, ata))
        B_vec(7) = -CONJG(gkl(j)) * (cavsig3(k, l, a, ef) + 0.5d0 * k2_cavsig3(k, l, a, ef))
        B_vec(8) = -gkl(k) * (cavsig3(j, l, ata, fe) + 0.5d0 * k2_cavsig3(j, l, ata, fe)) + &
                 & -gkl(l) * (cavsig3(j, k, ata, fe) + 0.5d0 * k2_cavsig3(j, k, ata, fe))
        ! Calculate k3
        k3_cavsig4(j, k, l, a, :) = dt * (MATMUL(Mat, (cavsig4(j, k, l, a, :) + 0.5d0 * k2_cavsig4(j, k, l, a, :))) + B_vec)

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j)) * (cavsig3(k, l, a, eg) + k3_cavsig3(k, l, a, eg)) + &
                 & -gkl(k) * (cavsig3(j, l, ata, ge) + k3_cavsig3(j, l, ata, ge)) + &
                 & -gkl(l) * (cavsig3(j, k, ata, ge) + k3_cavsig3(j, k, ata, ge))
        B_vec(2) = -CONJG(gkl(j)) * (cavsig3(k, l, a, ee) + k3_cavsig3(k, l, a, ee)) + &
                 & -gkl(k) * xi * (cavsig3(j, l, ata, gf) + k3_cavsig3(j, l, ata, gf)) + &
                 & -gkl(l) * xi * (cavsig3(j, k, ata, gf) + k3_cavsig3(j, k, ata, gf))
        B_vec(3) = -CONJG(gkl(j)) * xi * (cavsig3(k, l, a, fg) + k3_cavsig3(k, l, a, fg)) + &
                 & -gkl(k) * (cavsig3(j, l, ata, ee) + k3_cavsig3(j, l, ata, ee)) + &
                 & -gkl(l) * (cavsig3(j, k, ata, ee) + k3_cavsig3(j, k, ata, ee))
        B_vec(4) = Gamma * (xi ** 2) * (cav3(j, k, l, a) + k3_cav3(j, k, l, a)) + &
                 & -CONJG(gkl(j)) * xi * (cavsig3(k, l, a, fe) + k3_cavsig3(k, l, a, fe)) + &
                 & -gkl(k) * xi * (cavsig3(j, l, ata, ef) + k3_cavsig3(j, l, ata, ef)) + &
                 & -gkl(l) * xi * (cavsig3(j, k, ata, ef) + k3_cavsig3(j, k, ata, ef))
        B_vec(5) = i * xi * 0.5d0 * Omega * (cav3(j, k, l, a) + k3_cav3(j, k, l, a)) + &
                 & CONJG(gkl(j)) * xi * (cavsig3(k, l, a, gg) + k3_cavsig3(k, l, a, gg)) + &
                 & CONJG(gkl(j)) * xi * (cavsig3(k, l, a, ee) + k3_cavsig3(k, l, a, ee)) + &
                 & -CONJG(gkl(j)) * xi * (cav2(k, l, a) + k3_cav2(k, l, a))
        B_vec(6) = -i * xi * 0.5d0 * Omega * (cav3(j, k, l, a) + k3_cav3(j, k, l, a)) + &
                 & gkl(k) * xi * (cavsig3(j, l, ata, gg) + k3_cavsig3(j, l, ata, gg)) + &
                 & gkl(k) * xi * (cavsig3(j, l, ata, ee) + k3_cavsig3(j, l, ata, ee)) + &
                 & -gkl(k) * xi * (cav2(j, l, ata) + k3_cav2(j, l, ata)) + &
                 & gkl(l) * xi * (cavsig3(j, k, ata, gg) + k3_cavsig3(j, k, ata, gg)) + &
                 & gkl(l) * xi * (cavsig3(j, k, ata, ee) + k3_cavsig3(j, k, ata, ee)) + &
                 & -gkl(l) * xi * (cav2(j, k, ata) + k3_cav2(j, k, ata))
        B_vec(7) = -CONJG(gkl(j)) * (cavsig3(k, l, a, ef) + k3_cavsig3(k, l, a, ef))
        B_vec(8) = -gkl(k) * (cavsig3(j, l, ata, fe) + k3_cavsig3(j, l, ata, fe)) + &
                 & -gkl(l) * (cavsig3(j, k, ata, fe) + k3_cavsig3(j, k, ata, fe))
        ! Calculate k4
        k4_cavsig4(j, k, l, a, :) = dt * (MATMUL(Mat, (cavsig4(j, k, l, a, :) + k3_cavsig4(j, k, l, a, :))) + B_vec)

        !-------------------------------------------!
        ! < a^{\dagger}_{j} a^{\dagger}_{k} a_{l} > !
        !-------------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((3.0d0 * kappa) - i * (wl(j) + wl(k) - wl(l)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j)) * cavsig3(k, l, ata, eg) + &
                 & -CONJG(gkl(k)) * cavsig3(j, l, ata, eg) + &
                 & -gkl(l) * cavsig3(j, k, at, ge)
        B_vec(2) = -CONJG(gkl(j)) * cavsig3(k, l, ata, ee) + &
                 & -CONJG(gkl(k)) * cavsig3(j, l, ata, ee) + &
                 & -gkl(l) * xi * cavsig3(j, k, at, gf)
        B_vec(3) = -CONJG(gkl(j)) * xi * cavsig3(k, l, ata, fg) + &
                 & -CONJG(gkl(k)) * xi * cavsig3(j, l, ata, fg) + &
                 & -gkl(l) * cavsig3(j, k, at, ee)
        B_vec(4) = Gamma * (xi ** 2) * cav3(j, k, l, at) + &
                 & -CONJG(gkl(j)) * xi * cavsig3(k, l, ata, fe) + &
                 & -CONJG(gkl(k)) * xi * cavsig3(j, l, ata, fe) + &
                 & -gkl(l) * xi * cavsig3(j, k, at, ef)
        B_vec(5) = i * xi * 0.5d0 * Omega * cav3(j, k, l, at) + &
                 & CONJG(gkl(j)) * xi * cavsig3(k, l, ata, gg) + &
                 & CONJG(gkl(j)) * xi * cavsig3(k, l, ata, ee) + &
                 & -CONJG(gkl(j)) * xi * cav2(k, l, ata) + &
                 & CONJG(gkl(k)) * xi * cavsig3(j, l, ata, gg) + &
                 & CONJG(gkl(k)) * xi * cavsig3(j, l, ata, ee) + &
                 & -CONJG(gkl(k)) * xi * cav2(j, l, ata)
        B_vec(6) = -i * xi * 0.5d0 * Omega * cav3(j, k, l, at) + &
                 & gkl(l) * xi * cavsig3(j, k, at, gg) + &
                 & gkl(l) * xi * cavsig3(j, k, at, ee) + &
                 & -gkl(l) * xi * cav2(j, k, at)
        B_vec(7) = -CONJG(gkl(j)) * cavsig3(k, l, ata, ef) + &
                 & -CONJG(gkl(k)) * cavsig3(j, l, ata, ef)
        B_vec(8) = -gkl(l) * cavsig3(j, k, at, fe)
        ! Calculate k1
        k1_cavsig4(j, k, l, at, :) = dt * (MATMUL(Mat, cavsig4(j, k, l, at, :)) + B_vec)

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j)) * (cavsig3(k, l, ata, eg) + 0.5d0 * k1_cavsig3(k, l, ata, eg)) + &
                 & -CONJG(gkl(k)) * (cavsig3(j, l, ata, eg) + 0.5d0 * k1_cavsig3(j, l, ata, eg)) + &
                 & -gkl(l) * (cavsig3(j, k, at, ge) + 0.5d0 * k1_cavsig3(j, k, at, ge))
        B_vec(2) = -CONJG(gkl(j)) * (cavsig3(k, l, ata, ee) + 0.5d0 * k1_cavsig3(k, l, ata, ee)) + &
                 & -CONJG(gkl(k)) * (cavsig3(j, l, ata, ee) + 0.5d0 * k1_cavsig3(j, l, ata, ee)) + &
                 & -gkl(l) * xi * (cavsig3(j, k, at, gf) + 0.5d0 * k1_cavsig3(j, k, at, gf))
        B_vec(3) = -CONJG(gkl(j)) * xi * (cavsig3(k, l, ata, fg) + 0.5d0 * k1_cavsig3(k, l, ata, fg)) + &
                 & -CONJG(gkl(k)) * xi * (cavsig3(j, l, ata, fg) + 0.5d0 * k1_cavsig3(j, l, ata, fg)) + &
                 & -gkl(l) * (cavsig3(j, k, at, ee) + 0.5d0 * k1_cavsig3(j, k, at, ee))
        B_vec(4) = Gamma * (xi ** 2) * (cav3(j, k, l, at) + 0.5d0 * k1_cav3(j, k, l, at)) + &
                 & -CONJG(gkl(j)) * xi * (cavsig3(k, l, ata, fe) + 0.5d0 * k1_cavsig3(k, l, ata, fe)) + &
                 & -CONJG(gkl(k)) * xi * (cavsig3(j, l, ata, fe) + 0.5d0 * k1_cavsig3(j, l, ata, fe)) + &
                 & -gkl(l) * xi * (cavsig3(j, k, at, ef) + 0.5d0 * k1_cavsig3(j, k, at, ef))
        B_vec(5) = i * xi * 0.5d0 * Omega * (cav3(j, k, l, at) + 0.5d0 * k1_cav3(j, k, l, at)) + &
                 & CONJG(gkl(j)) * xi * (cavsig3(k, l, ata, gg) + 0.5d0 * k1_cavsig3(k, l, ata, gg)) + &
                 & CONJG(gkl(j)) * xi * (cavsig3(k, l, ata, ee) + 0.5d0 * k1_cavsig3(k, l, ata, ee)) + &
                 & -CONJG(gkl(j)) * xi * (cav2(k, l, ata) + 0.5d0 * k1_cav2(k, l, ata)) + &
                 & CONJG(gkl(k)) * xi * (cavsig3(j, l, ata, gg) + 0.5d0 * k1_cavsig3(j, l, ata, gg)) + &
                 & CONJG(gkl(k)) * xi * (cavsig3(j, l, ata, ee) + 0.5d0 * k1_cavsig3(j, l, ata, ee)) + &
                 & -CONJG(gkl(k)) * xi * (cav2(j, l, ata) + 0.5d0 * k1_cav2(j, l, ata))
        B_vec(6) = -i * xi * 0.5d0 * Omega * (cav3(j, k, l, at) + 0.5d0 * k1_cav3(j, k, l, at)) + &
                 & gkl(l) * xi * (cavsig3(j, k, at, gg) + 0.5d0 * k1_cavsig3(j, k, at, gg)) + &
                 & gkl(l) * xi * (cavsig3(j, k, at, ee) + 0.5d0 * k1_cavsig3(j, k, at, ee)) + &
                 & -gkl(l) * xi * (cav2(j, k, at) + 0.5d0 * k1_cav2(j, k, at))
        B_vec(7) = -CONJG(gkl(j)) * (cavsig3(k, l, ata, ef) + 0.5d0 * k1_cavsig3(k, l, ata, ef)) + &
                 & -CONJG(gkl(k)) * (cavsig3(j, l, ata, ef) + 0.5d0 * k1_cavsig3(j, l, ata, ef))
        B_vec(8) = -gkl(l) * (cavsig3(j, k, at, fe) + 0.5d0 * k1_cavsig3(j, k, at, fe))
        ! Calculate k2
        k2_cavsig4(j, k, l, at, :) = dt * (MATMUL(Mat, (cavsig4(j, k, l, at, :) + 0.5d0 * k1_cavsig4(j, k, l, at, :))) + B_vec)

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j)) * (cavsig3(k, l, ata, eg) + 0.5d0 * k2_cavsig3(k, l, ata, eg)) + &
                 & -CONJG(gkl(k)) * (cavsig3(j, l, ata, eg) + 0.5d0 * k2_cavsig3(j, l, ata, eg)) + &
                 & -gkl(l) * (cavsig3(j, k, at, ge) + 0.5d0 * k2_cavsig3(j, k, at, ge))
        B_vec(2) = -CONJG(gkl(j)) * (cavsig3(k, l, ata, ee) + 0.5d0 * k2_cavsig3(k, l, ata, ee)) + &
                 & -CONJG(gkl(k)) * (cavsig3(j, l, ata, ee) + 0.5d0 * k2_cavsig3(j, l, ata, ee)) + &
                 & -gkl(l) * xi * (cavsig3(j, k, at, gf) + 0.5d0 * k2_cavsig3(j, k, at, gf))
        B_vec(3) = -CONJG(gkl(j)) * xi * (cavsig3(k, l, ata, fg) + 0.5d0 * k2_cavsig3(k, l, ata, fg)) + &
                 & -CONJG(gkl(k)) * xi * (cavsig3(j, l, ata, fg) + 0.5d0 * k2_cavsig3(j, l, ata, fg)) + &
                 & -gkl(l) * (cavsig3(j, k, at, ee) + 0.5d0 * k2_cavsig3(j, k, at, ee))
        B_vec(4) = Gamma * (xi ** 2) * (cav3(j, k, l, at) + 0.5d0 * k2_cav3(j, k, l, at)) + &
                 & -CONJG(gkl(j)) * xi * (cavsig3(k, l, ata, fe) + 0.5d0 * k2_cavsig3(k, l, ata, fe)) + &
                 & -CONJG(gkl(k)) * xi * (cavsig3(j, l, ata, fe) + 0.5d0 * k2_cavsig3(j, l, ata, fe)) + &
                 & -gkl(l) * xi * (cavsig3(j, k, at, ef) + 0.5d0 * k2_cavsig3(j, k, at, ef))
        B_vec(5) = i * xi * 0.5d0 * Omega * (cav3(j, k, l, at) + 0.5d0 * k2_cav3(j, k, l, at)) + &
                 & CONJG(gkl(j)) * xi * (cavsig3(k, l, ata, gg) + 0.5d0 * k2_cavsig3(k, l, ata, gg)) + &
                 & CONJG(gkl(j)) * xi * (cavsig3(k, l, ata, ee) + 0.5d0 * k2_cavsig3(k, l, ata, ee)) + &
                 & -CONJG(gkl(j)) * xi * (cav2(k, l, ata) + 0.5d0 * k2_cav2(k, l, ata)) + &
                 & CONJG(gkl(k)) * xi * (cavsig3(j, l, ata, gg) + 0.5d0 * k2_cavsig3(j, l, ata, gg)) + &
                 & CONJG(gkl(k)) * xi * (cavsig3(j, l, ata, ee) + 0.5d0 * k2_cavsig3(j, l, ata, ee)) + &
                 & -CONJG(gkl(k)) * xi * (cav2(j, l, ata) + 0.5d0 * k2_cav2(j, l, ata))
        B_vec(6) = -i * xi * 0.5d0 * Omega * (cav3(j, k, l, at) + 0.5d0 * k2_cav3(j, k, l, at)) + &
                 & gkl(l) * xi * (cavsig3(j, k, at, gg) + 0.5d0 * k2_cavsig3(j, k, at, gg)) + &
                 & gkl(l) * xi * (cavsig3(j, k, at, ee) + 0.5d0 * k2_cavsig3(j, k, at, ee)) + &
                 & -gkl(l) * xi * (cav2(j, k, at) + 0.5d0 * k2_cav2(j, k, at))
        B_vec(7) = -CONJG(gkl(j)) * (cavsig3(k, l, ata, ef) + 0.5d0 * k2_cavsig3(k, l, ata, ef)) + &
             & -CONJG(gkl(k)) * (cavsig3(j, l, ata, ef) + 0.5d0 * k2_cavsig3(j, l, ata, ef))
        B_vec(8) = -gkl(l) * (cavsig3(j, k, at, fe) + 0.5d0 * k2_cavsig3(j, k, at, fe))
        ! Calculate k3
        k3_cavsig4(j, k, l, at, :) = dt * (MATMUL(Mat, (cavsig4(j, k, l, at, :) + 0.5d0 * k2_cavsig4(j, k, l, at, :))) + B_vec)

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -CONJG(gkl(j)) * (cavsig3(k, l, ata, eg) + k3_cavsig3(k, l, ata, eg)) + &
                 & -CONJG(gkl(k)) * (cavsig3(j, l, ata, eg) + k3_cavsig3(j, l, ata, eg)) + &
                 & -gkl(l) * (cavsig3(j, k, at, ge) + k3_cavsig3(j, k, at, ge))
        B_vec(2) = -CONJG(gkl(j)) * (cavsig3(k, l, ata, ee) + k3_cavsig3(k, l, ata, ee)) + &
                 & -CONJG(gkl(k)) * (cavsig3(j, l, ata, ee) + k3_cavsig3(j, l, ata, ee)) + &
                 & -gkl(l) * xi * (cavsig3(j, k, at, gf) + k3_cavsig3(j, k, at, gf))
        B_vec(3) = -CONJG(gkl(j)) * xi * (cavsig3(k, l, ata, fg) + k3_cavsig3(k, l, ata, fg)) + &
                 & -CONJG(gkl(k)) * xi * (cavsig3(j, l, ata, fg) + k3_cavsig3(j, l, ata, fg)) + &
                 & -gkl(l) * (cavsig3(j, k, at, ee) + k3_cavsig3(j, k, at, ee))
        B_vec(4) = Gamma * (xi ** 2) * (cav3(j, k, l, at) + k3_cav3(j, k, l, at)) + &
                 & -CONJG(gkl(j)) * xi * (cavsig3(k, l, ata, fe) + k3_cavsig3(k, l, ata, fe)) + &
                 & -CONJG(gkl(k)) * xi * (cavsig3(j, l, ata, fe) + k3_cavsig3(j, l, ata, fe)) + &
                 & -gkl(l) * xi * (cavsig3(j, k, at, ef) + k3_cavsig3(j, k, at, ef))
        B_vec(5) = i * xi * 0.5d0 * Omega * (cav3(j, k, l, at) + k3_cav3(j, k, l, at)) + &
                 & CONJG(gkl(j)) * xi * (cavsig3(k, l, ata, gg) + k3_cavsig3(k, l, ata, gg)) + &
                 & CONJG(gkl(j)) * xi * (cavsig3(k, l, ata, ee) + k3_cavsig3(k, l, ata, ee)) + &
                 & -CONJG(gkl(j)) * xi * (cav2(k, l, ata) + k3_cav2(k, l, ata)) + &
                 & CONJG(gkl(k)) * xi * (cavsig3(j, l, ata, gg) + k3_cavsig3(j, l, ata, gg)) + &
                 & CONJG(gkl(k)) * xi * (cavsig3(j, l, ata, ee) + k3_cavsig3(j, l, ata, ee)) + &
                 & -CONJG(gkl(k)) * xi * (cav2(j, l, ata) + k3_cav2(j, l, ata))
        B_vec(6) = -i * xi * 0.5d0 * Omega * (cav3(j, k, l, at) + k3_cav3(j, k, l, at)) + &
                 & gkl(l) * xi * (cavsig3(j, k, at, gg) + k3_cavsig3(j, k, at, gg)) + &
                 & gkl(l) * xi * (cavsig3(j, k, at, ee) + k3_cavsig3(j, k, at, ee)) + &
                 & -gkl(l) * xi * (cav2(j, k, at) + k3_cav2(j, k, at))
        B_vec(7) = -CONJG(gkl(j)) * (cavsig3(k, l, ata, ef) + k3_cavsig3(k, l, ata, ef)) + &
                 & -CONJG(gkl(k)) * (cavsig3(j, l, ata, ef) + k3_cavsig3(j, l, ata, ef))
        B_vec(8) = -gkl(l) * (cavsig3(j, k, at, fe) + k3_cavsig3(j, k, at, fe))
        ! Calculate k4
        k4_cavsig4(j, k, l, at, :) = dt * (MATMUL(Mat, (cavsig4(j, k, l, at, :) + k3_cavsig4(j, k, l, at, :))) + B_vec)

        ! Close j loop
      END DO
      ! CLose k loop
    END DO
    ! Close l loop
  END DO

  ! Cycle though modes
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
          k1_cav4(j, k, l, m) = -dt * (4.0d0 * kappa - i * (wl(j) + wl(k) - wl(l) - wl(m))) * cav4(j, k, l, m) + &
                              & -dt * CONJG(gkl(j)) * cavsig4(k, l, m, a, eg) + &
                              & -dt * CONJG(gkl(j)) * xi * cavsig4(k, l, m, a, fe) + &
                              & -dt * CONJG(gkl(k)) * cavsig4(j, l, m, a, eg) + &
                              & -dt * CONJG(gkl(k)) * xi * cavsig4(j, l, m, a, fe) + &
                              & -dt * gkl(l) * cavsig4(j, k, m, at, ge) + &
                              & -dt * gkl(l) * xi * cavsig4(j, k, m, at, ef) + &
                              & -dt * gkl(m) * cavsig4(j, k, l, at, ge) + &
                              & -dt * gkl(m) * xi * cavsig4(j, k, l, at, ef)

          k2_cav4(j, k, l, m) = -dt * (4.0d0 * kappa - i * (wl(j) + wl(k) - wl(l) - wl(m))) * (cav4(j, k, l, m) + 0.5d0 * k1_cav4(j, k, l, m)) + &
                              & -dt * CONJG(gkl(j)) * (cavsig4(k, l, m, a, eg) + 0.5d0 * k1_cavsig4(k, l, m, a, eg)) + &
                              & -dt * CONJG(gkl(j)) * xi * (cavsig4(k, l, m, a, fe) + 0.5d0 * k1_cavsig4(k, l, m, a, fe)) + &
                              & -dt * CONJG(gkl(k)) * (cavsig4(j, l, m, a, eg) + 0.5d0 * k1_cavsig4(j, l, m, a, eg)) + &
                              & -dt * CONJG(gkl(k)) * xi * (cavsig4(j, l, m, a, fe) + 0.5d0 * k1_cavsig4(j, l, m, a, fe)) + &
                              & -dt * gkl(l) * (cavsig4(j, k, m, at, ge) + 0.5d0 * k1_cavsig4(j, k, m, at, ge)) + &
                              & -dt * gkl(l) * xi * (cavsig4(j, k, m, at, ef) + 0.5d0 * k1_cavsig4(j, k, m, at, ef)) + &
                              & -dt * gkl(m) * (cavsig4(j, k, l, at, ge) + 0.5d0 * k1_cavsig4(j, k, l, at, ge)) + &
                              & -dt * gkl(m) * xi * (cavsig4(j, k, l, at, ef) + 0.5d0 * k1_cavsig4(j, k, l, at, ef))

          k3_cav4(j, k, l, m) = -dt * (4.0d0 * kappa - i * (wl(j) + wl(k) - wl(l) - wl(m))) * (cav4(j, k, l, m) + 0.5d0 * k2_cav4(j, k, l, m)) + &
                              & -dt * CONJG(gkl(j)) * (cavsig4(k, l, m, a, eg) + 0.5d0 * k2_cavsig4(k, l, m, a, eg)) + &
                              & -dt * CONJG(gkl(j)) * xi * (cavsig4(k, l, m, a, fe) + 0.5d0 * k2_cavsig4(k, l, m, a, fe)) + &
                              & -dt * CONJG(gkl(k)) * (cavsig4(j, l, m, a, eg) + 0.5d0 * k2_cavsig4(j, l, m, a, eg)) + &
                              & -dt * CONJG(gkl(k)) * xi * (cavsig4(j, l, m, a, fe) + 0.5d0 * k2_cavsig4(j, l, m, a, fe)) + &
                              & -dt * gkl(l) * (cavsig4(j, k, m, at, ge) + 0.5d0 * k2_cavsig4(j, k, m, at, ge)) + &
                              & -dt * gkl(l) * xi * (cavsig4(j, k, m, at, ef) + 0.5d0 * k2_cavsig4(j, k, m, at, ef)) + &
                              & -dt * gkl(m) * (cavsig4(j, k, l, at, ge) + 0.5d0 * k2_cavsig4(j, k, l, at, ge)) + &
                              & -dt * gkl(m) * xi * (cavsig4(j, k, l, at, ef) + 0.5d0 * k2_cavsig4(j, k, l, at, ef))

          k4_cav4(j, k, l, m) = -dt * (4.0d0 * kappa - i * (wl(j) + wl(k) - wl(l) - wl(m))) * (cav4(j, k, l, m) + k3_cav4(j, k, l, m)) + &
                              & -dt * CONJG(gkl(j)) * (cavsig4(k, l, m, a, eg) + k3_cavsig4(k, l, m, a, eg)) + &
                              & -dt * CONJG(gkl(j)) * xi * (cavsig4(k, l, m, a, fe) + k3_cavsig4(k, l, m, a, fe)) + &
                              & -dt * CONJG(gkl(k)) * (cavsig4(j, l, m, a, eg) + k3_cavsig4(j, l, m, a, eg)) + &
                              & -dt * CONJG(gkl(k)) * xi * (cavsig4(j, l, m, a, fe) + k3_cavsig4(j, l, m, a, fe)) + &
                              & -dt * gkl(l) * (cavsig4(j, k, m, at, ge) + k3_cavsig4(j, k, m, at, ge)) + &
                              & -dt * gkl(l) * xi * (cavsig4(j, k, m, at, ef) + k3_cavsig4(j, k, m, at, ef)) + &
                              & -dt * gkl(m) * (cavsig4(j, k, l, at, ge) + k3_cavsig4(j, k, l, at, ge)) + &
                              & -dt * gkl(m) * xi * (cavsig4(j, k, l, at, ef) + k3_cavsig4(j, k, l, at, ef))
        END DO
      END DO
    END DO
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
  ! Third-order
  cavsig3 = cavsig3 + xis * (k1_cavsig3 + 2.0d0 * (k2_cavsig3 + k3_cavsig3) + k4_cavsig3)
  cav3 = cav3 + xis * (k1_cav3 + 2.0d0 * (k2_cav3 + k3_cav3) + k4_cav3)
  ! Fourth-order
  cavsig4 = cavsig4 + xis * (k1_cavsig4 + 2.0d0 * (k2_cavsig4 + k3_cavsig4) + k4_cavsig4)
  cav4 = cav4 + xis * (k1_cav4 + 2.0d0 * (k2_cav4 + k3_cav4) + k4_cav4)
END DO

! Close file
CLOSE(2)
CLOSE(10)

! Save steady state moments
sigma_ss = sigma
cav1_ss = cav1
cavsig2_ss = cavsig2
cav2_ss = cav2
cavsig3_ss = cavsig3
cav3_ss = cav3
cavsig4_ss = cavsig4
cav4_ss = cav4
! and photon number
photon_ss = photon

!==============================================================================!
!                      TEST PRINT SOME STEADY STATE VALUES                     !
!==============================================================================!
! Set modes to output
IF (N .EQ. 0) THEN
  j = 0; k = 0; l = 0; m = 0
ELSE
  j = 1; k = -1; l = 0; m = 1
END IF

! WRITE(*, *) "=============================================="
! WRITE(*, *) "FIRST-ORDER: ATOM"
! WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< gg >_ss = ", REAL(sigma_ss(gg)), " + ", IMAG(sigma_ss(gg)), "i"
! WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< ee >_ss = ", REAL(sigma_ss(ee)), " + ", IMAG(sigma_ss(ee)), "i"
! WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< eg >_ss = ", REAL(sigma_ss(eg)), " + ", IMAG(sigma_ss(eg)), "i"
! WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "FIRST-ORDER: CAVITY"
! WRITE(*, '(A5,I2,A8,ES18.11E2 A3,ES18.11E2,A2)') " < a_", j, " >_ss = ", REAL(cav1_ss(j, a)), " + ", IMAG(cav1_ss(j, a)), "i"
! WRITE(*, '(A5,I2,A8,ES18.11E2 A3,ES18.11E2,A2)') "< at_", j, " >_ss = ", REAL(cav1_ss(j, at)), " + ", IMAG(cav1_ss(j, at)), "i"
! WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "SECOND-ORDER: CAVITY AND ATOM"
! WRITE(*, '(A5,I2,A11,ES18.11E2,A3,ES18.11E2,A2)') " < a_", j, " gg >_ss = ", REAL(cavsig2_ss(j, a, gg)), " + ", IMAG(cavsig2_ss(j, a, gg)), "i"
! WRITE(*, '(A5,I2,A11,ES18.11E2,A3,ES18.11E2,A2)') " < a_", j, " ee >_ss = ", REAL(cavsig2_ss(j, a, ee)), " + ", IMAG(cavsig2_ss(j, a, ee)), "i"
! WRITE(*, '(A5,I2,A11,ES18.11E2,A3,ES18.11E2,A2)') "< at_", j, " gg >_ss = ", REAL(cavsig2_ss(j, at, gg)), " + ", IMAG(cavsig2_ss(j, at, gg)), "i"
! WRITE(*, '(A5,I2,A11,ES18.11E2,A3,ES18.11E2,A2)') "< at_", j, " ee >_ss = ", REAL(cavsig2_ss(j, at, ee)), " + ", IMAG(cavsig2_ss(j, at, ee)), "i"
! WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "SECOND-ORDER: CAVITY"
! WRITE(*, '(A6,I2,A4,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') "  < a_", j, "  a_", k, " >_ss = ", REAL(cav2_ss(j, k, a)), " + ", AIMAG(cav2_ss(j, k, a)), "i"
! WRITE(*, '(A6,I2,A4,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < at_", j, " at_", k, " >_ss = ", REAL(cav2_ss(j, k, at)), " + ", AIMAG(cav2_ss(j, k, at)), "i"
! WRITE(*, '(A6,I2,A4,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < at_", j, "  a_", k, " >_ss = ", REAL(cav2_ss(j, k, ata)), " + ", AIMAG(cav2_ss(j, k, ata)), "i"
! WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "THIRD-ORDER: CAVITY AND ATOM"
! WRITE(*, '(A6,I2,A4,I2,A11,ES18.11E2,A3,ES18.11E2,A1)') "  < a_", j, "  a_", k, " fe >_ss = ", REAL(cavsig3_ss(j, k, a, fe)), " + ", AIMAG(cavsig3_ss(j, k, a, fe)), "i"
! WRITE(*, '(A6,I2,A4,I2,A11,ES18.11E2,A3,ES18.11E2,A1)') " < at_", j, " at_", k, " fe >_ss = ", REAL(cavsig3_ss(j, k, at, fe)), " + ", AIMAG(cavsig3_ss(j, k, at, fe)), "i"
! WRITE(*, '(A6,I2,A4,I2,A11,ES18.11E2,A3,ES18.11E2,A1)') " < at_", j, "  a_", k, " fe >_ss = ", REAL(cavsig3_ss(j, k, ata, fe)), " + ", AIMAG(cavsig3_ss(j, k, ata, fe)), "i"
! WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "THIRD-ORDER: CAVITY"
! WRITE(*, '(A6,I2,A4,I2,A3,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < at_", j, "  a_", k, " a_", l, " >_ss = ", REAL(cav3_ss(j, k, l, a)), " + ", AIMAG(cav3_ss(j, k, l, a)), "i"
! WRITE(*, '(A6,I2,A4,I2,A3,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < at_", j, " at_", k, " a_", l, " >_ss = ", REAL(cav3_ss(j, k, l, at)), " + ", AIMAG(cav3_ss(j, k, l, at)), "i"
! WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "FOURTH-ORDER: CAVITY AND ATOM"
WRITE(*, '(A26,ES18.11E2,A3,ES18.11E2,A1)') " < at_0 a_0 a_0 fe >_ss = ", REAL(cavsig4_ss(0, 0, 0, a, fe)), " + ", AIMAG(cavsig4_ss(0, 0, 0, a, fe)), "i"
WRITE(*, '(A26,ES18.11E2,A3,ES18.11E2,A1)') "< at_0 at_0 a_0 fe >_ss = ", REAL(cavsig4_ss(0, 0, 0, at, fe)), " + ", AIMAG(cavsig4_ss(0, 0, 0, at, fe)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "FOURTH-ORDER: CAVITY"
WRITE(*, '(A27,ES18.11E2,A3,ES18.11E2,A1)') "< at_0 at_0 a_0 a_0 >_ss = ", REAL(cav4_ss(0, 0, 0, 0)), " + ", AIMAG(cav4_ss(0, 0, 0, 0)), "i"
WRITE(*, *) "=============================================="

PRINT*, " "
PRINT*, "Mean photon number = ", photon_ss
PRINT*, " "

!==============================================================================!
!                                END OF PROGRAM                                !
!==============================================================================!

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

END PROGRAM THREE_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_RK4
