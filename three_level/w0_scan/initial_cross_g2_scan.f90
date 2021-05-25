! The system in this program is a resonantly driven two-level atom that is
! coupled into a multi-mode array of filter cavities, as described in some notes
! somewhere (or with this program if this goes to anyone).

! This program calculates the steady state photon number inside the filter
! cavity for a scanned range of central resonance frequency values (\omega_{0}).
! The matrix of omega values for both cavity filters are written to
! [filename_w0] which, by default, points to "./data_files/w0_cross.txt".
! The file is written contains a single column of the scanned \omega_{0} values,
! to be np.meshgrid-ed in Python:
!                              \omega_{0}
! The cross-correlation data is wrriten to [filename_cross_corr] which, by
! default, points to "./data_files/cross_corr.txt". The file is written as a
! matrix, where
!   M(i, j) = g^{(2)}_{cross}(0) for \omega_{0}^{(a)}(i) and \omega_{0}^{(b)}(j)

! The input parameters are taken from a NameList file [filename_ParamList] which,
! by default, points to "./ParamList.nml". The code can thus be compiled once,
! and parameters can be changed in the NameList file for subsequent runs.

! The parameters for each run are written to [filename_parameters] which, by
! default, points to "./data_files/parameters.txt".

! For the default filenames, the folder "./data_files/" and NameList file
! "./ParamList.nml" MUST EXIST IN THE WORKING DIRECTORY.

! To compiled the code, I use the Intel Parallel Studio compiler IFORT with the
! command
!       (LINUX): ifort -O3 -o g2 -mkl g2_RK4.f90
!     (WINDOWS): ifort /O3 /o g2 /Qmkl g2_RK4.f90
! where the -O3 (/O3) flag gives maximum optimisation, the -o (/o) g1 flag
! names the executable as "g2" ("g2.exe"), and the -mkl (/Qmkl) flag links
! the program to Intel's Math Kernel Library, to make use of the LAPACK routines.

PROGRAM TWO_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_INITIAL_CROSS_G2

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
REAL(KIND=8)                                           :: w0a, w0b
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                           :: kappaa, kappab
! Blackman window coefficient
REAL(KIND=8)                                           :: blackman
! Frequency spacing of modes
REAL(KIND=8)                                           :: dwa, dwb
! List of Delta values
REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE             :: wl
! List of mode dependent cascade coupling values
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: gkl
! Kappa values for both cavities
REAL(KIND=8), DIMENSION(2)                             :: kappa, w0, dw

! Scanning parameters
! Number of scans either side of centre
INTEGER                                                :: no_scans
! List of central frequencies
REAL(KIND=8), DIMENSION(:), ALLOCATABLE                :: w0l
! Scanning integer
INTEGER                                                :: scana, scanb
! Max/min value of D0a [Default = 30.0]
REAL(KIND=8)                                           :: scan_max
! Scan step size [Default = 0.5]
REAL(KIND=8)                                           :: scan_step

! ! Time stuff
! ! Time step
! REAL(KIND=8)                                           :: dt
! ! Maximum time to integrate for
! REAL(KIND=8)                                           :: t_max, tau1_max, tau2_max
! ! Maximum number of steps to integrate for
! INTEGER                                                :: t_steps, tau_steps
! ! Time step integer
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
! COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: cav1a_ss, cav1b_ss
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cav1_ss
! Second-order moments: Cavity and atom (< a \sigma >, < a^{\dagger} \sigma >
! COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cavsig2a_ss, cavsig2b_ss
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cavsig2_ss
! Second-order moments: Cavity (< a^{\dagger} a >)
! COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cav2_ss
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: cav2_ss
! Third-order moments: Cavity and atom (< a^{2} \sigma >, < a^{\dagger 2} \sigma >, < a^{\dagger} a \sigma >)
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cavsig3_ss
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: cavsig3_ss
! Third-order moments: Cavity (< a^{2} a^{\dagger} >, < a^{\dagger 2} a >)
! COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cav3_ss
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :, :), ALLOCATABLE :: cav3_ss
! Fourth-order moments: Cavity and atom ( < a^{\dagger} a^{2} \sigma >, < a^{\dagger 2} a \sigma >)
! COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: cavsig4_ss
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :, :, :), ALLOCATABLE :: cavsig4_ss
! Fourth-order moments: Cavity (< a^{\dagger 2} a^{2} >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cav4_ss
! COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :, :, :), ALLOCATABLE :: cav4_ss

! Integer indices for sigma operators
INTEGER, PARAMETER                                     :: gg = 1, ge = 2, eg = 3
INTEGER, PARAMETER                                     :: ee = 4, ef = 5, fe = 6
INTEGER, PARAMETER                                     :: gf = 7, fg = 8
! Integer indices for: a, b, a^{\dagger}, b^{\dagger}
INTEGER, PARAMETER                                     :: a = 1, at = 2, bt = 1, b = 2
! Integer indices for: ab, a^{\dagger} b^{\dagger}, a^{\dagger} a, b^{\dagger} b
INTEGER, PARAMETER                                     :: ab = 1, atbt = 2, ata = 3, btb = 4
! Integer indices for: a^{\dagger} b, b^{\dagger} a
INTEGER, PARAMETER                                     :: atb = 5, bta = 6
! Integer indices for: a^{\dagger} a b, b^{\dagger} b a
INTEGER, PARAMETER                                     :: atab = 1, btba = 2
! Integer indices for: b^{\dagger} a^{\dagger} a, a^{\dagger} b^{\dagger} b
INTEGER, PARAMETER                                     :: btata = 3, atbtb = 4

! Integer cavity counters
INTEGER                                                :: cj, ck, cl, cm

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
! Atomic population and photon number
! REAL(KIND=8)                                           :: popg, pope, popf, photon
! Steady state photon number
REAL(KIND=8), DIMENSION(2)                             :: photon_ss
! Complex data
COMPLEX(KIND=8)                                        :: moment_out, moment_out2
! Matrix of correlation values
REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE             :: corr_matrix

!------------------------!
!     FILENAME STUFF     !
!------------------------!
! Paramert Name List
CHARACTER(LEN=15), PARAMETER :: filename_ParamList = "./ParamList.nml"
! Filename of parameters
CHARACTER(LEN=33), PARAMETER :: filename_parameters = "./data_files/cross_parameters.txt"
! Filename for state population
CHARACTER(LEN=25), PARAMETER :: filename_w0 = "./data_files/w0_cross.txt"
! Filename for second-order correlation
CHARACTER(LEN=27), PARAMETER :: filename_cross_corr = "./data_files/cross_corr.txt"

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
NAMELIST /SCANLIST/ scan_max, scan_step
! NAMELIST /TIME/ dt, t_max, tau1_max, tau2_max

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

READ(IUNIT, NML=SCANLIST, IOSTAT=ISTAT)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') LINE
  CLOSE(IUNIT)
  PRINT *, "Invalid line in SCANLIST namelist: " // TRIM(line)
  CALL EXIT(1)
END IF

! Set system parameters
kappa(a) = kappaa; kappa(b) = kappab
w0(a) = w0a; w0(b) = w0b
dw(a) = dwa; dw(b) = dwb

! Max number of scans
no_scans = NINT(scan_max / scan_step)

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
ALLOCATE(wl(-N:N, 2))
wl = 0.0d0
ALLOCATE(gkl(-N:N, 2))
gkl = 0.0d0
DO j = -N, N
  IF (N == 0) THEN
    ! Cavity A
    wl(j, a) = w0a
    gkl(j, a) = DSQRT(0.5d0 * epsilon * Gamma * kappaa)
    ! Cavity B
    wl(j, b) = w0b
    gkl(j, b) = DSQRT(0.5d0 * epsilon * Gamma * kappab)
  ELSE
    ! Blackman window coefficient
    blackman = 1.0d0
    ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
    !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))
    ! Cavity A
    wl(j, a) = w0a + DBLE(j) * dwa
    ! Mode dependent phase difference
    gkl(j, a) = DSQRT((0.5d0 * epsilon / DBLE(2*N + 1)) * Gamma * kappaa) * blackman * EXP(i * DBLE(phase) * DBLE(j) * pi / DBLE(N))
    ! Cavity B
    wl(j, b) = w0b + DBLE(j) * dwb
    ! Mode dependent phase difference
    gkl(j, b) = DSQRT((0.5d0 * epsilon / DBLE(2*N + 1)) * Gamma * kappab) * blackman * EXP(i * DBLE(phase) * DBLE(j) * pi / DBLE(N))
  END IF
END DO

! Allocate central frequency list
ALLOCATE(w0l(-no_scans:no_scans))
w0l = 0.0d0
DO scana = -no_scans, no_scans
  w0l(scana) = 0.0d0 + DBLE(scana) * scan_step
END DO

! Allocate data matrix
ALLOCATE(corr_matrix(-no_scans:no_scans, -no_scans:no_scans))
corr_matrix = 0.0d0

!------------------------------------------!
!     INITALISE OPERATOR MOMENT ARRAYS     !
!------------------------------------------!
! Steady states
! First-order: Cavity
! ALLOCATE(cav1a_ss(-N:N, 2)); cav1a_ss = 0.0d0
! ALLOCATE(cav1b_ss(-N:N, 2)); cav1b_ss = 0.0d0
ALLOCATE(cav1_ss(-N:N, 2, 2)); cav1_ss = 0.0d0
! Second-order: Cavity and Atom
! ALLOCATE(cavsig2a_ss(-N:N, 2, N_mat)); cavsig2a_ss = 0.0d0
! ALLOCATE(cavsig2b_ss(-N:N, 2, N_mat)); cavsig2b_ss = 0.0d0
ALLOCATE(cavsig2_ss(-N:N, 2, N_mat, 2)); cavsig2_ss = 0.0d0
! Second-order: Cavity
! ALLOCATE(cav2_ss(-N:N, -N:N, 6)); cav2_ss = 0.0d0
ALLOCATE(cav2_ss(-N:N, -N:N, 3, 2, 2)); cav2_ss = 0.0d0
! Third-order: Cavity and Atom
! ALLOCATE(cavsig3_ss(-N:N, -N:N, 6, N_mat)); cavsig3_ss = 0.0d0
ALLOCATE(cavsig3_ss(-N:N, -N:N, 3, N_mat, 2, 2)); cavsig3_ss = 0.0d0
! Third-order: Cavity
! ALLOCATE(cav3_ss(-N:N, -N:N, -N:N, 4)); cav3_ss = 0.0d0
ALLOCATE(cav3_ss(-N:N, -N:N, -N:N, 2, 2, 2, 2)); cav3_ss = 0.0d0
! Fourth-order: Cavity and atom
! ALLOCATE(cavsig4_ss(-N:N, -N:N, -N:N, 4, N_mat)); cavsig4_ss = 0.0d0
ALLOCATE(cavsig4_ss(-N:N, -N:N, -N:N, 2, N_mat, 2, 2, 2)); cavsig4_ss = 0.0d0
! Fourth-order: Cavity
ALLOCATE(cav4_ss(-N:N, -N:N, -N:N, -N:N)); cav4_ss = 0.0d0
! ALLOCATE(cav4_ss(-N:N, -N:N, -N:N, -N:N, 2, 2, 2, 2)); cav4_ss = 0.0d0

! !----------------------------!
! !     INITIAL CONDITIONS     !
! !----------------------------!
! sigma = 0.0d0
! ! Atom in ground state, cjity in vacuum.
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
! Cycle through each different \omega_{0} value and calculate steady states

! Open file to write time and data to
OPEN(UNIT=2, FILE=filename_w0, STATUS='REPLACE', ACTION='WRITE', RECL=4000)
OPEN(UNIT=3, FILE=filename_cross_corr, STATUS='REPLACE', ACTION='WRITE', RECL=16000)

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

! Start scan loop
DO scana = -no_scans, no_scans
  DO scanb = -no_scans, no_scans
    ! Set central resonance frequency
    w0a = w0l(scana)
    w0b = w0l(scanb)
    ! Initialise mode resonance frequency list
    wl = 0.0d0
    DO j = -N, N
      IF (N == 0) THEN
        ! Cavity A
        wl(j, a) = w0a
        ! Cavity B
        wl(j, b) = w0b
      ELSE
      ! Cavity A
      wl(j, a) = w0a + DBLE(j) * dwa
      ! Cavity B
      wl(j, b) = w0b + DBLE(j) * dwb
      END IF
    END DO

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
        ! Invert matrix
        CALL SquareMatrixInverse(N_mat, Mat_inv)

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
        ! Invert matrix
        CALL SquareMatrixInverse(N_mat, Mat_inv)

        ! Calculate steady state
        cavsig2_ss(j, at, :, cj) = -MATMUL(Mat_inv, B_vec)

        ! Close j loop
      END DO
      ! Close cj loop
    END DO

    moment_out = 0.0d0; moment_out2 = 0.0d0
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

            ! Update phootn numbers
            IF (cj .EQ. 1 .AND. ck .EQ. 1) THEN
              moment_out = moment_out + cav2_ss(j, k, ata, cj, ck)
            ELSE IF (cj .EQ. 2 .AND. ck .EQ. 2) THEN
              moment_out2 = moment_out2 + cav2_ss(j, k, ata, cj, ck)
            END IF

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
            ! Invert matrix
            CALL SquareMatrixInverse(N_mat, Mat_inv)

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
            ! Invert matrix
            CALL SquareMatrixInverse(N_mat, Mat_inv)

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
            ! Invert matrix
            CALL SquareMatrixInverse(N_mat, Mat_inv)

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

    ! Update photon number
    photon_ss = 0.0d0
    photon_ss(a) = moment_out
    photon_ss(b) = moment_out2

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
                ! Invert matrix
                CALL SquareMatrixInverse(N_mat, Mat_inv)

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
                ! Invert matrix
                CALL SquareMatrixInverse(N_mat, Mat_inv)

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

    moment_out = 0.0d0
    ! Cycle through modes
    DO m = -N, N
      DO l = -N, N
        DO k = -N, N
          DO j = -N, N
            !------------------------------!
            !     FOURTH-ORDER: CAVITY     !
            !------------------------------!
            !-------------------------------------------------!
            ! < a^{\dagger}_{j} b^{\dagger}_{k} b_{l} a_{m} > !
            !-------------------------------------------------!
            cav4_ss(j, k, l, m) = -CONJG(gkl(j, a)) * cavsig4_ss(k, l, m, a, eg, 2, 2, 1) + &
                                & -CONJG(gkl(j, a)) * xi * cavsig4_ss(k, l, m, a, fe, 2, 2, 1) + &
                                & -CONJG(gkl(k, b)) * cavsig4_ss(j, l, m, a, eg, 1, 2, 1) + &
                                & -CONJG(gkl(k, b)) * xi * cavsig4_ss(j, l, m, a, fe, 1, 2, 1) + &
                                & -gkl(l, b) * cavsig4_ss(j, k, m, at, ge, 1, 2, 1) + &
                                & -gkl(l, b) * xi * cavsig4_ss(j, k, m, at, ef, 1, 2, 1) + &
                                & -gkl(m, a) * cavsig4_ss(j, k, l, at, ge, 1, 2, 2) + &
                                & -gkl(m, a) * xi * cavsig4_ss(j, k, l, at, ef, 1, 2, 2)
            cav4_ss(j, k, l, m) = cav4_ss(j, k, l, m) / &
                                & ((2.0d0 * kappa(a)) + (2.0d0 * kappa(b)) - i * (wl(j, a) + wl(k, b)) + i * (wl(l, b) + wl(m, a)))

            ! Update initial correlation value
            moment_out = moment_out + cav4_ss(j, k, l, m)
            ! Close j loop
          END DO
          ! Close k loop
        END DO
        ! Close l loop
      END DO
      ! Close m loop
    END DO

    ! Normalise correlation value
    moment_out = moment_out / (photon_ss(a) * photon_ss(b))

    ! Write value to data matrix
    corr_matrix(scana, scanb) = REAL(moment_out)

    ! Close scanb loop
  END DO
  !-----------------------------!
  !     WRITE DATA TO FILES     !
  !-----------------------------!
  ! Write \omega_{0} values
  WRITE(2, *) w0l(scana)
  ! Write correlation data
  ! WRITE(3, *) corr_matrix(:, scana)
  WRITE(3, *) corr_matrix(scana, :)

  ! Close scana loop
END DO

! Close file
CLOSE(2)
CLOSE(3)

!==============================================================================!
!                                END OF PROGRAM                                !
!==============================================================================!

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

END PROGRAM TWO_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_INITIAL_CROSS_G2

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
