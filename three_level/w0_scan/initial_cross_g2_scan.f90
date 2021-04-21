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
REAL(KIND=8)                                           :: w0a, w0b
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                           :: kappaa, kappab
! Percentage of fluorecence aimed at cavity
REAL(KIND=8)                                           :: epsilon
! Number of mode either side of w0, 2N + 1 total mode
INTEGER                                                :: N
! Frequency spacing of modes
REAL(KIND=8)                                           :: dwa, dwb
! Phase modulation of mode coupling
INTEGER                                                :: phase
! List of Delta values
REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE              :: wl
! List of mode dependent cascade coupling values
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: gkl
! Blackman window coefficient
REAL(KIND=8)                                           :: blackman
! kappa, w0, and dw values
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
INTEGER, PARAMETER                                     :: N_mat = 3
! M matrix (filled as transpose)
COMPLEX(KIND=8), DIMENSION(N_mat, N_mat)               :: Mat, Mat_OG, Mat_inv
! Non-homogeneous vector
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: B_vec, B_OG

! Steady state arrays
! First-order moments: Atomic equations (< \sigma >)
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: sigma_ss
! First-order moments: Cavity (< a >, < a^{\dagger} >)
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: cav1_ss, cav1b_ss
! Second-order moments: Cavity and atom (< a \sigma >, < a^{\dagger} \sigma >
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cavsig2_ss, cavsig2b_ss
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
INTEGER, PARAMETER                                     :: sm = 1, sp = 2, sz = 3
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
! REAL(KIND=8)                                           :: popg, pope, popf, photon
! Steady state photon number
REAL(KIND=8), DIMENSION(2)                             :: photon_ss
! Complex data
COMPLEX(KIND=8)                                        :: moment_out, moment_out2
! Matrix of correlation values
REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE             :: corr_matrix

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

! READ(IUNIT, NML=TIME, IOSTAT=ISTAT)
! IF (ISTAT .NE. 0) THEN
!   BACKSPACE(IUNIT)
!   READ(IUNIT, FMT='(A)') LINE
!   CLOSE(IUNIT)
!   PRINT *, "Invalid line in TIME namelist: " // TRIM(line)
!   CALL EXIT(1)
! END IF
! CLOSE(IUNIT)

! Number of time-steps
! t_steps = NINT(t_max / dt)

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
! Allocate array of Delta and gka values
ALLOCATE(wl(-N:N, 2))
wl = 0.0d0
ALLOCATE(gkl(-N:N, 2))
gkl = 0.0d0
DO j = -N, N
  IF (N == 0) THEN
    ! Cavity A
    wl(j, a) = w0(a)
    gkl(j, a) = DSQRT(epsilon * Gamma * kappa(a))
    ! Cavity B
    wl(j, b) = w0(b)
    gkl(j, b) = DSQRT(epsilon * Gamma * kappa(b))
  ELSE
    ! Blackman window coefficient
    blackman = 1.0d0
    ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N))) + &
    !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + j) / (2.0d0 * DBLE(N)))
    ! Cavity A
    wl(j, a) = w0(a) + DBLE(j) * dw(a)
    ! Mode dependent phase difference
    gkl(j, a) = DSQRT((epsilon / DBLE(2*N + 1)) * Gamma * kappa(a)) * blackman * EXP(i * DBLE(phase) * DBLE(j) * pi / DBLE(N))
    ! Cavity B
    wl(j, b) = w0(b) + DBLE(j) * dw(b)
    ! Mode dependent phase difference
    gkl(j, b) = DSQRT((epsilon / DBLE(2*N + 1)) * Gamma * kappa(b)) * blackman * EXP(i * DBLE(phase) * DBLE(j) * pi / DBLE(N))
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
ALLOCATE(cav1_ss(-N:N, 2)); cav1_ss = 0.0d0
ALLOCATE(cav1b_ss(-N:N, 2)); cav1b_ss = 0.0d0
! Second-order: Cavity and Atom
ALLOCATE(cavsig2_ss(-N:N, 2, N_mat)); cavsig2_ss = 0.0d0
ALLOCATE(cavsig2b_ss(-N:N, 2, N_mat)); cavsig2b_ss = 0.0d0
! Second-order: Cavity
ALLOCATE(cav2_ss(-N:N, -N:N, 6)); cav2_ss = 0.0d0
! Third-order: Cavity and Atom
ALLOCATE(cavsig3_ss(-N:N, -N:N, 6, N_mat)); cavsig3_ss = 0.0d0
! Third-order: Cavity
ALLOCATE(cav3_ss(-N:N, -N:N, -N:N, 4)); cav3_ss = 0.0d0
! Fourth-order: Cavity and atom
ALLOCATE(cavsig4_ss(-N:N, -N:N, -N:N, 4, N_mat)); cavsig4_ss = 0.0d0
! Fourth-order: Cavity
ALLOCATE(cav4_ss(-N:N, -N:N, -N:N, -N:N)); cav4_ss = 0.0d0

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

    ! Cycle through modes
    DO j = -N, N
      !-----------------------------!
      !     FIRST-ORDER: CAVITY     !
      !-----------------------------!
      !-----------!
      ! < a_{j} > !
      !-----------!
      cav1_ss(j, a) = -gkl(j, a) * sigma_ss(sm)
      cav1_ss(j, a) = cav1_ss(j, a) / &
                    & (kappa(a) + i * wl(j, a))
      !---------------------!
      ! < a^{\dagger}_{j} > !
      !---------------------!
      cav1_ss(j, at) = -CONJG(gkl(j, a)) * sigma_ss(sp)
      cav1_ss(j, at) = cav1_ss(j, at) / &
                     & (kappa(a) - i * wl(j, a))

      !-----------!
      ! < b_{j} > !
      !-----------!
      cav1b_ss(j, b) = -gkl(j, b) * sigma_ss(sm)
      cav1b_ss(j, b) = cav1b_ss(j, b) / &
                     & (kappa(b) + i * wl(j, b))
      !---------------------!
      ! < b^{\dagger}_{j} > !
      !---------------------!
      cav1b_ss(j, bt) = -CONJG(gkl(j, b)) * sigma_ss(sp)
      cav1b_ss(j, bt) = cav1b_ss(j, bt) / &
                      & (kappa(b) - i * wl(j, b))
      !---------------------------------------!
      !     SECOND-ORDER: CAVITY AND ATOM     !
      !---------------------------------------!
      !------------------!
      ! < a_{j} \sigma > !
      !------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - (kappa(a) + i * wl(j, a))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = 0.0d0
      B_vec(2) = -0.5d0 * gkl(j, a) * (sigma_ss(sz) + 1.0d0)
      B_vec(3) = -gamma * cav1_ss(j, a) + &
               & gkl(j, a) * sigma_ss(sm)

      ! Set inverse matrix
      Mat_inv = Mat
      ! Perform LU-factorization of matrix
      CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
      ! Invert Matrix (Optimal LWORK = 3)
      LWORK = 3
      CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

      ! Calculate steady state
      cavsig2_ss(j, a, :) = -MATMUL(Mat_inv, B_vec)

      !----------------------------!
      ! < a^{\dagger}_{j} \sigma > !
      !----------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - (kappa(a) - i * wl(j, a))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -0.5d0 * CONJG(gkl(j, a)) * (sigma_ss(sz) + 1.0d0)
      B_vec(2) = 0.0d0
      B_vec(3) = -gamma * cav1_ss(j, at) + &
               & CONJG(gkl(j, a)) * sigma_ss(sp)

      ! Set inverse matrix
      Mat_inv = Mat
      ! Perform LU-factorization of matrix
      CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
      ! Invert Matrix (Optimal LWORK = 3)
      LWORK = 3
      CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

      ! Calculate steady state
      cavsig2_ss(j, at, :) = -MATMUL(Mat_inv, B_vec)

      !------------------!
      ! < b_{j} \sigma > !
      !------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - (kappa(b) + i * wl(j, b))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = 0.0d0
      B_vec(2) = -0.5d0 * gkl(j, b) * (sigma_ss(sz) + 1.0d0)
      B_vec(3) = -gamma * cav1b_ss(j, b) + &
               & gkl(j, b) * sigma_ss(sm)

      ! Set inverse matrix
      Mat_inv = Mat
      ! Perform LU-factorization of matrix
      CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
      ! Invert Matrix (Optimal LWORK = 3)
      LWORK = 3
      CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

      ! Calculate steady state
      cavsig2b_ss(j, b, :) = -MATMUL(Mat_inv, B_vec)

      !----------------------------!
      ! < b^{\dagger}_{j} \sigma > !
      !----------------------------!
      ! Set the diagonal matrix elements for M
      Mat = Mat_OG
      DO x = 1, N_mat
        Mat(x, x) = Mat(x, x) - (kappa(b) - i * wl(j, b))
      END DO

      ! Set the non-homogeneous vector
      B_vec = 0.0d0
      B_vec(1) = -0.5d0 * CONJG(gkl(j, b)) * (sigma_ss(sz) + 1.0d0)
      B_vec(2) = 0.0d0
      B_vec(3) = -gamma * cav1b_ss(j, bt) + &
               & CONJG(gkl(j, b)) * sigma_ss(sp)

      ! Set inverse matrix
      Mat_inv = Mat
      ! Perform LU-factorization of matrix
      CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
      ! Invert Matrix (Optimal LWORK = 3)
      LWORK = 3
      CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

      ! Calculate steady state
      cavsig2b_ss(j, bt, :) = -MATMUL(Mat_inv, B_vec)
      ! Close j loop
    END DO

    moment_out = 0.0d0
    moment_out2 = 0.0d0
    ! Cycle through modes
    DO k = -N, N
      DO j = -N, N
        !------------------------------!
        !     SECOND-ORDER: CAVITY     !
        !------------------------------!
        !-----------------!
        ! < a_{j} b_{k} > !
        !-----------------!
        cav2_ss(j, k, ab) = -gkl(j, a) * cavsig2b_ss(k, b, sm) + &
                          & -gkl(k, b) * cavsig2_ss(j, a, sm)
        cav2_ss(j, k, ab) = cav2_ss(j, k, ab) / &
                          & (kappa(a) + kappa(b) + i * (wl(j, a) + wl(k, b)))

        !-------------------------------------!
        ! < a^{\dagger}_{j} b^{\dagger}_{k} > !
        !-------------------------------------!
        cav2_ss(j, k, atbt) = -CONJG(gkl(j, a)) * cavsig2b_ss(k, bt, sp) + &
                           & -CONJG(gkl(k, b)) * cavsig2_ss(j, at, sp)
        cav2_ss(j, k, atbt) = cav2_ss(j, k, atbt) / &
                            & (kappa(a) + kappa(b) - i * (wl(j, a) + wl(k, b)))

        !---------------------------!
        ! < a^{\dagger}_{j} a_{k} > !
        !---------------------------!
        cav2_ss(j, k, ata) = -CONJG(gkl(j, a)) * cavsig2_ss(k, a, sp) + &
                           & -gkl(k, a) * cavsig2_ss(j, at, sm)
        cav2_ss(j, k, ata) = cav2_ss(j, k, ata) / &
                           & ((2.0d0 * kappa(a)) - i * (wl(j, a) - wl(k, a)))
        ! Update photon number
        moment_out = moment_out + cav2_ss(j, k, ata)

        !---------------------------!
        ! < b^{\dagger}_{j} b_{k} > !
        !---------------------------!
        cav2_ss(j, k, btb) = -CONJG(gkl(j, b)) * cavsig2b_ss(k, b, sp) + &
                           & -gkl(k, b) * cavsig2b_ss(j, bt, sm)
        cav2_ss(j, k, btb) = cav2_ss(j, k, btb) / &
                           & ((2.0d0 * kappa(b)) - i * (wl(j, b) - wl(k, b)))
        ! Update photon number
        moment_out2 = moment_out2 + cav2_ss(j, k, btb)

        !---------------------------!
        ! < a^{\dagger}_{j} b_{k} > !
        !---------------------------!
        cav2_ss(j, k, atb) = -CONJG(gkl(j, a)) * cavsig2b_ss(k, b, sp) + &
                           & -gkl(k, b) * cavsig2_ss(j, at, sm)
        cav2_ss(j, k, atb) = cav2_ss(j, k, atb) / &
                           & (kappa(a) + kappa(b) - i * (wl(j, a) - wl(k, b)))

        !---------------------------!
        ! < b^{\dagger}_{j} a_{k} > !
        !---------------------------!
        cav2_ss(j, k, bta) = -CONJG(gkl(j, b)) * cavsig2_ss(k, a, sp) + &
                           & -gkl(k, a) * cavsig2b_ss(j, bt, sm)
        cav2_ss(j, k, bta) = cav2_ss(j, k, bta) / &
                           & (kappa(a) + kappa(b) - i * (wl(j, b) - wl(k, a)))

        !--------------------------------------!
        !     THIRD-ORDER: CAVITY AND ATOM     !
        !--------------------------------------!
        !------------------------!
        ! < a_{j} b_{k} \sigma > !
        !------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - (kappa(a) + kappa(b) + i * (wl(j, a) + wl(k, b)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = 0.0d0
        B_vec(2) = -0.5d0 * gkl(j, a) * cavsig2b_ss(k, b, sz) + &
                 & -0.5d0 * gkl(j, a) * cav1b_ss(k, b) + &
                 & -0.5d0 * gkl(k, b) * cavsig2_ss(j, a, sz) + &
                 & -0.5d0 * gkl(k, b) * cav1_ss(j, a)
        B_vec(3) = -gamma * cav2_ss(j, k, ab) + &
                 & gkl(j, a) * cavsig2b_ss(k, b, sm) + &
                 & gkl(k, b) * cavsig2_ss(j, a, sm)

        ! Set inverse matrix
        Mat_inv = Mat
        ! Perform LU-factorization of matrix
        CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
        ! Invert Matrix (Optimal LWORK = 3)
        LWORK = 3
        CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

        ! Calculate steady state
        cavsig3_ss(j, k, ab, :) = -MATMUL(Mat_inv, B_vec)

        !--------------------------------------------!
        ! < a^{\dagger}_{j} b^{\dagger}_{k} \sigma > !
        !--------------------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - (kappa(a) + kappa(b) - i * (wl(j, a) + wl(k, b)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -0.5d0 * CONJG(gkl(j, a)) * cavsig2b_ss(k, bt, sz) + &
                 & -0.5d0 * CONJG(gkl(j, a)) * cav1b_ss(k, bt) + &
                 & -0.5d0 * CONJG(gkl(k, b)) * cavsig2_ss(j, at, sz) + &
                 & -0.5d0 * CONJG(gkl(k, b)) * cav1_ss(j, at)
        B_vec(2) = 0.0d0
        B_vec(3) = -gamma * cav2_ss(j, k, atbt) + &
                 & CONJG(gkl(j, a)) * cavsig2b_ss(k, bt, sp) + &
                 & CONJG(gkl(k, b)) * cavsig2_ss(j, at, sp)

        ! Set inverse matrix
        Mat_inv = Mat
        ! Perform LU-factorization of matrix
        CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
        ! Invert Matrix (Optimal LWORK = 3)
        LWORK = 3
        CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

        ! Calculate steady state
        cavsig3_ss(j, k, atbt, :) = -MATMUL(Mat_inv, B_vec)

        !----------------------------------!
        ! < a^{\dagger}_{j} a_{k} \sigma > !
        !----------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa(a)) - i * (wl(j, a) - wl(k, a)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -0.5d0 * CONJG(gkl(j, a)) * cavsig2_ss(k, a, sz) + &
                 & -0.5d0 * CONJG(gkl(j, a)) * cav1_ss(k, a)
        B_vec(2) = -0.5d0 * gkl(k, a) * cavsig2_ss(j, at, sz) + &
                 & -0.5d0 * gkl(k, a) * cav1_ss(j, at)
        B_vec(3) = -gamma * cav2_ss(j, k, ata) + &
                 & CONJG(gkl(j, a)) * cavsig2_ss(k, a, sp) + &
                 & gkl(k, a) * cavsig2_ss(j, at, sm)

        ! Set inverse matrix
        Mat_inv = Mat
        ! Perform LU-factorization of matrix
        CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
        ! Invert Matrix (Optimal LWORK = 3)
        LWORK = 3
        CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

        ! Calculate steady state
        cavsig3_ss(j, k, ata, :) = -MATMUL(Mat_inv, B_vec)

        !----------------------------------!
        ! < b^{\dagger}_{j} b_{k} \sigma > !
        !----------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa(b)) - i * (wl(j, b) - wl(k, b)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -0.5d0 * CONJG(gkl(j, b)) * cavsig2b_ss(k, b, sz) + &
                 & -0.5d0 * CONJG(gkl(j, b)) * cav1b_ss(k, b)
        B_vec(2) = -0.5d0 * gkl(k, b) * cavsig2b_ss(j, bt, sz) + &
                 & -0.5d0 * gkl(k, b) * cav1b_ss(j, bt)
        B_vec(3) = -gamma * cav2_ss(j, k, btb) + &
                 & CONJG(gkl(j, b)) * cavsig2b_ss(k, b, sp) + &
                 & gkl(k, b) * cavsig2b_ss(j, bt, sm)

        ! Set inverse matrix
        Mat_inv = Mat
        ! Perform LU-factorization of matrix
        CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
        ! Invert Matrix (Optimal LWORK = 3)
        LWORK = 3
        CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

        ! Calculate steady state
        cavsig3_ss(j, k, btb, :) = -MATMUL(Mat_inv, B_vec)

        !----------------------------------!
        ! < a^{\dagger}_{j} b_{k} \sigma > !
        !----------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - (kappa(a) + kappa(b) - i * (wl(j, a) - wl(k, b)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -0.5d0 * CONJG(gkl(j, a)) * cavsig2b_ss(k, b, sz) + &
                 & -0.5d0 * CONJG(gkl(j, a)) * cav1b_ss(k, b)
        B_vec(2) = -0.5d0 * gkl(k, b) * cavsig2_ss(j, at, sz) + &
                 & -0.5d0 * gkl(k, b) * cav1_ss(j, at)
        B_vec(3) = -gamma * cav2_ss(j, k, atb) + &
                 & CONJG(gkl(j, a)) * cavsig2b_ss(k, b, sp) + &
                 & gkl(k, b) * cavsig2_ss(j, at, sm)

        ! Set inverse matrix
        Mat_inv = Mat
        ! Perform LU-factorization of matrix
        CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
        ! Invert Matrix (Optimal LWORK = 3)
        LWORK = 3
        CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

        ! Calculate steady state
        cavsig3_ss(j, k, atb, :) = -MATMUL(Mat_inv, B_vec)

        !----------------------------------!
        ! < b^{\dagger}_{j} a_{k} \sigma > !
        !----------------------------------!
        ! Set the diagonal matrix elements for M
        Mat = Mat_OG
        DO x = 1, N_mat
          Mat(x, x) = Mat(x, x) - (kappa(a) + kappa(b) - i * (wl(j, b) - wl(k, a)))
        END DO

        ! Set the non-homogeneous vector
        B_vec = 0.0d0
        B_vec(1) = -0.5d0 * CONJG(gkl(j, b)) * cavsig2_ss(k, a, sz) + &
                 & -0.5d0 * CONJG(gkl(j, b)) * cav1_ss(k, a)
        B_vec(2) = -0.5d0 * gkl(k, a) * cavsig2b_ss(j, bt, sz) + &
                 & -0.5d0 * gkl(k, a) * cav1b_ss(j, bt)
        B_vec(3) = -gamma * cav2_ss(j, k, bta) + &
                 & CONJG(gkl(j, b)) * cavsig2_ss(k, a, sp) + &
                 & gkl(k, a) * cavsig2b_ss(j, bt, sm)

        ! Set inverse matrix
        Mat_inv = Mat
        ! Perform LU-factorization of matrix
        CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
        ! Invert Matrix (Optimal LWORK = 3)
        LWORK = 3
        CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

        ! Calculate steady state
        cavsig3_ss(j, k, bta, :) = -MATMUL(Mat_inv, B_vec)
        ! Close j loop
      END DO
      ! Close k loop
    END DO

    ! Update steady state photon number
    photon_ss = 0.0d0
    photon_ss(a) = REAL(moment_out)
    photon_ss(b) = REAL(moment_out2)

    ! Cycle through modes
    DO l = -N, N
      DO k = -N, N
        DO j = -N, N
          !-----------------------------!
          !     THIRD-ORDER: CAVITY     !
          !-----------------------------!
          !---------------------------------!
          ! < a^{\dagger}_{j} a_{k} b_{l} > !
          !---------------------------------!
          cav3_ss(j, k, l, atab) = -CONJG(gkl(j, a)) * cavsig3_ss(k, l, ab, sp) + &
                                 & -gkl(k, a) * cavsig3_ss(j, l, atb, sm) + &
                                 & -gkl(l, b) * cavsig3_ss(j, k, ata, sm)
          cav3_ss(j, k, l, atab) = cav3_ss(j, k, l, atab) / &
                                 & ((2.0d0 * kappa(a)) + kappa(b) - i * (wl(j, a) - wl(k, a) - wl(l, b)))

          !---------------------------------!
          ! < b^{\dagger}_{j} b_{k} a_{l} > !
          !---------------------------------!
          cav3_ss(j, k, l, btba) = -CONJG(gkl(j, b)) * cavsig3_ss(l, k, ab, sp) + &
                                 & -gkl(k, b) * cavsig3_ss(j, l, bta, sm) + &
                                 & -gkl(l, a) * cavsig3_ss(j, k, btb, sm)
          cav3_ss(j, k, l, btba) = cav3_ss(j, k, l, btba) / &
                                 & (kappa(a) + (2.0d0 * kappa(b)) - i * (wl(j, b) - wl(k, b) - wl(l, a)))

          !-------------------------------------------!
          ! < b^{\dagger}_{j} a^{\dagger}_{k} a_{l} > !
          !-------------------------------------------!
          cav3_ss(j, k, l, btata) = -CONJG(gkl(j, b)) * cavsig3_ss(k, l, ata, sp) + &
                                  & -CONJG(gkl(k, a)) * cavsig3_ss(j, l, bta, sp) + &
                                  & -gkl(l, a) * cavsig3_ss(k, j, atbt, sm)
          cav3_ss(j, k, l, btata) = cav3_ss(j, k, l, btata) / &
                                  & ((2.0d0 * kappa(a)) + kappa(b) - i * (wl(j, b) + wl(k, a) - wl(l, a)))

          !-------------------------------------------!
          ! < a^{\dagger}_{j} b^{\dagger}_{k} b{l} > !
          !-------------------------------------------!
          cav3_ss(j, k, l, atbtb) = -CONJG(gkl(j, a)) * cavsig3_ss(k, l, btb, sp) + &
                                  & -CONJG(gkl(k, b)) * cavsig3_ss(j, l, atb, sp) + &
                                  & -gkl(l, b) * cavsig3_ss(j, k, atbt, sm)
          cav3_ss(j, k, l, atbtb) = cav3_ss(j, k, l, atbtb) / &
                                  & (kappa(a) + (2.0d0 * kappa(b)) - i * (wl(j, a) + wl(k, b) - wl(l, b)))

          !--------------------------------------!
          !     FOURTH-ORDER: CAVITY AND ATOM    !
          !--------------------------------------!
          !----------------------------------------!
          ! < a^{\dagger}_{j} a_{k} b_{l} \sigma > !
          !----------------------------------------!
          ! Set the diagonal matrix elements for M
          Mat = Mat_OG
          DO x = 1, N_mat
            Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa(a)) + kappa(b) - i * (wl(j, a) - wl(k, a) - wl(l, b)))
          END DO

          ! Set the non-homogeneous vector
          B_vec = 0.0d0
          B_vec(1) = -0.5d0 * CONJG(gkl(j, a)) * cavsig3_ss(k, l, ab, sz) + &
                   & -0.5d0 * CONJG(gkl(j, a)) * cav2_ss(k, l, ab)
          B_vec(2) = -0.5d0 * gkl(k, a) * cavsig3_ss(j, l, atb, sz) + &
                   & -0.5d0 * gkl(k, a) * cav2_ss(j, l, atb) + &
                   & -0.5d0 * gkl(l, b) * cavsig3_ss(j, k, ata, sz) + &
                   & -0.5d0 * gkl(l, b) * cav2_ss(j, k, ata)
          B_vec(3) = -gamma * cav3_ss(j, k, l, atab) + &
                   & CONJG(gkl(j, a)) * cavsig3_ss(k, l, ab, sp) + &
                   & gkl(k, a) * cavsig3_ss(j, l, atb, sm) + &
                   & gkl(l, b) * cavsig3_ss(j, k, ata, sm)

          ! Set inverse matrix
          Mat_inv = Mat
          ! Perform LU-factorization of matrix
          CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
          ! Invert Matrix (Optimal LWORK = 3)
          LWORK = 3
          CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

          ! Calculate steady state
          cavsig4_ss(j, k, l, atab, :) = -MATMUL(Mat_inv, B_vec)

          !----------------------------------------!
          ! < b^{\dagger}_{j} b_{k} a_{l} \sigma > !
          !----------------------------------------!
          ! Set the diagonal matrix elements for M
          Mat = Mat_OG
          DO x = 1, N_mat
            Mat(x, x) = Mat(x, x) - (kappa(a) + (2.0d0 * kappa(b)) - i * (wl(j, b) - wl(k, b) - wl(l, a)))
          END DO

          ! Set the non-homogeneous vector
          B_vec = 0.0d0
          B_vec(1) = -0.5d0 * CONJG(gkl(j, b)) * cavsig3_ss(l, k, ab, sz) + &
                   & -0.5d0 * CONJG(gkl(j, b)) * cav2_ss(l, k, ab)
          B_vec(2) = -0.5d0 * gkl(k, b) * cavsig3_ss(j, l, bta, sz) + &
                   & -0.5d0 * gkl(k, b) * cav2_ss(j, l, bta) + &
                   & -0.5d0 * gkl(l, a) * cavsig3_ss(j, k, btb, sz) + &
                   & -0.5d0 * gkl(l, a) * cav2_ss(j, k, btb)
          B_vec(3) = -gamma * cav3_ss(j, k, l, btba) + &
                   & CONJG(gkl(j, b)) * cavsig3_ss(l, k, ab, sp) + &
                   & gkl(k, b) * cavsig3_ss(j, l, bta, sm) + &
                   & gkl(l, a) * cavsig3_ss(j, k, btb, sm)

          ! Set inverse matrix
          Mat_inv = Mat
          ! Perform LU-factorization of matrix
          CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
          ! Invert Matrix (Optimal LWORK = 3)
          LWORK = 3
          CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

          ! Calculate steady state
          cavsig4_ss(j, k, l, btba, :) = -MATMUL(Mat_inv, B_vec)

          !-------------------------------------------!
          ! < b^{\dagger}_{j} a^{\dagger}_{k} a_{l} > !
          !-------------------------------------------!
          ! Set the diagonal matrix elements for M
          Mat = Mat_OG
          DO x = 1, N_mat
            Mat(x, x) = Mat(x, x) - ((2.0d0 * kappa(a)) + kappa(b) - i * (wl(j, b) + wl(k, a) - wl(l, a)))
          END DO

          ! Set the non-homogeneous vector
          B_vec = 0.0d0
          B_vec(1) = -0.5d0 * CONJG(gkl(j, b)) * cavsig3_ss(k, l, ata, sz) + &
                   & -0.5d0 * CONJG(gkl(j, b)) * cav2_ss(k, l, ata) + &
                   & -0.5d0 * CONJG(gkl(k, a)) * cavsig3_ss(j, l, bta, sz) + &
                   & -0.5d0 * CONJG(gkl(k, a)) * cav2_ss(j, l, bta)
          B_vec(2) = -0.5d0 * gkl(l, a) * cavsig3_ss(k, j, atbt, sz) + &
                   & -0.5d0 * gkl(l, a) * cav2_ss(k, j, atbt)
          B_vec(3) = -gamma * cav3_ss(j, k, l, btata) + &
                   & CONJG(gkl(j, b)) * cavsig3_ss(k, l, ata, sp) + &
                   & CONJG(gkl(k, a)) * cavsig3_ss(j, l, bta, sp) + &
                   & gkl(l, a) * cavsig3_ss(k, j, atbt, sm)

          ! Set inverse matrix
          Mat_inv = Mat
          ! Perform LU-factorization of matrix
          CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
          ! Invert Matrix (Optimal LWORK = 3)
          LWORK = 3
          CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

          ! Calculate steady state
          cavsig4_ss(j, k, l, btata, :) = -MATMUL(Mat_inv, B_vec)

          !-------------------------------------------!
          ! < a^{\dagger}_{j} b^{\dagger}_{k} b_{l} > !
          !-------------------------------------------!
          ! Set the diagonal matrix elements for M
          Mat = Mat_OG
          DO x = 1, N_mat
            Mat(x, x) = Mat(x, x) - (kappa(a) + (2.0d0 * kappa(b)) - i * (wl(j, a) + wl(k, b) - wl(l, b)))
          END DO

          ! Set the non-homogeneous vector
          B_vec = 0.0d0
          B_vec(1) = -0.5d0 * CONJG(gkl(j, a)) * cavsig3_ss(k, l, btb, sz) + &
                   & -0.5d0 * CONJG(gkl(j, a)) * cav2_ss(k, l, btb) + &
                   & -0.5d0 * CONJG(gkl(k, b)) * cavsig3_ss(j, l, atb, sz) + &
                   & -0.5d0 * CONJG(gkl(k, b)) * cav2_ss(j, l, atb)
          B_vec(2) = -0.5d0 * gkl(l, b) * cavsig3_ss(j, k, atbt, sz) + &
                   & -0.5d0 * gkl(l, b) * cav2_ss(j, k, atbt)
          B_vec(3) = -gamma * cav3_ss(j, k, l, atbtb) + &
                   & CONJG(gkl(j, a)) * cavsig3_ss(k, l, btb, sp) + &
                   & CONJG(gkl(k, b)) * cavsig3_ss(j, l, atb, sp) + &
                   & gkl(l, b) * cavsig3_ss(j, k, atbt, sm)

          ! Set inverse matrix
          Mat_inv = Mat
          ! Perform LU-factorization of matrix
          CALL zGETRF(N_mat, N_mat, Mat_inv, N_mat, IPIV, INFO)
          ! Invert Matrix (Optimal LWORK = 3)
          LWORK = 3
          CALL zGETRI(N_mat, Mat_inv, N_mat, IPIV, WORK, LWORK, INFO)

          ! Calculate steady state
          cavsig4_ss(j, k, l, atbtb, :) = -MATMUL(Mat_inv, B_vec)

          ! Close j loop
        END DO
        ! Close k loop
      END DO
      ! Close l loop
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
            cav4_ss(j, k, l, m) = -CONJG(gkl(j, a)) * cavsig4_ss(k, l, m, btba, sp) + &
                                & -CONJG(gkl(k, b)) * cavsig4_ss(j, m, l, atab, sp) + &
                                & -gkl(l, b) * cavsig4_ss(k, j, m, btata, sm) + &
                                & -gkl(m, a) * cavsig4_ss(j, k, l, atbtb, sm)
            cav4_ss(j, k, l, m) = cav4_ss(j, k, l, m) / &
                                & (((2.0d0 * kappa(a)) + (2.0d0 * kappa(b))) - i * (wl(j, a) + wl(k, b)) + i * (wl(l, b) + wl(m, a)))

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
