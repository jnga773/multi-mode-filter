! The system in this program is a resonantly driven two-level atom that is
! coupled into a multi-mode array of filter cavities, as described in some notes
! somewhere (or with this program if this goes to anyone).

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

PROGRAM TWO_LEVEL_ATOM_MULTI_MODE_FILTER_STEADY_STATES

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
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                           :: kappaa
! Percentage of fluorecence aimed at cavity
REAL(KIND=8)                                           :: epsilon
! Number of mode either side of w0, 2N + 1 total mode
INTEGER                                                :: N
! Frequency spacing of modes
REAL(KIND=8)                                           :: dwa
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
INTEGER, PARAMETER                                     :: N_mat = 3
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
NAMELIST /ATOM/ Gamma, Omega
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
! sigma = 0.0d0
! ! Atom in ground state, cavity in vacuum.
! sigma(sz) = -1.0d0

!==============================================================================!
!                           WRITE PARAMETERS TO FILE                           !
!==============================================================================!

! ! Open file to write time to
! OPEN(UNIT=1, FILE=filename_parameters, STATUS='REPLACE', ACTION='WRITE')
! ! Write parameter
! WRITE(1,*) "Parameters are in the following order:"
! WRITE(1,"(A10,F25.15)") "gamma =", gamma
! WRITE(1,"(A10,F25.15)") "Omega =", Omega
! WRITE(1,"(A10,F25.15)") "w0 =", w0
! WRITE(1,"(A10,F25.15)") "kappa =", kappa
! WRITE(1,"(A11,F25.15)") "dw = ", dw
! WRITE(1,"(A10,F25.15)") "epsilon =", epsilon
! WRITE(1,"(A11,I9)") "N = ", N
! WRITE(1,"(A11,I9)") "phase = ", phase
! WRITE(1,"(A10,F25.15)") "dt =", dt
! WRITE(1,"(A10,F25.15)") "Max time =", t_max
! ! Close file
! CLOSE(1)

!==============================================================================!
!                        CALCULATE STEADY-STATE MOMENTS                        !
!==============================================================================!
!---------------------------!
!     FIRST-ORDER: ATOM     !
!---------------------------!
! Set matrix to be inverted
Mat_inv = Mat_OG
! Invert matrix
CALL SquareMatrixInverse(N_mat, Mat_inv)

! Calculate steady states
sigma_ss = 0.0d0
sigma_ss = -MATMUL(Mat_inv, B_OG)

! Cycle through modes
DO j = -N, N
  !-----------------------------!
  !     FIRST-ORDER: CAVITY     !
  !-----------------------------!
  !-----------!
  ! < a_{j} > !
  !-----------!
  cav1_ss(j, a) = -gkl(j) * sigma_ss(sm)
  cav1_ss(j, a) = cav1_ss(j, a) / &
                & (kappa + i * wl(j))
  !---------------------!
  ! < a^{\dagger}_{j} > !
  !---------------------!
  cav1_ss(j, at) = -CONJG(gkl(j)) * sigma_ss(sp)
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
  B_vec(1) = 0.0d0
  B_vec(2) = -0.5d0 * gkl(j) * (sigma_ss(sz) + 1.0d0)
  B_vec(3) = -gamma * cav1_ss(j, a) + &
           & gkl(j) * sigma_ss(sm)

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
  B_vec(1) = -0.5d0 * CONJG(gkl(j)) * (sigma_ss(sz) + 1.0d0)
  B_vec(2) = 0.0d0
  B_vec(3) = -gamma * cav1_ss(j, at) + &
           & CONJG(gkl(j)) * sigma_ss(sp)

  ! Set inverse matrix
  Mat_inv = Mat
  ! Invert matrix
  CALL SquareMatrixInverse(N_mat, Mat_inv)

  ! Calculate steady state
  cavsig2_ss(j, at, :) = -MATMUL(Mat_inv, B_vec)

  ! Close j loop
END DO

photon_ss = 0.0d0
! Cycle through modes
DO k = -N, N
  DO j = -N, N
    !------------------------------!
    !     SECOND-ORDER: CAVITY     !
    !------------------------------!
    !-----------------!
    ! < a_{j} a_{k} > !
    !-----------------!
    cav2_ss(j, k, a) = -gkl(j) * cavsig2_ss(k, a, sm) + &
                     & -gkl(k) * cavsig2_ss(j, a, sm)
    cav2_ss(j, k, a) = cav2_ss(j, k, a) / &
                     & (2.0d0 * kappa + i * (wl(j) + wl(k)))

    !-------------------------------------!
    ! < a^{\dagger}_{j} a^{\dagger}_{k} > !
    !-------------------------------------!
    cav2_ss(j, k, at) = -CONJG(gkl(j)) * cavsig2_ss(k, at, sp) + &
                      & -CONJG(gkl(k)) * cavsig2_ss(j, at, sp)
    cav2_ss(j, k, at) = cav2_ss(j, k, at) / &
                      & (2.0d0 * kappa - i * (wl(j) + wl(k)))

    !---------------------------!
    ! < a^{\dagger}_{j} a_{k} > !
    !---------------------------!
    cav2_ss(j, k, ata) = -CONJG(gkl(j)) * cavsig2_ss(k, a, sp) + &
                       & -gkl(k) * cavsig2_ss(j, at, sm)
    cav2_ss(j, k, ata) = cav2_ss(j, k, ata) / &
                       & (2.0d0 * kappa - i * (wl(j) - wl(k)))

    ! Update photon number
    photon_ss = photon_ss + cav2_ss(j, k, ata)

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
    B_vec(1) = 0.0d0
    B_vec(2) = -0.5d0 * gkl(j) * cavsig2_ss(k, a, sz) + &
             & -0.5d0 * gkl(j) * cav1_ss(k, a) + &
             & -0.5d0 * gkl(k) * cavsig2_ss(j, a, sz) + &
             & -0.5d0 * gkl(k) * cav1_ss(j, a)
    B_vec(3) = -gamma * cav2_ss(j, k, a) + &
             & gkl(j) * cavsig2_ss(k, a, sm) + &
             & gkl(k) * cavsig2_ss(j, a, sm)

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
    B_vec(1) = -0.5d0 * CONJG(gkl(j)) * cavsig2_ss(k, at, sz) + &
             & -0.5d0 * CONJG(gkl(j)) * cav1_ss(k, at) + &
             & -0.5d0 * CONJG(gkl(k)) * cavsig2_ss(j, at, sz) + &
             & -0.5d0 * CONJG(gkl(k)) * cav1_ss(j, at)
    B_vec(2) = 0.0d0
    B_vec(3) = -gamma * cav2_ss(j, k, at) + &
             & CONJG(gkl(j)) * cavsig2_ss(k, at, sp) + &
             & CONJG(gkl(k)) * cavsig2_ss(j, at, sp)

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
    B_vec(1) = -0.5d0 * CONJG(gkl(j)) * cavsig2_ss(k, a, sz) + &
             & -0.5d0 * CONJG(gkl(j)) * cav1_ss(k, a)
    B_vec(2) = -0.5d0 * gkl(k) * cavsig2_ss(j, at, sz) + &
             & -0.5d0 * gkl(k) * cav1_ss(j, at)
    B_vec(3) = -gamma * cav2_ss(j, k, ata) + &
             & CONJG(gkl(j)) * cavsig2_ss(k, a, sp) + &
             & gkl(k) * cavsig2_ss(j, at, sm)

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
      cav3_ss(j, k, l, a) = -CONJG(gkl(j)) * cavsig3_ss(k, l, a, sp) + &
                          & -gkl(k) * cavsig3_ss(j, l, ata, sm) + &
                          & -gkl(l) * cavsig3_ss(j, k, ata, sm)
      cav3_ss(j, k, l, a) = cav3_ss(j, k, l, a) / &
                          & (3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l)))

      !-------------------------------------------!
      ! < a^{\dagger}_{j} a^{\dagger}_{k} a_{l} > !
      !-------------------------------------------!
      cav3_ss(j, k, l, at) = -CONJG(gkl(j)) * cavsig3_ss(k, l, ata, sp) + &
                           & -CONJG(gkl(k)) * cavsig3_ss(j, l, ata, sp) + &
                           & -gkl(l) * cavsig3_ss(j, k, at, sm)
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
      B_vec(1) = -0.5d0 * CONJG(gkl(j)) * cavsig3_ss(k, l, a, sz) + &
               & -0.5d0 * CONJG(gkl(j)) * cav2_ss(k, l, a)
      B_vec(2) = -0.5d0 * gkl(k) * cavsig3_ss(j, l, ata, sz) + &
               & -0.5d0 * gkl(k) * cav2_ss(j, l, ata) + &
               & -0.5d0 * gkl(l) * cavsig3_ss(j, k, ata, sz) + &
               & -0.5d0 * gkl(l) * cav2_ss(j, k, ata)
      B_vec(3) = -gamma * cav3_ss(j, k, l, a) + &
               & CONJG(gkl(j)) * cavsig3_ss(k, l, a, sp) + &
               & gkl(k) * cavsig3_ss(j, l, ata, sm) + &
               & gkl(l) * cavsig3_ss(j, k, ata, sm)

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
      B_vec(1) = -0.5d0 * CONJG(gkl(j)) * cavsig3_ss(k, l, ata, sz) + &
               & -0.5d0 * CONJG(gkl(j)) * cav2_ss(k, l, ata) + &
               & -0.5d0 * CONJG(gkl(k)) * cavsig3_ss(j, l, ata, sz) + &
               & -0.5d0 * CONJG(gkl(k)) * cav2_ss(j, l, ata)
      B_vec(2) = -0.5d0 * gkl(l) * cavsig3_ss(j, k, at, sz) + &
               & -0.5d0 * gkl(l) * cav2_ss(j, k, at)
      B_vec(3) = -gamma * cav3_ss(j, k, l, at) + &
               & CONJG(gkl(j)) * cavsig3_ss(k, l, ata, sp) + &
               & CONJG(gkl(k)) * cavsig3_ss(j, l, ata, sp) + &
               & gkl(l) * cavsig3_ss(j, k, at, sm)

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
        cav4_ss(j, k, l, m) = -CONJG(gkl(j)) * cavsig4_ss(k, l, m, a, sp) + &
                            & -CONJG(gkl(k)) * cavsig4_ss(j, l, m, a, sp) + &
                            & -gkl(l) * cavsig4_ss(j, k, m, at, sm) + &
                            & -gkl(m) * cavsig4_ss(j, k, l, at, sm)
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
!                      TEST PRINT SOME STEADY STATE VALUES                     !
!==============================================================================!
! WRITE(*, *) "=============================================="
! WRITE(*, *) "FIRST-ORDER: ATOM"
! WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< sm >_ss = ", REAL(sigma_ss(sm)), " + ", IMAG(sigma_ss(sm)), "i"
! WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< sp >_ss = ", REAL(sigma_ss(sp)), " + ", IMAG(sigma_ss(sp)), "i"
! WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< sz >_ss = ", REAL(sigma_ss(sz)), " + ", IMAG(sigma_ss(sz)), "i"
! WRITE(*, *) "=============================================="
!
WRITE(*, *) "=============================================="
WRITE(*, *) "FIRST-ORDER: CAVITY"
WRITE(*, '(A14,ES18.11E2 A3,ES18.11E2,A2)') " < a_0 >_ss = ", REAL(cav1_ss(0, a)), " + ", IMAG(cav1_ss(0, a)), "i"
WRITE(*, '(A14,ES18.11E2 A3,ES18.11E2,A2)') "< at_0 >_ss = ", REAL(cav1_ss(0, at)), " + ", IMAG(cav1_ss(0, at)), "i"
WRITE(*, *) "=============================================="
!
! WRITE(*, *) "=============================================="
! WRITE(*, *) "SECOND-ORDER: CAVITY AND ATOM"
! WRITE(*, '(A17,ES18.11E2,A3,ES18.11E2,A2)') " < a_0 sm >_ss = ", REAL(cavsig2_ss(0, a, sm)), " + ", IMAG(cavsig2_ss(0, a, sm)), "i"
! WRITE(*, '(A17,ES18.11E2,A3,ES18.11E2,A2)') " < a_0 sp >_ss = ", REAL(cavsig2_ss(0, a, sp)), " + ", IMAG(cavsig2_ss(0, a, sp)), "i"
! WRITE(*, '(A17,ES18.11E2,A3,ES18.11E2,A2)') "< at_0 sm >_ss = ", REAL(cavsig2_ss(0, at, sm)), " + ", IMAG(cavsig2_ss(0, at, sm)), "i"
! WRITE(*, '(A17,ES18.11E2,A3,ES18.11E2,A2)') "< at_0 sp >_ss = ", REAL(cavsig2_ss(0, at, sp)), " + ", IMAG(cavsig2_ss(0, at, sp)), "i"
! WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "SECOND-ORDER: CAVITY"
! WRITE(*, '(A20,ES18.11E2,A3,ES18.11E2,A1)') "   < a_0 a_0 >_ss = ", REAL(cav2_ss(0, 0, a)), " + ", AIMAG(cav2_ss(0, 0, a)), "i"
! WRITE(*, '(A20,ES18.11E2,A3,ES18.11E2,A1)') " < at_0 at_0 >_ss = ", REAL(cav2_ss(0, 0, at)), " + ", AIMAG(cav2_ss(0, 0, at)), "i"
! WRITE(*, '(A20,ES18.11E2,A3,ES18.11E2,A1)') "  < at_0 a_0 >_ss = ", REAL(cav2_ss(0, 0, ata)), " + ", AIMAG(cav2_ss(0, 0, ata)), "i"
! WRITE(*, *) "=============================================="
!
! WRITE(*, *) "=============================================="
! WRITE(*, *) "THIRD-ORDER: CAVITY AND ATOM"
! WRITE(*, '(A23,ES18.11E2,A3,ES18.11E2,A1)') "   < a_0 a_0 sp >_ss = ", REAL(cavsig3_ss(0, 0, a, sp)), " + ", AIMAG(cavsig3_ss(0, 0, a, sp)), "i"
! WRITE(*, '(A23,ES18.11E2,A3,ES18.11E2,A1)') " < at_0 at_0 sp >_ss = ", REAL(cavsig3_ss(0, 0, at, sp)), " + ", AIMAG(cavsig3_ss(0, 0, at, sp)), "i"
! WRITE(*, '(A23,ES18.11E2,A3,ES18.11E2,A1)') "  < at_0 a_0 sp >_ss = ", REAL(cavsig3_ss(0, 0, ata, sp)), " + ", AIMAG(cavsig3_ss(0, 0, ata, sp)), "i"
! WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "THIRD-ORDER: CAVITY"
! WRITE(*, '(A24,ES18.11E2,A3,ES18.11E2,A1)') "  < at_0 a_0 a_0 >_ss = ", REAL(cav3_ss(0, 0, 0, a)), " + ", AIMAG(cav3_ss(0, 0, 0, a)), "i"
! WRITE(*, '(A24,ES18.11E2,A3,ES18.11E2,A1)') " < at_0 at_0 a_0 >_ss = ", REAL(cav3_ss(0, 0, 0, at)), " + ", AIMAG(cav3_ss(0, 0, 0, at)), "i"
! WRITE(*, *) "=============================================="
!
! WRITE(*, *) "=============================================="
! WRITE(*, *) "FOURTH-ORDER: CAVITY AND ATOM"
! WRITE(*, '(A26,ES18.11E2,A3,ES18.11E2,A1)') " < at_0 a_0 a_0 sm >_ss = ", REAL(cavsig4_ss(0, 0, 0, a, sm)), " + ", AIMAG(cavsig4_ss(0, 0, 0, a, sm)), "i"
! WRITE(*, '(A26,ES18.11E2,A3,ES18.11E2,A1)') "< at_0 at_0 a_0 sm >_ss = ", REAL(cavsig4_ss(0, 0, 0, at, sm)), " + ", AIMAG(cavsig4_ss(0, 0, 0, at, sm)), "i"
! WRITE(*, *) "=============================================="
!
WRITE(*, *) "=============================================="
WRITE(*, *) "FOURTH-ORDER: CAVITY"
WRITE(*, '(A27,ES18.11E2,A3,ES18.11E2,A1)') "< at_0 at_0 a_0 a_0 >_ss = ", REAL(cav4_ss(0, 0, 0, 0)), " + ", AIMAG(cav4_ss(0, 0, 0, 0)), "i"
WRITE(*, *) "=============================================="

PRINT*, " "
PRINT*, "Mean photon number =", photon_ss
PRINT*, " "

!==============================================================================!
!                                END OF PROGRAM                                !
!==============================================================================!

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

END PROGRAM TWO_LEVEL_ATOM_MULTI_MODE_FILTER_STEADY_STATES

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
