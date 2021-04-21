! This program will test solve for the moment equations of the two-photon
! resonantly driven three-level ladder-type atom, with equation:
!     d/dt A = M A + B,
! where A, M, and B are defined in a TeX document.
! This equation has general solution:
!     A(t) = e^{M t} A(0) + (1 - e^{M t}) A_{ss},
! which can be found with the eigenvalues and eigenvectors of matrix M.

! Must be compiled with the -mkl (/Qmkl if Windows) flag.

PROGRAM THREE_LEVEL_ATOM_FILTER_ARRAY_MOMENT_EQUATIONS_RK4

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
REAL(KIND=8)                                           :: w0
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                           :: kappa
! Percentage of fluorecence aimed at cavity
REAL(KIND=8)                                           :: epsilon
! Number of mode either side of w0, 2N + 1 total mode
INTEGER                                                :: N
! Frequency spacing of modes
REAL(KIND=8)                                           :: dw
! Phase modulation of mode coupling
INTEGER                                                :: phase
! List of Delta values
REAL(KIND=8), DIMENSION(:), ALLOCATABLE                :: wl
! List of mode dependent cascade coupling values
COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE             :: gkl
! Blackman window coefficient
REAL(KIND=8)                                           :: blackman

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
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: B, B_OG

! First-order moments: Atomic equations (< \sigma >)
COMPLEX(KIND=8), DIMENSION(N_mat)                      :: sigma, sigma_ss
! First-order moments: Cavity (< a >, < a^{\dagger} >)
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: cav1, cav1_ss
! Second-order moments: Cavity (< a^{\dagger} a >)
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cav2, cav2_ss
! Second-order moments: Cavity and atom (< a \sigma >, < a^{\dagger} \sigma >
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: cavsig2, cavsig2_ss
! Third-order moments: Cavity and atom (< a^{2} \sigma >, < a^{\dagger 2} \sigma >, < a^{\dagger} a \sigma >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cavsig3, cavsig3_ss
! Third-order moments: Cavity (< a^{2} a^{\dagger} >, < a^{\dagger 2} a >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cav3, cav3_ss
! Fourth-order moments: Cavity and atom ( < a^{\dagger} a^{2} \sigma >, < a^{\dagger 2} a \sigma >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: cavsig4, cavsig4_ss
! Fourth-order moments: Cavity (< a^{\dagger 2} a^{2} >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: cav4, cav4_ss

! Integer indices for sigma operators
INTEGER, PARAMETER                                     :: sm = 1, sp = 2, sz = 3
! Integer indices for: a, a^{\dagger}, a^{\dagger} a
INTEGER, PARAMETER                                     :: a = 1, at = 2, ata = 3

!----------------------------!
!     OTHER USEFUL STUFF     !
!----------------------------!
! Imaginary i
COMPLEX(KIND=8), PARAMETER                             :: i = CMPLX(0.0d0, 1.0d0, 8)
! pi
REAL(KIND=8), PARAMETER                                :: pi = 3.1415926535897932384d0
! 1 / 6
REAL(KIND=8), PARAMETER                                :: xis = 1.0d0 / 6.0d0
! Data
REAL(KIND=8)                                           :: popg, pope, popf, photon
! Complex data
COMPLEX(KIND=8)                                        :: moment_out
! Integer counter
INTEGER                                                :: j, k, l, m, x, y
! Temporal complex values
COMPLEX(KIND=8)                                        :: alpha

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
INTEGER :: ISTAT, IUNIT
! Line to be read from file
CHARACTER(LEN=512) :: LINE
! Namelist parameters
NAMELIST /ATOM/ gamma, Omega
NAMELIST /FILTER/ w0, kappa, dw, epsilon, N, phase
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
! First-order: Cavity
ALLOCATE(cav1(-N:N, 2)); cav1 = 0.0d0
ALLOCATE(cav1_ss(-N:N, 2)); cav1_ss = 0.0d0

! Second-order: Cavity and Atom
ALLOCATE(cavsig2(-N:N, 2, N_mat)); cavsig2 = 0.0d0
ALLOCATE(cavsig2_ss(-N:N, 2, N_mat)); cavsig2_ss = 0.0d0

! Second-order: Cavity
ALLOCATE(cav2(-N:N, -N:N, 3)); cav2 = 0.0d0
ALLOCATE(cav2_ss(-N:N, -N:N, 3)); cav2_ss = 0.0d0

! Third-order: Cavity and Atom
ALLOCATE(cavsig3(-N:N, -N:N, 3, N_mat)); cavsig3 = 0.0d0
ALLOCATE(cavsig3_ss(-N:N, -N:N, 3, N_mat)); cavsig3_ss = 0.0d0

! Third-order: Cavity
ALLOCATE(cav3(-N:N, -N:N, -N:N, 2)); cav3 = 0.0d0
ALLOCATE(cav3_ss(-N:N, -N:N, -N:N, 2)); cav3_ss = 0.0d0

! Fourth-order: Cavity and atom
ALLOCATE(cavsig4(-N:N, -N:N, -N:N, 2, N_mat)); cavsig4 = 0.0d0
ALLOCATE(cavsig4_ss(-N:N, -N:N, -N:N, 2, N_mat)); cavsig4_ss = 0.0d0

! Fourth-order: Cavity
ALLOCATE(cav4(-N:N, -N:N, -N:N, -N:N)); cav4 = 0.0d0
ALLOCATE(cav4_ss(-N:N, -N:N, -N:N, -N:N)); cav4_ss = 0.0d0

!----------------------------!
!     INITIAL CONDITIONS     !
!----------------------------!
sigma = 0.0d0
! Atom in ground state, cavity in vacuum.
sigma(sz) = -1.0d0

!==============================================================================!
!                           WRITE PARAMETERS TO FILE                           !
!==============================================================================!

! Open file to write time to
OPEN(UNIT=1, FILE=filename_parameters, STATUS='REPLACE', ACTION='WRITE')
! Write parameter
WRITE(1,*) "Parameters are in the following order:"
WRITE(1,"(A10,F25.15)") "gamma =", gamma
WRITE(1,"(A10,F25.15)") "Omega =", Omega
WRITE(1,"(A10,F25.15)") "w0 =", w0
WRITE(1,"(A10,F25.15)") "kappa =", kappa
WRITE(1,"(A11,F25.15)") "dw = ", dw
WRITE(1,"(A10,F25.15)") "epsilon =", epsilon
WRITE(1,"(A11,I9)") "N = ", N
WRITE(1,"(A11,I9)") "phase = ", phase
WRITE(1,"(A10,F25.15)") "dt =", dt
WRITE(1,"(A10,F25.15)") "Max time =", t_max
! Close file
CLOSE(1)

!==============================================================================!
!                        CALCULATE STEADY-STATE MOMENTS                        !
!==============================================================================!
!---------------------------!
!     FIRST-ORDER: ATOM     !
!---------------------------!
!---------------------------!
!     FIRST-ORDER: ATOM     !
!---------------------------!
! Calculate steady state expressions
! First-order: Bloch equations
sigma_ss(sm) = -i * gamma * Omega / ((2.0d0 * Omega ** 2) + (gamma ** 2))
sigma_ss(sp) = i * gamma * Omega / ((2.0d0 * Omega ** 2) + (gamma ** 2))
sigma_ss(sz) = -(gamma ** 2) / ((2.0d0 * Omega ** 2) + (gamma ** 2))
print*, "sz_ss = ", sigma_ss(sz)
! Cycle through modes
DO j = -N, N
  !-----------------------------!
  !     First-order: Cavity     !
  !-----------------------------!
  cav1_ss(j, a) = - gkl(j) * sigma_ss(sm) / (kappa + i * wl(j))
  cav1_ss(j, at) = - CONJG(gkl(j)) * sigma_ss(sp) / (kappa - i * wl(j))
  !---------------------------------------!
  !     Second-order: Cavity and atom     !
  !---------------------------------------!
  ! <a \sigma>
  alpha = -(0.5d0 * gamma + kappa + i * wl(j))
  ! Set matrix elements
  Mat_inv = INVERT_MATRIX(gamma, Omega, alpha)
  B = 0.0d0
  B(1) = 0.0d0
  B(2) = -0.5 * gkl(j) * (sigma_ss(sz) + 1.0d0)
  B(3) = -gamma * cav1_ss(j, a) + &
       & gkl(j) * sigma_ss(sm)
  ! Matrix multiplication
  cavsig2_ss(j, a, :) = -MATMUL(Mat_inv, B)

  ! <at \sigma>
  alpha = -(0.5d0 * gamma + kappa - i * wl(j))
  ! Set matrix elements
  Mat_inv = INVERT_MATRIX(gamma, Omega, alpha)
  B = 0.0d0
  B(1) = -0.5 * CONJG(gkl(j)) * (sigma_ss(sz) + 1.0d0)
  B(2) = 0.0d0
  B(3) = -gamma * cav1_ss(j, at) + &
       & CONJG(gkl(j)) * sigma_ss(sp)
  ! Matrix multiplication
  cavsig2_ss(j, at, :) = -MATMUL(Mat_inv, B)
END DO

! Cycle through modes
DO j = -N, N
  DO k = -N, N
    !------------------------------!
    !     Second-order: Cavity     !
    !------------------------------!
    ! <a a>
    cav2_ss(j, k, a) = -(gkl(j) * cavsig2_ss(k, a, sm) + gkl(k) * cavsig2_ss(j, a, sm)) / &
                     & (2.0d0 * kappa + i * (wl(j) + wl(k)))
    ! <at at>
    cav2_ss(j, k, at) = -(CONJG(gkl(j)) * cavsig2_ss(k, at, sp) + CONJG(gkl(k)) * cavsig2_ss(j, at, sp)) / &
                      & (2.0d0 * kappa - i * (wl(j) + wl(k)))
    ! <at a>
    cav2_ss(j, k, ata) = -(CONJG(gkl(j)) * cavsig2_ss(k, a, sp) + gkl(k) * cavsig2_ss(j, at, sm)) / &
                       & (2.0d0 * kappa - i * (wl(j) - wl(k)))

    !--------------------------------------!
    !     Third-Order: Cavity and Atom     !
    !--------------------------------------!
    ! <a a \sigma>
    alpha = -(0.5d0 * gamma + 2.0d0 * kappa + i * (wl(j) + wl(k)))
    ! Set Matrix elements
    Mat_inv = INVERT_MATRIX(gamma, Omega, alpha)
    B(1) = 0.0d0
    B(2) = -0.5d0 * gkl(j) * (cavsig2_ss(k, a, sz) + cav1_ss(k, a)) + &
            & -0.5d0 * gkl(k) * (cavsig2_ss(j, a, sz) + cav1_ss(j, a))
    B(3) = -gamma * cav2_ss(j, k, a) + &
         & gkl(j) * cavsig2_ss(k, a, sm) + &
         & gkl(k) * cavsig2_ss(j, a, sm)
    ! Matrix multiplication
    cavsig3_ss(j, k, a, :) = -MATMUL(Mat_inv, B)

    ! <at at \sigma>
    alpha = -(0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) + wl(k)))
    ! Set Matrix elements
    Mat_inv = INVERT_MATRIX(gamma, Omega, alpha)
    B(1) = -0.5d0 * CONJG(gkl(j)) * (cavsig2_ss(k, at, sz) + cav1_ss(k, at)) + &
         & -0.5d0 * CONJG(gkl(k)) * (cavsig2_ss(j, at, sz) + cav1_ss(j, at))
    B(2) = 0.0d0
    B(3) = -gamma * cav2_ss(j, k, at) + &
         & CONJG(gkl(j)) * cavsig2_ss(k, at, sp) + &
         & CONJG(gkl(k)) * cavsig2_ss(j, at, sp)
    ! Matrix multiplication
    cavsig3_ss(j, k, at, :) = -MATMUL(Mat_inv, B)

    ! <at a \sigma>
    alpha = -(0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) - wl(k)))
    ! Set Matrix elements
    Mat_inv = INVERT_MATRIX(gamma, Omega, alpha)
    B(1) = -0.5d0 * CONJG(gkl(j)) * (cavsig2_ss(k, a, sz) + cav1_ss(k, a))
    B(2) = -0.5d0 * gkl(k) * (cavsig2_ss(j, at, sz) + cav1_ss(j, at))
    B(3) = -gamma * cav2_ss(j, k, ata) + &
         & CONJG(gkl(j)) * cavsig2_ss(k, a, sp) + &
         & gkl(k) * cavsig2_ss(j, at, sm)
    ! Matrix multiplication
    cavsig3_ss(j, k, ata, :) = -MATMUL(Mat_inv, B)
  END DO
END DO

!==============================================================================!
!                      TEST PRINT SOME STEADY STATE VALUES                     !
!==============================================================================!
WRITE(*, *) "=============================================="
WRITE(*, *) "FIRST-ORDER: ATOM"
WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< sm >_ss = ", REAL(sigma_ss(sm)), " + ", IMAG(sigma_ss(sm)), "i"
WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< sp >_ss = ", REAL(sigma_ss(sp)), " + ", IMAG(sigma_ss(sp)), "i"
WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< sz >_ss = ", REAL(sigma_ss(sz)), " + ", IMAG(sigma_ss(sz)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "FIRST-ORDER: CAVITY"
WRITE(*, '(A14,ES18.11E2 A3,ES18.11E2,A2)') " < a_0 >_ss = ", REAL(cav1_ss(0, a)), " + ", IMAG(cav1_ss(0, a)), "i"
WRITE(*, '(A14,ES18.11E2 A3,ES18.11E2,A2)') "< at_0 >_ss = ", REAL(cav1_ss(0, at)), " + ", IMAG(cav1_ss(0, at)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "SECOND-ORDER: CAVITY AND ATOM"
WRITE(*, '(A17,ES18.11E2,A3,ES18.11E2,A2)') " < a_0 sm >_ss = ", REAL(cavsig2_ss(0, a, sm)), " + ", IMAG(cavsig2_ss(0, a, sm)), "i"
WRITE(*, '(A17,ES18.11E2,A3,ES18.11E2,A2)') " < a_0 sp >_ss = ", REAL(cavsig2_ss(0, a, sp)), " + ", IMAG(cavsig2_ss(0, a, sp)), "i"
WRITE(*, '(A17,ES18.11E2,A3,ES18.11E2,A2)') "< at_0 sm >_ss = ", REAL(cavsig2_ss(0, at, sm)), " + ", IMAG(cavsig2_ss(0, at, sm)), "i"
WRITE(*, '(A17,ES18.11E2,A3,ES18.11E2,A2)') "< at_0 sp >_ss = ", REAL(cavsig2_ss(0, at, sp)), " + ", IMAG(cavsig2_ss(0, at, sp)), "i"
WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "SECOND-ORDER: CAVITY"
! WRITE(*, '(A20,ES18.11E2,A3,ES18.11E2,A1)') "   < a_0 a_1 >_ss = ", REAL(cav2_ss(0, 1, a)), " + ", AIMAG(cav2_ss(0, 1, a)), "i"
! WRITE(*, '(A20,ES18.11E2,A3,ES18.11E2,A1)') " < at_0 at_1 >_ss = ", REAL(cav2_ss(0, 1, at)), " + ", AIMAG(cav2_ss(0, 1, at)), "i"
! WRITE(*, '(A20,ES18.11E2,A3,ES18.11E2,A1)') "  < at_0 a_1 >_ss = ", REAL(cav2_ss(0, 1, ata)), " + ", AIMAG(cav2_ss(0, 1, ata)), "i"
! WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "THIRD-ORDER: CAVITY AND ATOM"
! WRITE(*, '(A27,ES18.11E2,A3,ES18.11E2,A1)') "   < a_0 a_0 |e><e| >_ss = ", REAL(cavsig3_ss(0, 0, a, ee)), " + ", AIMAG(cavsig3(0, 0, a, ee)), "i"
! WRITE(*, '(A27,ES18.11E2,A3,ES18.11E2,A1)') " < at_0 at_0 |e><e| >_ss = ", REAL(cavsig3_ss(0, 0, at, ee)), " + ", AIMAG(cavsig3(0, 0, at, ee)), "i"
! WRITE(*, '(A27,ES18.11E2,A3,ES18.11E2,A1)') "  < at_0 a_0 |e><e| >_ss = ", REAL(cavsig3_ss(0, 0, ata, ee)), " + ", AIMAG(cavsig3(0, 0, ata, ee)), "i"
! WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "THIRD-ORDER: CAVITY"
! WRITE(*, '(A24,ES18.11E2,A3,ES18.11E2,A1)') "  < at_0 a_1 a_0 >_ss = ", REAL(cav3_ss(0, 1, 0, a)), " + ", AIMAG(cav3_ss(0, 1, 0, a)), "i"
! WRITE(*, '(A24,ES18.11E2,A3,ES18.11E2,A1)') " < at_0 at_1 a_0 >_ss = ", REAL(cav3_ss(0, 1, 0, at)), " + ", AIMAG(cav3_ss(0, 1, 0, at)), "i"
! WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "FOURTH-ORDER: CAVITY AND ATOM"
! WRITE(*, '(A30,ES18.11E2,A3,ES18.11E2,A1)') " < at_0 a_1 a_0 |g><e| >_ss = ", REAL(cavsig4(0, 1, 0, a, ge)), " + ", AIMAG(cavsig4(0, 1, 0, a, ge)), "i"
! WRITE(*, '(A30,ES18.11E2,A3,ES18.11E2,A1)') "< at_0 at_1 a_0 |g><e| >_ss = ", REAL(cavsig4(0, 1, 0, at, ge)), " + ", AIMAG(cavsig4(0, 1, 0, at, ge)), "i"
! WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "FOURTH-ORDER: CAVITY"
! WRITE(*, '(A27,ES18.11E2,A3,ES18.11E2,A1)') "< at_0 at_1 a_0 a_1 >_ss = ", REAL(cav4(0, 0, 0, 0)), " + ", AIMAG(cav4(0, 0, 0, 0)), "i"
! WRITE(*, *) "=============================================="

PRINT*, " "
PRINT*, "Mean photon number in cavity A =", photon
PRINT*, " "

!==============================================================================!
!                                END OF PROGRAM                                !
!==============================================================================!

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

CONTAINS

FUNCTION INVERT_MATRIX(gamma_in, Omega_in, alpha_in)
  IMPLICIT NONE
  ! Input arguments
  REAL(KIND=8), INTENT(IN) :: gamma_in, Omega_in
  COMPLEX(KIND=8), INTENT(IN) :: alpha_in
  ! Output matrix
  COMPLEX(KIND=8), DIMENSION(3, 3) :: INVERT_MATRIX
  ! Temporal denominator
  COMPLEX(KIND=8) :: denominator
  denominator = 1.0d0 / (alpha_in * (2.0d0 * (Omega_in ** 2) - gamma_in * alpha_in + 2.0d0 * (alpha_in ** 2)))
  ! Fill matrix
  !----------------------------------------------------------------------------!
  INVERT_MATRIX(1, 1) = denominator * ((Omega_in ** 2) - gamma_in * alpha_in + 2.0d0 * (alpha_in ** 2))
  INVERT_MATRIX(1, 2) = denominator * (Omega_in ** 2)
  INVERT_MATRIX(1, 3) = denominator * (CMPLX(0.0d0, -1.0d0, 8) * Omega_in * alpha_in)
  !----------------------------------------------------------------------------!
  INVERT_MATRIX(2, 1) = denominator * (Omega_in ** 2)
  INVERT_MATRIX(2, 2) = denominator * ((Omega_in ** 2) - gamma_in * alpha_in + 2.0d0 * (alpha_in ** 2))
  INVERT_MATRIX(2, 3) = denominator * (CMPLX(0.0d0, 1.0d0, 8) * Omega_in * alpha_in)
  !----------------------------------------------------------------------------!
  INVERT_MATRIX(3, 1) = denominator * (CMPLX(0.0d0, -2.0d0, 8) * Omega_in * alpha_in)
  INVERT_MATRIX(3, 2) = denominator * (CMPLX(0.0d0, 2.0d0, 8) * Omega * alpha_in)
  INVERT_MATRIX(3, 3) = denominator * (2.0d0 * (alpha_in ** 2))
  !----------------------------------------------------------------------------!
END FUNCTION INVERT_MATRIX

END PROGRAM THREE_LEVEL_ATOM_FILTER_ARRAY_MOMENT_EQUATIONS_RK4
