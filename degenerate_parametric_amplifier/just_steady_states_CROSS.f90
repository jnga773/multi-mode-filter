! The system in this program is a parametric oscillator in the parametric
! approximation that is coupled into a multi-mode array of filter cavities,
! as described in some notes somewhere (or with this program if this goes
! to anyone).

! After verifying the moment equations with "all_order_moments_RK4.f90" with
! QuTiP, this program calculates the steady state values for each moment by
! inverting the Lindblad matrix for the atomic and cavity-atom coupled moments.
! The steady state values are printed to the console and can be compared with
! similar outouts from the RK4 program.

! The input parameters are taken from a NameList file [filename_ParamList] which,
! by default, points to "./ParamList.nml". The code can thus be compiled once,
! and parameters can be changed in the NameList file for subsequent runs.

! For the default filenames, the folder "./data_files/" and NameList file
! "./ParamList.nml" MUST EXIST IN THE WORKING DIRECTORY.

! To compiled the code, I use the Intel Parallel Studio compiler IFORT with the
! command:
!       (LINUX): ifort -O3 -o RK4 -qmkl 
!                     ./just_steady_states_cross.f90 ./MODULE_single_filter.f90
!     (WINDOWS): ifort /O3 /o RK4 /Qmkl
!                     ./just_steady_states_cross.f90 ./MODULE_single_filter.f90
! where the -O3 (/O3) flag gives maximum optimisation, and the -o (/o) RK4 flag
! names the executable as "RK4" ("RK4.exe")

PROGRAM PARAMETRIC_OSCILLATOR_STEADY_STATES

!==============================================================================!
!                    DEFINING AND DECLARING VARIABLES/ARRAYS                   !
!==============================================================================!

! Import subroutines from module filter
USE TWO_FILTER_SUBROUTINES

IMPLICIT NONE

!---------------------------------!
!     SYSTEM PARAMETERS STUFF     !
!---------------------------------!
! Parametric cavity decay rate
REAL(KIND=8)                                             :: kappa_p
! Drive detuning from cavity resonance
REAL(KIND=8)                                             :: Delta
! Driving amplitude
REAL(KIND=8)                                             :: lambda

! Filter parameter stuff
! Percentage of fluorecence aimed at cavity
REAL(KIND=8)                                             :: epsilon
! Number of mode either side of w0, 2N + 1 total mode
INTEGER                                                  :: N
! Phase modulation of mode coupling
INTEGER                                                  :: phase

! Central mode frequency of the filter cavity, with N mode frequencies either side
REAL(KIND=8)                                             :: w0a, w0b
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                             :: kappaa, kappab
! Frequency spacing of modes
REAL(KIND=8)                                             :: dwa, dwb

! Runtime variables
REAL(KIND=8)                                             :: start_time, end_time

!------------------------------------!
!     MOMENT EQUATION ARRAY STUFF    !
!------------------------------------!
! Parametric Oscillator: Evolution matrices
! First order matrix
INTEGER, PARAMETER                                       :: N1 = 2
! Second-order matrix
INTEGER, PARAMETER                                       :: N2 = 3
! Third-order matrix
INTEGER, PARAMETER                                       :: N3 = 4
! Second-order matrix
INTEGER, PARAMETER                                       :: N4 = 5

! Steady state arrays
! Parametric oscillator moments
! Second-order: < a^{2} >, < a^{\dagger} a >, < a^{\dagger}^{2} >
COMPLEX(KIND=8), DIMENSION(N2)                           :: p2_ss
! Fourth-order: < a^{4} >, < a^{\dagger} a^{3} >, < a^{\dagger}^{2} a^{2} >,
!               < a^{\dagger}^{3} a >< a^{\dagger}^{4} >
COMPLEX(KIND=8), DIMENSION(N4)                           :: p4_ss

! First-order: Filter / First-order: Parametric
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE         :: f1p1_ss
! First-order: Filter / Third-order: Parametric
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE         :: f1p3_ss

! Second-order: Filter operators
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE         :: f2_ss
! Second-order: Filter / Second-order: Parametric
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE      :: f2p2_ss

! Third-order: Filter / First-order: Parametric
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE   :: f3p1_ss

! Fourth-order: Filter
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE      :: f4_ss

! Integer indices for sigma operators! Integer indices for sigma operators
INTEGER, PARAMETER                                       :: a = 1, at = 2
INTEGER, PARAMETER                                       :: a2 = 1, ata = 2, at2 = 3
INTEGER, PARAMETER                                       :: a3 = 1, ata2 = 2, at2a = 3, at3 = 4
INTEGER, PARAMETER                                       :: a4 = 1, ata3 = 2, at2a2 = 3, at3a = 4, at4 = 5
! Integer indices for: a, f^{\dagger}, f^{\dagger} a
INTEGER, PARAMETER                                       :: f = 1, ft = 3, g = 2, gt = 4
INTEGER, PARAMETER                                       :: fg = 1, ftgt = 2, ftf = 3, gtg = 4, ftg = 5, gtf = 6
INTEGER, PARAMETER                                       :: ftfg = 1, gtgf = 2, gtftf = 3, ftgtg = 4

! INTEGER, PARAMETER                                       :: f = 1, ft = 2,  g = 2, gt = 4
! INTEGER, PARAMETER                                       :: ff = 1, ftf = 2, ft2 = 3
! INTEGER, PARAMETER                                       :: ftf2 = 1, ft2f = 2

!----------------------------!
!     OTHER USEFUL STUFF     !
!----------------------------!
! Integer counter
INTEGER                                                  :: j, k, l, m
! Parametric: Photon number
REAL(KIND=8)                                             :: photon_p_ss
! Filter: Photon number
REAL(KIND=8), DIMENSION(2)                               :: photon_f_ss
  ! kappa_f, w0, and dw values
  REAL(KIND=8), DIMENSION(2)                             :: kappa_f, w0, dw

!------------------------!
!     FILENAME STUFF     !
!------------------------!
! Paramert Name List
CHARACTER(LEN=99), PARAMETER :: filename_ParamList = "./ParamList.nml"
! Filename of parameters
CHARACTER(LEN=99), PARAMETER :: filename_parameters = "./data_files/cross_parameters.txt"

!==============================================================================!
!                 NAMELIST AND PARAMETERS TO BE READ FROM FILE                 !
!==============================================================================!
! NameList things
! Status and unit integers
INTEGER :: ISTAT, IUNIT
! Line to be read from file
CHARACTER(LEN=512) :: LINE
! Namelist parameters
NAMELIST /PARAMETRIC/ kappa_p, Delta, lambda
NAMELIST /FILTER/ epsilon, N, phase
NAMELIST /CAVITYA/ kappaa, w0a, dwa
NAMELIST /CAVITYB/ kappab, w0b, dwb
! NAMELIST /TIME/ dt, t_max, tau1_max, tau2_max

! Call start time from CPU_TIME
CALL CPU_TIME(start_time)

! Read the parameters from the NAMELIST file
IUNIT = 420
OPEN(IUNIT, FILE=filename_ParamList, STATUS="OLD", DELIM="QUOTE")

READ(IUNIT, NML=PARAMETRIC, IOSTAT=ISTAT)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') LINE
  CLOSE(IUNIT)
  PRINT *, "Invalid line in PARAMETRIC namelist: " // TRIM(line)
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

! READ(IUNIT, NML=TIME, IOSTAT=ISTAT)
! IF (ISTAT .NE. 0) THEN
!   BACKSPACE(IUNIT)
!   READ(IUNIT, FMT='(A)') LINE
!   CLOSE(IUNIT)
!   PRINT *, "Invalid line in TIME namelist: " // TRIM(line)
!   CALL EXIT(1)
! END IF

CLOSE(IUNIT)

! Set system parameters
kappa_f = kappaa
w0 = w0a
dw = dwa

kappa_f(f) = kappaa; kappa_f(g) = kappab
w0(f) = w0a; w0(g) = w0b
dw(f) = dwa; dw(g) = dwb

!==============================================================================!
!                DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS                !
!==============================================================================!
!-------------------------------------------------------!
!     INITALISE STEADY STATE OPERATOR MOMENT ARRAYS     !
!-------------------------------------------------------!
!-----------------------------!
!     First-Order: Filter     !
!-----------------------------!
! First-order: Filter / First-order: Parametric
ALLOCATE(f1p1_ss(-N:N, 4, N1)); f1p1_ss = 0.0d0
! First-order: Filter / Third-order: Parametric
ALLOCATE(f1p3_ss(-N:N, 4, N3)); f1p3_ss = 0.0d0

!------------------------------!
!     Second-Order: Filter     !
!------------------------------!
! Second-order: Filter
ALLOCATE(f2_ss(-N:N, -N:N, 6)); f2_ss = 0.0d0
! Second-order: Filter / Second-order Parametric
ALLOCATE(f2p2_ss(-N:N, -N:N, 6, N2)); f2p2_ss = 0.0d0

!-----------------------------!
!     Third-Order: Filter     !
!-----------------------------!
! Third-order: Filter / First-order Parametric
ALLOCATE(f3p1_ss(-N:N, -N:N, -N:N, 4, N1)); f3p1_ss = 0.0d0

!------------------------------!
!     Fourth-Order: Filter     !
!------------------------------!
! Fourth-order: Filter
ALLOCATE(f4_ss(-N:N, -N:N, -N:N, -N:N)); f4_ss = 0.0d0

!==============================================================================!
!                           WRITE PARAMETERS TO FILE                           !
!==============================================================================!
! Open file to write time to
OPEN(UNIT=1, FILE=filename_parameters, STATUS='REPLACE', ACTION='WRITE')

! Write parameters
WRITE(1,"(A10,F25.15)") "kappa_p =", kappa_p
WRITE(1,"(A10,F25.15)") "Delta =", Delta
WRITE(1,"(A10,F25.15)") "lambda =", lambda

! WRITE(1,"(A10,F25.15)") "epsilon =", epsilon
WRITE(1,"(A11,I9)") "N = ", N
! WRITE(1,"(A11,I9)") "phase = ", phase

WRITE(1,"(A10,F25.15)") "w0a =", w0a
WRITE(1,"(A10,F25.15)") "kappaa =", kappaa
WRITE(1,"(A11,F25.15)") "dwa = ", dwa

WRITE(1,"(A10,F25.15)") "w0b =", w0b
WRITE(1,"(A10,F25.15)") "kappab =", kappab
WRITE(1,"(A11,F25.15)") "dwb = ", dwb

! WRITE(1,"(A10,F25.15)") "dt =", dt
! WRITE(1,"(A10,F25.15)") "Max t =", t_max
! WRITE(1,"(A10,F25.15)") "Max tau1 =", tau1_max
! WRITE(1,"(A10,F25.15)") "Max tau2 =", tau2_max

! Close file
CLOSE(1)

!==============================================================================!
!                        CALCULATE STEADY-STATE MOMENTS                        !
!==============================================================================!
CALL SteadyStateMoments(kappa_p, Delta, lambda, &
                      & epsilon, N, phase, &
                      & w0(f), kappa_f(f), dw(f), &
                      & w0(g), kappa_f(g), dw(g), &
                      & photon_f_ss, &
                      & p2_ss, p4_ss, &
                      & f1p1_ss, f1p3_ss, &
                      & f2_ss, f2p2_ss, &
                      & f3p1_ss, &
                      & f4_ss)

!==============================================================================!
!                      TEST PRINT SOME STEADY STATE VALUES                     !
!==============================================================================!
! Set some modes ow
IF (N .EQ. 0) THEN
  j = 0; k = 0; l = 0; m = 0
ELSE
  j = 1; k = -1; l = 0; m = 1
END IF

WRITE(*, *) "=============================================="
WRITE(*, '(A23, ES18.11E2)') "< F^{\dagger} F >_ss = ", photon_f_ss(f)
WRITE(*, '(A23, ES18.11E2)') "< G^{\dagger} G >_ss = ", photon_f_ss(g)
WRITE(*, *) "=============================================="

! WRITE(*, *) "=============================================="
! WRITE(*, *) "PARAMETRIC"
! WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< a2 >_ss = ", REAL(p2_ss(a2)), " + ", IMAG(p2_ss(a2)), "i"
! WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< a4 >_ss = ", REAL(p4_ss(a4)), " + ", IMAG(p4_ss(a4)), "i"
! WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "FIRST-ORDER: FILTER / FIRST-ORDER: PARAMETRIC"
WRITE(*, '(A5,I2,A11,ES18.11E2,A3,ES18.11E2,A2)') " < f_", j, " a >_ss = ", &
        & REAL(f1p1_ss(j, f, a)), " + ", IMAG(f1p1_ss(j, f, a)), "i"
WRITE(*, '(A5,I2,A11,ES18.11E2,A3,ES18.11E2,A2)') " < ft_", j, " at >_ss = ", &
        & REAL(f1p1_ss(j, ft, at)), " + ", IMAG(f1p1_ss(j, ft, at)), "i"
WRITE(*, '(A5,I2,A11,ES18.11E2,A3,ES18.11E2,A2)') "< g_", j, " a >_ss = ", &
        & REAL(f1p1_ss(j, g, a)), " + ", IMAG(f1p1_ss(j, g, a)), "i"
WRITE(*, '(A5,I2,A11,ES18.11E2,A3,ES18.11E2,A2)') "< gt_", j, " at >_ss = ", &
        & REAL(f1p1_ss(j, gt, at)), " + ", IMAG(f1p1_ss(j, gt, at)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "FIRST-ORDER: FILTER / THIRD-ORDER: PARAMETRIC"
WRITE(*, '(A6,I2,A15,ES18.11E2,A3,ES18.11E2,A2)') " < f_", j, " a3 >_ss = ", &
        & REAL(f1p3_ss(j, f, a3)), " + ", IMAG(f1p3_ss(j, f, a3)), "i"
WRITE(*, '(A6,I2,A15,ES18.11E2,A3,ES18.11E2,A2)') " < ft_", j, " ata2 >_ss = ", &
        & REAL(f1p3_ss(j, ft, ata2)), " + ", IMAG(f1p3_ss(j, ft, ata2)), "i"
WRITE(*, '(A6,I2,A15,ES18.11E2,A3,ES18.11E2,A2)') "< g_", j, " at2a >_ss = ", &
        & REAL(f1p3_ss(j, g, at2a)), " + ", IMAG(f1p3_ss(j, g, at2a)), "i"
WRITE(*, '(A6,I2,A15,ES18.11E2,A3,ES18.11E2,A2)') "< gt_", j, " at3 >_ss = ", &
        & REAL(f1p3_ss(j, gt, at3)), " + ", IMAG(f1p3_ss(j, gt, at3)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "SECOND-ORDER: FILTER"
WRITE(*, '(A6,I2,A6,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') "  < f_", j, "  g_", k, " >_ss = ", &
        & REAL(f2_ss(j, k, fg)), " + ", AIMAG(f2_ss(j, k, fg)), "i"
WRITE(*, '(A6,I2,A6,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, " gt_", k, " >_ss = ", &
        & REAL(f2_ss(j, k, ftgt)), " + ", AIMAG(f2_ss(j, k, ftgt)), "i"
WRITE(*, '(A6,I2,A6,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, " f_", k, " >_ss = ", &
        & REAL(f2_ss(j, k, ftf)), " + ", AIMAG(f2_ss(j, k, ftf)), "i"
WRITE(*, '(A6,I2,A6,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < gt_", j, " g_", k, " >_ss = ", &
        & REAL(f2_ss(j, k, gtg)), " + ", AIMAG(f2_ss(j, k, gtg)), "i"
WRITE(*, '(A6,I2,A6,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, " g_", k, " >_ss = ", &
        & REAL(f2_ss(j, k, ftg)), " + ", AIMAG(f2_ss(j, k, ftg)), "i"
WRITE(*, '(A6,I2,A6,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < gt_", j, " f_", k, " >_ss = ", &
        & REAL(f2_ss(j, k, gtf)), " + ", AIMAG(f2_ss(j, k, gtf)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "SECOND-ORDER: FILTER / SECOND-ORDER: PARAMETRIC"
WRITE(*, '(A6,I2,A4,I2,A14,ES18.11E2,A3,ES18.11E2,A1)') " < f_", j, " g_", k, " a2 >_ss = ", &
  & REAL(f2p2_ss(j, k, fg, a2)), " + ", AIMAG(f2p2_ss(j, k, fg, a2)), "i"
WRITE(*, '(A6,I2,A4,I2,A14,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, " gt_", k, " ata >_ss = ", &
  & REAL(f2p2_ss(j, k, ftgt, ata)), " + ", AIMAG(f2p2_ss(j, k, ftgt, ata)), "i"
WRITE(*, '(A6,I2,A4,I2,A14,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, "  f_", k, " at2 >_ss = ", &
  & REAL(f2p2_ss(j, k, ftf, at2)), " + ", AIMAG(f2p2_ss(j, k, ftf, at2)), "i"
WRITE(*, '(A6,I2,A4,I2,A14,ES18.11E2,A3,ES18.11E2,A1)') " < gt_", j, "  g_", k, " at2 >_ss = ", &
  & REAL(f2p2_ss(j, k, gtg, at2)), " + ", AIMAG(f2p2_ss(j, k, gtg, at2)), "i"
WRITE(*, '(A6,I2,A4,I2,A14,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, "  g_", k, " at2 >_ss = ", &
  & REAL(f2p2_ss(j, k, ftg, at2)), " + ", AIMAG(f2p2_ss(j, k, ftg, at2)), "i"
WRITE(*, '(A6,I2,A4,I2,A14,ES18.11E2,A3,ES18.11E2,A1)') " < gt_", j, "  f_", k, " at2 >_ss = ", &
  & REAL(f2p2_ss(j, k, gtf, at2)), " + ", AIMAG(f2p2_ss(j, k, gtf, at2)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "THIRD-ORDER: FILTER / FIRST-ORDER: PARAMETRIC"
WRITE(*, '(A6,I2,A4,I2,A3,I2,A11,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, "  f_", k, " g_", l,  " a >_ss = ", &
  & REAL(f3p1_ss(j, k, l, ftfg, a)), " + ", AIMAG(f3p1_ss(j, k, l, ftfg, a)), "i"
WRITE(*, '(A6,I2,A4,I2,A3,I2,A11,ES18.11E2,A3,ES18.11E2,A1)') " < gt_", j, "  g_", k, " f_", l,  " at >_ss = ", &
  & REAL(f3p1_ss(j, k, l, gtgf, at)), " + ", AIMAG(f3p1_ss(j, k, l, gtgf, at)), "i"
WRITE(*, '(A6,I2,A4,I2,A3,I2,A11,ES18.11E2,A3,ES18.11E2,A1)') " < gt_", j, " ft_", k, " f_", l, " a >_ss = ", &
  & REAL(f3p1_ss(j, k, l, gtftf, a)), " + ", AIMAG(f3p1_ss(j, k, l, gtftf, a)), "i"
WRITE(*, '(A6,I2,A4,I2,A3,I2,A11,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, " gt_", k, " g_", l, " at >_ss = ", &
  & REAL(f3p1_ss(j, k, l, ftgtg, at)), " + ", AIMAG(f3p1_ss(j, k, l, ftgtg, at)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "FOURTH-ORDER: FILTER"
WRITE(*, '(A6,I2,A4,I2,A3,I2,A3,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, " gt_", k, " g_", l, " f_", m, " >_ss = ", &
  & REAL(f4_ss(j, k, l, m)), " + ", AIMAG(f4_ss(j, k, l, m)), "i"
WRITE(*, *) "=============================================="

!==============================================================================!
!                                END OF PROGRAM                                !
!==============================================================================!

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

END PROGRAM PARAMETRIC_OSCILLATOR_STEADY_STATES
