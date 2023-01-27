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

! To compiled the code, I use the Intel oneAPI IFORT compiler with:
!        (UNIX): ifort -O3 -qmkl -heap-arrays ./MODULE_single_filter.f90
!                  ./just_steady_states.f90 -o [NAME]
!     (WINDOWS): ifort /O3 /Qmkl /heap-arrays ./MODULE_single_filter.f90
!                  ./just_steady_states.f90 -o [NAME]
! where the -O3 (/O3) flag gives maximum optimisation, the -o (/o) [NAME] flag
! names the executable as "[NAME]" ("[NAME].exe"), the -qmkl (/Qmkl) flag links
! the program to Intel's Math Kernel Library, to make use of the LAPACK routines,
! and the -heap-arrays (/heap-arrays) lets arrays to be dumped to memory, allowing
! for larger values of N.

! You can also compile it with GFORTRAN (provided you have LAPACK and BLAS
! installed correctly) with:
!    gfortran -O3 ./MODULE_single_filter.f90 ./just_steady_states.f90 -o [NAME]
!               -I/path/to/LAPACK -L/path/to/LAPACK -llapack -lblas

! In order for the program to compile, the module file
! [./MODULE_single_filter.f90] must be added to the compilation BEFORE this code.

PROGRAM TWO_LEVEL_ATOM_MULTI_MODE_FILTER_STEADY_STATES

!==============================================================================!
!                    DEFINING AND DECLARING VARIABLES/ARRAYS                   !
!==============================================================================!

! Import subroutines from the module file
USE SINGLE_FILTER_SUBROUTINES

IMPLICIT NONE

!---------------------------------!
!     SYSTEM PARAMETERS STUFF     !
!---------------------------------!
! Atomic decay rate
REAL(KIND=8)                                           :: gamma
! Driving amplitude
REAL(KIND=8)                                           :: Omega

! Filter parameter stuff
! Number of mode either side of w0, 2N + 1 total mode
INTEGER                                                :: N
! Halfwidth of filter (\kappa if N = 0, N \delta\omega otherwise)
REAL(KIND=8)                                           :: halfwidth
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                           :: kappa
! Frequency spacing of modes
REAL(KIND=8)                                           :: dw
! Phase modulation of mode coupling
REAL(KIND=8)                                           :: phase

! Central mode frequency of the filter cavity, with N mode frequencies either side
REAL(KIND=8)                                           :: w0a

! Percentage of fluorecence aimed at cavity
REAL(KIND=8), PARAMETER                                :: epsilon = 1.0d0
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
! First-order moments: Cavity (< a >, < f^{\dagger} >)
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE          :: f1_ss
! Second-order moments: Cavity and atom (< a \sigma >, < f^{\dagger} \sigma >
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: f1sig_ss
! Second-order moments: Cavity (< f^{\dagger} a >)
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE       :: f2_ss
! Third-order moments: Cavity and atom (< a^{2} \sigma >, < a^{\dagger 2} \sigma >, < f^{\dagger} a \sigma >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: f2sig_ss
! Third-order moments: Cavity (< a^{2} f^{\dagger} >, < a^{\dagger 2} a >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: f3_ss
! Fourth-order moments: Cavity and atom ( < f^{\dagger} a^{2} \sigma >, < a^{\dagger 2} a \sigma >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: f3sig_ss
! Fourth-order moments: Cavity (< a^{\dagger 2} a^{2} >)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE    :: f4_ss

! Integer indices for sigma operators
INTEGER, PARAMETER                                     :: sm = 1, sp = 2, sz = 3
! Integer indices for: f, f^{\dagger}, f^{2}, f^{\dagger}^{2} ... etc
INTEGER, PARAMETER                                     :: f = 1, ft = 2
INTEGER, PARAMETER                                     :: ff = 1, ftf = 2, ft2 = 3
INTEGER, PARAMETER                                     :: ftf2 = 1, ft2f = 2

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
! ! Filename of parameters
! CHARACTER(LEN=27), PARAMETER :: filename_parameters = "./data_files/parameters.txt"
! ! Filename for state population
! CHARACTER(LEN=23), PARAMETER :: filename_states = "./data_files/states.txt"
! ! Filename for operators
! CHARACTER(LEN=26), PARAMETER :: filename_cavity = "./data_files/operators.txt"

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
NAMELIST /FILTER/ N, halfwidth, kappa, phase
NAMELIST /CAVITYA/ w0a
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
IF (N .EQ. 0) THEN
  ! Single-mode
  ! \kappa = the halfwidth
  kappa = halfwidth
  ! Set dw to zero
  dw = 0.0d0
ELSE
  ! Multi-mode
  ! Set dw = halfwidth / N
  dw = halfwidth / DBLE(N)
END IF

!==============================================================================!
!                DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS                !
!==============================================================================!
!------------------------------------------!
!     INITALISE OPERATOR MOMENT ARRAYS     !
!------------------------------------------!
! Steady states
! First-order: Cavity
ALLOCATE(f1_ss(-N:N, 2)); f1_ss = 0.0d0
! Second-order: Cavity and Atom
ALLOCATE(f1sig_ss(-N:N, 2, N_mat)); f1sig_ss = 0.0d0
! Second-order: Cavity
ALLOCATE(f2_ss(-N:N, -N:N, 3)); f2_ss = 0.0d0
! Third-order: Cavity and Atom
ALLOCATE(f2sig_ss(-N:N, -N:N, 3, N_mat)); f2sig_ss = 0.0d0
! Third-order: Cavity
ALLOCATE(f3_ss(-N:N, -N:N, -N:N, 2)); f3_ss = 0.0d0
! Fourth-order: Cavity and atom
ALLOCATE(f3sig_ss(-N:N, -N:N, -N:N, 2, N_mat)); f3sig_ss = 0.0d0
! Fourth-order: Cavity
ALLOCATE(f4_ss(-N:N, -N:N, -N:N, -N:N)); f4_ss = 0.0d0

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

! ! Write parameters
! WRITE(1,"(A15,F25.15)") "gamma =", Gamma
! WRITE(1,"(A15,F25.15)") "Omega =", Omega

! WRITE(1,"(A15,I9)") "N = ", N
! WRITE(1,"(A15,F25.15)") "halfwidth =", halfwidth
! WRITE(1,"(A15,F25.15)") "kappa =", kappa
! WRITE(1,"(A15,F25.15)") "dw =", dw
! WRITE(1,"(A15,F25.15)") "m =", phase

! WRITE(1,"(A15,F25.15)") "w0a =", w0a

! ! WRITE(1,"(A15,F25.15)") "dt =", dt
! ! WRITE(1,"(A15,F25.15)") "Max t =", t_max
! ! WRITE(1,"(A15,F25.15)") "Max tau1 =", tau1_max
! ! WRITE(1,"(A15,F25.15)") "Max tau2 =", tau2_max

! ! Close file
! CLOSE(1)

!==============================================================================!
!                        CALCULATE STEADY-STATE MOMENTS                        !
!==============================================================================!
CALL SteadyStateMoments(Gamma, Omega, &
                      & epsilon, N, phase, &
                      & w0a, kappa, dw, &
                      & photon_ss, sigma_ss, .TRUE., &
                      & f1_ss, f1sig_ss, &
                      & f2_ss, f2sig_ss, &
                      & f3_ss, f3sig_ss, &
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
WRITE(*, *) "ATOM"
WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< sm >_ss = ", REAL(sigma_ss(sm)), " + ", IMAG(sigma_ss(sm)), "i"
WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< sp >_ss = ", REAL(sigma_ss(sp)), " + ", IMAG(sigma_ss(sp)), "i"
WRITE(*, '(A12,ES18.11E2,A3,ES18.11E2,A2)') "< sz >_ss = ", REAL(sigma_ss(sz)), " + ", IMAG(sigma_ss(sz)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "FIRST-ORDER: FILTER"
WRITE(*, '(A5,I2,A8,ES18.11E2 A3,ES18.11E2,A2)') " < f_", j, " >_ss = ", REAL(f1_ss(j, f)), " + ", IMAG(f1_ss(j, f)), "i"
WRITE(*, '(A5,I2,A8,ES18.11E2 A3,ES18.11E2,A2)') "< ft_", j, " >_ss = ", REAL(f1_ss(j, ft)), " + ", IMAG(f1_ss(j, ft)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "FIRST-ORDER: FILTER / FIRST-ORDER: PARAMETRIC"
WRITE(*, '(A5,I2,A11,ES18.11E2,A3,ES18.11E2,A2)') " < f_", j, " sm >_ss = ", &
        & REAL(f1sig_ss(j, f, sm)), " + ", IMAG(f1sig_ss(j, f, sm)), "i"
WRITE(*, '(A5,I2,A11,ES18.11E2,A3,ES18.11E2,A2)') " < f_", j, " sp >_ss = ", &
        & REAL(f1sig_ss(j, f, sp)), " + ", IMAG(f1sig_ss(j, f, sp)), "i"
WRITE(*, '(A5,I2,A11,ES18.11E2,A3,ES18.11E2,A2)') "< ft_", j, " sm >_ss = ", &
        & REAL(f1sig_ss(j, ft, sm)), " + ", IMAG(f1sig_ss(j, ft, sm)), "i"
WRITE(*, '(A5,I2,A11,ES18.11E2,A3,ES18.11E2,A2)') "< ft_", j, " sp >_ss = ", &
        & REAL(f1sig_ss(j, ft, sp)), " + ", IMAG(f1sig_ss(j, ft, sp)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "SECOND-ORDER: FILTER"
WRITE(*, '(A6,I2,A6,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') "  < f2_", j, "  f2_", k, " >_ss = ", &
        & REAL(f2_ss(j, k, ff)), " + ", AIMAG(f2_ss(j, k, ff)), "i"
WRITE(*, '(A6,I2,A6,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < ft2_", j, " ft2_", k, " >_ss = ", &
        & REAL(f2_ss(j, k, ft2)), " + ", AIMAG(f2_ss(j, k, ft2)), "i"
WRITE(*, '(A6,I2,A6,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < ftf_", j, " ftf_", k, " >_ss = ", &
        & REAL(f2_ss(j, k, ftf)), " + ", AIMAG(f2_ss(j, k, ftf)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "SECOND-ORDER: FILTER / ATOM"
WRITE(*, '(A7,I2,A5,I2,A11,ES18.11E2,A3,ES18.11E2,A1)') "  < f2_", j, "  f2_", k, " sm >_ss = ", &
  & REAL(f2sig_ss(j, k, ff, sm)), " + ", AIMAG(f2sig_ss(j, k, ff, sm)), "i"
WRITE(*, '(A7,I2,A5,I2,A11,ES18.11E2,A3,ES18.11E2,A1)') " < ft2_", j, " ft2_", k, " sp >_ss = ", &
  & REAL(f2sig_ss(j, k, ft2, sp)), " + ", AIMAG(f2sig_ss(j, k, ft2, sp)), "i"
WRITE(*, '(A7,I2,A5,I2,A11,ES18.11E2,A3,ES18.11E2,A1)') " < ftf_", j, " ftf_", k, " sz >_ss = ", &
  & REAL(f2sig_ss(j, k, ftf, sz)), " + ", AIMAG(f2sig_ss(j, k, ftf, sz)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "THIRD-ORDER: FILTER"
WRITE(*, '(A6,I2,A4,I2,A3,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, "  f_", k, " f_", l, " >_ss = ", &
  & REAL(f3_ss(j, k, l, ftf2)), " + ", AIMAG(f3_ss(j, k, l, ftf2)), "i"
WRITE(*, '(A6,I2,A4,I2,A3,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, " ft_", k, " f_", l, " >_ss = ", &
  & REAL(f3_ss(j, k, l, ft2f)), " + ", AIMAG(f3_ss(j, k, l, ft2f)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "THIRD-ORDER: FILTER / ATOM"
WRITE(*, '(A6,I2,A4,I2,A3,I2,A11,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, "  f_", k, " f_", l,  " sm >_ss = ", &
  & REAL(f3sig_ss(j, k, l, ftf2, sm)), " + ", AIMAG(f3sig_ss(j, k, l, ftf2, sm)), "i"
WRITE(*, '(A6,I2,A4,I2,A3,I2,A11,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, " ft_", k, " f_", l, " sz >_ss = ", &
  & REAL(f3sig_ss(j, k, l, ft2f, sz)), " + ", AIMAG(f3sig_ss(j, k, l, ft2f, sz)), "i"
WRITE(*, *) "=============================================="

WRITE(*, *) "=============================================="
WRITE(*, *) "FOURTH-ORDER: FILTER"
WRITE(*, '(A6,I2,A4,I2,A3,I2,A3,I2,A8,ES18.11E2,A3,ES18.11E2,A1)') " < ft_", j, " ft_", k, " f_", l, " f_", m, " >_ss = ", &
  & REAL(f4_ss(j, k, l, m)), " + ", AIMAG(f4_ss(j, k, l, m)), "i"
WRITE(*, *) "=============================================="

! PRINT*, " "
PRINT*, "Mean photon number in cavity A =", photon_ss
! PRINT*, " "

!==============================================================================!
!                                END OF PROGRAM                                !
!==============================================================================!

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

END PROGRAM TWO_LEVEL_ATOM_MULTI_MODE_FILTER_STEADY_STATES
