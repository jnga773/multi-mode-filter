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

! To compiled the code, I use the Intel oneAPI IFORT compiler with:
!        (UNIX): ifort -O3 -qmkl -heap-arrays ./MODULE_two_filters.f90
!                  ./g2_cross_RK4.f90 -o [NAME]
!     (WINDOWS): ifort /O3 /Qmkl /heap-arrays ./MODULE_two_filters.f90
!                  ./g2_cross_RK4.f90 -o [NAME]
! where the -O3 (/O3) flag gives maximum optimisation, the -o (/o) [NAME] flag
! names the executable as "[NAME]" ("[NAME].exe"), the -qmkl (/Qmkl) flag links
! the program to Intel's Math Kernel Library, to make use of the LAPACK routines,
! and the -heap-arrays (/heap-arrays) lets arrays to be dumped to memory, allowing
! for larger values of N.

! You can also compile it with GFORTRAN (provided you have LAPACK and BLAS
! installed correctly) with:
!    gfortran -O3 ./MODULE_two_filters.f90 ./g2_cross_RK4.f90 -o [NAME]
!               -I/path/to/LAPACK -L/path/to/LAPACK -llapack -lblas

! In order for the program to compile, the module file
! [./MODULE_two_filters.f90] must be added to the compilation BEFORE this code.

PROGRAM PARAMETRIC_OSCILLATOR_TWO_FILTER_CROSS_CORRELATION

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
! Number of mode either side of w0, 2N + 1 total mode
INTEGER                                                  :: N
! Halfwidth of filter (\kappa if N = 0, N \delta\omega otherwise)
REAL(KIND=8)                                             :: halfwidth
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                             :: kappa_f
! Frequency spacing of modes
REAL(KIND=8)                                             :: dw
! Phase modulation of mode coupling
REAL(KIND=8)                                             :: phase

! Central mode frequency of the filter cavity, with N mode frequencies either side
REAL(KIND=8)                                             :: w0a, w0b

! Percentage of fluorecence aimed at cavity
REAL(KIND=8), PARAMETER                                  :: epsilon = 1.0d0

! Time stuff
! Time step
REAL(KIND=8)                                             :: dt
! Maximum time to integrate for
REAL(KIND=8)                                             :: t_max, tau1_max, tau2_max
! Maximum number of steps to integrate for
INTEGER                                                  :: tau_steps
! Runtime variables
REAL(KIND=8)                                             :: start_time, end_time

!----------------------------!
!     OTHER USEFUL STUFF     !
!----------------------------!
! Correlation data array
REAL(KIND=8), DIMENSION(:), ALLOCATABLE                  :: g2_positive
! Integer indices for the filter modes F and G
INTEGER, PARAMETER                                       :: f = 1, g = 2

!------------------------!
!     FILENAME STUFF     !
!------------------------!
! Paramert Name List
CHARACTER(LEN=99), PARAMETER :: filename_ParamList = "./ParamList.nml"
! Filename of parameters
CHARACTER(LEN=99), PARAMETER :: filename_parameters = "./data_files/g2_cross_parameters.txt"
! Filename for state population
CHARACTER(LEN=99), PARAMETER :: filename_g2 = "./data_files/g2_cross_corr.txt"

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
NAMELIST /FILTER/ N, halfwidth, kappa_f, phase
NAMELIST /CAVITYA/ w0a
NAMELIST /CAVITYB/ w0b
NAMELIST /TIME/ dt, t_max, tau1_max, tau2_max

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

READ(IUNIT, NML=TIME, IOSTAT=ISTAT)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') LINE
  CLOSE(IUNIT)
  PRINT *, "Invalid line in TIME namelist: " // TRIM(line)
  CALL EXIT(1)
END IF

CLOSE(IUNIT)

! Set system parameters
IF (N .EQ. 0) THEN
  ! Single-mode
  ! \kappa = the halfwidth
  kappa_f = halfwidth
  ! Set dw to zero
  dw = 0.0d0
ELSE
  ! Multi-mode
  ! Set dw = halfwidth / N
  dw = halfwidth / DBLE(N)
END IF

!==============================================================================!
!                           WRITE PARAMETERS TO FILE                           !
!==============================================================================!

! Open file to write time to
OPEN(UNIT=1, FILE=filename_parameters, STATUS='REPLACE', ACTION='WRITE')

! Write parameters
WRITE(1,"(A10,F25.15)") "kappa_p =", kappa_p
WRITE(1,"(A10,F25.15)") "Delta =", Delta
WRITE(1,"(A10,F25.15)") "lambda =", lambda

WRITE(1,"(A15,I9)") "N = ", N
WRITE(1,"(A15,F25.15)") "halfwidth =", halfwidth
WRITE(1,"(A15,F25.15)") "kappa_f =", kappa_f
WRITE(1,"(A15,F25.15)") "dw =", dw
! WRITE(1,"(A15,F25.15)") "m =", phase

WRITE(1,"(A15,F25.15)") "w0a =", w0a
WRITE(1,"(A15,F25.15)") "w0b =", w0b

! WRITE(1,"(A10,F25.15)") "dt =", dt
! WRITE(1,"(A10,F25.15)") "Max t =", t_max
! WRITE(1,"(A10,F25.15)") "Max tau1 =", tau1_max
! WRITE(1,"(A10,F25.15)") "Max tau2 =", tau2_max

! Close file
CLOSE(1)

!==============================================================================!
!                  CALCULATE SECOND-ORDER CORRELATION FUNCTION                 !
!==============================================================================!
! Number of time-steps
tau_steps = NINT(tau2_max / dt)
CALL G2_CalculateRK4(kappa_p, Delta, lambda, &
                   & epsilon, N, phase, &
                   & w0a, kappa_f, dw, &     ! Cavity F
                   & w0b, kappa_f, dw, &     ! Cavity G (with same kappa
                   & dt, tau_steps, &          !           and dw)
                   & g2_positive, .TRUE., filename_g2)

!==============================================================================!
!                                END OF PROGRAM                                !
!==============================================================================!

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

END PROGRAM PARAMETRIC_OSCILLATOR_TWO_FILTER_CROSS_CORRELATION