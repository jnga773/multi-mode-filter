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

! The normalised second-order correlation function is written to [filename_g2]
! which, by default, points to "./data_files/g2_corr.txt". The file has two
! columns:
!                            t     REAL(g2_cross)

! For the default filenames, the folder "./data_files/" and NameList file
! "./ParamList.nml" MUST EXIST IN THE WORKING DIRECTORY.

! To compiled the code, I use the Intel oneAPI IFORT compiler with:
!        (UNIX): ifort -O3 -qmkl -heap-arrays ./MODULE_two_filters.f90
!                  ./g2_cross_twotime_RK4.f90 -o [NAME]
!     (WINDOWS): ifort /O3 /Qmkl /heap-arrays ./MODULE_two_filters.f90
!                  ./g2_cross_twotime_RK4.f90 -o [NAME]
! where the -O3 (/O3) flag gives maximum optimisation, the -o (/o) [NAME] flag
! names the executable as "[NAME]" ("[NAME].exe"), the -qmkl (/Qmkl) flag links
! the program to Intel's Math Kernel Library, to make use of the LAPACK routines,
! and the -heap-arrays (/heap-arrays) lets arrays to be dumped to memory, allowing
! for larger values of N.

! You can also compile it with GFORTRAN (provided you have LAPACK and BLAS
! installed correctly) with:
!    gfortran -O3 ./MODULE_two_filters.f90 ./g2_cross_twotime_RK4.f90 -o [NAME]
!               -I/path/to/LAPACK -L/path/to/LAPACK -llapack -lblas

! In order for the program to compile, the module file
! [./MODULE_two_filters.f90] must be added to the compilation BEFORE this code.

PROGRAM THREE_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_G2_CROSS_TWOTIME

!==============================================================================!
!                    DEFINING AND DECLARING VARIABLES/ARRAYS                   !
!==============================================================================!

! Import subroutines from module_cross_subroutines.f90
USE TWO_FILTER_SUBROUTINES

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
REAL(KIND=8)                                           :: w0a, w0b

! Percentage of fluorecence aimed at cavity
REAL(KIND=8), PARAMETER                                :: epsilon = 1.0d0

! Time stuff
! Time step integer
INTEGER                                                :: t
! Time step
REAL(KIND=8)                                           :: dt
! Maximum time to integrate for
REAL(KIND=8)                                           :: t_max, tau1_max, tau2_max
! Maximum number of steps to integrate for
INTEGER                                                :: tau_steps
! Runtime variables
REAL(KIND=8)                                           :: start_time, end_time

!----------------------------!
!     OTHER USEFUL STUFF     !
!----------------------------!
! Data array
REAL(KIND=8), DIMENSION(:), ALLOCATABLE                :: g2_positive, g2_negative

!------------------------!
!     FILENAME STUFF     !
!------------------------!
! Paramert Name List
CHARACTER(LEN=15), PARAMETER :: filename_ParamList = "./ParamList.nml"
! Filename of parameters
CHARACTER(LEN=34), PARAMETER :: filename_parameters = "./data_files/g2_cross_parameters.txt"
! Filename for second-order correlation
CHARACTER(LEN=99), PARAMETER :: filename_g2 = "./data_files/g2_cross.txt"

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
NAMELIST /FILTER/ N, halfwidth, kappa, phase
NAMELIST /CAVITYA/ w0a
NAMELIST /CAVITYB/ w0b
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
tau_steps = NINT(tau2_max / dt)

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
! Data array
ALLOCATE(g2_positive(0:tau_steps)); g2_positive = 0.0d0
ALLOCATE(g2_negative(0:tau_steps)); g2_negative = 0.0d0

!==============================================================================!
!                           WRITE PARAMETERS TO FILE                           !
!==============================================================================!

! Open file to write time to
OPEN(UNIT=1, FILE=filename_parameters, STATUS='REPLACE', ACTION='WRITE')

! Write parameters
WRITE(1,"(A15,F25.15)") "Gamma =", Gamma
WRITE(1,"(A15,F25.15)") "Omega =", Omega
WRITE(1,"(A15,F25.15)") "alpha =", alpha
WRITE(1,"(A15,F25.15)") "delta =", delta
WRITE(1,"(A15,F25.15)") "xi =", xi

WRITE(1,"(A15,I9)") "N = ", N
WRITE(1,"(A15,F25.15)") "halfwidth =", halfwidth
WRITE(1,"(A15,F25.15)") "kappa =", kappa
WRITE(1,"(A15,F25.15)") "dw =", dw
! WRITE(1,"(A15,F25.15)") "m =", phase

WRITE(1,"(A15,F25.15)") "w0a =", w0a
WRITE(1,"(A15,F25.15)") "w0b =", w0b

! WRITE(1,"(A15,F25.15)") "dt =", dt
! WRITE(1,"(A15,F25.15)") "Max t =", t_max
! WRITE(1,"(A15,F25.15)") "Max tau1 =", tau1_max
! WRITE(1,"(A15,F25.15)") "Max tau2 =", tau2_max

! Close file
CLOSE(1)

!==============================================================================!
!                  CALCULATE SECOND-ORDER CORRELATION FUNCTION                 !
!==============================================================================!
! Calculate from subroutine and write data to file as it's calculated
CALL G2_CalculateRK4(Gamma, Omega, alpha, delta, xi, &
                   & epsilon, N, phase, &
                   & w0a, kappa, dw, &
                   & w0b, kappa, dw, &
                   & dt, tau_steps, &
                   & g2_positive, .TRUE., "./data_files/g2_cross_positive.txt")

! Swap w0a and w0b, calculate the negative side
CALL G2_CalculateRK4(Gamma, Omega, alpha, delta, xi, &
                   & epsilon, N, phase, &
                   & w0b, kappa, dw, &
                   & w0a, kappa, dw, &
                   & dt, tau_steps, &
                   & g2_negative, .TRUE., "./data_files/g2_cross_negative.txt")                     

!==============================================================================!
!                             WRITE DATA TO FILE                               !
!==============================================================================!

! Open file to write time and data to
OPEN(UNIT=4, FILE=filename_g2, STATUS='REPLACE', ACTION='WRITE', RECL=4000)

DO t = -tau_steps, tau_steps
  ! WRITE(4, *) DBLE(t) * dt

  IF (t .LT. 0) THEN
    ! Write the negative (flipped) g2
    WRITE(4, *) DBLE(t) * dt, g2_negative(-t)
    ! WRITE(4, *) DBLE(t) * dt, tau_steps - t

  ELSE IF (t .GE. 0) THEN
    ! Write the positive g2
    WRITE(4, *) DBLE(t) * dt, g2_positive(t)

  END IF

END DO

! Close file
CLOSE(4)

!==============================================================================!
!                                END OF PROGRAM                                !
!==============================================================================!

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

END PROGRAM THREE_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_G2_CROSS_TWOTIME