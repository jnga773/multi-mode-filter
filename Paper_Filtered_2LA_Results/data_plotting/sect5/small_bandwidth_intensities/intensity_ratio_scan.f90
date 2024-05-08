! The system in this program is a resonantly driven two-level atom that is
! coupled into a multi-mode array of filter cavities, as described in some notes
! somewhere (or with this program if this goes to anyone).

! The operator moment steady states are calculated using the inverse matrix/
! analytical method. Using quantum regression theorem, we then use Runge-Kutta
! fourth-order to numerically integrate the moment equations, suing the steady
! states and initial conditions, to calculate the normalised second-order
! correlation function.

! The input parameters are taken from a NameList file [filename_ParamList] which,
! by default, points to "./ParamList.nml". The code can thus be compiled once,
! and parameters can be changed in the NameList file for subsequent runs.

! The parameters for each run are written to [filename_parameters] which, by
! default, points to "./data_files/g2_initial_phase_scan/parameters.txt".

! The normalised second-order correlation function is written to [filename_g2]
! which, by default, points to "./data_files/g2_initial_phase_scan/g2_corr.txt". The file has two
! columns:
!                            t     REAL(g2)

! For the default filenames, the folder "./data_files/g2_initial_phase_scan/" and NameList file
! "./ParamList.nml" MUST EXIST IN THE WORKING DIRECTORY.

! To compiled the code, I use the Intel Parallel Studio compiler IFORT with the
! command
!       (LINUX): ifort -O3 -o g2 -qmkl -heap-arrays
!                    ./MODULE_single_filter.f90 g2_RK4.f90
!     (WINDOWS): ifort /O3 /o g2 /Qmkl /heap-arrays
!                    ./MODULE_single_filter.f90 g2_RK4.f90
! where the -O3 (/O3) flag gives maximum optimisation, the -o (/o) g2 flag
! names the executable as "g2" ("g2.exe"), and the -mkl (/Qmkl) flag links
! the program to Intel's Math Kernel Library, to make use of the LAPACK routines.

PROGRAM TWO_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_G2

! Import subroutines from the module file
USE SINGLE_FILTER_SUBROUTINES

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
! Runtime variables
REAL(KIND=8)                                           :: start_time, end_time

! Scan parameters
! Starting scan value
REAL(KIND=8)                                           :: scan_start
! Final scan value
REAL(KIND=8)                                           :: scan_end
! Scan step size
INTEGER                                                :: scan_steps
! Number of scan steps
INTEGER                                                :: number_of_scans

! Scan stuff
! Array of width values to scan over
REAL(KIND=8), DIMENSION(:), ALLOCATABLE                :: halfwidth_array
! Boolean for logarithmic scale
LOGICAL, PARAMETER                                     :: halfwidth_log_scale = .TRUE.
! LOGICAL, PARAMETER                                     :: halfwidth_log_scale = .FALSE.
! Number of log-scale "steps"
INTEGER                                                :: log_scale_steps

!----------------------------!
!     OTHER USEFUL STUFF     !
!----------------------------!
! Correlation data
REAL(KIND=8)                                           :: I_ratio
REAL(KIND=8), DIMENSION(:), ALLOCATABLE                :: I_ratio_array
! Index integer
INTEGER                                                :: index
! Temporal halfwidth value
REAL(KIND=8)                                           :: halfwidth_temp, d_halfwidth
! Index integers
INTEGER                                                :: j, k, l, m
! Scan variables
REAL(KIND=8)                                           :: halfwidth_scan, kappa_scan, dw_scan
! Run counter
INTEGER                                                :: run_counter

!------------------------!
!     FILENAME STUFF     !
!------------------------!
! Paramert Name List
CHARACTER(LEN=15), PARAMETER :: filename_ParamList = "./ParamList.nml"
! Filename of parameters
CHARACTER(LEN=99)            :: filename_parameters
! Filename for second-order correlation
CHARACTER(LEN=99)            :: filename_g2

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
NAMELIST /SCANPARAMS/ scan_start, scan_end, scan_steps

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

READ(IUNIT, NML=SCANPARAMS, IOSTAT=ISTAT)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') LINE
  CLOSE(IUNIT)
  PRINT *, "Invalid line in CAVITYA namelist: " // TRIM(line)
  CALL EXIT(1)
END IF

CLOSE(IUNIT)

! ! Set system parameters
! IF (N .EQ. 0) THEN
!   ! Single-mode
!   ! \kappa = the halfwidth
!   kappa = halfwidth
!   ! Set dw to zero
!   dw = 0.0d0
! ELSE
!   ! Multi-mode
!   ! Set dw = halfwidth / N
!   dw = halfwidth / DBLE(N)

!   ! kappa = 2.5 dw
!   kappa = 2.5d0 * dw
! END IF

! Number of scan steps
! number_of_scans = NINT((scan_end - scan_start) / scan_steps)
number_of_scans = scan_steps

!==============================================================================!
!                         ALLOCATING ARRAYS AND STUFF                          !
!==============================================================================!
! Set filenames
IF (N .EQ. 0) THEN
  ! Single-mode
  IF (w0a .EQ. 0.0d0) THEN
    ! Central peak
    filename_parameters = "./data_files/centre/intensity_ratio_parameters_single.txt"
    filename_g2 = "./data_files/centre/intensity_ratio_single.txt"
  ELSE IF (w0a .EQ. Omega) THEN
    ! Right peak
    filename_parameters = "./data_files/right/intensity_ratio_parameters_single.txt"
    filename_g2 = "./data_files/right/intensity_ratio_single.txt"
  END IF
ELSE
  ! Multi-mode
  IF (w0a .EQ. 0.0d0) THEN
    ! Central peak
    filename_parameters = "./data_files/centre/intensity_ratio_parameters_multi.txt"
    filename_g2 = "./data_files/centre/intensity_ratio_multi.txt"
  ELSE IF (w0a .EQ. Omega) THEN
    ! Right peak
    filename_parameters = "./data_files/right/intensity_ratio_parameters_multi.txt"
    filename_g2 = "./data_files/right/intensity_ratio_multi.txt"
  END IF
END IF

!------------------------------------------!
!     INITALISE BANDWIDTH VALUES ARRAY     !
!------------------------------------------!
IF (halfwidth_log_scale .EQV. .TRUE.) THEN
  ! Calculate the number of log scale "steps", i.e., number of powers of 10
  ! difference between bandwidth_min and bandwidth_max
  log_scale_steps = NINT(LOG10(scan_end) - LOG10(scan_start))

  ! Number of steps
  number_of_scans = INT(log_scale_steps * DBLE(scan_steps - 1))

  ! Allocate the array
  ALLOCATE(halfwidth_array(INT(log_scale_steps * DBLE(scan_steps - 1))))
  halfwidth_array = 0.0d0
  
  ! Allocate the data array
  ALLOCATE(I_ratio_array(INT(log_scale_steps * DBLE(scan_steps - 1))))
  I_ratio_array = 0.0d0

  DO j = 0, log_scale_steps-1
    ! Calculate the bandwidth step size
    d_halfwidth = scan_start * (10.0d0 ** DBLE(j+1)) / scan_steps

    ! Cycle through number of steps for each scale
    m = scan_steps - INT(scan_steps / 10)
    DO k = 1, m
      ! Calculate bandwidth of this index
      halfwidth_temp = scan_start * (10.0d0 ** DBLE(j)) + (d_halfwidth * DBLE(k-1))
      
      ! Calculate the index of the array
      ! l = NINT(DBLE(j) * DBLE(bandwidth_steps - 1) + DBLE(k))
      l = NINT(DBLE(j) * DBLE(scan_steps) * 0.9d0 + DBLE(k))

      ! Update array
      halfwidth_array(l) = halfwidth_temp

    END DO
  END DO

ELSE IF (halfwidth_log_scale .EQV. .FALSE.) THEN
  ! Calculate the bandwidth step size
  d_halfwidth = (scan_end - scan_start) / DBLE(number_of_scans)

  ! Allocate the array
  ALLOCATE(halfwidth_array(number_of_scans))
  halfwidth_array = 0.0d0
  
  ! Allocate the data array
  ALLOCATE(I_ratio_array(number_of_scans))
  I_ratio_array = 0.0d0

  ! Set the values
  DO j = 1, number_of_scans
    halfwidth_array(j) = scan_start + d_halfwidth * DBLE(j-1)
  END DO

END IF
!==============================================================================!
!                           WRITE PARAMETERS TO FILE                           !
!==============================================================================!
! Open file to write parameters to
OPEN(UNIT=1, FILE=filename_parameters, STATUS='REPLACE', ACTION='WRITE')

! Write parameters
WRITE(1,"(A15, F25.15)") "gamma =", Gamma
WRITE(1,"(A15, F25.15)") "Omega =", Omega

WRITE(1,"(A15,I9)") "N = ", N
! WRITE(1,"(A15,F25.15)") "halfwidth =", halfwidth
! WRITE(1,"(A15,F25.15)") "kappa =", kappa
! WRITE(1,"(A15,F25.15)") "dw =", dw
! WRITE(1,"(A15,F25.15)") "m =", phase

WRITE(1,"(A15,F25.15)") "w0a =", w0a

! WRITE(1,"(A15, F25.15)") "dt =", dt
! WRITE(1,"(A15, F25.15)") "Max tau2 =", tau2_max

! Close file
CLOSE(1)

!===============================================================================!
!                  CALCULATE SECOND-ORDER CORRELATION FUNCTION                  !
!===============================================================================!
! Reset run_counter
run_counter = 0

!$OMP PARALLEL DO PRIVATE(index, halfwidth_scan, kappa_scan, dw_scan, I_ratio) SHARED(run_counter)
DO index = 0, number_of_scans
  IF (N .EQ. 0) THEN
    ! Single-mode (N = 0)
    ! Grab halfwidth value
    halfwidth_scan = halfwidth_array(index)

    ! Set width of single mode to the scan value
    kappa_scan = halfwidth_scan
  
  ELSE
    ! Multi-mode (N > 0)
    ! Grab halfwidth value
    halfwidth_scan = halfwidth_array(index)

    ! Set mode spacing to scan value / N
    dw_scan = halfwidth_scan / DBLE(N)

    ! Set kappa
    kappa_scan = kappa * dw_scan
    ! kappa_scan = 2.5d0 * dw_scan

  END IF

  ! Calculate g2
  CALL CalcIntensityRatio(gamma, Omega, &
                        & epsilon, N, phase, &
                        & w0a, kappa_scan, dw_scan, &
                        & I_ratio)

  ! Update array
  I_ratio_array(index) = I_ratio

  ! Print completion
  WRITE(*, "(I4, A3, I4, A15)") run_counter+1, " / ", number_of_scans+1, " scans complete"
  run_counter = run_counter + 1

  ! Close scan loop
END DO
!$OMP END PARALLEL DO

!----------------------------!
!     Write Data to File     !
!----------------------------!
! Open file to write data to
OPEN(UNIT=2, FILE=filename_g2, STATUS='REPLACE', ACTION='WRITE', RECL=12800000)

! Cycle through halfwidth values
DO index = 0, number_of_scans
  ! Write to file
  WRITE(2, *) halfwidth_array(index), I_ratio_array(index)

  ! Close scan loop
END DO
! Close files
CLOSE(2)

!==============================================================================!
!                                END OF PROGRAM                                !
!==============================================================================!

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

END PROGRAM TWO_LEVEL_ATOM_MULTI_MODE_FILTER_MOMENTS_G2
