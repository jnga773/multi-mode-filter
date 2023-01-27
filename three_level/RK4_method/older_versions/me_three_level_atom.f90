! This program solves the Lindblad master equation for a driven three-level atom
! (with ground state |g>, intermediate state |e>, and final state |f>) with
! spontaneous emission into a single channel. We numerically solve the master
! equation using the Runge-Kutta 4th Order method.
!
! The Hamiltonian, in a frame rotating at the drive frequency \omega_{d}, is:
! - H = -(\alpha/2 + \delta) |e><e| - 2 \delta |f><f|
! -         + \Omega/2 (\Sigma + \Sigma^{\dagger}),
! where
! - \alpha = \omega_{fe} - \omega_{eg},
! is the anharmonicity of the atom (\omega_{ij} = \omega_{i} - \omega_{j}),
! - \delta = \omega_{d} - \omega_{fg}/2,
! is the drive detuning from the two-photon resonance frequency, \xi is the
! dipole moment ratio, and
! - \Sigma = |g><e| + \xi |e><f|
! is the atomic lowering operator. The Lindblad master equation is then
! - d/dt \rho = 1/i\hbar [H, \rho] +
!      \Gamma/2 (2 \Sigma \rho \Sigma^{\dagger} - \Sigma^{\dagger} \Sigma \rho -
!                \rho \Sigma^{\dagger} \Sigma),
! where \Gamma is the atomic decay rate.

! The input parameters (gamma, alpha, Omega, delta, xi, dt, t_max, tau1_max, tau2_max)
! are taken from a NameList file (ParamList.nml).

! The program outputs the atomic population of the ground, intermediate, and
! final states, the first-order correlation function, the second-order correlation
! function, and a parameters file.
! - ./data_files/me_parameters.txt,
! - ./data_files/me_states.txt (time, |g><g|, |e><e|, |f><f|),
! - ./data_files/g1_corr.txt (time, Real(g1), Imag(g1)),
! - ./data_files/g2_corr.txt (time, g2).

PROGRAM Master_Equation

IMPLICIT NONE

! Parameters in terms of decay rate gamma
! Decay Rate
REAL(KIND=8) :: gamma
! Atomic anharmonicity
REAL(KIND=8) :: alpha
! Drive strength
REAL(KIND=8) :: Omega
! Drive detuning from two-photon resonance
REAL(KIND=8) :: delta
! Dipole moment ratio
REAL(KIND=8) :: xi

! Quantum object stuff
! Density Operators rho and steady state density operator
COMPLEX(KIND=8), DIMENSION(9) :: rho, rho_ss, rho_corr
! Runge-kutta vectors
COMPLEX(KIND=8), DIMENSION(9) :: k1, k2, k3, k4
! Atomic Density Operator integers (ij => |i><j|)
INTEGER, PARAMETER :: gg = 1
INTEGER, PARAMETER :: ge = 2
INTEGER, PARAMETER :: gf = 3
INTEGER, PARAMETER :: eg = 4
INTEGER, PARAMETER :: ee = 5
INTEGER, PARAMETER :: ef = 6
INTEGER, PARAMETER :: fg = 7
INTEGER, PARAMETER :: fe = 8
INTEGER, PARAMETER :: ff = 9

! Time stuff
! Time step
REAL(KIND=8) :: dt
! Maximum time to integrate for
REAL(KIND=8) :: t_max
! Maximum tau time to calculate correlation function for
REAL(KIND=8) :: tau1_max, tau2_max
! Maximum number of steps to integrate for
INTEGER :: t_steps, tau_steps
! Runtime variables
REAL :: start_time, end_time, loop_start_time, loop_check_time
REAL :: loop_run_time, loop_remaining_time
! Ten percent of time steps for loop
INTEGER :: ten_percent
! Percentage of time steps completed
INTEGER :: percentage

! Useful stuff
! General integer counters
INTEGER :: t
! Complex i = SQRT(-1)
COMPLEX(KIND=8), PARAMETER :: i = CMPLX(0.d0, 1.0d0, 8)
! 1.0 / 6.0
REAL(KIND=8), PARAMETER :: xis = 1.0d0 / 6.0d0
! List of xi multiplications for easy calculations
REAL(KIND=8), DIMENSION(3) :: xi_list

! Data stuff
! State population
REAL(KIND=8) :: popg, pope, popf, trace_ss
! Bloch equation stuff
COMPLEX(KIND=8) :: sm, sp, sz
! First(second)-order correlation
COMPLEX(KIND=8) :: corr

! Filename stuff
! I/O Format for percentage
CHARACTER(LEN=28), PARAMETER :: FMT_ss = "(T2,I3,A28,F9.2,A19,F9.2,A1)"
CHARACTER(LEN=28), PARAMETER :: FMT_corr = "(T2,I3,A30,F9.2,A19,F9.2,A1)"
! Filename of parameters
CHARACTER(LEN=32), PARAMETER :: filename_parameters = "./data_files/atom_parameters.txt"
! Filename for state population
CHARACTER(LEN=28), PARAMETER :: filename_state = "./data_files/atom_states.txt"
! Filename for state population
! CHARACTER(LEN=27), PARAMETER :: filename_bloch = "./data_files/atom_bloch.txt"
! Filename for first-order correlation
CHARACTER(LEN=29), PARAMETER :: filename_g1 = "./data_files/atom_g1_corr.txt"
! Filename for second-order correlation
CHARACTER(LEN=29), PARAMETER :: filename_g2 = "./data_files/atom_g2_corr.txt"
! Paramert Name List
CHARACTER(LEN=15), PARAMETER :: filename_ParamList = "./ParamList.nml"

!----------------------------------------------------------------------------!
!                   Calling parameters from ParamList_me.nml                 !
!----------------------------------------------------------------------------!
! NameList things
! Status and unit integers
INTEGER :: ISTAT, IUNIT
! Line to be read from file
CHARACTER(LEN=512) :: LINE
! Namelist parameters
NAMELIST /ATOM/ Gamma, Omega, alpha, delta, xi
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

!----------------------------------------------------------------------------!
!                    END parameters from ParamList_me.nml                    !
!----------------------------------------------------------------------------!

! Initialising matrices
!Density operator p. Using the initial condition
!p = |n,1><n,1| for n initial photons in chamber
rho = 0.0d0
rho(gg) = 1.0d0

! Initialise Runge-Kutta Vectors
k1 = 0.0d0
k2 = 0.0d0
k3 = 0.0d0
k4 = 0.0d0

! Open file to write time to
OPEN(UNIT=1, FILE=filename_parameters, STATUS='REPLACE', ACTION='WRITE')
! Write parameter
WRITE(1,*) "Parameters are in the following order:"
WRITE(1,"(A10,F25.15)") "Gamma =", Gamma
WRITE(1,"(A10,F25.15)") "Omega =", Omega
WRITE(1,"(A10,F25.15)") "alpha =", alpha
WRITE(1,"(A10,F25.15)") "delta =", delta
WRITE(1,"(A10,F25.15)") "xi =", xi
WRITE(1,"(A10,F25.15)") "dt =", dt
WRITE(1,"(A10,F25.15)") "Max time =", t_max
WRITE(1,"(A10,F25.15)") "Max tau1 =", tau1_max
WRITE(1,"(A10,F25.15)") "Max tau2 =", tau2_max
! Close file
CLOSE(1)

!##############################################################################!
!                                 STATE SOLVER                                 !
!##############################################################################!

! Open file to write time and data to
OPEN(UNIT=2, FILE=filename_state, STATUS='REPLACE', ACTION='WRITE', RECL=4000)
! OPEN(UNIT=6, FILE=filename_bloch, STATUS='REPLACE', ACTION='WRITE', RECL=4000)

! Ten percent of time steps
ten_percent = NINT((1.0 * t_steps / 10.0))

! Call CPU clock time
CALL CPU_TIME(loop_start_time)

PRINT*, "Calculate steady state ..."
PRINT*, " "

DO t = 0, t_steps
  ! Calculate state probabilities
  popg = 0.0d0
  popg = REAL(rho(gg))
  pope = 0.0d0
  pope = REAL(rho(ee))
  popf = 0.0d0
  popf = REAL(rho(ff))

  ! Calculate Bloch equation stuff
  sm = 0.0d0
  sm = rho(eg) + xi * rho(fe)
  sp = 0.0d0
  sp = rho(ge) + xi * rho(ef)
  sz = 0.0d0
  sz = (xi ** 2) * rho(ff) + (1.0d0 - (xi ** 2)) * rho(ee) - rho(gg)

  ! Write data to file
  WRITE(2, *) DBLE(t) * dt, popg, pope, popf
  ! WRITE(6, *) REAL(sm), IMAG(sm), REAL(sp), IMAG(sp), REAL(sz), IMAG(sz)

  ! Calculate k1
  k1 = 0.0d0
  ! Three-level equations of motion
  k1(gg) = dt * i * 0.5d0 * Omega * rho(ge) - &
         & dt * i * 0.5d0 * Omega * rho(eg) + &
         & dt * gamma * rho(ee)
  k1(ge) = dt * i * 0.5d0 * Omega * rho(gg) - &
         & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * rho(ge) + &
         & dt * i * 0.5d0 * xi * Omega * rho(gf) - &
         & dt * i * 0.5d0 * Omega * rho(ee) + &
         & dt * gamma * xi * rho(ef)
  k1(gf) = dt * i * 0.5d0 * xi * Omega * rho(ge) - &
         & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * rho(gf) - &
         & dt * i * 0.5d0 * Omega * rho(ef)

  k1(eg) = dt * -i * 0.5d0 * Omega * rho(gg) + &
         & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * rho(eg) + &
         & dt * i * 0.5d0 * Omega * rho(ee) - &
         & dt * i * 0.5d0 * xi * Omega * rho(fg) + &
         & dt * gamma * xi * rho(fe)
  k1(ee) = dt * -i * 0.5d0 * Omega * rho(ge) + &
         & dt * i * 0.5d0 * Omega * rho(eg) - &
         & dt * gamma * rho(ee) + &
         & dt * i * 0.5d0 * xi * Omega * rho(ef) - &
         & dt * i * 0.5d0 * xi * Omega * rho(fe) + &
         & dt * gamma * (xi ** 2) * rho(ff)
  k1(ef) = dt * -i * 0.5d0 * Omega * rho(gf) + &
         & dt * i * 0.5d0 * xi * Omega * rho(ee) + &
         & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * rho(ef) - &
         & dt * i * 0.5d0 * xi * Omega * rho(ff)

  k1(fg) = dt * -i * 0.5d0 * xi * Omega * rho(eg) + &
         & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * rho(fg) + &
         & dt * i * 0.5d0 * Omega * rho(fe)
  k1(fe) = -dt * i * 0.5d0 * xi * Omega * rho(ee) + &
         & dt * i * 0.5d0 * Omega * rho(fg) - &
         & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * rho(fe) + &
         & dt * i * 0.5d0 * xi * Omega * rho(ff)
  k1(ff) = -dt * i * 0.5d0 * xi * Omega * rho(ef) + &
         & dt * i * 0.5d0 * xi * Omega * rho(fe) - &
         & dt * gamma * (xi ** 2) * rho(ff)

  ! Calculate k2
  k2 = 0.0d0
  ! Three-level equations of motion
  k2(gg) = dt * i * 0.5d0 * Omega * (rho(ge) + 0.5d0 * k1(ge)) - &
         & dt * i * 0.5d0 * Omega * (rho(eg) + 0.5d0 * k1(eg)) + &
         & dt * gamma * (rho(ee) + 0.5d0 * k1(ee))
  k2(ge) = dt * i * 0.5d0 * Omega * (rho(gg) + 0.5d0 * k1(gg)) - &
         & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(ge) + 0.5d0 * k1(ge)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(gf) + 0.5d0 * k1(gf)) - &
         & dt * i * 0.5d0 * Omega * (rho(ee) + 0.5d0 * k1(ee)) + &
         & dt * gamma * xi * (rho(ef) + 0.5d0 * k1(ef))
  k2(gf) = dt * i * 0.5d0 * xi * Omega * (rho(ge) + 0.5d0 * k1(ge)) - &
         & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(gf) + 0.5d0 * k1(gf)) - &
         & dt * i * 0.5d0 * Omega * (rho(ef) + 0.5d0 * k1(ef))
  k2(eg) = dt * -i * 0.5d0 * Omega * (rho(gg) + 0.5d0 * k1(gg)) + &
         & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(eg) + 0.5d0 * k1(eg)) + &
         & dt * i * 0.5d0 * Omega * (rho(ee) + 0.5d0 * k1(ee)) - &
         & dt * i * 0.5d0 * xi * Omega * (rho(fg) + 0.5d0 * k1(fg)) + &
         & dt * gamma * xi * (rho(fe) + 0.5d0 * k1(fe))
  k2(ee) = dt * -i * 0.5d0 * Omega * (rho(ge) + 0.5d0 * k1(ge)) + &
         & dt * i * 0.5d0 * Omega * (rho(eg) + 0.5d0 * k1(eg)) - &
         & dt * gamma * (rho(ee) + 0.5d0 * k1(ee)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(ef) + 0.5d0 * k1(ef)) - &
         & dt * i * 0.5d0 * xi * Omega * (rho(fe) + 0.5d0 * k1(fe)) + &
         & dt * gamma * (xi ** 2) * (rho(ff) + 0.5d0 * k1(ff))
  k2(ef) = dt * -i * 0.5d0 * Omega * (rho(gf) + 0.5d0 * k1(gf)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(ee) + 0.5d0 * k1(ee)) + &
         & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(ef) + 0.5d0 * k1(ef)) - &
         & dt * i * 0.5d0 * xi * Omega * (rho(ff) + 0.5d0 * k1(ff))
  k2(fg) = dt * -i * 0.5d0 * xi * Omega * (rho(eg) + 0.5d0 * k1(eg)) + &
         & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(fg) + 0.5d0 * k1(fg)) + &
         & dt * i * 0.5d0 * Omega * (rho(fe) + 0.5d0 * k1(fe))
  k2(fe) = -dt * i * 0.5d0 * xi * Omega * (rho(ee) + 0.5d0 * k1(ee)) + &
         & dt * i * 0.5d0 * Omega * (rho(fg) + 0.5d0 * k1(fg)) - &
         & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(fe) + 0.5d0 * k1(fe)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(ff) + 0.5d0 * k1(ff))
  k2(ff) = -dt * i * 0.5d0 * xi * Omega * (rho(ef) + 0.5d0 * k1(ef)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(fe) + 0.5d0 * k1(fe)) - &
         & dt * gamma * (xi ** 2) * (rho(ff) + 0.5d0 * k1(ff))

  ! Calculate k3
  k3 = 0.0d0
  ! Three-level equations of motion
  k3(gg) = dt * i * 0.5d0 * Omega * (rho(ge) + 0.5d0 * k2(ge)) - &
         & dt * i * 0.5d0 * Omega * (rho(eg) + 0.5d0 * k2(eg)) + &
         & dt * gamma * (rho(ee) + 0.5d0 * k2(ee))
  k3(ge) = dt * i * 0.5d0 * Omega * (rho(gg) + 0.5d0 * k2(gg)) - &
         & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(ge) + 0.5d0 * k2(ge)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(gf) + 0.5d0 * k2(gf)) - &
         & dt * i * 0.5d0 * Omega * (rho(ee) + 0.5d0 * k2(ee)) + &
         & dt * gamma * xi * (rho(ef) + 0.5d0 * k2(ef))
  k3(gf) = dt * i * 0.5d0 * xi * Omega * (rho(ge) + 0.5d0 * k2(ge)) - &
         & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(gf) + 0.5d0 * k2(gf)) - &
         & dt * i * 0.5d0 * Omega * (rho(ef) + 0.5d0 * k2(ef))
  k3(eg) = dt * -i * 0.5d0 * Omega * (rho(gg) + 0.5d0 * k2(gg)) + &
         & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(eg) + 0.5d0 * k2(eg)) + &
         & dt * i * 0.5d0 * Omega * (rho(ee) + 0.5d0 * k2(ee)) - &
         & dt * i * 0.5d0 * xi * Omega * (rho(fg) + 0.5d0 * k2(fg)) + &
         & dt * gamma * xi * (rho(fe) + 0.5d0 * k2(fe))
  k3(ee) = dt * -i * 0.5d0 * Omega * (rho(ge) + 0.5d0 * k2(ge)) + &
         & dt * i * 0.5d0 * Omega * (rho(eg) + 0.5d0 * k2(eg)) - &
         & dt * gamma * (rho(ee) + 0.5d0 * k2(ee)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(ef) + 0.5d0 * k2(ef)) - &
         & dt * i * 0.5d0 * xi * Omega * (rho(fe) + 0.5d0 * k2(fe)) + &
         & dt * gamma * (xi ** 2) * (rho(ff) + 0.5d0 * k2(ff))
  k3(ef) = dt * -i * 0.5d0 * Omega * (rho(gf) + 0.5d0 * k2(gf)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(ee) + 0.5d0 * k2(ee)) + &
         & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(ef) + 0.5d0 * k2(ef)) - &
         & dt * i * 0.5d0 * xi * Omega * (rho(ff) + 0.5d0 * k2(ff))
  k3(fg) = dt * -i * 0.5d0 * xi * Omega * (rho(eg) + 0.5d0 * k2(eg)) + &
         & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(fg) + 0.5d0 * k2(fg)) + &
         & dt * i * 0.5d0 * Omega * (rho(fe) + 0.5d0 * k2(fe))
  k3(fe) = -dt * i * 0.5d0 * xi * Omega * (rho(ee) + 0.5d0 * k2(ee)) + &
         & dt * i * 0.5d0 * Omega * (rho(fg) + 0.5d0 * k2(fg)) - &
         & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(fe) + 0.5d0 * k2(fe)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(ff) + 0.5d0 * k2(ff))
  k3(ff) = -dt * i * 0.5d0 * xi * Omega * (rho(ef) + 0.5d0 * k2(ef)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(fe) + 0.5d0 * k2(fe)) - &
         & dt * gamma * (xi ** 2) * (rho(ff) + 0.5d0 * k2(ff))

  ! Calculate k4
  k4 = 0.0d0
  ! Three-level equations of motion
  k4(gg) = dt * i * 0.5d0 * Omega * (rho(ge) + k3(ge)) - &
         & dt * i * 0.5d0 * Omega * (rho(eg) + k3(eg)) + &
         & dt * gamma * (rho(ee) + k3(ee))
  k4(ge) = dt * i * 0.5d0 * Omega * (rho(gg) + k3(gg)) - &
         & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(ge) + k3(ge)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(gf) + k3(gf)) - &
         & dt * i * 0.5d0 * Omega * (rho(ee) + k3(ee)) + &
         & dt * gamma * xi * (rho(ef) + k3(ef))
  k4(gf) = dt * i * 0.5d0 * xi * Omega * (rho(ge) + k3(ge)) - &
         & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(gf) + k3(gf)) - &
         & dt * i * 0.5d0 * Omega * (rho(ef) + k3(ef))
  k4(eg) = dt * -i * 0.5d0 * Omega * (rho(gg) + k3(gg)) + &
         & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(eg) + k3(eg)) + &
         & dt * i * 0.5d0 * Omega * (rho(ee) + k3(ee)) - &
         & dt * i * 0.5d0 * xi * Omega * (rho(fg) + k3(fg)) + &
         & dt * gamma * xi * (rho(fe) + k3(fe))
  k4(ee) = dt * -i * 0.5d0 * Omega * (rho(ge) + k3(ge)) + &
         & dt * i * 0.5d0 * Omega * (rho(eg) + k3(eg)) - &
         & dt * gamma * (rho(ee) + k3(ee)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(ef) + k3(ef)) - &
         & dt * i * 0.5d0 * xi * Omega * (rho(fe) + k3(fe)) + &
         & dt * gamma * (xi ** 2) * (rho(ff) + k3(ff))
  k4(ef) = dt * -i * 0.5d0 * Omega * (rho(gf) + k3(gf)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(ee) + k3(ee)) + &
         & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(ef) + k3(ef)) - &
         & dt * i * 0.5d0 * xi * Omega * (rho(ff) + k3(ff))
  k4(fg) = dt * -i * 0.5d0 * xi * Omega * (rho(eg) + k3(eg)) + &
         & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(fg) + k3(fg)) + &
         & dt * i * 0.5d0 * Omega * (rho(fe) + k3(fe))
  k4(fe) = -dt * i * 0.5d0 * xi * Omega * (rho(ee) + k3(ee)) + &
         & dt * i * 0.5d0 * Omega * (rho(fg) + k3(fg)) - &
         & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(fe) + k3(fe)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(ff) + k3(ff))
  k4(ff) = -dt * i * 0.5d0 * xi * Omega * (rho(ef) + k3(ef)) + &
         & dt * i * 0.5d0 * xi * Omega * (rho(fe) + k3(fe)) - &
         & dt * gamma * (xi ** 2) * (rho(ff) + k3(ff))

  ! Update
  rho = rho + xis * (k1 + 2.0d0 * (k2 + k3) + k4)

  ! Check percentage
  IF (ten_percent /= 0) THEN
    IF (MOD(t, ten_percent) == 0 .AND. t /= 0) THEN
      CALL CPU_TIME(loop_check_time)
      percentage = NINT((100.0 * t) / (1.0 * t_steps))
      loop_run_time = loop_check_time - loop_start_time
      loop_remaining_time = ((100.0 * (loop_check_time - loop_start_time)) / (1.0 * percentage)) - loop_run_time
      WRITE(*, FMT_ss) percentage, "%. Run time (steady state): ", loop_run_time, "s. Est. time left: ", loop_remaining_time, "s"
    END IF
  END IF

  ! Close time loop
END DO

CLOSE(2)
CLOSE(6)

PRINT*, "Steady state found"

! Set steady state rho
rho_ss = rho
trace_ss = pope + (xi ** 2) * popf
PRINT*, " "
PRINT*, "Steady state found!"
PRINT*, "< \Sigma^{\dagger} \Sigma >_{ss} =", trace_ss

!##############################################################################!
!                        FIRST-ORDER CORRELATION FUNCTION                      !
!##############################################################################!
IF (tau1_max /= 0) THEN
  ! Propagate \rho(ss) \Sigma^{\dagger} (lower atom states on <bra|)
  rho = 0.0d0
  ! \rho_{ss} |e><g|
  rho(gg) = rho_ss(ge)
  rho(eg) = rho_ss(ee)
  rho(fg) = rho_ss(fe)
  ! \rho_{ss} \xi |f><e|
  rho(ge) = xi * rho_ss(gf)
  rho(ee) = xi * rho_ss(ef)
  rho(fe) = xi * rho_ss(ff)

  ! Re-initialising the 4-th order vectors
  k1 = 0.0d0
  k2 = 0.0d0
  k3 = 0.0d0
  k4 = 0.0d0
  ! Initialise rho_corr
  rho_corr = 0.0d0

  ! Change tau_steps for g1
  tau_steps = NINT(tau1_max / dt)
  ! Ten percent of time steps
  ten_percent = NINT((1.0 * tau_steps / 10.0))

  ! Open file to write data to
  OPEN(UNIT=3, FILE=filename_g1, STATUS='REPLACE', ACTION='WRITE')

  PRINT*, " "
  PRINT*, "Calculating first-order correlation ..."
  PRINT*, " "

  CALL CPU_TIME(loop_start_time)

  DO t = 0, tau_steps
    ! \Sigma \rho(ss) (lower atom states on |ket>)
    rho_corr = 0.0d0
    ! |g><e| \rho
    rho_corr(gg) = rho(eg)
    rho_corr(ge) = rho(ee)
    rho_corr(gf) = rho(ef)
    ! \xi |e><f| \rho
    rho_corr(eg) = xi * rho(fg)
    rho_corr(ee) = xi * rho(fe)
    rho_corr(ef) = xi * rho(ff)

    ! Calculate the trace
    corr = 0.0d0
    corr = rho_corr(gg) + rho_corr(ee) + rho_corr(ff)

    IF (trace_ss /= 0.0) THEN
      corr = corr / trace_ss
    END IF

    ! Write the results to file
    WRITE(3,*) DBLE(t) * dt, REAL(corr), IMAG(corr)

    ! Calculate k1
    k1 = 0.0d0
    ! Three-level equations of motion
    k1(gg) = dt * i * 0.5d0 * Omega * rho(ge) - &
           & dt * i * 0.5d0 * Omega * rho(eg) + &
           & dt * gamma * rho(ee)
    k1(ge) = dt * i * 0.5d0 * Omega * rho(gg) - &
           & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * rho(ge) + &
           & dt * i * 0.5d0 * xi * Omega * rho(gf) - &
           & dt * i * 0.5d0 * Omega * rho(ee) + &
           & dt * gamma * xi * rho(ef)
    k1(gf) = dt * i * 0.5d0 * xi * Omega * rho(ge) - &
           & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * rho(gf) - &
           & dt * i * 0.5d0 * Omega * rho(ef)
    k1(eg) = dt * -i * 0.5d0 * Omega * rho(gg) + &
           & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * rho(eg) + &
           & dt * i * 0.5d0 * Omega * rho(ee) - &
           & dt * i * 0.5d0 * xi * Omega * rho(fg) + &
           & dt * gamma * xi * rho(fe)
    k1(ee) = dt * -i * 0.5d0 * Omega * rho(ge) + &
           & dt * i * 0.5d0 * Omega * rho(eg) - &
           & dt * gamma * rho(ee) + &
           & dt * i * 0.5d0 * xi * Omega * rho(ef) - &
           & dt * i * 0.5d0 * xi * Omega * rho(fe) + &
           & dt * gamma * (xi ** 2) * rho(ff)
    k1(ef) = dt * -i * 0.5d0 * Omega * rho(gf) + &
           & dt * i * 0.5d0 * xi * Omega * rho(ee) + &
           & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * rho(ef) - &
           & dt * i * 0.5d0 * xi * Omega * rho(ff)
    k1(fg) = dt * -i * 0.5d0 * xi * Omega * rho(eg) + &
           & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * rho(fg) + &
           & dt * i * 0.5d0 * Omega * rho(fe)
    k1(fe) = -dt * i * 0.5d0 * xi * Omega * rho(ee) + &
           & dt * i * 0.5d0 * Omega * rho(fg) - &
           & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * rho(fe) + &
           & dt * i * 0.5d0 * xi * Omega * rho(ff)
    k1(ff) = -dt * i * 0.5d0 * xi * Omega * rho(ef) + &
           & dt * i * 0.5d0 * xi * Omega * rho(fe) - &
           & dt * gamma * (xi ** 2) * rho(ff)

    ! Calculate k2
    k2 = 0.0d0
    ! Three-level equations of motion
    k2(gg) = dt * i * 0.5d0 * Omega * (rho(ge) + 0.5d0 * k1(ge)) - &
           & dt * i * 0.5d0 * Omega * (rho(eg) + 0.5d0 * k1(eg)) + &
           & dt * gamma * (rho(ee) + 0.5d0 * k1(ee))
    k2(ge) = dt * i * 0.5d0 * Omega * (rho(gg) + 0.5d0 * k1(gg)) - &
           & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(ge) + 0.5d0 * k1(ge)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(gf) + 0.5d0 * k1(gf)) - &
           & dt * i * 0.5d0 * Omega * (rho(ee) + 0.5d0 * k1(ee)) + &
           & dt * gamma * xi * (rho(ef) + 0.5d0 * k1(ef))
    k2(gf) = dt * i * 0.5d0 * xi * Omega * (rho(ge) + 0.5d0 * k1(ge)) - &
           & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(gf) + 0.5d0 * k1(gf)) - &
           & dt * i * 0.5d0 * Omega * (rho(ef) + 0.5d0 * k1(ef))
    k2(eg) = dt * -i * 0.5d0 * Omega * (rho(gg) + 0.5d0 * k1(gg)) + &
           & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(eg) + 0.5d0 * k1(eg)) + &
           & dt * i * 0.5d0 * Omega * (rho(ee) + 0.5d0 * k1(ee)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(fg) + 0.5d0 * k1(fg)) + &
           & dt * gamma * xi * (rho(fe) + 0.5d0 * k1(fe))
    k2(ee) = dt * -i * 0.5d0 * Omega * (rho(ge) + 0.5d0 * k1(ge)) + &
           & dt * i * 0.5d0 * Omega * (rho(eg) + 0.5d0 * k1(eg)) - &
           & dt * gamma * (rho(ee) + 0.5d0 * k1(ee)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ef) + 0.5d0 * k1(ef)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(fe) + 0.5d0 * k1(fe)) + &
           & dt * gamma * (xi ** 2) * (rho(ff) + 0.5d0 * k1(ff))
    k2(ef) = dt * -i * 0.5d0 * Omega * (rho(gf) + 0.5d0 * k1(gf)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ee) + 0.5d0 * k1(ee)) + &
           & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(ef) + 0.5d0 * k1(ef)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(ff) + 0.5d0 * k1(ff))
    k2(fg) = dt * -i * 0.5d0 * xi * Omega * (rho(eg) + 0.5d0 * k1(eg)) + &
           & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(fg) + 0.5d0 * k1(fg)) + &
           & dt * i * 0.5d0 * Omega * (rho(fe) + 0.5d0 * k1(fe))
    k2(fe) = -dt * i * 0.5d0 * xi * Omega * (rho(ee) + 0.5d0 * k1(ee)) + &
           & dt * i * 0.5d0 * Omega * (rho(fg) + 0.5d0 * k1(fg)) - &
           & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(fe) + 0.5d0 * k1(fe)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ff) + 0.5d0 * k1(ff))
    k2(ff) = -dt * i * 0.5d0 * xi * Omega * (rho(ef) + 0.5d0 * k1(ef)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(fe) + 0.5d0 * k1(fe)) - &
           & dt * gamma * (xi ** 2) * (rho(ff) + 0.5d0 * k1(ff))

    ! Calculate k3
    k3 = 0.0d0
    ! Three-level equations of motion
    k3(gg) = dt * i * 0.5d0 * Omega * (rho(ge) + 0.5d0 * k2(ge)) - &
           & dt * i * 0.5d0 * Omega * (rho(eg) + 0.5d0 * k2(eg)) + &
           & dt * gamma * (rho(ee) + 0.5d0 * k2(ee))
    k3(ge) = dt * i * 0.5d0 * Omega * (rho(gg) + 0.5d0 * k2(gg)) - &
           & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(ge) + 0.5d0 * k2(ge)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(gf) + 0.5d0 * k2(gf)) - &
           & dt * i * 0.5d0 * Omega * (rho(ee) + 0.5d0 * k2(ee)) + &
           & dt * gamma * xi * (rho(ef) + 0.5d0 * k2(ef))
    k3(gf) = dt * i * 0.5d0 * xi * Omega * (rho(ge) + 0.5d0 * k2(ge)) - &
           & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(gf) + 0.5d0 * k2(gf)) - &
           & dt * i * 0.5d0 * Omega * (rho(ef) + 0.5d0 * k2(ef))
    k3(eg) = dt * -i * 0.5d0 * Omega * (rho(gg) + 0.5d0 * k2(gg)) + &
           & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(eg) + 0.5d0 * k2(eg)) + &
           & dt * i * 0.5d0 * Omega * (rho(ee) + 0.5d0 * k2(ee)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(fg) + 0.5d0 * k2(fg)) + &
           & dt * gamma * xi * (rho(fe) + 0.5d0 * k2(fe))
    k3(ee) = dt * -i * 0.5d0 * Omega * (rho(ge) + 0.5d0 * k2(ge)) + &
           & dt * i * 0.5d0 * Omega * (rho(eg) + 0.5d0 * k2(eg)) - &
           & dt * gamma * (rho(ee) + 0.5d0 * k2(ee)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ef) + 0.5d0 * k2(ef)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(fe) + 0.5d0 * k2(fe)) + &
           & dt * gamma * (xi ** 2) * (rho(ff) + 0.5d0 * k2(ff))
    k3(ef) = dt * -i * 0.5d0 * Omega * (rho(gf) + 0.5d0 * k2(gf)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ee) + 0.5d0 * k2(ee)) + &
           & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(ef) + 0.5d0 * k2(ef)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(ff) + 0.5d0 * k2(ff))
    k3(fg) = dt * -i * 0.5d0 * xi * Omega * (rho(eg) + 0.5d0 * k2(eg)) + &
           & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(fg) + 0.5d0 * k2(fg)) + &
           & dt * i * 0.5d0 * Omega * (rho(fe) + 0.5d0 * k2(fe))
    k3(fe) = -dt * i * 0.5d0 * xi * Omega * (rho(ee) + 0.5d0 * k2(ee)) + &
           & dt * i * 0.5d0 * Omega * (rho(fg) + 0.5d0 * k2(fg)) - &
           & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(fe) + 0.5d0 * k2(fe)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ff) + 0.5d0 * k2(ff))
    k3(ff) = -dt * i * 0.5d0 * xi * Omega * (rho(ef) + 0.5d0 * k2(ef)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(fe) + 0.5d0 * k2(fe)) - &
           & dt * gamma * (xi ** 2) * (rho(ff) + 0.5d0 * k2(ff))

    ! Calculate k4
    k4 = 0.0d0
    ! Three-level equations of motion
    k4(gg) = dt * i * 0.5d0 * Omega * (rho(ge) + k3(ge)) - &
           & dt * i * 0.5d0 * Omega * (rho(eg) + k3(eg)) + &
           & dt * gamma * (rho(ee) + k3(ee))
    k4(ge) = dt * i * 0.5d0 * Omega * (rho(gg) + k3(gg)) - &
           & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(ge) + k3(ge)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(gf) + k3(gf)) - &
           & dt * i * 0.5d0 * Omega * (rho(ee) + k3(ee)) + &
           & dt * gamma * xi * (rho(ef) + k3(ef))
    k4(gf) = dt * i * 0.5d0 * xi * Omega * (rho(ge) + k3(ge)) - &
           & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(gf) + k3(gf)) - &
           & dt * i * 0.5d0 * Omega * (rho(ef) + k3(ef))
    k4(eg) = dt * -i * 0.5d0 * Omega * (rho(gg) + k3(gg)) + &
           & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(eg) + k3(eg)) + &
           & dt * i * 0.5d0 * Omega * (rho(ee) + k3(ee)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(fg) + k3(fg)) + &
           & dt * gamma * xi * (rho(fe) + k3(fe))
    k4(ee) = dt * -i * 0.5d0 * Omega * (rho(ge) + k3(ge)) + &
           & dt * i * 0.5d0 * Omega * (rho(eg) + k3(eg)) - &
           & dt * gamma * (rho(ee) + k3(ee)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ef) + k3(ef)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(fe) + k3(fe)) + &
           & dt * gamma * (xi ** 2) * (rho(ff) + k3(ff))
    k4(ef) = dt * -i * 0.5d0 * Omega * (rho(gf) + k3(gf)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ee) + k3(ee)) + &
           & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(ef) + k3(ef)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(ff) + k3(ff))
    k4(fg) = dt * -i * 0.5d0 * xi * Omega * (rho(eg) + k3(eg)) + &
           & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(fg) + k3(fg)) + &
           & dt * i * 0.5d0 * Omega * (rho(fe) + k3(fe))
    k4(fe) = -dt * i * 0.5d0 * xi * Omega * (rho(ee) + k3(ee)) + &
           & dt * i * 0.5d0 * Omega * (rho(fg) + k3(fg)) - &
           & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(fe) + k3(fe)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ff) + k3(ff))
    k4(ff) = -dt * i * 0.5d0 * xi * Omega * (rho(ef) + k3(ef)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(fe) + k3(fe)) - &
           & dt * gamma * (xi ** 2) * (rho(ff) + k3(ff))

    ! Update
    rho = rho + xis * (k1 + 2.0d0 * (k2 + k3) + k4)

    ! Check percentage
    IF (ten_percent /= 0) THEN
      IF (MOD(t, ten_percent) == 0 .AND. t /= 0) THEN
        CALL CPU_TIME(loop_check_time)
        percentage = NINT((100.0 * t) / (1.0 * tau_steps))
        loop_run_time = loop_check_time - loop_start_time
        loop_remaining_time = ((100.0 * (loop_check_time - loop_start_time)) / (1.0 * percentage)) - loop_run_time
        WRITE(*, FMT_corr) percentage, "%. Run time (g1 correlation): ", loop_run_time, "s. Est. time left: ", loop_remaining_time, "s"
      END IF
    END IF

    ! Close time loop
  END DO

  ! Close file
  CLOSE(3)
END IF

!##############################################################################!
!                       SECOND-ORDER CORRELATION FUNCTION                      !
!##############################################################################!
IF (tau2_max /= 0) THEN
  ! Propagate \Sigma \rho_{ss} \Sigma^{\dagger} (lower atom on both sides)
  rho = 0.0d0
  rho(gg) = rho_ss(ee)
  rho(ge) = xi * rho_ss(ef)
  rho(eg) = xi * rho_ss(fe)
  rho(ee) = (xi ** 2) * rho_ss(ff)

  PRINT*, " "
  PRINT*, "Calculating second-order correlation ..."
  PRINT*, " "

  ! Change tau_steps for g2
  tau_steps = NINT(tau2_max / dt)
  ! Ten percent of time steps
  ten_percent = NINT((1.0 * tau_steps / 10.0))

  ! Open file to write time and correlation to
  OPEN(UNIT=4, FILE=filename_g2, STATUS='REPLACE', ACTION='WRITE', RECL=4000)

  ! Call CPU clock time
  CALL CPU_TIME(loop_start_time)

  DO t = 0, tau_steps
    ! \Sigma^{\dagger} \Sigma \rho
    rho_corr = 0.0d0
    ! |e><e| + \xi^{2} |f><f|
    rho_corr(gg) = rho(ee)
    rho_corr(ge) = xi * rho(ef)
    rho_corr(eg) = xi * rho(fe)
    rho_corr(ee) = (xi ** 2) * rho(ff)

    ! Calculate correlation
    corr = 0.0d0
    corr = REAL(rho_corr(gg) + rho_corr(ee) + rho_corr(ff))
    ! Normalise by steady-state
    IF (trace_ss /= 0.0) THEN
      corr = REAL(corr) / (trace_ss ** 2)
    END IF

    ! Write the results to file
    WRITE(4,*) DBLE(t) * dt, REAL(corr)

    ! Calculate k1
    k1 = 0.0d0
    ! Three-level equations of motion
    k1(gg) = dt * i * 0.5d0 * Omega * rho(ge) - &
           & dt * i * 0.5d0 * Omega * rho(eg) + &
           & dt * gamma * rho(ee)
    k1(ge) = dt * i * 0.5d0 * Omega * rho(gg) - &
           & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * rho(ge) + &
           & dt * i * 0.5d0 * xi * Omega * rho(gf) - &
           & dt * i * 0.5d0 * Omega * rho(ee) + &
           & dt * gamma * xi * rho(ef)
    k1(gf) = dt * i * 0.5d0 * xi * Omega * rho(ge) - &
           & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * rho(gf) - &
           & dt * i * 0.5d0 * Omega * rho(ef)
    k1(eg) = dt * -i * 0.5d0 * Omega * rho(gg) + &
           & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * rho(eg) + &
           & dt * i * 0.5d0 * Omega * rho(ee) - &
           & dt * i * 0.5d0 * xi * Omega * rho(fg) + &
           & dt * gamma * xi * rho(fe)
    k1(ee) = dt * -i * 0.5d0 * Omega * rho(ge) + &
           & dt * i * 0.5d0 * Omega * rho(eg) - &
           & dt * gamma * rho(ee) + &
           & dt * i * 0.5d0 * xi * Omega * rho(ef) - &
           & dt * i * 0.5d0 * xi * Omega * rho(fe) + &
           & dt * gamma * (xi ** 2) * rho(ff)
    k1(ef) = dt * -i * 0.5d0 * Omega * rho(gf) + &
           & dt * i * 0.5d0 * xi * Omega * rho(ee) + &
           & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * rho(ef) - &
           & dt * i * 0.5d0 * xi * Omega * rho(ff)
    k1(fg) = dt * -i * 0.5d0 * xi * Omega * rho(eg) + &
           & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * rho(fg) + &
           & dt * i * 0.5d0 * Omega * rho(fe)
    k1(fe) = -dt * i * 0.5d0 * xi * Omega * rho(ee) + &
           & dt * i * 0.5d0 * Omega * rho(fg) - &
           & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * rho(fe) + &
           & dt * i * 0.5d0 * xi * Omega * rho(ff)
    k1(ff) = -dt * i * 0.5d0 * xi * Omega * rho(ef) + &
           & dt * i * 0.5d0 * xi * Omega * rho(fe) - &
           & dt * gamma * (xi ** 2) * rho(ff)

    ! Calculate k2
    k2 = 0.0d0
    ! Three-level equations of motion
    k2(gg) = dt * i * 0.5d0 * Omega * (rho(ge) + 0.5d0 * k1(ge)) - &
           & dt * i * 0.5d0 * Omega * (rho(eg) + 0.5d0 * k1(eg)) + &
           & dt * gamma * (rho(ee) + 0.5d0 * k1(ee))
    k2(ge) = dt * i * 0.5d0 * Omega * (rho(gg) + 0.5d0 * k1(gg)) - &
           & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(ge) + 0.5d0 * k1(ge)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(gf) + 0.5d0 * k1(gf)) - &
           & dt * i * 0.5d0 * Omega * (rho(ee) + 0.5d0 * k1(ee)) + &
           & dt * gamma * xi * (rho(ef) + 0.5d0 * k1(ef))
    k2(gf) = dt * i * 0.5d0 * xi * Omega * (rho(ge) + 0.5d0 * k1(ge)) - &
           & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(gf) + 0.5d0 * k1(gf)) - &
           & dt * i * 0.5d0 * Omega * (rho(ef) + 0.5d0 * k1(ef))
    k2(eg) = dt * -i * 0.5d0 * Omega * (rho(gg) + 0.5d0 * k1(gg)) + &
           & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(eg) + 0.5d0 * k1(eg)) + &
           & dt * i * 0.5d0 * Omega * (rho(ee) + 0.5d0 * k1(ee)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(fg) + 0.5d0 * k1(fg)) + &
           & dt * gamma * xi * (rho(fe) + 0.5d0 * k1(fe))
    k2(ee) = dt * -i * 0.5d0 * Omega * (rho(ge) + 0.5d0 * k1(ge)) + &
           & dt * i * 0.5d0 * Omega * (rho(eg) + 0.5d0 * k1(eg)) - &
           & dt * gamma * (rho(ee) + 0.5d0 * k1(ee)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ef) + 0.5d0 * k1(ef)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(fe) + 0.5d0 * k1(fe)) + &
           & dt * gamma * (xi ** 2) * (rho(ff) + 0.5d0 * k1(ff))
    k2(ef) = dt * -i * 0.5d0 * Omega * (rho(gf) + 0.5d0 * k1(gf)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ee) + 0.5d0 * k1(ee)) + &
           & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(ef) + 0.5d0 * k1(ef)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(ff) + 0.5d0 * k1(ff))
    k2(fg) = dt * -i * 0.5d0 * xi * Omega * (rho(eg) + 0.5d0 * k1(eg)) + &
           & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(fg) + 0.5d0 * k1(fg)) + &
           & dt * i * 0.5d0 * Omega * (rho(fe) + 0.5d0 * k1(fe))
    k2(fe) = -dt * i * 0.5d0 * xi * Omega * (rho(ee) + 0.5d0 * k1(ee)) + &
           & dt * i * 0.5d0 * Omega * (rho(fg) + 0.5d0 * k1(fg)) - &
           & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(fe) + 0.5d0 * k1(fe)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ff) + 0.5d0 * k1(ff))
    k2(ff) = -dt * i * 0.5d0 * xi * Omega * (rho(ef) + 0.5d0 * k1(ef)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(fe) + 0.5d0 * k1(fe)) - &
           & dt * gamma * (xi ** 2) * (rho(ff) + 0.5d0 * k1(ff))

    ! Calculate k3
    k3 = 0.0d0
    ! Three-level equations of motion
    k3(gg) = dt * i * 0.5d0 * Omega * (rho(ge) + 0.5d0 * k2(ge)) - &
           & dt * i * 0.5d0 * Omega * (rho(eg) + 0.5d0 * k2(eg)) + &
           & dt * gamma * (rho(ee) + 0.5d0 * k2(ee))
    k3(ge) = dt * i * 0.5d0 * Omega * (rho(gg) + 0.5d0 * k2(gg)) - &
           & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(ge) + 0.5d0 * k2(ge)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(gf) + 0.5d0 * k2(gf)) - &
           & dt * i * 0.5d0 * Omega * (rho(ee) + 0.5d0 * k2(ee)) + &
           & dt * gamma * xi * (rho(ef) + 0.5d0 * k2(ef))
    k3(gf) = dt * i * 0.5d0 * xi * Omega * (rho(ge) + 0.5d0 * k2(ge)) - &
           & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(gf) + 0.5d0 * k2(gf)) - &
           & dt * i * 0.5d0 * Omega * (rho(ef) + 0.5d0 * k2(ef))
    k3(eg) = dt * -i * 0.5d0 * Omega * (rho(gg) + 0.5d0 * k2(gg)) + &
           & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(eg) + 0.5d0 * k2(eg)) + &
           & dt * i * 0.5d0 * Omega * (rho(ee) + 0.5d0 * k2(ee)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(fg) + 0.5d0 * k2(fg)) + &
           & dt * gamma * xi * (rho(fe) + 0.5d0 * k2(fe))
    k3(ee) = dt * -i * 0.5d0 * Omega * (rho(ge) + 0.5d0 * k2(ge)) + &
           & dt * i * 0.5d0 * Omega * (rho(eg) + 0.5d0 * k2(eg)) - &
           & dt * gamma * (rho(ee) + 0.5d0 * k2(ee)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ef) + 0.5d0 * k2(ef)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(fe) + 0.5d0 * k2(fe)) + &
           & dt * gamma * (xi ** 2) * (rho(ff) + 0.5d0 * k2(ff))
    k3(ef) = dt * -i * 0.5d0 * Omega * (rho(gf) + 0.5d0 * k2(gf)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ee) + 0.5d0 * k2(ee)) + &
           & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(ef) + 0.5d0 * k2(ef)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(ff) + 0.5d0 * k2(ff))
    k3(fg) = dt * -i * 0.5d0 * xi * Omega * (rho(eg) + 0.5d0 * k2(eg)) + &
           & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(fg) + 0.5d0 * k2(fg)) + &
           & dt * i * 0.5d0 * Omega * (rho(fe) + 0.5d0 * k2(fe))
    k3(fe) = -dt * i * 0.5d0 * xi * Omega * (rho(ee) + 0.5d0 * k2(ee)) + &
           & dt * i * 0.5d0 * Omega * (rho(fg) + 0.5d0 * k2(fg)) - &
           & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(fe) + 0.5d0 * k2(fe)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ff) + 0.5d0 * k2(ff))
    k3(ff) = -dt * i * 0.5d0 * xi * Omega * (rho(ef) + 0.5d0 * k2(ef)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(fe) + 0.5d0 * k2(fe)) - &
           & dt * gamma * (xi ** 2) * (rho(ff) + 0.5d0 * k2(ff))

    ! Calculate k4
    k4 = 0.0d0
    ! Three-level equations of motion
    k4(gg) = dt * i * 0.5d0 * Omega * (rho(ge) + k3(ge)) - &
           & dt * i * 0.5d0 * Omega * (rho(eg) + k3(eg)) + &
           & dt * gamma * (rho(ee) + k3(ee))
    k4(ge) = dt * i * 0.5d0 * Omega * (rho(gg) + k3(gg)) - &
           & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(ge) + k3(ge)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(gf) + k3(gf)) - &
           & dt * i * 0.5d0 * Omega * (rho(ee) + k3(ee)) + &
           & dt * gamma * xi * (rho(ef) + k3(ef))
    k4(gf) = dt * i * 0.5d0 * xi * Omega * (rho(ge) + k3(ge)) - &
           & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(gf) + k3(gf)) - &
           & dt * i * 0.5d0 * Omega * (rho(ef) + k3(ef))
    k4(eg) = dt * -i * 0.5d0 * Omega * (rho(gg) + k3(gg)) + &
           & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(eg) + k3(eg)) + &
           & dt * i * 0.5d0 * Omega * (rho(ee) + k3(ee)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(fg) + k3(fg)) + &
           & dt * gamma * xi * (rho(fe) + k3(fe))
    k4(ee) = dt * -i * 0.5d0 * Omega * (rho(ge) + k3(ge)) + &
           & dt * i * 0.5d0 * Omega * (rho(eg) + k3(eg)) - &
           & dt * gamma * (rho(ee) + k3(ee)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ef) + k3(ef)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(fe) + k3(fe)) + &
           & dt * gamma * (xi ** 2) * (rho(ff) + k3(ff))
    k4(ef) = dt * -i * 0.5d0 * Omega * (rho(gf) + k3(gf)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ee) + k3(ee)) + &
           & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(ef) + k3(ef)) - &
           & dt * i * 0.5d0 * xi * Omega * (rho(ff) + k3(ff))
    k4(fg) = dt * -i * 0.5d0 * xi * Omega * (rho(eg) + k3(eg)) + &
           & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(fg) + k3(fg)) + &
           & dt * i * 0.5d0 * Omega * (rho(fe) + k3(fe))
    k4(fe) = -dt * i * 0.5d0 * xi * Omega * (rho(ee) + k3(ee)) + &
           & dt * i * 0.5d0 * Omega * (rho(fg) + k3(fg)) - &
           & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(fe) + k3(fe)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(ff) + k3(ff))
    k4(ff) = -dt * i * 0.5d0 * xi * Omega * (rho(ef) + k3(ef)) + &
           & dt * i * 0.5d0 * xi * Omega * (rho(fe) + k3(fe)) - &
           & dt * gamma * (xi ** 2) * (rho(ff) + k3(ff))

    ! Update
    rho = rho + xis * (k1 + 2.0d0 * (k2 + k3) + k4)

    ! Check percentage
    IF (ten_percent /= 0) THEN
      IF (MOD(t, ten_percent) == 0 .AND. t /= 0) THEN
        CALL CPU_TIME(loop_check_time)
        percentage = NINT((100.0 * t) / (1.0 * tau_steps))
        loop_run_time = loop_check_time - loop_start_time
        loop_remaining_time = ((100.0 * (loop_check_time - loop_start_time)) / (1.0 * percentage)) - loop_run_time
        WRITE(*, FMT_corr) percentage, "%. Run time (g2 correlation): ", loop_run_time, "s. Est. time left: ", loop_remaining_time, "s"
      END IF
    END IF

    ! Close time loop
  END DO
END IF



! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

END PROGRAM Master_Equation
