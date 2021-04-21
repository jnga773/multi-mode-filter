! This program solves the Lindblad master equation for the state of an atom
! coupled into a multi-mode filter cavity.

! The coupling is achieved via a cascaded systems approach.

! The Hamiltonian is in the form H = H_a + H_c + H_I, where
! -  H_a = -(\alpha/2 + \delta) |e><e| - 2 \delta |f><f|
! -         + \Omega/2 (\Sigma + \Sigma^{\dagger}),
! -  H_c = \sum_{ja=-N}^{N} \Delta_{ja} a^{\dagger}_{N} a_{N},
! and
! -  H_I = -0.5 i \sqrt{\epsilon \Gamma \kappa / 2N+1} *
! \sum_{ja=-N}^{N} (e^{i \phi_{j}} a^{\dagger}_{ja} \Sigma - e^{-i \phi_{j}} a_{ja} \Sigma^{\dagger})
! with
! - \Sigma = |g><e| + \xi |e><f|,
! -  \omega_{j} = \omega_{0} + (ja \delta\omega),
! -  \phi_{ja} = \ja \pi / N

! The master equation is in the Lindblad form:
! - d/dt \rho = 1/(i\hbar) [H, \rho] +
!    \Gamma(1 - \epsilon) /2 (2 \Sigma \rho \Sigma^{\dagger} - \Sigma^{\dagger} \Sigma \rho -
!                \rho \Sigma^{\dagger} \Sigma),
!    \frac{\kappa}{2} \sum_{j=-N}^{N} (2 a_{j} \rho a_{j}^{\dagger} -
!                             a_{j}^{\dagger} a_{j} \rho - \rho a_{j}^{\dagger} a_{j}) +
!    \sum_{j} 1/2 (2 C_{j} \rho C^{\dagger}_{j} - C^{\dagger}_{j} C_{j} \rho - \rho C^{\dagger}_{j} C_{j}),
! where
! - C_{j} = \sqrt{\epsilon \Gamma / 2N+1} e^{i \phi_{j}} \Sigma + \sqrt{2\kappa} a_{j}.

! The program truncates the Fock basis to one photon, allowing for
! - the vacuum: |0>, and
! - a photon in one mode: |1_{j}>
! - two photons in one modes: |2_{j}>, and
! - one photon in two different modes: |1_{j}>|1_{k}> j =/ k.

! The input parameters (gamma, Omega, alpha, delta, xi, w0, epsilon, phase,
! kappa, dw, N, dt, t_max) are taken from a NameList file (ParamList.nml).

! The program outputs the atomic population of the ground and excited states,
! the mean photon number in the multi-mode cavity (<A^{\dagger} A>), and a
! parameter output file:
! - ./data_files/parameters.txt,
! - ./data_files/states.txt (time, rho_gg, rho_ee, rho_ff, <A^{\dagger} A>),
! - ./data_files/g1_corr.txt (time, REAL(g1), IMAG(g1)),
! - ./data_files/g2_corr.txt (time, g2).

PROGRAM THREE_LEVEL_BANDPASS_MASTER_EQUATION_TWO_PHOTON

IMPLICIT NONE

! Parameters in terms of decay rate gamma
! Atomic ecay Rate
REAL(KIND=8) :: gamma
! Drive strength
REAL(KIND=8) :: Omega
! Atomic anharmonicity
REAL(KIND=8) :: alpha
! Drive detuning from two-photon resonance
REAL(KIND=8) :: delta
! Dipole moment ratio
REAL(KIND=8) :: xi

! Filter parameter stuff
! Percentage of fluorecence aimed at cavity
REAL(KIND=8)                                           :: epsilon
! Number of mode either side of w0, 2N + 1 total mode
INTEGER                                                :: N
! Phase modulation of mode coupling
INTEGER                                                :: phase
! Central mode frequency of the filter cavity, with N mode frequencies either side
REAL(KIND=8)                                           :: w0a
! REAL(KIND=8)                                           :: w0b
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8)                                           :: kappaa
! REAL(KIND=8)                                           :: kappab
! Frequency spacing of modes
REAL(KIND=8)                                           :: dwa
! REAL(KIND=8)                                           :: dwb
! kappa, w0, and dw values
REAL(KIND=8)                                           :: kappa, w0, dw

! Quantum Object stuff
! Fock basis truncation
INTEGER :: Fock = 2
! Density operator
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE :: rho, rho_ss, rho_corr
! Runge-kutta vectors
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE :: k1, k2, k3, k4
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

! Useful things
! Time step integer
INTEGER :: t
! Photon number integers
INTEGER :: nja, njac, nka, nkac
! Filter number integer
INTEGER :: ja, jac, ka, kac
! Index variable for state vector
INTEGER :: pa, pac, pn, pnc
! General purpose integer counter
INTEGER :: z
! Max number of indice for array of state
INTEGER :: no_states, no_states_two_photon, no_states_one_photon
! Mapping matrix (nja, nka, ja, ka) --> i
INTEGER, DIMENSION(:, :, :, :), ALLOCATABLE :: MAP_n2p
! Mapping array i --> ;(nja, nka, ja, ka)
INTEGER, DIMENSION(:, :), ALLOCATABLE :: MAP_p2n
! List of square roots
REAL(KIND=8), DIMENSION(0:20) :: sqrtl
! List of Delta values
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: wl
! List of mode dependent cascade coupling values
COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: gkal
! Blackman window coefficient
REAL(KIND=8) :: blackman
! Complex i = SQRT(-1)
COMPLEX(KIND=8), PARAMETER :: i = CMPLX(0.0d0, 1.0d0, 8)
! 1.0 / 6.0
REAL(KIND=8), PARAMETER :: xis = 1.0d0 / 6.0d0
! pi
REAL(KIND=8), PARAMETER :: pi = 3.14159265358979323846d0

! Data stuff
! Sample rate for state populations
INTEGER :: sample_rate
! State population
REAL(KIND=8) :: rho_gg, rho_ee, rho_ff, trace_ss
! Mean photon number in cavity
REAL(KIND=8) :: photona, photona_ss
! First(second)-order correlation
COMPLEX(KIND=8) :: corr

! Filename stuff
! I/O Format for percentage
CHARACTER(LEN=28), PARAMETER :: FMT_ss = "(T2,I3,A28,F9.2,A19,F9.2,A1)"
CHARACTER(LEN=28), PARAMETER :: FMT_corr = "(T2,I3,A30,F9.2,A19,F9.2,A1)"
! Filename of parameters
CHARACTER(LEN=32), PARAMETER :: filename_parameters = "./data_files/me_parameters.txt"
! Filename for state population
CHARACTER(LEN=28), PARAMETER :: filename_state = "./data_files/me_states.txt"
! Filename for first-order correlation
CHARACTER(LEN=29), PARAMETER :: filename_g1 = "./data_files/me_g1_corr.txt"
! Filename for second-order correlation
CHARACTER(LEN=29), PARAMETER :: filename_g2 = "./data_files/me_g2_corr.txt"
! Paramert Name List
CHARACTER(LEN=15), PARAMETER :: filename_ParamList = "./ParamList.nml"

!------------------------------------------------------------------------------!
!                    Calling parameters from ParamList_me.nml                  !
!------------------------------------------------------------------------------!
! NameList things
! Status and unit integers
INTEGER            :: ISTAT, IUNIT
! Line to be read from file
CHARACTER(LEN=512) :: LINE
! Namelist parameters
NAMELIST /ATOM/ Gamma, Omega, alpha, delta, xi
NAMELIST /FILTER/ epsilon, N, phase
NAMELIST /CAVITYA/ kappaa, w0a, dwa
! NAMELIST /CAVITYB/ kappab, w0b, dwb
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

! READ(IUNIT, NML=CAVITYB, IOSTAT=ISTAT)
! IF (ISTAT .NE. 0) THEN
!   BACKSPACE(IUNIT)
!   READ(IUNIT, FMT='(A)') LINE
!   CLOSE(IUNIT)
!   PRINT *, "Invalid line in CAVITYB namelist: " // TRIM(line)
!   CALL EXIT(1)
! END IF

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
! Max number of state
no_states_two_photon = NINT(1.0 + (2.0 * (2.0 * N + 1)) + ((2.0 * (N ** 2) + N)))
no_states_one_photon = NINT(1.0 + (2.0 * N + 1))

! Howard uses \kappa as the halfwidth of a single-mode cavity, for a rate of
! \sqrt{2 \kappa}. If you want to use that definition, un-comment the next line.
! kappa = 2.0d0 * kappa

!------------------------------------------------------------------------------!
!                           Mapping of variable                                !
!------------------------------------------------------------------------------!
no_states = no_states_two_photon
Fock = 2

! Mapping (ja, ka, nja, nka) --> index number i
ALLOCATE(MAP_n2p(-N:N, -N:N, 0:Fock, 0:Fock))
MAP_n2p = 0
DO nka = 0, Fock
  DO nja = 0, Fock
    DO ja = -N, N
      DO ka = -N, N
        ! Vacuum state |0_{ja}>|0_{ka}>
        IF (nja == 0 .AND. nka == 0) THEN
          MAP_n2p(ja, ka, nja, nka) = 1
        END IF

        ! Single-photon state |1_{ja}>|0_{ka}>
        IF (nja == 1 .AND. nka == 0) THEN
          MAP_n2p(ja, ka, nja, nka) = 1 + (ja + (N + 1))
          MAP_n2p(ka, ja, nka, nja) = 1 + (ja + (N + 1))
        END IF

        ! Two-Photon state |2_{ja}>|0_{ka}>
        IF (nja == 2 .AND. nka == 0) THEN
          MAP_n2p(ja, ka, nja, nka) = ((2 * N) + 2) + (ja + (N + 1))
          MAP_n2p(ka, ja, nka, nja) = ((2 * N) + 2) + (ja + (N + 1))
        END IF
        ! if ja == ka and nka = ka = 1
        IF (nja == 1 .AND. nka == 1 .AND. ja == ka) THEN
          MAP_n2p(ja, ka, nja, nka) = ((2 * N) + 2) + (ja + (N + 1))
        END IF

        ! One photon in each mode |1_{ja}>|1_{ka}>
        IF (nja == 1 .AND. nka == 1 .AND. ja /= ka) THEN
          MAP_n2p(ja, ka, nja, nka) = ((4 * N) + 3) + NINT(-0.5 * ((N + ja) ** 2)) + &
                                    & NINT(((2*N) + 0.5) * ABS(N + ja)) + ABS(ja - ka)
        END IF
      END DO
    END DO
  END DO
END DO

! Mapping index number i --> (ja, ka, nja, nka)
ALLOCATE(MAP_p2n(no_states_two_photon, 4))
MAP_p2n = 0

DO pa = 1, no_states_two_photon
  DO nja = 0, Fock
    DO nka = 0, Fock
      DO ja = -N, N
        DO ka = ja, N
          IF (MAP_n2p(ja, ka, nja, nka) == pa) THEN
            ! If nja > 0 and nka = 0 then have ja /= ka
            IF (N == 0) THEN
              MAP_p2n(pa, 1) = ja
              MAP_p2n(pa, 2) = ka
              MAP_p2n(pa, 3) = nja
              MAP_p2n(pa, 4) = nka
            ELSE
              IF (nja == 0 .AND. nka == 0 .AND. ja == ka) THEN
                MAP_p2n(pa, 1) = ja
                MAP_p2n(pa, 2) = ka-1
                MAP_p2n(pa, 3) = nja
                MAP_p2n(pa, 4) = nka
              ELSE IF (nja > 0 .AND. nka == 0 .AND. ja == ka) THEN
                MAP_p2n(pa, 1) = ja
                MAP_p2n(pa, 2) = ka-1
                MAP_p2n(pa, 3) = nja
                MAP_p2n(pa, 4) = nka
              ELSE
                MAP_p2n(pa, 1) = ja
                MAP_p2n(pa, 2) = ka
                MAP_p2n(pa, 3) = nja
                MAP_p2n(pa, 4) = nka
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO
END DO

! Cycle through and flip indices of MAP_n2p
DO nja = 0, Fock
  DO nka = 0, Fock
    DO ja = -N, N
      DO ka = ja, N
        MAP_n2p(ka, ja, nka, nja) = MAP_n2p(ja, ka, nja, nka)
      END DO
    END DO
  END DO
END DO

!------------------------------------------------------------------------------!
!                          Allocate and initialise                             !
!------------------------------------------------------------------------------!
! Allocate list of square root value
sqrtl = 0.0d0
DO ja = 0, 20
  sqrtl(ja) = DSQRT(DBLE(ja))
END DO

! Allocate array of Delta and gka values
ALLOCATE(wl(-N:N))
wl = 0.0d0
ALLOCATE(gkal(-N:N))
gkal = 0.0d0
DO ja = -N, N
  IF (N == 0) THEN
    wl(ja) = w0
    gkal(ja) = DSQRT(epsilon * gamma * kappa)
  ELSE
    wl(ja) = w0 + DBLE(ja) * dw
    ! Blackman window coefficient
    blackman = 1.0d0
    ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + ja) / (2.0d0 * DBLE(N))) + &
    !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + ja) / (2.0d0 * DBLE(N)))
    ! Mode dependent phase difference
    gkal(ja) = DSQRT((epsilon / DBLE(2*N + 1)) * gamma * kappa) * blackman * EXP(i * DBLE(phase) * DBLE(ja) * pi / DBLE(N))
  END IF
END DO

! If not calculating second-order correlation, set to single-photon truncation.
! Otherwise, set to two photon truncation
IF (tau2_max > 0.0) THEN
  no_states = no_states_two_photon
  Fock = 2
ELSE
  no_states = no_states_one_photon
  Fock = 1
END IF

! Set up density operator (n_ja, n_ja', n_ka, n_ka', ja, ja', ka, ka', +/-, +/-')'
ALLOCATE(rho(no_states, no_states, 9))
rho = 0.0d0
! State initially in ground state and zero photons
rho(1, 1, gg) = 1.0d0

! Allocate Runge-kutta vectors
ALLOCATE(k1(no_states, no_states, 9))
k1 = 0.0d0
ALLOCATE(k2(no_states, no_states, 9))
k2 = 0.0d0
ALLOCATE(k3(no_states, no_states, 9))
k3 = 0.0d0
ALLOCATE(k4(no_states, no_states, 9))
k4 = 0.0d0

rho_gg = 0.0d0
rho_ee = 0.0d0
photona = 0.0d0

! Open file to write time to
OPEN(UNIT=1, FILE=filename_parameters, STATUS='REPLACE', ACTION='WRITE')
! Write parameter
WRITE(1,*) "Parameters are in the following order:"
WRITE(1,"(A10,F25.15)") "Gamma =", Gamma
WRITE(1,"(A10,F25.15)") "Omega =", Omega
WRITE(1,"(A10,F25.15)") "alpha =", alpha
WRITE(1,"(A10,F25.15)") "delta =", delta
WRITE(1,"(A10,F25.15)") "xi =", xi
WRITE(1,"(A10,F25.15)") "w0 =", w0
WRITE(1,"(A10,F25.15)") "kappa =", kappa
WRITE(1,"(A11,F25.15)") "dw = ", dw
WRITE(1,"(A10,F25.15)") "epsilon =", epsilon
WRITE(1,"(A11,I9)") "N = ", N
WRITE(1,"(A11,I9)") "phase = ", phase
WRITE(1,"(A10,F25.15)") "dt =", dt
WRITE(1,"(A10,F25.15)") "Max time =", t_max
WRITE(1,"(A11,F25.15)") "Max tau1 = ", tau1_max
! Close file
CLOSE(1)

!##############################################################################!
!                                 STATE SOLVER                                 !
!##############################################################################!

! Open file to write time and data to
OPEN(UNIT=2, FILE=filename_state, STATUS='REPLACE', ACTION='WRITE', RECL=4000)

! Ten percent of time steps
ten_percent = NINT((1.0 * t_steps / 10.0))

! Let max number of data points be 100,000 ~ 10mb file.
sample_rate = NINT(DBLE(t_steps) / 1d5)

! Call CPU clock time
CALL CPU_TIME(loop_start_time)

PRINT*, "Calculate steady state ..."
PRINT*, " "

! Start time integration to find steady state
DO t = 0, t_steps
  ! Calculate state probabilities and mean photon number
  rho_gg = 0.0d0
  rho_ee = 0.0d0
  rho_ff = 0.0d0
  photona = 0.0d0
  DO pac = 1, no_states
    ! Atomic state
    rho_gg = rho_gg + REAL(rho(pac, pac, gg))
    rho_ee = rho_ee + REAL(rho(pac, pac, ee))
    rho_ff = rho_ff + REAL(rho(pac, pac, ff))
    DO pa = 1, no_states
      ja = MAP_p2n(pa, 1)
      jac = MAP_p2n(pac, 1)
      ka = MAP_p2n(pa, 2)
      kac = MAP_p2n(pac, 2)
      nja = MAP_p2n(pa, 3)
      njac = MAP_p2n(pac, 3)
      nka = MAP_p2n(pa, 4)
      nkac = MAP_p2n(pac, 4)
      ! To calculate < A^{\dagger} A > = <p' | A^{\dagger} A | p>, remove a
      ! photon from either side the, if the indice are the same, update the
      ! photona variable.
      IF (nja > 0 .AND. njac > 0) THEN
        pn = MAP_n2p(ja, ka, nja-1, nka)
        pnc = MAP_n2p(jac, kac, njac-1, nkac)
        IF (pn == pnc) THEN
          photona = photona + sqrtl(nja) * sqrtl(njac) * REAL(rho(pa, pac, gg) + rho(pa, pac, ee) + rho(pa, pac, ff))
        END IF
      END IF
      IF (nja > 0 .AND. nkac > 0) THEN
        pn = MAP_n2p(ja, ka, nja-1, nka)
        pnc = MAP_n2p(jac, kac, njac, nkac-1)
        IF (pn == pnc) THEN
          photona = photona + sqrtl(nja) * sqrtl(nkac) * REAL(rho(pa, pac, gg) + rho(pa, pac, ee) + rho(pa, pac, ff))
        END IF
      END IF
      IF (nka > 0 .AND. njac > 0) THEN
        pn = MAP_n2p(ja, ka, nja, nka-1)
        pnc = MAP_n2p(jac, kac, njac-1, nkac)
        IF (pn == pnc) THEN
          photona = photona + sqrtl(nka) * sqrtl(njac) * REAL(rho(pa, pac, gg) + rho(pa, pac, ee) + rho(pa, pac, ff))
        END IF
      END IF
      IF (nka > 0 .AND. nkac > 0) THEN
        pn = MAP_n2p(ja, ka, nja, nka-1)
        pnc = MAP_n2p(jac, kac, njac, nkac-1)
        IF (pn == pnc) THEN
          photona = photona + sqrtl(nka) * sqrtl(nkac) * REAL(rho(pa, pac, gg) + rho(pa, pac, ee) + rho(pa, pac, ff))
        END IF
      END IF
    END DO
  END DO

  ! If t_max is really big, only take a sample of results to write to file
  ! so file size isn't huge-mongous.
  IF (t_max > 100.0) THEN
    ! Check if mod
    IF (MOD(t, sample_rate) == 0) THEN
      ! Write data
      WRITE(2,*) DBLE(t) * dt, rho_gg, rho_ee, rho_ff, photona
    END IF
  ELSE
    ! Write data as normal
    WRITE(2,*) DBLE(t) * dt, rho_gg, rho_ee, rho_ff, photona
  END IF

  ! Calculate k1
  k1 = 0.0d0
  DO pac = 1, no_states
    DO pa = 1, no_states
      ! Get mode and photon number
      ja = MAP_p2n(pa, 1)
      jac = MAP_p2n(pac, 1)
      ka = MAP_p2n(pa, 2)
      kac = MAP_p2n(pac, 2)
      nja = MAP_p2n(pa, 3)
      njac = MAP_p2n(pac, 3)
      nka = MAP_p2n(pa, 4)
      nkac = MAP_p2n(pac, 4)

      ! Three-level equations of motion
      k1(pa, pac, gg) = k1(pa, pac, gg) + dt * i * 0.5d0 * Omega * rho(pa, pac, ge) - &
                      & dt * i * 0.5d0 * Omega * rho(pa, pac, eg) + &
                      & dt * gamma * rho(pa, pac, ee)
      k1(pa, pac, ge) = k1(pa, pac, ge) + dt * i * 0.5d0 * Omega * rho(pa, pac, gg) - &
                      & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * rho(pa, pac, ge) + &
                      & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, gf) - &
                      & dt * i * 0.5d0 * Omega * rho(pa, pac, ee) + &
                      & dt * gamma * xi * rho(pa, pac, ef)
      k1(pa, pac, gf) = k1(pa, pac, gf) + dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ge) - &
                      & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * rho(pa, pac, gf) - &
                      & dt * i * 0.5d0 * Omega * rho(pa, pac, ef)
      k1(pa, pac, eg) = k1(pa, pac, eg) - dt * i * 0.5d0 * Omega * rho(pa, pac, gg) + &
                      & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * rho(pa, pac, eg) + &
                      & dt * i * 0.5d0 * Omega * rho(pa, pac, ee) - &
                      & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, fg) + &
                      & dt * gamma * xi * rho(pa, pac, fe)
      k1(pa, pac, ee) = k1(pa, pac, ee) - dt * i * 0.5d0 * Omega * rho(pa, pac, ge) + &
                      & dt * i * 0.5d0 * Omega * rho(pa, pac, eg) - &
                      & dt * gamma * rho(pa, pac, ee) + &
                      & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ef) - &
                      & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, fe) + &
                      & dt * gamma * (xi ** 2) * rho(pa, pac, ff)
      k1(pa, pac, ef) = k1(pa, pac, ef) - dt * i * 0.5d0 * Omega * rho(pa, pac, gf) + &
                      & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ee) + &
                      & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * rho(pa, pac, ef) - &
                      & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ff)
      k1(pa, pac, fg) = k1(pa, pac, fg) - dt * i * 0.5d0 * xi * Omega * rho(pa, pac, eg) + &
                      & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * rho(pa, pac, fg) + &
                      & dt * i * 0.5d0 * Omega * rho(pa, pac, fe)
      k1(pa, pac, fe) = k1(pa, pac, fe) - dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ee) + &
                      & dt * i * 0.5d0 * Omega * rho(pa, pac, fg) - &
                      & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * rho(pa, pac, fe) + &
                      & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ff)
      k1(pa, pac, ff) = k1(pa, pac, ff) - dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ef) + &
                      & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, fe) - &
                      & dt * gamma * (xi ** 2) * rho(pa, pac, ff)

      ! Mode Hamiltonians
      k1(pa, pac, :) = k1(pa, pac, :) - i * dt * (&
                     & (wl(ja) * DBLE(nja)) - (wl(jac) * DBLE(njac)) + &
                     & (wl(ka) * DBLE(nka)) - (wl(kac) * DBLE(nkac))) * &
                     & rho(pa, pac, :)

      !-----------------------!
      !  Commutator -i[H, p]  !
      !-----------------------!
      ! -i * H * rho
      ! Cascade Hamiltonian
      ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{i \phi_{j}} a^{\dagger}_{j} \sigma_{-} \rho
      IF (nja < Fock .AND. nka == 0) THEN
        ! Cycle through modes and create a photon
        DO z = -N, N
          IF (z == ja) THEN
            ! Add a photon to ja
            pn = MAP_n2p(ja, ka, nja+1, nka)
            ! |g><e| \rho
            k1(pn, pac, gg) = k1(pn, pac, gg) - dt * gkal(ja) * sqrtl(nja+1) * &
                            & rho(pa, pac, eg)
            k1(pn, pac, ge) = k1(pn, pac, ge) - dt * gkal(ja) * sqrtl(nja+1) * &
                            & rho(pa, pac, ee)
            k1(pn, pac, gf) = k1(pn, pac, gf) - dt * gkal(ja) * sqrtl(nja+1) * &
                            & rho(pa, pac, ef)
            ! \xi |e><f| \rho
            k1(pn, pac, eg) = k1(pn, pac, eg) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                            & rho(pa, pac, fg)
            k1(pn, pac, ee) = k1(pn, pac, ee) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                            & rho(pa, pac, fe)
            k1(pn, pac, ef) = k1(pn, pac, ef) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                            & rho(pa, pac, ff)
          ELSE
            ! Add a photon to other modes
            pn = MAP_n2p(ja, z, nja, nka+1)
            ! |g><e| \rho
            k1(pn, pac, gg) = k1(pn, pac, gg) - dt * gkal(z) * &
                            & rho(pa, pac, eg)
            k1(pn, pac, ge) = k1(pn, pac, ge) - dt * gkal(z) * &
                            & rho(pa, pac, ee)
            k1(pn, pac, gf) = k1(pn, pac, gf) - dt * gkal(z) * &
                            & rho(pa, pac, ef)
            ! \xi |e><f| \rho
            k1(pn, pac, eg) = k1(pn, pac, eg) - dt * gkal(z) * xi * &
                            & rho(pa, pac, fg)
            k1(pn, pac, ee) = k1(pn, pac, ee) - dt * gkal(z) * xi * &
                            & rho(pa, pac, fe)
            k1(pn, pac, ef) = k1(pn, pac, ef) - dt * gkal(z) * xi * &
                            & rho(pa, pac, ff)
          END IF
        END DO
      END IF

      ! i * rho * H
      ! Cascade Hamiltonian
      ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{-i \phi_{j}} \rho \sigma_{+} a_{j}
      IF (njac < Fock .AND. nkac == 0) THEN
        ! Cycle through modes and create a photon
        DO z = -N, N
          IF (z == jac) THEN
            ! Add a photon to mode ja
            pn = MAP_n2p(jac, kac, njac+1, nkac)
            ! \rho |e><g|
            k1(pa, pn, gg) = k1(pa, pn, gg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                           & rho(pa, pac, ge)
            k1(pa, pn, eg) = k1(pa, pn, eg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                           & rho(pa, pac, ee)
            k1(pa, pn, fg) = k1(pa, pn, fg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                           & rho(pa, pac, fe)
            ! \xi \rho |f><e|
            k1(pa, pn, ge) = k1(pa, pn, ge) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                           & rho(pa, pac, gf)
            k1(pa, pn, ee) = k1(pa, pn, ee) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                           & rho(pa, pac, ef)
            k1(pa, pn, fe) = k1(pa, pn, fe) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                           & rho(pa, pac, ff)
          ELSE
            ! Add a photon to other modes
            pn = MAP_n2p(jac, z, njac, nkac+1)
            ! \rho |e><g|
            k1(pa, pn, gg) = k1(pa, pn, gg) - dt * CONJG(gkal(z)) * &
                           & rho(pa, pac, ge)
            k1(pa, pn, eg) = k1(pa, pn, eg) - dt * CONJG(gkal(z)) * &
                           & rho(pa, pac, ee)
            k1(pa, pn, fg) = k1(pa, pn, fg) - dt * CONJG(gkal(z)) * &
                           & rho(pa, pac, fe)
            ! \xi \rho |f><e|
            k1(pa, pn, ge) = k1(pa, pn, ge) - dt * CONJG(gkal(z)) * xi * &
                           & rho(pa, pac, gf)
            k1(pa, pn, ee) = k1(pa, pn, ee) - dt * CONJG(gkal(z)) * xi * &
                           & rho(pa, pac, ef)
            k1(pa, pn, fe) = k1(pa, pn, fe) - dt * CONJG(gkal(z)) * xi * &
                           & rho(pa, pac, ff)
          END IF
        END DO
      END IF

      !-----------------------------!
      !      LINDBLAD DECAY         !
      !-----------------------------!
      ! The atom fluorecence is split by a beam splitter into two paths:
      ! one goes towards the cavity at a fraction \epsilon, and the other
      ! goes to a separate detector at a fraction 1 - \epsilon.
      ! We split the cascaded systems decay operator:
      ! ja = \sqrt{\gamma} \sigma_{-} + \sqrt{\kappa} A,
      ! into separate components correponding to atomic decay, cavity decay
      ! and two cross terms. We also include a separate decay term to account
      ! for the detector placed outside the ring cavity.

      !--------------------!
      !    Cavity Decay    !
      !--------------------!
      ! Lindblad is (accounting for both cascaded and cavity decay)
      ! 2.0 * 0.5 \kappa (2 A \rho A^{\dagger} - A^{\dagger} A \rho
      !                                        - \rho A^{\dagger} A_),
      ! where A is a sum of all mode annihilation operators.

      ! \kappa * A * \rho * A^{\dagger} term
      ! For |nja>_{ja}|0>_{ka} <nja|_{ja'} <0|_{ka'}
      IF (nja > 0 .AND. njac > 0 .AND. ja == jac) THEN
        pn = MAP_n2p(ja, ka, nja-1, nka)
        pnc = MAP_n2p(jac, kac, njac-1, nkac)
        k1(pn, pnc, :) = k1(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * sqrtl(njac) * &
                       & rho(pa, pac, :)
      END IF
      IF (nja > 0 .AND. nkac > 0 .AND. ja == kac) THEN
        pn = MAP_n2p(ja, ka, nja-1, nka)
        pnc = MAP_n2p(jac, kac, njac, nkac-1)
        k1(pn, pnc, :) = k1(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * &
                       & rho(pa, pac, :)
      END IF
      IF (nka > 0 .AND. njac > 0 .AND. ka == jac) THEN
        pn = MAP_n2p(ja, ka, nja, nka-1)
        pnc = MAP_n2p(jac, kac, njac-1, nkac)
        k1(pn, pnc, :) = k1(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(njac) * &
                       & rho(pa, pac, :)
      END IF
      IF (nka > 0 .AND. nkac > 0 .AND. ka == kac) THEN
        pn = MAP_n2p(ja, ka, nja, nka-1)
        pnc = MAP_n2p(jac, kac, njac, nkac-1)
        k1(pn, pnc, :) = k1(pn, pnc, :) + dt * 2.0d0 * kappa * &
                       & rho(pa, pac, :)
      END IF

      ! -a^{\dagger}_{j} a_{j} \rho - \rho a^{\dagger}_{j} a_{j}
      k1(pa, pac, :) = k1(pa, pac, :) - dt * kappa * (DBLE(nja+nka) + DBLE(njac+nkac)) * &
                     & rho(pa, pac, :)

      !--------------------!
      !   Cascade Decay    !
      !--------------------!
      ! \sqrt(2 epsilon gamma kappa) e^{i \phi_{j}} \sigma_{-} \rho a^{\dagger}_{j}
      ! Lower atom |\pm> and annihilate photon from mode jac/kac
      IF (njac > 0) THEN
        ! |1_{j},0_{k}> --> |0>, |2_{j},0_{k}> -> |1_{j},0_{k}>
        ! |1_{j},1_{k}> --> |0_{j},1_{k}>
        ! Annihilate a photon from mode jac
        pn = MAP_n2p(jac, kac, njac-1, nkac)
        ! |g><e| \rho
        k1(pa, pn, gg) = k1(pa, pn, gg) + dt * gkal(jac) * sqrtl(njac) * &
                       & rho(pa, pac, eg)
        k1(pa, pn, ge) = k1(pa, pn, ge) + dt * gkal(jac) * sqrtl(njac) * &
                       & rho(pa, pac, ee)
        k1(pa, pn, gf) = k1(pa, pn, gf) + dt * gkal(jac) * sqrtl(njac) * &
                       & rho(pa, pac, ef)
        ! \xi |e><f| \rho
        k1(pa, pn, eg) = k1(pa, pn, eg) + dt * gkal(jac) * xi * sqrtl(njac) * &
                       & rho(pa, pac, fg)
        k1(pa, pn, ee) = k1(pa, pn, ee) + dt * gkal(jac) * xi * sqrtl(njac) * &
                       & rho(pa, pac, fe)
        k1(pa, pn, ef) = k1(pa, pn, ef) + dt * gkal(jac) * xi * sqrtl(njac) * &
                       & rho(pa, pac, ff)
        IF (nkac > 0) THEN
          ! |1_{j},1_{k}> --> |1_{j},0_{k}>
          ! Annihilate a photon from mode kac
          pn = MAP_n2p(jac, kac, njac, nkac-1)
          ! |g><e| \rho
          k1(pa, pn, gg) = k1(pa, pn, gg) + dt * gkal(kac) * &
                         & rho(pa, pac, eg)
          k1(pa, pn, ge) = k1(pa, pn, ge) + dt * gkal(kac) * &
                         & rho(pa, pac, ee)
          k1(pa, pn, gf) = k1(pa, pn, gf) + dt * gkal(kac) * &
                         & rho(pa, pac, ef)
          ! \xi |e><f| \rho
          k1(pa, pn, eg) = k1(pa, pn, eg) + dt * gkal(kac) * xi * &
                         & rho(pa, pac, fg)
          k1(pa, pn, ee) = k1(pa, pn, ee) + dt * gkal(kac) * xi * &
                         & rho(pa, pac, fe)
          k1(pa, pn, ef) = k1(pa, pn, ef) + dt * gkal(kac) * xi * &
                         & rho(pa, pac, ff)
        END IF
      END IF

      ! \sqrt(2 epsilon gamma kappa) e^{-i \phi_{j}} a_{j} \rho \sigma_{+}
      ! Lower atom <\pm| and annihilate photon from mode ja/ka
      IF (nja > 0) THEN
        ! <1_{j}|<0_{k}| or <2_{j}|<0_{k}|
        ! Annihilate a photon from mode ja
        pn = MAP_n2p(ja, ka, nja-1, nka)
        ! \rho |e><g|
        k1(pn, pac, gg) = k1(pn, pac, gg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                        & rho(pa, pac, ge)
        k1(pn, pac, eg) = k1(pn, pac, eg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                        & rho(pa, pac, ee)
        k1(pn, pac, fg) = k1(pn, pac, fg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                        & rho(pa, pac, fe)
        ! \xi \rho |f><e|
        k1(pn, pac, ge) = k1(pn, pac, ge) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                        & rho(pa, pac, gf)
        k1(pn, pac, ee) = k1(pn, pac, ee) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                        & rho(pa, pac, ef)
        k1(pn, pac, fe) = k1(pn, pac, fe) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                        & rho(pa, pac, ff)
        IF (nka > 0) THEN
          ! |1_{j}>|1_{k}>
          ! Annihilate a photon from mode ka
          pn = MAP_n2p(ja, ka, nja, nka-1)
          ! \rho |e><g|
          k1(pn, pac, gg) = k1(pn, pac, gg) + dt * CONJG(gkal(ka)) * &
                          & rho(pa, pac, ge)
          k1(pn, pac, eg) = k1(pn, pac, eg) + dt * CONJG(gkal(ka)) * &
                          & rho(pa, pac, ee)
          k1(pn, pac, fg) = k1(pn, pac, fg) + dt * CONJG(gkal(ka)) * &
                          & rho(pa, pac, fe)
          ! \xi \rho |f><e|
          k1(pn, pac, ge) = k1(pn, pac, ge) + dt * CONJG(gkal(ka)) * xi * &
                          & rho(pa, pac, gf)
          k1(pn, pac, ee) = k1(pn, pac, ee) + dt * CONJG(gkal(ka)) * xi * &
                          & rho(pa, pac, ef)
          k1(pn, pac, fe) = k1(pn, pac, fe) + dt * CONJG(gkal(ka)) * xi * &
                          & rho(pa, pac, ff)
        END IF
      END IF

      ! Close p loop
    END DO
    ! Close pc loop
  END DO

  ! Calculate k2
  k2 = 0.0d0
  DO pac = 1, no_states
    DO pa = 1, no_states
      ! Get mode and photon number
      ja = MAP_p2n(pa, 1)
      jac = MAP_p2n(pac, 1)
      ka = MAP_p2n(pa, 2)
      kac = MAP_p2n(pac, 2)
      nja = MAP_p2n(pa, 3)
      njac = MAP_p2n(pac, 3)
      nka = MAP_p2n(pa, 4)
      nkac = MAP_p2n(pac, 4)

      ! Three-level equations of motion
      k2(pa, pac, gg) = k2(pa, pac, gg) + dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge)) - &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg)) + &
                      & dt * gamma * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
      k2(pa, pac, ge) = k2(pa, pac, ge) + dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + 0.5d0 * k1(pa, pac, gg)) - &
                      & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf)) - &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) + &
                      & dt * gamma * xi * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
      k2(pa, pac, gf) = k2(pa, pac, gf) + dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge)) - &
                      & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf)) - &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
      k2(pa, pac, eg) = k2(pa, pac, eg) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + 0.5d0 * k1(pa, pac, gg)) + &
                      & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg)) + &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) - &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg)) + &
                      & dt * gamma * xi * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
      k2(pa, pac, ee) = k2(pa, pac, ee) - dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge)) + &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg)) - &
                      & dt * gamma * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef)) - &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe)) + &
                      & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
      k2(pa, pac, ef) = k2(pa, pac, ef) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) + &
                      & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef)) - &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
      k2(pa, pac, fg) = k2(pa, pac, fg) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg)) + &
                      & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg)) + &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
      k2(pa, pac, fe) = k2(pa, pac, fe) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) + &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg)) - &
                      & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
      k2(pa, pac, ff) = k2(pa, pac, ff) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe)) - &
                      & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))

      ! Mode Hamiltonians
      k2(pa, pac, :) = k2(pa, pac, :) - i * dt * (&
                     & (wl(ja) * DBLE(nja)) - (wl(jac) * DBLE(njac)) + &
                     & (wl(ka) * DBLE(nka)) - (wl(kac) * DBLE(nkac))) * &
                     & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))

      !-----------------------!
      !  Commutator -i[H, p]  !
      !-----------------------!
      ! -i * H * rho
      ! Cascade Hamiltonian
      ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{i \phi_{j}} a^{\dagger}_{j} \sigma_{-} \rho
      IF (nja < Fock .AND. nka == 0) THEN
        ! Cycle through modes and create a photon
        DO z = -N, N
          IF (z == ja) THEN
            ! Add a photon to ja
            pn = MAP_n2p(ja, ka, nja+1, nka)
            ! |g><e| \rho
            k2(pn, pac, gg) = k2(pn, pac, gg) - dt * gkal(ja) * sqrtl(nja+1) * &
                            & (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg))
            k2(pn, pac, ge) = k2(pn, pac, ge) - dt * gkal(ja) * sqrtl(nja+1) * &
                            & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
            k2(pn, pac, gf) = k2(pn, pac, gf) - dt * gkal(ja) * sqrtl(nja+1) * &
                            & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
            ! \xi |e><f| \rho
            k2(pn, pac, eg) = k2(pn, pac, eg) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                            & (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg))
            k2(pn, pac, ee) = k2(pn, pac, ee) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                            & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
            k2(pn, pac, ef) = k2(pn, pac, ef) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                            & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
          ELSE
            ! Add a photon to other modes
            pn = MAP_n2p(ja, z, nja, nka+1)
            ! |g><e| \rho
            k2(pn, pac, gg) = k2(pn, pac, gg) - dt * gkal(z) * &
                            & (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg))
            k2(pn, pac, ge) = k2(pn, pac, ge) - dt * gkal(z) * &
                            & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
            k2(pn, pac, gf) = k2(pn, pac, gf) - dt * gkal(z) * &
                            & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
            ! \xi |e><f| \rho
            k2(pn, pac, eg) = k2(pn, pac, eg) - dt * gkal(z) * xi * &
                            & (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg))
            k2(pn, pac, ee) = k2(pn, pac, ee) - dt * gkal(z) * xi * &
                            & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
            k2(pn, pac, ef) = k2(pn, pac, ef) - dt * gkal(z) * xi * &
                            & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
          END IF
        END DO
      END IF

      ! i * rho * H
      ! Cascade Hamiltonian
      ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{-i \phi_{j}} \rho \sigma_{+} a_{j}
      IF (njac < Fock .AND. nkac == 0) THEN
        ! Cycle through modes and create a photon
        DO z = -N, N
          IF (z == jac) THEN
            ! Add a photon to mode ja
            pn = MAP_n2p(jac, kac, njac+1, nkac)
            ! \rho |e><g|
            k2(pa, pn, gg) = k2(pa, pn, gg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                           & (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge))
            k2(pa, pn, eg) = k2(pa, pn, eg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                           & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
            k2(pa, pn, fg) = k2(pa, pn, fg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                           & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
            ! \xi \rho |f><e|
            k2(pa, pn, ge) = k2(pa, pn, ge) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                           & (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf))
            k2(pa, pn, ee) = k2(pa, pn, ee) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                           & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
            k2(pa, pn, fe) = k2(pa, pn, fe) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                           & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
          ELSE
            ! Add a photon to other modes
            pn = MAP_n2p(jac, z, njac, nkac+1)
            ! \rho |e><g|
            k2(pa, pn, gg) = k2(pa, pn, gg) - dt * CONJG(gkal(z)) * &
                           & (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge))
            k2(pa, pn, eg) = k2(pa, pn, eg) - dt * CONJG(gkal(z)) * &
                           & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
            k2(pa, pn, fg) = k2(pa, pn, fg) - dt * CONJG(gkal(z)) * &
                           & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
            ! \xi \rho |f><e|
            k2(pa, pn, ge) = k2(pa, pn, ge) - dt * CONJG(gkal(z)) * xi * &
                           & (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf))
            k2(pa, pn, ee) = k2(pa, pn, ee) - dt * CONJG(gkal(z)) * xi * &
                           & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
            k2(pa, pn, fe) = k2(pa, pn, fe) - dt * CONJG(gkal(z)) * xi * &
                           & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
          END IF
        END DO
      END IF

      !-----------------------------!
      !      LINDBLAD DECAY         !
      !-----------------------------!
      ! The atom fluorecence is split by a beam splitter into two paths:
      ! one goes towards the cavity at a fraction \epsilon, and the other
      ! goes to a separate detector at a fraction 1 - \epsilon.
      ! We split the cascaded systems decay operator:
      ! ja = \sqrt{\gamma} \sigma_{-} + \sqrt{\kappa} A,
      ! into separate components correponding to atomic decay, cavity decay
      ! and two cross terms. We also include a separate decay term to account
      ! for the detector placed outside the ring cavity.

      !--------------------!
      !    Cavity Decay    !
      !--------------------!
      ! Lindblad is (accounting for both cascaded and cavity decay)
      ! 2.0 * 0.5 \kappa (2 A \rho A^{\dagger} - A^{\dagger} A \rho
      !                                        - \rho A^{\dagger} A_),
      ! where A is a sum of all mode annihilation operators.

      ! \kappa * A * \rho * A^{\dagger} term
      ! For |nja>_{ja}|0>_{ka} <nja|_{ja'} <0|_{ka'}
      IF (nja > 0 .AND. njac > 0 .AND. ja == jac) THEN
        pn = MAP_n2p(ja, ka, nja-1, nka)
        pnc = MAP_n2p(jac, kac, njac-1, nkac)
        k2(pn, pnc, :) = k2(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * sqrtl(njac) * &
                       & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))
      END IF
      IF (nja > 0 .AND. nkac > 0 .AND. ja == kac) THEN
        pn = MAP_n2p(ja, ka, nja-1, nka)
        pnc = MAP_n2p(jac, kac, njac, nkac-1)
        k2(pn, pnc, :) = k2(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * &
                       & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))
      END IF
      IF (nka > 0 .AND. njac > 0 .AND. ka == jac) THEN
        pn = MAP_n2p(ja, ka, nja, nka-1)
        pnc = MAP_n2p(jac, kac, njac-1, nkac)
        k2(pn, pnc, :) = k2(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(njac) * &
                       & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))
      END IF
      IF (nka > 0 .AND. nkac > 0 .AND. ka == kac) THEN
        pn = MAP_n2p(ja, ka, nja, nka-1)
        pnc = MAP_n2p(jac, kac, njac, nkac-1)
        k2(pn, pnc, :) = k2(pn, pnc, :) + dt * 2.0d0 * kappa * &
                       & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))
      END IF

      ! -a^{\dagger}_{j} a_{j} \rho - \rho a^{\dagger}_{j} a_{j}
      k2(pa, pac, :) = k2(pa, pac, :) - dt * kappa * (DBLE(nja+nka) + DBLE(njac+nkac)) * &
                     & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))

      !--------------------!
      !   Cascade Decay    !
      !--------------------!
      ! \sqrt(2 epsilon gamma kappa) e^{i \phi_{j}} \sigma_{-} \rho a^{\dagger}_{j}
      ! Lower atom |\pm> and annihilate photon from mode jac/kac
      IF (njac > 0) THEN
        ! |1_{j},0_{k}> --> |0>, |2_{j},0_{k}> -> |1_{j},0_{k}>
        ! |1_{j},1_{k}> --> |0_{j},1_{k}>
        ! Annihilate a photon from mode jac
        pn = MAP_n2p(jac, kac, njac-1, nkac)
        ! |g><e| \rho
        k2(pa, pn, gg) = k2(pa, pn, gg) + dt * gkal(jac) * sqrtl(njac) * &
                       & (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg))
        k2(pa, pn, ge) = k2(pa, pn, ge) + dt * gkal(jac) * sqrtl(njac) * &
                       & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
        k2(pa, pn, gf) = k2(pa, pn, gf) + dt * gkal(jac) * sqrtl(njac) * &
                       & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
        ! \xi |e><f| \rho
        k2(pa, pn, eg) = k2(pa, pn, eg) + dt * gkal(jac) * xi * sqrtl(njac) * &
                       & (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg))
        k2(pa, pn, ee) = k2(pa, pn, ee) + dt * gkal(jac) * xi * sqrtl(njac) * &
                       & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
        k2(pa, pn, ef) = k2(pa, pn, ef) + dt * gkal(jac) * xi * sqrtl(njac) * &
                       & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
        IF (nkac > 0) THEN
          ! |1_{j},1_{k}> --> |1_{j},0_{k}>
          ! Annihilate a photon from mode kac
          pn = MAP_n2p(jac, kac, njac, nkac-1)
          ! |g><e| \rho
          k2(pa, pn, gg) = k2(pa, pn, gg) + dt * gkal(kac) * &
                         & (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg))
          k2(pa, pn, ge) = k2(pa, pn, ge) + dt * gkal(kac) * &
                         & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
          k2(pa, pn, gf) = k2(pa, pn, gf) + dt * gkal(kac) * &
                         & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
          ! \xi |e><f| \rho
          k2(pa, pn, eg) = k2(pa, pn, eg) + dt * gkal(kac) * xi * &
                         & (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg))
          k2(pa, pn, ee) = k2(pa, pn, ee) + dt * gkal(kac) * xi * &
                         & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
          k2(pa, pn, ef) = k2(pa, pn, ef) + dt * gkal(kac) * xi * &
                         & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
        END IF
      END IF

      ! \sqrt(2 epsilon gamma kappa) e^{-i \phi_{j}} a_{j} \rho \sigma_{+}
      ! Lower atom <\pm| and annihilate photon from mode ja/ka
      IF (nja > 0) THEN
        ! <1_{j}|<0_{k}| or <2_{j}|<0_{k}|
        ! Annihilate a photon from mode ja
        pn = MAP_n2p(ja, ka, nja-1, nka)
        ! \rho |e><g|
        k2(pn, pac, gg) = k2(pn, pac, gg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                        & (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge))
        k2(pn, pac, eg) = k2(pn, pac, eg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                        & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
        k2(pn, pac, fg) = k2(pn, pac, fg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                        & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
        ! \xi \rho |f><e|
        k2(pn, pac, ge) = k2(pn, pac, ge) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                        & (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf))
        k2(pn, pac, ee) = k2(pn, pac, ee) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                        & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
        k2(pn, pac, fe) = k2(pn, pac, fe) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                        & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
        IF (nka > 0) THEN
          ! |1_{j}>|1_{k}>
          ! Annihilate a photon from mode ka
          pn = MAP_n2p(ja, ka, nja, nka-1)
          ! \rho |e><g|
          k2(pn, pac, gg) = k2(pn, pac, gg) + dt * CONJG(gkal(ka)) * &
                          & (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge))
          k2(pn, pac, eg) = k2(pn, pac, eg) + dt * CONJG(gkal(ka)) * &
                          & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
          k2(pn, pac, fg) = k2(pn, pac, fg) + dt * CONJG(gkal(ka)) * &
                          & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
          ! \xi \rho |f><e|
          k2(pn, pac, ge) = k2(pn, pac, ge) + dt * CONJG(gkal(ka)) * xi * &
                          & (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf))
          k2(pn, pac, ee) = k2(pn, pac, ee) + dt * CONJG(gkal(ka)) * xi * &
                          & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
          k2(pn, pac, fe) = k2(pn, pac, fe) + dt * CONJG(gkal(ka)) * xi * &
                          & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
        END IF
      END IF

      ! Close p loop
    END DO
    ! Close pc loop
  END DO

  ! Calculate k3
  k3 = 0.0d0
  DO pac = 1, no_states
    DO pa = 1, no_states
      ! Get mode and photon number
      ja = MAP_p2n(pa, 1)
      jac = MAP_p2n(pac, 1)
      ka = MAP_p2n(pa, 2)
      kac = MAP_p2n(pac, 2)
      nja = MAP_p2n(pa, 3)
      njac = MAP_p2n(pac, 3)
      nka = MAP_p2n(pa, 4)
      nkac = MAP_p2n(pac, 4)

      ! Three-level equations of motion
      k3(pa, pac, gg) = k3(pa, pac, gg) + dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge)) - &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg)) + &
                      & dt * gamma * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
      k3(pa, pac, ge) = k3(pa, pac, ge) + dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + 0.5d0 * k2(pa, pac, gg)) - &
                      & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf)) - &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) + &
                      & dt * gamma * xi * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
      k3(pa, pac, gf) = k3(pa, pac, gf) + dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge)) - &
                      & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf)) - &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
      k3(pa, pac, eg) = k3(pa, pac, eg) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + 0.5d0 * k2(pa, pac, gg)) + &
                      & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg)) + &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) - &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg)) + &
                      & dt * gamma * xi * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
      k3(pa, pac, ee) = k3(pa, pac, ee) - dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge)) + &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg)) - &
                      & dt * gamma * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef)) - &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe)) + &
                      & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
      k3(pa, pac, ef) = k3(pa, pac, ef) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) + &
                      & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef)) - &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
      k3(pa, pac, fg) = k3(pa, pac, fg) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg)) + &
                      & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg)) + &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
      k3(pa, pac, fe) = k3(pa, pac, fe) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) + &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg)) - &
                      & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
      k3(pa, pac, ff) = k3(pa, pac, ff) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe)) - &
                      & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))

      ! Mode Hamiltonians
      k3(pa, pac, :) = k3(pa, pac, :) - i * dt * (&
                     & (wl(ja) * DBLE(nja)) - (wl(jac) * DBLE(njac)) + &
                     & (wl(ka) * DBLE(nka)) - (wl(kac) * DBLE(nkac))) * &
                     & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))

      !-----------------------!
      !  Commutator -i[H, p]  !
      !-----------------------!
      ! -i * H * rho
      ! Cascade Hamiltonian
      ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{i \phi_{j}} a^{\dagger}_{j} \sigma_{-} \rho
      IF (nja < Fock .AND. nka == 0) THEN
        ! Cycle through modes and create a photon
        DO z = -N, N
          IF (z == ja) THEN
            ! Add a photon to ja
            pn = MAP_n2p(ja, ka, nja+1, nka)
            ! |g><e| \rho
            k3(pn, pac, gg) = k3(pn, pac, gg) - dt * gkal(ja) * sqrtl(nja+1) * &
                            & (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg))
            k3(pn, pac, ge) = k3(pn, pac, ge) - dt * gkal(ja) * sqrtl(nja+1) * &
                            & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
            k3(pn, pac, gf) = k3(pn, pac, gf) - dt * gkal(ja) * sqrtl(nja+1) * &
                            & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
            ! \xi |e><f| \rho
            k3(pn, pac, eg) = k3(pn, pac, eg) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                            & (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg))
            k3(pn, pac, ee) = k3(pn, pac, ee) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                            & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
            k3(pn, pac, ef) = k3(pn, pac, ef) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                            & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
          ELSE
            ! Add a photon to other modes
            pn = MAP_n2p(ja, z, nja, nka+1)
            ! |g><e| \rho
            k3(pn, pac, gg) = k3(pn, pac, gg) - dt * gkal(z) * &
                            & (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg))
            k3(pn, pac, ge) = k3(pn, pac, ge) - dt * gkal(z) * &
                            & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
            k3(pn, pac, gf) = k3(pn, pac, gf) - dt * gkal(z) * &
                            & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
            ! \xi |e><f| \rho
            k3(pn, pac, eg) = k3(pn, pac, eg) - dt * gkal(z) * xi * &
                            & (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg))
            k3(pn, pac, ee) = k3(pn, pac, ee) - dt * gkal(z) * xi * &
                            & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
            k3(pn, pac, ef) = k3(pn, pac, ef) - dt * gkal(z) * xi * &
                            & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
          END IF
        END DO
      END IF

      ! i * rho * H
      ! Cascade Hamiltonian
      ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{-i \phi_{j}} \rho \sigma_{+} a_{j}
      IF (njac < Fock .AND. nkac == 0) THEN
        ! Cycle through modes and create a photon
        DO z = -N, N
          IF (z == jac) THEN
            ! Add a photon to mode ja
            pn = MAP_n2p(jac, kac, njac+1, nkac)
            ! \rho |e><g|
            k3(pa, pn, gg) = k3(pa, pn, gg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                           & (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge))
            k3(pa, pn, eg) = k3(pa, pn, eg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                           & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
            k3(pa, pn, fg) = k3(pa, pn, fg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                           & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
            ! \xi \rho |f><e|
            k3(pa, pn, ge) = k3(pa, pn, ge) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                           & (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf))
            k3(pa, pn, ee) = k3(pa, pn, ee) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                           & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
            k3(pa, pn, fe) = k3(pa, pn, fe) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                           & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
          ELSE
            ! Add a photon to other modes
            pn = MAP_n2p(jac, z, njac, nkac+1)
            ! \rho |e><g|
            k3(pa, pn, gg) = k3(pa, pn, gg) - dt * CONJG(gkal(z)) * &
                           & (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge))
            k3(pa, pn, eg) = k3(pa, pn, eg) - dt * CONJG(gkal(z)) * &
                           & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
            k3(pa, pn, fg) = k3(pa, pn, fg) - dt * CONJG(gkal(z)) * &
                           & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
            ! \xi \rho |f><e|
            k3(pa, pn, ge) = k3(pa, pn, ge) - dt * CONJG(gkal(z)) * xi * &
                           & (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf))
            k3(pa, pn, ee) = k3(pa, pn, ee) - dt * CONJG(gkal(z)) * xi * &
                           & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
            k3(pa, pn, fe) = k3(pa, pn, fe) - dt * CONJG(gkal(z)) * xi * &
                           & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
          END IF
        END DO
      END IF

      !-----------------------------!
      !      LINDBLAD DECAY         !
      !-----------------------------!
      ! The atom fluorecence is split by a beam splitter into two paths:
      ! one goes towards the cavity at a fraction \epsilon, and the other
      ! goes to a separate detector at a fraction 1 - \epsilon.
      ! We split the cascaded systems decay operator:
      ! ja = \sqrt{\gamma} \sigma_{-} + \sqrt{\kappa} A,
      ! into separate components correponding to atomic decay, cavity decay
      ! and two cross terms. We also include a separate decay term to account
      ! for the detector placed outside the ring cavity.

      !--------------------!
      !    Cavity Decay    !
      !--------------------!
      ! Lindblad is (accounting for both cascaded and cavity decay)
      ! 2.0 * 0.5 \kappa (2 A \rho A^{\dagger} - A^{\dagger} A \rho
      !                                        - \rho A^{\dagger} A_),
      ! where A is a sum of all mode annihilation operators.

      ! \kappa * A * \rho * A^{\dagger} term
      ! For |nja>_{ja}|0>_{ka} <nja|_{ja'} <0|_{ka'}
      IF (nja > 0 .AND. njac > 0 .AND. ja == jac) THEN
        pn = MAP_n2p(ja, ka, nja-1, nka)
        pnc = MAP_n2p(jac, kac, njac-1, nkac)
        k3(pn, pnc, :) = k3(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * sqrtl(njac) * &
                       & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))
      END IF
      IF (nja > 0 .AND. nkac > 0 .AND. ja == kac) THEN
        pn = MAP_n2p(ja, ka, nja-1, nka)
        pnc = MAP_n2p(jac, kac, njac, nkac-1)
        k3(pn, pnc, :) = k3(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * &
                       & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))
      END IF
      IF (nka > 0 .AND. njac > 0 .AND. ka == jac) THEN
        pn = MAP_n2p(ja, ka, nja, nka-1)
        pnc = MAP_n2p(jac, kac, njac-1, nkac)
        k3(pn, pnc, :) = k3(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(njac) * &
                       & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))
      END IF
      IF (nka > 0 .AND. nkac > 0 .AND. ka == kac) THEN
        pn = MAP_n2p(ja, ka, nja, nka-1)
        pnc = MAP_n2p(jac, kac, njac, nkac-1)
        k3(pn, pnc, :) = k3(pn, pnc, :) + dt * 2.0d0 * kappa * &
                       & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))
      END IF

      ! -a^{\dagger}_{j} a_{j} \rho - \rho a^{\dagger}_{j} a_{j}
      k3(pa, pac, :) = k3(pa, pac, :) - dt * kappa * (DBLE(nja+nka) + DBLE(njac+nkac)) * &
                     & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))

      !--------------------!
      !   Cascade Decay    !
      !--------------------!
      ! \sqrt(2 epsilon gamma kappa) e^{i \phi_{j}} \sigma_{-} \rho a^{\dagger}_{j}
      ! Lower atom |\pm> and annihilate photon from mode jac/kac
      IF (njac > 0) THEN
        ! |1_{j},0_{k}> --> |0>, |2_{j},0_{k}> -> |1_{j},0_{k}>
        ! |1_{j},1_{k}> --> |0_{j},1_{k}>
        ! Annihilate a photon from mode jac
        pn = MAP_n2p(jac, kac, njac-1, nkac)
        ! |g><e| \rho
        k3(pa, pn, gg) = k3(pa, pn, gg) + dt * gkal(jac) * sqrtl(njac) * &
                       & (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg))
        k3(pa, pn, ge) = k3(pa, pn, ge) + dt * gkal(jac) * sqrtl(njac) * &
                       & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
        k3(pa, pn, gf) = k3(pa, pn, gf) + dt * gkal(jac) * sqrtl(njac) * &
                       & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
        ! \xi |e><f| \rho
        k3(pa, pn, eg) = k3(pa, pn, eg) + dt * gkal(jac) * xi * sqrtl(njac) * &
                       & (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg))
        k3(pa, pn, ee) = k3(pa, pn, ee) + dt * gkal(jac) * xi * sqrtl(njac) * &
                       & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
        k3(pa, pn, ef) = k3(pa, pn, ef) + dt * gkal(jac) * xi * sqrtl(njac) * &
                       & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
        IF (nkac > 0) THEN
          ! |1_{j},1_{k}> --> |1_{j},0_{k}>
          ! Annihilate a photon from mode kac
          pn = MAP_n2p(jac, kac, njac, nkac-1)
          ! |g><e| \rho
          k3(pa, pn, gg) = k3(pa, pn, gg) + dt * gkal(kac) * &
                         & (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg))
          k3(pa, pn, ge) = k3(pa, pn, ge) + dt * gkal(kac) * &
                         & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
          k3(pa, pn, gf) = k3(pa, pn, gf) + dt * gkal(kac) * &
                         & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
          ! \xi |e><f| \rho
          k3(pa, pn, eg) = k3(pa, pn, eg) + dt * gkal(kac) * xi * &
                         & (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg))
          k3(pa, pn, ee) = k3(pa, pn, ee) + dt * gkal(kac) * xi * &
                         & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
          k3(pa, pn, ef) = k3(pa, pn, ef) + dt * gkal(kac) * xi * &
                         & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
        END IF
      END IF

      ! \sqrt(2 epsilon gamma kappa) e^{-i \phi_{j}} a_{j} \rho \sigma_{+}
      ! Lower atom <\pm| and annihilate photon from mode ja/ka
      IF (nja > 0) THEN
        ! <1_{j}|<0_{k}| or <2_{j}|<0_{k}|
        ! Annihilate a photon from mode ja
        pn = MAP_n2p(ja, ka, nja-1, nka)
        ! \rho |e><g|
        k3(pn, pac, gg) = k3(pn, pac, gg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                        & (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge))
        k3(pn, pac, eg) = k3(pn, pac, eg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                        & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
        k3(pn, pac, fg) = k3(pn, pac, fg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                        & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
        ! \xi \rho |f><e|
        k3(pn, pac, ge) = k3(pn, pac, ge) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                        & (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf))
        k3(pn, pac, ee) = k3(pn, pac, ee) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                        & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
        k3(pn, pac, fe) = k3(pn, pac, fe) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                        & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
        IF (nka > 0) THEN
          ! |1_{j}>|1_{k}>
          ! Annihilate a photon from mode ka
          pn = MAP_n2p(ja, ka, nja, nka-1)
          ! \rho |e><g|
          k3(pn, pac, gg) = k3(pn, pac, gg) + dt * CONJG(gkal(ka)) * &
                          & (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge))
          k3(pn, pac, eg) = k3(pn, pac, eg) + dt * CONJG(gkal(ka)) * &
                          & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
          k3(pn, pac, fg) = k3(pn, pac, fg) + dt * CONJG(gkal(ka)) * &
                          & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
          ! \xi \rho |f><e|
          k3(pn, pac, ge) = k3(pn, pac, ge) + dt * CONJG(gkal(ka)) * xi * &
                          & (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf))
          k3(pn, pac, ee) = k3(pn, pac, ee) + dt * CONJG(gkal(ka)) * xi * &
                          & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
          k3(pn, pac, fe) = k3(pn, pac, fe) + dt * CONJG(gkal(ka)) * xi * &
                          & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
        END IF
      END IF

      ! Close p loop
    END DO
    ! Close pc loop
  END DO

  ! Calculate k4
  k4 = 0.0d0
  DO pac = 1, no_states
    DO pa = 1, no_states
      ! Get mode and photon number
      ja = MAP_p2n(pa, 1)
      jac = MAP_p2n(pac, 1)
      ka = MAP_p2n(pa, 2)
      kac = MAP_p2n(pac, 2)
      nja = MAP_p2n(pa, 3)
      njac = MAP_p2n(pac, 3)
      nka = MAP_p2n(pa, 4)
      nkac = MAP_p2n(pac, 4)

      ! Three-level equations of motion
      k4(pa, pac, gg) = k4(pa, pac, gg) + dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + k3(pa, pac, ge)) - &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + k3(pa, pac, eg)) + &
                      & dt * gamma * (rho(pa, pac, ee) + k3(pa, pac, ee))
      k4(pa, pac, ge) = k4(pa, pac, ge) + dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + k3(pa, pac, gg)) - &
                      & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(pa, pac, ge) + k3(pa, pac, ge)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, gf) + k3(pa, pac, gf)) - &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + k3(pa, pac, ee)) + &
                      & dt * gamma * xi * (rho(pa, pac, ef) + k3(pa, pac, ef))
      k4(pa, pac, gf) = k4(pa, pac, gf) + dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ge) + k3(pa, pac, ge)) - &
                      & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, gf) + k3(pa, pac, gf)) - &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, ef) + k3(pa, pac, ef))
      k4(pa, pac, eg) = k4(pa, pac, eg) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + k3(pa, pac, gg)) + &
                      & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(pa, pac, eg) + k3(pa, pac, eg)) + &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + k3(pa, pac, ee)) - &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fg) + k3(pa, pac, fg)) + &
                      & dt * gamma * xi * (rho(pa, pac, fe) + k3(pa, pac, fe))
      k4(pa, pac, ee) = k4(pa, pac, ee) - dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + k3(pa, pac, ge)) + &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + k3(pa, pac, eg)) - &
                      & dt * gamma * (rho(pa, pac, ee) + k3(pa, pac, ee)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + k3(pa, pac, ef)) - &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + k3(pa, pac, fe)) + &
                      & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + k3(pa, pac, ff))
      k4(pa, pac, ef) = k4(pa, pac, ef) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gf) + k3(pa, pac, gf)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + k3(pa, pac, ee)) + &
                      & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, ef) + k3(pa, pac, ef)) - &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + k3(pa, pac, ff))
      k4(pa, pac, fg) = k4(pa, pac, fg) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, eg) + k3(pa, pac, eg)) + &
                      & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, fg) + k3(pa, pac, fg)) + &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, fe) + k3(pa, pac, fe))
      k4(pa, pac, fe) = k4(pa, pac, fe) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + k3(pa, pac, ee)) + &
                      & dt * i * 0.5d0 * Omega * (rho(pa, pac, fg) + k3(pa, pac, fg)) - &
                      & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, fe) + k3(pa, pac, fe)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + k3(pa, pac, ff))
      k4(pa, pac, ff) = k4(pa, pac, ff) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + k3(pa, pac, ef)) + &
                      & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + k3(pa, pac, fe)) - &
                      & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + k3(pa, pac, ff))

      ! Mode Hamiltonians
      k4(pa, pac, :) = k4(pa, pac, :) - i * dt * (&
                     & (wl(ja) * DBLE(nja)) - (wl(jac) * DBLE(njac)) + &
                     & (wl(ka) * DBLE(nka)) - (wl(kac) * DBLE(nkac))) * &
                     & (rho(pa, pac, :) + k3(pa, pac, :))

      !-----------------------!
      !  Commutator -i[H, p]  !
      !-----------------------!
      ! -i * H * rho
      ! Cascade Hamiltonian
      ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{i \phi_{j}} a^{\dagger}_{j} \sigma_{-} \rho
      IF (nja < Fock .AND. nka == 0) THEN
        ! Cycle through modes and create a photon
        DO z = -N, N
          IF (z == ja) THEN
            ! Add a photon to ja
            pn = MAP_n2p(ja, ka, nja+1, nka)
            ! |g><e| \rho
            k4(pn, pac, gg) = k4(pn, pac, gg) - dt * gkal(ja) * sqrtl(nja+1) * &
                            & (rho(pa, pac, eg) + k3(pa, pac, eg))
            k4(pn, pac, ge) = k4(pn, pac, ge) - dt * gkal(ja) * sqrtl(nja+1) * &
                            & (rho(pa, pac, ee) + k3(pa, pac, ee))
            k4(pn, pac, gf) = k4(pn, pac, gf) - dt * gkal(ja) * sqrtl(nja+1) * &
                            & (rho(pa, pac, ef) + k3(pa, pac, ef))
            ! \xi |e><f| \rho
            k4(pn, pac, eg) = k4(pn, pac, eg) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                            & (rho(pa, pac, fg) + k3(pa, pac, fg))
            k4(pn, pac, ee) = k4(pn, pac, ee) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                            & (rho(pa, pac, fe) + k3(pa, pac, fe))
            k4(pn, pac, ef) = k4(pn, pac, ef) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                            & (rho(pa, pac, ff) + k3(pa, pac, ff))
          ELSE
            ! Add a photon to other modes
            pn = MAP_n2p(ja, z, nja, nka+1)
            ! |g><e| \rho
            k4(pn, pac, gg) = k4(pn, pac, gg) - dt * gkal(z) * &
                            & (rho(pa, pac, eg) + k3(pa, pac, eg))
            k4(pn, pac, ge) = k4(pn, pac, ge) - dt * gkal(z) * &
                            & (rho(pa, pac, ee) + k3(pa, pac, ee))
            k4(pn, pac, gf) = k4(pn, pac, gf) - dt * gkal(z) * &
                            & (rho(pa, pac, ef) + k3(pa, pac, ef))
            ! \xi |e><f| \rho
            k4(pn, pac, eg) = k4(pn, pac, eg) - dt * gkal(z) * xi * &
                            & (rho(pa, pac, fg) + k3(pa, pac, fg))
            k4(pn, pac, ee) = k4(pn, pac, ee) - dt * gkal(z) * xi * &
                            & (rho(pa, pac, fe) + k3(pa, pac, fe))
            k4(pn, pac, ef) = k4(pn, pac, ef) - dt * gkal(z) * xi * &
                            & (rho(pa, pac, ff) + k3(pa, pac, ff))
          END IF
        END DO
      END IF

      ! i * rho * H
      ! Cascade Hamiltonian
      ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{-i \phi_{j}} \rho \sigma_{+} a_{j}
      IF (njac < Fock .AND. nkac == 0) THEN
        ! Cycle through modes and create a photon
        DO z = -N, N
          IF (z == jac) THEN
            ! Add a photon to mode ja
            pn = MAP_n2p(jac, kac, njac+1, nkac)
            ! \rho |e><g|
            k4(pa, pn, gg) = k4(pa, pn, gg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                           & (rho(pa, pac, ge) + k3(pa, pac, ge))
            k4(pa, pn, eg) = k4(pa, pn, eg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                           & (rho(pa, pac, ee) + k3(pa, pac, ee))
            k4(pa, pn, fg) = k4(pa, pn, fg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                           & (rho(pa, pac, fe) + k3(pa, pac, fe))
            ! \xi \rho |f><e|
            k4(pa, pn, ge) = k4(pa, pn, ge) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                           & (rho(pa, pac, gf) + k3(pa, pac, gf))
            k4(pa, pn, ee) = k4(pa, pn, ee) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                           & (rho(pa, pac, ef) + k3(pa, pac, ef))
            k4(pa, pn, fe) = k4(pa, pn, fe) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                           & (rho(pa, pac, ff) + k3(pa, pac, ff))
          ELSE
            ! Add a photon to other modes
            pn = MAP_n2p(jac, z, njac, nkac+1)
            ! \rho |e><g|
            k4(pa, pn, gg) = k4(pa, pn, gg) - dt * CONJG(gkal(z)) * &
                           & (rho(pa, pac, ge) + k3(pa, pac, ge))
            k4(pa, pn, eg) = k4(pa, pn, eg) - dt * CONJG(gkal(z)) * &
                           & (rho(pa, pac, ee) + k3(pa, pac, ee))
            k4(pa, pn, fg) = k4(pa, pn, fg) - dt * CONJG(gkal(z)) * &
                           & (rho(pa, pac, fe) + k3(pa, pac, fe))
            ! \xi \rho |f><e|
            k4(pa, pn, ge) = k4(pa, pn, ge) - dt * CONJG(gkal(z)) * xi * &
                           & (rho(pa, pac, gf) + k3(pa, pac, gf))
            k4(pa, pn, ee) = k4(pa, pn, ee) - dt * CONJG(gkal(z)) * xi * &
                           & (rho(pa, pac, ef) + k3(pa, pac, ef))
            k4(pa, pn, fe) = k4(pa, pn, fe) - dt * CONJG(gkal(z)) * xi * &
                           & (rho(pa, pac, ff) + k3(pa, pac, ff))
          END IF
        END DO
      END IF

      !-----------------------------!
      !      LINDBLAD DECAY         !
      !-----------------------------!
      ! The atom fluorecence is split by a beam splitter into two paths:
      ! one goes towards the cavity at a fraction \epsilon, and the other
      ! goes to a separate detector at a fraction 1 - \epsilon.
      ! We split the cascaded systems decay operator:
      ! ja = \sqrt{\gamma} \sigma_{-} + \sqrt{\kappa} A,
      ! into separate components correponding to atomic decay, cavity decay
      ! and two cross terms. We also include a separate decay term to account
      ! for the detector placed outside the ring cavity.

      !--------------------!
      !    Cavity Decay    !
      !--------------------!
      ! Lindblad is (accounting for both cascaded and cavity decay)
      ! 2.0 * 0.5 \kappa (2 A \rho A^{\dagger} - A^{\dagger} A \rho
      !                                        - \rho A^{\dagger} A_),
      ! where A is a sum of all mode annihilation operators.

      ! \kappa * A * \rho * A^{\dagger} term
      ! For |nja>_{ja}|0>_{ka} <nja|_{ja'} <0|_{ka'}
      IF (nja > 0 .AND. njac > 0 .AND. ja == jac) THEN
        pn = MAP_n2p(ja, ka, nja-1, nka)
        pnc = MAP_n2p(jac, kac, njac-1, nkac)
        k4(pn, pnc, :) = k4(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * sqrtl(njac) * &
                       & (rho(pa, pac, :) + k3(pa, pac, :))
      END IF
      IF (nja > 0 .AND. nkac > 0 .AND. ja == kac) THEN
        pn = MAP_n2p(ja, ka, nja-1, nka)
        pnc = MAP_n2p(jac, kac, njac, nkac-1)
        k4(pn, pnc, :) = k4(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * &
                       & (rho(pa, pac, :) + k3(pa, pac, :))
      END IF
      IF (nka > 0 .AND. njac > 0 .AND. ka == jac) THEN
        pn = MAP_n2p(ja, ka, nja, nka-1)
        pnc = MAP_n2p(jac, kac, njac-1, nkac)
        k4(pn, pnc, :) = k4(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(njac) * &
                       & (rho(pa, pac, :) + k3(pa, pac, :))
      END IF
      IF (nka > 0 .AND. nkac > 0 .AND. ka == kac) THEN
        pn = MAP_n2p(ja, ka, nja, nka-1)
        pnc = MAP_n2p(jac, kac, njac, nkac-1)
        k4(pn, pnc, :) = k4(pn, pnc, :) + dt * 2.0d0 * kappa * &
                       & (rho(pa, pac, :) + k3(pa, pac, :))
      END IF

      ! -a^{\dagger}_{j} a_{j} \rho - \rho a^{\dagger}_{j} a_{j}
      k4(pa, pac, :) = k4(pa, pac, :) - dt * kappa * (DBLE(nja+nka) + DBLE(njac+nkac)) * &
                     & (rho(pa, pac, :) + k3(pa, pac, :))

      !--------------------!
      !   Cascade Decay    !
      !--------------------!
      ! \sqrt(2 epsilon gamma kappa) e^{i \phi_{j}} \sigma_{-} \rho a^{\dagger}_{j}
      ! Lower atom |\pm> and annihilate photon from mode jac/kac
      IF (njac > 0) THEN
        ! |1_{j},0_{k}> --> |0>, |2_{j},0_{k}> -> |1_{j},0_{k}>
        ! |1_{j},1_{k}> --> |0_{j},1_{k}>
        ! Annihilate a photon from mode jac
        pn = MAP_n2p(jac, kac, njac-1, nkac)
        ! |g><e| \rho
        k4(pa, pn, gg) = k4(pa, pn, gg) + dt * gkal(jac) * sqrtl(njac) * &
                       & (rho(pa, pac, eg) + k3(pa, pac, eg))
        k4(pa, pn, ge) = k4(pa, pn, ge) + dt * gkal(jac) * sqrtl(njac) * &
                       & (rho(pa, pac, ee) + k3(pa, pac, ee))
        k4(pa, pn, gf) = k4(pa, pn, gf) + dt * gkal(jac) * sqrtl(njac) * &
                       & (rho(pa, pac, ef) + k3(pa, pac, ef))
        ! \xi |e><f| \rho
        k4(pa, pn, eg) = k4(pa, pn, eg) + dt * gkal(jac) * xi * sqrtl(njac) * &
                       & (rho(pa, pac, fg) + k3(pa, pac, fg))
        k4(pa, pn, ee) = k4(pa, pn, ee) + dt * gkal(jac) * xi * sqrtl(njac) * &
                       & (rho(pa, pac, fe) + k3(pa, pac, fe))
        k4(pa, pn, ef) = k4(pa, pn, ef) + dt * gkal(jac) * xi * sqrtl(njac) * &
                       & (rho(pa, pac, ff) + k3(pa, pac, ff))
        IF (nkac > 0) THEN
          ! |1_{j},1_{k}> --> |1_{j},0_{k}>
          ! Annihilate a photon from mode kac
          pn = MAP_n2p(jac, kac, njac, nkac-1)
          ! |g><e| \rho
          k4(pa, pn, gg) = k4(pa, pn, gg) + dt * gkal(kac) * &
                         & (rho(pa, pac, eg) + k3(pa, pac, eg))
          k4(pa, pn, ge) = k4(pa, pn, ge) + dt * gkal(kac) * &
                         & (rho(pa, pac, ee) + k3(pa, pac, ee))
          k4(pa, pn, gf) = k4(pa, pn, gf) + dt * gkal(kac) * &
                         & (rho(pa, pac, ef) + k3(pa, pac, ef))
          ! \xi |e><f| \rho
          k4(pa, pn, eg) = k4(pa, pn, eg) + dt * gkal(kac) * xi * &
                         & (rho(pa, pac, fg) + k3(pa, pac, fg))
          k4(pa, pn, ee) = k4(pa, pn, ee) + dt * gkal(kac) * xi * &
                         & (rho(pa, pac, fe) + k3(pa, pac, fe))
          k4(pa, pn, ef) = k4(pa, pn, ef) + dt * gkal(kac) * xi * &
                         & (rho(pa, pac, ff) + k3(pa, pac, ff))
        END IF
      END IF

      ! \sqrt(2 epsilon gamma kappa) e^{-i \phi_{j}} a_{j} \rho \sigma_{+}
      ! Lower atom <\pm| and annihilate photon from mode ja/ka
      IF (nja > 0) THEN
        ! <1_{j}|<0_{k}| or <2_{j}|<0_{k}|
        ! Annihilate a photon from mode ja
        pn = MAP_n2p(ja, ka, nja-1, nka)
        ! \rho |e><g|
        k4(pn, pac, gg) = k4(pn, pac, gg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                        & (rho(pa, pac, ge) + k3(pa, pac, ge))
        k4(pn, pac, eg) = k4(pn, pac, eg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                        & (rho(pa, pac, ee) + k3(pa, pac, ee))
        k4(pn, pac, fg) = k4(pn, pac, fg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                        & (rho(pa, pac, fe) + k3(pa, pac, fe))
        ! \xi \rho |f><e|
        k4(pn, pac, ge) = k4(pn, pac, ge) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                        & (rho(pa, pac, gf) + k3(pa, pac, gf))
        k4(pn, pac, ee) = k4(pn, pac, ee) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                        & (rho(pa, pac, ef) + k3(pa, pac, ef))
        k4(pn, pac, fe) = k4(pn, pac, fe) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                        & (rho(pa, pac, ff) + k3(pa, pac, ff))
        IF (nka > 0) THEN
          ! |1_{j}>|1_{k}>
          ! Annihilate a photon from mode ka
          pn = MAP_n2p(ja, ka, nja, nka-1)
          ! \rho |e><g|
          k4(pn, pac, gg) = k4(pn, pac, gg) + dt * CONJG(gkal(ka)) * &
                          & (rho(pa, pac, ge) + k3(pa, pac, ge))
          k4(pn, pac, eg) = k4(pn, pac, eg) + dt * CONJG(gkal(ka)) * &
                          & (rho(pa, pac, ee) + k3(pa, pac, ee))
          k4(pn, pac, fg) = k4(pn, pac, fg) + dt * CONJG(gkal(ka)) * &
                          & (rho(pa, pac, fe) + k3(pa, pac, fe))
          ! \xi \rho |f><e|
          k4(pn, pac, ge) = k4(pn, pac, ge) + dt * CONJG(gkal(ka)) * xi * &
                          & (rho(pa, pac, gf) + k3(pa, pac, gf))
          k4(pn, pac, ee) = k4(pn, pac, ee) + dt * CONJG(gkal(ka)) * xi * &
                          & (rho(pa, pac, ef) + k3(pa, pac, ef))
          k4(pn, pac, fe) = k4(pn, pac, fe) + dt * CONJG(gkal(ka)) * xi * &
                          & (rho(pa, pac, ff) + k3(pa, pac, ff))
        END IF
      END IF

      ! Close p loop
    END DO
    ! Close pc loop
  END DO

  ! Update rho
  rho = rho + xis * (k1 + 2.0d0 * (k2 + k3) + k4)

  ! Check percentage
  IF (ten_percent /= 0) THEN
    IF (MOD(t, ten_percent) == 0 .AND. t /= 0) THEN
      CALL CPU_TIME(loop_check_time)
      percentage = NINT((100.0 * t) / (1.0 * t_steps))
      loop_run_time = loop_check_time - loop_start_time
      loop_remaining_time = ((100.0 * (loop_check_time - loop_start_time)) / (1.0 * percentage)) - loop_run_time
      WRITE(*, FMT_ss) percentage, "%. Run time (steady state): ", loop_run_time, &
                     & "s. Est. time left: ", loop_remaining_time, "s"
    END IF
  END IF

  ! Close time integration DO loop
END DO

! Close state file
CLOSE(2)

photona_ss = photona
PRINT*, " "
PRINT*, "Steady state found!"
PRINT*, "Mean photon number in cavity A =", photona_ss

! Save steady state density operator
ALLOCATE(rho_ss(no_states, no_states, 9))
rho_ss = rho

! Change into single-photon truncation to make calculations way faster
IF (tau2_max > 0.0) THEN
  ! Reallocate density operator
  DEALLOCATE(rho)
  ALLOCATE(rho(no_states_one_photon, no_states_one_photon, 9))
  rho = 0.0
  ! Reallocate Runge-Kutta vectors
  DEALLOCATE(k1)
  ALLOCATE(k1(no_states_one_photon, no_states_one_photon, 9))
  k1 = 0.0
  DEALLOCATE(k2)
  ALLOCATE(k2(no_states_one_photon, no_states_one_photon, 9))
  k2 = 0.0
  DEALLOCATE(k3)
  ALLOCATE(k3(no_states_one_photon, no_states_one_photon, 9))
  k3 = 0.0
  DEALLOCATE(k4)
  ALLOCATE(k4(no_states_one_photon, no_states_one_photon, 9))
  k4 = 0.0
END IF

!##############################################################################!
!                        FIRST-ORDER CORRELATION FUNCTION                      !
!##############################################################################!
IF (tau1_max > 0.0) THEN
  ! Propagate \rho(ss) * A^{\dagger} (remove a photon from <bra|)
  rho = 0.0d0
  DO pac = 1, no_states
    DO pa = 1, no_states_one_photon
      jac = MAP_p2n(pac, 1)
      kac = MAP_p2n(pac, 2)
      njac = MAP_p2n(pac, 3)
      nkac = MAP_p2n(pac, 4)
      ! If njac > 0, annihilate a photon
      IF (njac > 0) THEN
        pn = MAP_n2p(jac, kac, njac-1, nkac)
        rho(pa, pn, :) = rho(pa, pn, :) + sqrtl(njac) * rho_ss(pa, pac, :)
      END IF
      IF (nkac > 0) THEN
        pn = MAP_n2p(jac, kac, njac, nkac-1)
        rho(pa, pn, :) = rho(pa, pn, :) + sqrtl(nkac) * rho_ss(pa, pac, :)
      END IF
    END DO
  END DO

  ! Allocate rho_corr
  ALLOCATE(rho_corr(no_states_one_photon, no_states_one_photon, 9))
  rho_corr = 0.0

  ! Open file to write time and correlation to
  OPEN(UNIT=3, FILE=filename_g1, STATUS='REPLACE', ACTION='WRITE', RECL=4000)

  PRINT*, " "
  PRINT*, "Calculating first-order correlation ..."
  PRINT*, " "

  ! Change tau_steps for g1
  tau_steps = NINT(tau1_max / dt)
  ! Ten percent of time steps
  ten_percent = NINT((1.0 * tau_steps / 10.0))
  ! Call CPU clock time
  CALL CPU_TIME(loop_start_time)

  ! Change number of states to single_photon
  no_states = no_states_one_photon
  Fock = 1

  ! Start time integration to find first-order correlation
  DO t = 0, tau_steps
    ! Update rho_corr and leave rho intact
    rho_corr = 0.0d0
    ! rho_corr = rho
    ! Propagate <A \rho>
    DO pac = 1, no_states
      DO pa = 1, no_states
        ja = MAP_p2n(pa, 1)
        ka = MAP_p2n(pa, 2)
        nja = MAP_p2n(pa, 3)
        nka = MAP_p2n(pa, 4)
        ! If nja > 0, annihilate a photon
        IF (nja > 0) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          rho_corr(pn, pac, :) = rho_corr(pn, pac, :) + sqrtl(nja) * rho(pa, pac, :)
        END IF
        IF (nka > 0) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          rho_corr(pn, pac, :) = rho_corr(pn, pac, :) + sqrtl(nka) * rho(pa, pac, :)
        END IF
      END DO
    END DO
    ! Calculate correlation
    corr = 0.0d0
    DO pa = 1, no_states
      corr = corr + (rho_corr(pa, pa, gg) + rho_corr(pa, pa, ee) + rho_corr(pa, pa, ff))
    END DO

    ! Normalise by steady-state mean-photon number
    IF (photona_ss /= 0.0) THEN
      corr = corr / (photona_ss)
    END IF

    ! Write data
    WRITE(3,*) DBLE(t) * dt, REAL(corr), -IMAG(corr)

    ! Calculate k1
    k1 = 0.0d0
    DO pac = 1, no_states
      DO pa = 1, no_states
        ! Get mode and photon number
        ja = MAP_p2n(pa, 1)
        jac = MAP_p2n(pac, 1)
        ka = MAP_p2n(pa, 2)
        kac = MAP_p2n(pac, 2)
        nja = MAP_p2n(pa, 3)
        njac = MAP_p2n(pac, 3)
        nka = MAP_p2n(pa, 4)
        nkac = MAP_p2n(pac, 4)

        ! Three-level equations of motion
        k1(pa, pac, gg) = k1(pa, pac, gg) + dt * i * 0.5d0 * Omega * rho(pa, pac, ge) - &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, eg) + &
                        & dt * gamma * rho(pa, pac, ee)
        k1(pa, pac, ge) = k1(pa, pac, ge) + dt * i * 0.5d0 * Omega * rho(pa, pac, gg) - &
                        & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * rho(pa, pac, ge) + &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, gf) - &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, ee) + &
                        & dt * gamma * xi * rho(pa, pac, ef)
        k1(pa, pac, gf) = k1(pa, pac, gf) + dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ge) - &
                        & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * rho(pa, pac, gf) - &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, ef)
        k1(pa, pac, eg) = k1(pa, pac, eg) - dt * i * 0.5d0 * Omega * rho(pa, pac, gg) + &
                        & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * rho(pa, pac, eg) + &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, ee) - &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, fg) + &
                        & dt * gamma * xi * rho(pa, pac, fe)
        k1(pa, pac, ee) = k1(pa, pac, ee) - dt * i * 0.5d0 * Omega * rho(pa, pac, ge) + &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, eg) - &
                        & dt * gamma * rho(pa, pac, ee) + &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ef) - &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, fe) + &
                        & dt * gamma * (xi ** 2) * rho(pa, pac, ff)
        k1(pa, pac, ef) = k1(pa, pac, ef) - dt * i * 0.5d0 * Omega * rho(pa, pac, gf) + &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ee) + &
                        & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * rho(pa, pac, ef) - &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ff)
        k1(pa, pac, fg) = k1(pa, pac, fg) - dt * i * 0.5d0 * xi * Omega * rho(pa, pac, eg) + &
                        & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * rho(pa, pac, fg) + &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, fe)
        k1(pa, pac, fe) = k1(pa, pac, fe) - dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ee) + &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, fg) - &
                        & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * rho(pa, pac, fe) + &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ff)
        k1(pa, pac, ff) = k1(pa, pac, ff) - dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ef) + &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, fe) - &
                        & dt * gamma * (xi ** 2) * rho(pa, pac, ff)

        ! Mode Hamiltonians
        k1(pa, pac, :) = k1(pa, pac, :) - i * dt * (&
                       & (wl(ja) * DBLE(nja)) - (wl(jac) * DBLE(njac)) + &
                       & (wl(ka) * DBLE(nka)) - (wl(kac) * DBLE(nkac))) * &
                       & rho(pa, pac, :)

        !-----------------------!
        !  Commutator -i[H, p]  !
        !-----------------------!
        ! -i * H * rho
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{i \phi_{j}} a^{\dagger}_{j} \sigma_{-} \rho
        IF (nja < Fock .AND. nka == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == ja) THEN
              ! Add a photon to ja
              pn = MAP_n2p(ja, ka, nja+1, nka)
              ! |g><e| \rho
              k1(pn, pac, gg) = k1(pn, pac, gg) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & rho(pa, pac, eg)
              k1(pn, pac, ge) = k1(pn, pac, ge) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & rho(pa, pac, ee)
              k1(pn, pac, gf) = k1(pn, pac, gf) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & rho(pa, pac, ef)
              ! \xi |e><f| \rho
              k1(pn, pac, eg) = k1(pn, pac, eg) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & rho(pa, pac, fg)
              k1(pn, pac, ee) = k1(pn, pac, ee) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & rho(pa, pac, fe)
              k1(pn, pac, ef) = k1(pn, pac, ef) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & rho(pa, pac, ff)
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(ja, z, nja, nka+1)
              ! |g><e| \rho
              k1(pn, pac, gg) = k1(pn, pac, gg) - dt * gkal(z) * &
                              & rho(pa, pac, eg)
              k1(pn, pac, ge) = k1(pn, pac, ge) - dt * gkal(z) * &
                              & rho(pa, pac, ee)
              k1(pn, pac, gf) = k1(pn, pac, gf) - dt * gkal(z) * &
                              & rho(pa, pac, ef)
              ! \xi |e><f| \rho
              k1(pn, pac, eg) = k1(pn, pac, eg) - dt * gkal(z) * xi * &
                              & rho(pa, pac, fg)
              k1(pn, pac, ee) = k1(pn, pac, ee) - dt * gkal(z) * xi * &
                              & rho(pa, pac, fe)
              k1(pn, pac, ef) = k1(pn, pac, ef) - dt * gkal(z) * xi * &
                              & rho(pa, pac, ff)
            END IF
          END DO
        END IF

        ! i * rho * H
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{-i \phi_{j}} \rho \sigma_{+} a_{j}
        IF (njac < Fock .AND. nkac == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == jac) THEN
              ! Add a photon to mode ja
              pn = MAP_n2p(jac, kac, njac+1, nkac)
              ! \rho |e><g|
              k1(pa, pn, gg) = k1(pa, pn, gg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & rho(pa, pac, ge)
              k1(pa, pn, eg) = k1(pa, pn, eg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & rho(pa, pac, ee)
              k1(pa, pn, fg) = k1(pa, pn, fg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & rho(pa, pac, fe)
              ! \xi \rho |f><e|
              k1(pa, pn, ge) = k1(pa, pn, ge) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & rho(pa, pac, gf)
              k1(pa, pn, ee) = k1(pa, pn, ee) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & rho(pa, pac, ef)
              k1(pa, pn, fe) = k1(pa, pn, fe) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & rho(pa, pac, ff)
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(jac, z, njac, nkac+1)
              ! \rho |e><g|
              k1(pa, pn, gg) = k1(pa, pn, gg) - dt * CONJG(gkal(z)) * &
                             & rho(pa, pac, ge)
              k1(pa, pn, eg) = k1(pa, pn, eg) - dt * CONJG(gkal(z)) * &
                             & rho(pa, pac, ee)
              k1(pa, pn, fg) = k1(pa, pn, fg) - dt * CONJG(gkal(z)) * &
                             & rho(pa, pac, fe)
              ! \xi \rho |f><e|
              k1(pa, pn, ge) = k1(pa, pn, ge) - dt * CONJG(gkal(z)) * xi * &
                             & rho(pa, pac, gf)
              k1(pa, pn, ee) = k1(pa, pn, ee) - dt * CONJG(gkal(z)) * xi * &
                             & rho(pa, pac, ef)
              k1(pa, pn, fe) = k1(pa, pn, fe) - dt * CONJG(gkal(z)) * xi * &
                             & rho(pa, pac, ff)
            END IF
          END DO
        END IF

        !-----------------------------!
        !      LINDBLAD DECAY         !
        !-----------------------------!
        ! The atom fluorecence is split by a beam splitter into two paths:
        ! one goes towards the cavity at a fraction \epsilon, and the other
        ! goes to a separate detector at a fraction 1 - \epsilon.
        ! We split the cascaded systems decay operator:
        ! ja = \sqrt{\gamma} \sigma_{-} + \sqrt{\kappa} A,
        ! into separate components correponding to atomic decay, cavity decay
        ! and two cross terms. We also include a separate decay term to account
        ! for the detector placed outside the ring cavity.

        !--------------------!
        !    Cavity Decay    !
        !--------------------!
        ! Lindblad is (accounting for both cascaded and cavity decay)
        ! 2.0 * 0.5 \kappa (2 A \rho A^{\dagger} - A^{\dagger} A \rho
        !                                        - \rho A^{\dagger} A_),
        ! where A is a sum of all mode annihilation operators.

        ! \kappa * A * \rho * A^{\dagger} term
        ! For |nja>_{ja}|0>_{ka} <nja|_{ja'} <0|_{ka'}
        IF (nja > 0 .AND. njac > 0 .AND. ja == jac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k1(pn, pnc, :) = k1(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * sqrtl(njac) * &
                         & rho(pa, pac, :)
        END IF
        IF (nja > 0 .AND. nkac > 0 .AND. ja == kac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k1(pn, pnc, :) = k1(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * &
                         & rho(pa, pac, :)
        END IF
        IF (nka > 0 .AND. njac > 0 .AND. ka == jac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k1(pn, pnc, :) = k1(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(njac) * &
                         & rho(pa, pac, :)
        END IF
        IF (nka > 0 .AND. nkac > 0 .AND. ka == kac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k1(pn, pnc, :) = k1(pn, pnc, :) + dt * 2.0d0 * kappa * &
                         & rho(pa, pac, :)
        END IF

        ! -a^{\dagger}_{j} a_{j} \rho - \rho a^{\dagger}_{j} a_{j}
        k1(pa, pac, :) = k1(pa, pac, :) - dt * kappa * (DBLE(nja+nka) + DBLE(njac+nkac)) * &
                       & rho(pa, pac, :)

        !--------------------!
        !   Cascade Decay    !
        !--------------------!
        ! \sqrt(2 epsilon gamma kappa) e^{i \phi_{j}} \sigma_{-} \rho a^{\dagger}_{j}
        ! Lower atom |\pm> and annihilate photon from mode jac/kac
        IF (njac > 0) THEN
          ! |1_{j},0_{k}> --> |0>, |2_{j},0_{k}> -> |1_{j},0_{k}>
          ! |1_{j},1_{k}> --> |0_{j},1_{k}>
          ! Annihilate a photon from mode jac
          pn = MAP_n2p(jac, kac, njac-1, nkac)
          ! |g><e| \rho
          k1(pa, pn, gg) = k1(pa, pn, gg) + dt * gkal(jac) * sqrtl(njac) * &
                         & rho(pa, pac, eg)
          k1(pa, pn, ge) = k1(pa, pn, ge) + dt * gkal(jac) * sqrtl(njac) * &
                         & rho(pa, pac, ee)
          k1(pa, pn, gf) = k1(pa, pn, gf) + dt * gkal(jac) * sqrtl(njac) * &
                         & rho(pa, pac, ef)
          ! \xi |e><f| \rho
          k1(pa, pn, eg) = k1(pa, pn, eg) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & rho(pa, pac, fg)
          k1(pa, pn, ee) = k1(pa, pn, ee) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & rho(pa, pac, fe)
          k1(pa, pn, ef) = k1(pa, pn, ef) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & rho(pa, pac, ff)
          IF (nkac > 0) THEN
            ! |1_{j},1_{k}> --> |1_{j},0_{k}>
            ! Annihilate a photon from mode kac
            pn = MAP_n2p(jac, kac, njac, nkac-1)
            ! |g><e| \rho
            k1(pa, pn, gg) = k1(pa, pn, gg) + dt * gkal(kac) * &
                           & rho(pa, pac, eg)
            k1(pa, pn, ge) = k1(pa, pn, ge) + dt * gkal(kac) * &
                           & rho(pa, pac, ee)
            k1(pa, pn, gf) = k1(pa, pn, gf) + dt * gkal(kac) * &
                           & rho(pa, pac, ef)
            ! \xi |e><f| \rho
            k1(pa, pn, eg) = k1(pa, pn, eg) + dt * gkal(kac) * xi * &
                           & rho(pa, pac, fg)
            k1(pa, pn, ee) = k1(pa, pn, ee) + dt * gkal(kac) * xi * &
                           & rho(pa, pac, fe)
            k1(pa, pn, ef) = k1(pa, pn, ef) + dt * gkal(kac) * xi * &
                           & rho(pa, pac, ff)
          END IF
        END IF

        ! \sqrt(2 epsilon gamma kappa) e^{-i \phi_{j}} a_{j} \rho \sigma_{+}
        ! Lower atom <\pm| and annihilate photon from mode ja/ka
        IF (nja > 0) THEN
          ! <1_{j}|<0_{k}| or <2_{j}|<0_{k}|
          ! Annihilate a photon from mode ja
          pn = MAP_n2p(ja, ka, nja-1, nka)
          ! \rho |e><g|
          k1(pn, pac, gg) = k1(pn, pac, gg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & rho(pa, pac, ge)
          k1(pn, pac, eg) = k1(pn, pac, eg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & rho(pa, pac, ee)
          k1(pn, pac, fg) = k1(pn, pac, fg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & rho(pa, pac, fe)
          ! \xi \rho |f><e|
          k1(pn, pac, ge) = k1(pn, pac, ge) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & rho(pa, pac, gf)
          k1(pn, pac, ee) = k1(pn, pac, ee) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & rho(pa, pac, ef)
          k1(pn, pac, fe) = k1(pn, pac, fe) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & rho(pa, pac, ff)
          IF (nka > 0) THEN
            ! |1_{j}>|1_{k}>
            ! Annihilate a photon from mode ka
            pn = MAP_n2p(ja, ka, nja, nka-1)
            ! \rho |e><g|
            k1(pn, pac, gg) = k1(pn, pac, gg) + dt * CONJG(gkal(ka)) * &
                            & rho(pa, pac, ge)
            k1(pn, pac, eg) = k1(pn, pac, eg) + dt * CONJG(gkal(ka)) * &
                            & rho(pa, pac, ee)
            k1(pn, pac, fg) = k1(pn, pac, fg) + dt * CONJG(gkal(ka)) * &
                            & rho(pa, pac, fe)
            ! \xi \rho |f><e|
            k1(pn, pac, ge) = k1(pn, pac, ge) + dt * CONJG(gkal(ka)) * xi * &
                            & rho(pa, pac, gf)
            k1(pn, pac, ee) = k1(pn, pac, ee) + dt * CONJG(gkal(ka)) * xi * &
                            & rho(pa, pac, ef)
            k1(pn, pac, fe) = k1(pn, pac, fe) + dt * CONJG(gkal(ka)) * xi * &
                            & rho(pa, pac, ff)
          END IF
        END IF

        ! Close p loop
      END DO
      ! Close pc loop
    END DO

    ! Calculate k2
    k2 = 0.0d0
    DO pac = 1, no_states
      DO pa = 1, no_states
        ! Get mode and photon number
        ja = MAP_p2n(pa, 1)
        jac = MAP_p2n(pac, 1)
        ka = MAP_p2n(pa, 2)
        kac = MAP_p2n(pac, 2)
        nja = MAP_p2n(pa, 3)
        njac = MAP_p2n(pac, 3)
        nka = MAP_p2n(pa, 4)
        nkac = MAP_p2n(pac, 4)

        ! Three-level equations of motion
        k2(pa, pac, gg) = k2(pa, pac, gg) + dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg)) + &
                        & dt * gamma * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
        k2(pa, pac, ge) = k2(pa, pac, ge) + dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + 0.5d0 * k1(pa, pac, gg)) - &
                        & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) + &
                        & dt * gamma * xi * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
        k2(pa, pac, gf) = k2(pa, pac, gf) + dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge)) - &
                        & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
        k2(pa, pac, eg) = k2(pa, pac, eg) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + 0.5d0 * k1(pa, pac, gg)) + &
                        & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg)) + &
                        & dt * gamma * xi * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
        k2(pa, pac, ee) = k2(pa, pac, ee) - dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg)) - &
                        & dt * gamma * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe)) + &
                        & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
        k2(pa, pac, ef) = k2(pa, pac, ef) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) + &
                        & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
        k2(pa, pac, fg) = k2(pa, pac, fg) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg)) + &
                        & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
        k2(pa, pac, fe) = k2(pa, pac, fe) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg)) - &
                        & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
        k2(pa, pac, ff) = k2(pa, pac, ff) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe)) - &
                        & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))

        ! Mode Hamiltonians
        k2(pa, pac, :) = k2(pa, pac, :) - i * dt * (&
                       & (wl(ja) * DBLE(nja)) - (wl(jac) * DBLE(njac)) + &
                       & (wl(ka) * DBLE(nka)) - (wl(kac) * DBLE(nkac))) * &
                       & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))

        !-----------------------!
        !  Commutator -i[H, p]  !
        !-----------------------!
        ! -i * H * rho
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{i \phi_{j}} a^{\dagger}_{j} \sigma_{-} \rho
        IF (nja < Fock .AND. nka == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == ja) THEN
              ! Add a photon to ja
              pn = MAP_n2p(ja, ka, nja+1, nka)
              ! |g><e| \rho
              k2(pn, pac, gg) = k2(pn, pac, gg) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg))
              k2(pn, pac, ge) = k2(pn, pac, ge) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
              k2(pn, pac, gf) = k2(pn, pac, gf) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
              ! \xi |e><f| \rho
              k2(pn, pac, eg) = k2(pn, pac, eg) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg))
              k2(pn, pac, ee) = k2(pn, pac, ee) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
              k2(pn, pac, ef) = k2(pn, pac, ef) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(ja, z, nja, nka+1)
              ! |g><e| \rho
              k2(pn, pac, gg) = k2(pn, pac, gg) - dt * gkal(z) * &
                              & (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg))
              k2(pn, pac, ge) = k2(pn, pac, ge) - dt * gkal(z) * &
                              & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
              k2(pn, pac, gf) = k2(pn, pac, gf) - dt * gkal(z) * &
                              & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
              ! \xi |e><f| \rho
              k2(pn, pac, eg) = k2(pn, pac, eg) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg))
              k2(pn, pac, ee) = k2(pn, pac, ee) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
              k2(pn, pac, ef) = k2(pn, pac, ef) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
            END IF
          END DO
        END IF

        ! i * rho * H
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{-i \phi_{j}} \rho \sigma_{+} a_{j}
        IF (njac < Fock .AND. nkac == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == jac) THEN
              ! Add a photon to mode ja
              pn = MAP_n2p(jac, kac, njac+1, nkac)
              ! \rho |e><g|
              k2(pa, pn, gg) = k2(pa, pn, gg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge))
              k2(pa, pn, eg) = k2(pa, pn, eg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
              k2(pa, pn, fg) = k2(pa, pn, fg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
              ! \xi \rho |f><e|
              k2(pa, pn, ge) = k2(pa, pn, ge) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf))
              k2(pa, pn, ee) = k2(pa, pn, ee) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
              k2(pa, pn, fe) = k2(pa, pn, fe) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(jac, z, njac, nkac+1)
              ! \rho |e><g|
              k2(pa, pn, gg) = k2(pa, pn, gg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge))
              k2(pa, pn, eg) = k2(pa, pn, eg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
              k2(pa, pn, fg) = k2(pa, pn, fg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
              ! \xi \rho |f><e|
              k2(pa, pn, ge) = k2(pa, pn, ge) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf))
              k2(pa, pn, ee) = k2(pa, pn, ee) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
              k2(pa, pn, fe) = k2(pa, pn, fe) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
            END IF
          END DO
        END IF

        !-----------------------------!
        !      LINDBLAD DECAY         !
        !-----------------------------!
        ! The atom fluorecence is split by a beam splitter into two paths:
        ! one goes towards the cavity at a fraction \epsilon, and the other
        ! goes to a separate detector at a fraction 1 - \epsilon.
        ! We split the cascaded systems decay operator:
        ! ja = \sqrt{\gamma} \sigma_{-} + \sqrt{\kappa} A,
        ! into separate components correponding to atomic decay, cavity decay
        ! and two cross terms. We also include a separate decay term to account
        ! for the detector placed outside the ring cavity.

        !--------------------!
        !    Cavity Decay    !
        !--------------------!
        ! Lindblad is (accounting for both cascaded and cavity decay)
        ! 2.0 * 0.5 \kappa (2 A \rho A^{\dagger} - A^{\dagger} A \rho
        !                                        - \rho A^{\dagger} A_),
        ! where A is a sum of all mode annihilation operators.

        ! \kappa * A * \rho * A^{\dagger} term
        ! For |nja>_{ja}|0>_{ka} <nja|_{ja'} <0|_{ka'}
        IF (nja > 0 .AND. njac > 0 .AND. ja == jac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k2(pn, pnc, :) = k2(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * sqrtl(njac) * &
                         & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))
        END IF
        IF (nja > 0 .AND. nkac > 0 .AND. ja == kac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k2(pn, pnc, :) = k2(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * &
                         & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))
        END IF
        IF (nka > 0 .AND. njac > 0 .AND. ka == jac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k2(pn, pnc, :) = k2(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(njac) * &
                         & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))
        END IF
        IF (nka > 0 .AND. nkac > 0 .AND. ka == kac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k2(pn, pnc, :) = k2(pn, pnc, :) + dt * 2.0d0 * kappa * &
                         & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))
        END IF

        ! -a^{\dagger}_{j} a_{j} \rho - \rho a^{\dagger}_{j} a_{j}
        k2(pa, pac, :) = k2(pa, pac, :) - dt * kappa * (DBLE(nja+nka) + DBLE(njac+nkac)) * &
                       & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))

        !--------------------!
        !   Cascade Decay    !
        !--------------------!
        ! \sqrt(2 epsilon gamma kappa) e^{i \phi_{j}} \sigma_{-} \rho a^{\dagger}_{j}
        ! Lower atom |\pm> and annihilate photon from mode jac/kac
        IF (njac > 0) THEN
          ! |1_{j},0_{k}> --> |0>, |2_{j},0_{k}> -> |1_{j},0_{k}>
          ! |1_{j},1_{k}> --> |0_{j},1_{k}>
          ! Annihilate a photon from mode jac
          pn = MAP_n2p(jac, kac, njac-1, nkac)
          ! |g><e| \rho
          k2(pa, pn, gg) = k2(pa, pn, gg) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg))
          k2(pa, pn, ge) = k2(pa, pn, ge) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
          k2(pa, pn, gf) = k2(pa, pn, gf) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
          ! \xi |e><f| \rho
          k2(pa, pn, eg) = k2(pa, pn, eg) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg))
          k2(pa, pn, ee) = k2(pa, pn, ee) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
          k2(pa, pn, ef) = k2(pa, pn, ef) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
          IF (nkac > 0) THEN
            ! |1_{j},1_{k}> --> |1_{j},0_{k}>
            ! Annihilate a photon from mode kac
            pn = MAP_n2p(jac, kac, njac, nkac-1)
            ! |g><e| \rho
            k2(pa, pn, gg) = k2(pa, pn, gg) + dt * gkal(kac) * &
                           & (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg))
            k2(pa, pn, ge) = k2(pa, pn, ge) + dt * gkal(kac) * &
                           & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
            k2(pa, pn, gf) = k2(pa, pn, gf) + dt * gkal(kac) * &
                           & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
            ! \xi |e><f| \rho
            k2(pa, pn, eg) = k2(pa, pn, eg) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg))
            k2(pa, pn, ee) = k2(pa, pn, ee) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
            k2(pa, pn, ef) = k2(pa, pn, ef) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
          END IF
        END IF

        ! \sqrt(2 epsilon gamma kappa) e^{-i \phi_{j}} a_{j} \rho \sigma_{+}
        ! Lower atom <\pm| and annihilate photon from mode ja/ka
        IF (nja > 0) THEN
          ! <1_{j}|<0_{k}| or <2_{j}|<0_{k}|
          ! Annihilate a photon from mode ja
          pn = MAP_n2p(ja, ka, nja-1, nka)
          ! \rho |e><g|
          k2(pn, pac, gg) = k2(pn, pac, gg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge))
          k2(pn, pac, eg) = k2(pn, pac, eg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
          k2(pn, pac, fg) = k2(pn, pac, fg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
          ! \xi \rho |f><e|
          k2(pn, pac, ge) = k2(pn, pac, ge) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf))
          k2(pn, pac, ee) = k2(pn, pac, ee) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
          k2(pn, pac, fe) = k2(pn, pac, fe) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
          IF (nka > 0) THEN
            ! |1_{j}>|1_{k}>
            ! Annihilate a photon from mode ka
            pn = MAP_n2p(ja, ka, nja, nka-1)
            ! \rho |e><g|
            k2(pn, pac, gg) = k2(pn, pac, gg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge))
            k2(pn, pac, eg) = k2(pn, pac, eg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
            k2(pn, pac, fg) = k2(pn, pac, fg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
            ! \xi \rho |f><e|
            k2(pn, pac, ge) = k2(pn, pac, ge) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf))
            k2(pn, pac, ee) = k2(pn, pac, ee) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
            k2(pn, pac, fe) = k2(pn, pac, fe) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
          END IF
        END IF

        ! Close p loop
      END DO
      ! Close pc loop
    END DO

    ! Calculate k3
    k3 = 0.0d0
    DO pac = 1, no_states
      DO pa = 1, no_states
        ! Get mode and photon number
        ja = MAP_p2n(pa, 1)
        jac = MAP_p2n(pac, 1)
        ka = MAP_p2n(pa, 2)
        kac = MAP_p2n(pac, 2)
        nja = MAP_p2n(pa, 3)
        njac = MAP_p2n(pac, 3)
        nka = MAP_p2n(pa, 4)
        nkac = MAP_p2n(pac, 4)

        ! Three-level equations of motion
        k3(pa, pac, gg) = k3(pa, pac, gg) + dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg)) + &
                        & dt * gamma * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
        k3(pa, pac, ge) = k3(pa, pac, ge) + dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + 0.5d0 * k2(pa, pac, gg)) - &
                        & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) + &
                        & dt * gamma * xi * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
        k3(pa, pac, gf) = k3(pa, pac, gf) + dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge)) - &
                        & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
        k3(pa, pac, eg) = k3(pa, pac, eg) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + 0.5d0 * k2(pa, pac, gg)) + &
                        & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg)) + &
                        & dt * gamma * xi * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
        k3(pa, pac, ee) = k3(pa, pac, ee) - dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg)) - &
                        & dt * gamma * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe)) + &
                        & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
        k3(pa, pac, ef) = k3(pa, pac, ef) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) + &
                        & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
        k3(pa, pac, fg) = k3(pa, pac, fg) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg)) + &
                        & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
        k3(pa, pac, fe) = k3(pa, pac, fe) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg)) - &
                        & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
        k3(pa, pac, ff) = k3(pa, pac, ff) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe)) - &
                        & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))

        ! Mode Hamiltonians
        k3(pa, pac, :) = k3(pa, pac, :) - i * dt * (&
                       & (wl(ja) * DBLE(nja)) - (wl(jac) * DBLE(njac)) + &
                       & (wl(ka) * DBLE(nka)) - (wl(kac) * DBLE(nkac))) * &
                       & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))

        !-----------------------!
        !  Commutator -i[H, p]  !
        !-----------------------!
        ! -i * H * rho
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{i \phi_{j}} a^{\dagger}_{j} \sigma_{-} \rho
        IF (nja < Fock .AND. nka == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == ja) THEN
              ! Add a photon to ja
              pn = MAP_n2p(ja, ka, nja+1, nka)
              ! |g><e| \rho
              k3(pn, pac, gg) = k3(pn, pac, gg) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg))
              k3(pn, pac, ge) = k3(pn, pac, ge) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
              k3(pn, pac, gf) = k3(pn, pac, gf) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
              ! \xi |e><f| \rho
              k3(pn, pac, eg) = k3(pn, pac, eg) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg))
              k3(pn, pac, ee) = k3(pn, pac, ee) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
              k3(pn, pac, ef) = k3(pn, pac, ef) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(ja, z, nja, nka+1)
              ! |g><e| \rho
              k3(pn, pac, gg) = k3(pn, pac, gg) - dt * gkal(z) * &
                              & (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg))
              k3(pn, pac, ge) = k3(pn, pac, ge) - dt * gkal(z) * &
                              & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
              k3(pn, pac, gf) = k3(pn, pac, gf) - dt * gkal(z) * &
                              & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
              ! \xi |e><f| \rho
              k3(pn, pac, eg) = k3(pn, pac, eg) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg))
              k3(pn, pac, ee) = k3(pn, pac, ee) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
              k3(pn, pac, ef) = k3(pn, pac, ef) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
            END IF
          END DO
        END IF

        ! i * rho * H
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{-i \phi_{j}} \rho \sigma_{+} a_{j}
        IF (njac < Fock .AND. nkac == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == jac) THEN
              ! Add a photon to mode ja
              pn = MAP_n2p(jac, kac, njac+1, nkac)
              ! \rho |e><g|
              k3(pa, pn, gg) = k3(pa, pn, gg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge))
              k3(pa, pn, eg) = k3(pa, pn, eg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
              k3(pa, pn, fg) = k3(pa, pn, fg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
              ! \xi \rho |f><e|
              k3(pa, pn, ge) = k3(pa, pn, ge) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf))
              k3(pa, pn, ee) = k3(pa, pn, ee) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
              k3(pa, pn, fe) = k3(pa, pn, fe) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(jac, z, njac, nkac+1)
              ! \rho |e><g|
              k3(pa, pn, gg) = k3(pa, pn, gg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge))
              k3(pa, pn, eg) = k3(pa, pn, eg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
              k3(pa, pn, fg) = k3(pa, pn, fg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
              ! \xi \rho |f><e|
              k3(pa, pn, ge) = k3(pa, pn, ge) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf))
              k3(pa, pn, ee) = k3(pa, pn, ee) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
              k3(pa, pn, fe) = k3(pa, pn, fe) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
            END IF
          END DO
        END IF

        !-----------------------------!
        !      LINDBLAD DECAY         !
        !-----------------------------!
        ! The atom fluorecence is split by a beam splitter into two paths:
        ! one goes towards the cavity at a fraction \epsilon, and the other
        ! goes to a separate detector at a fraction 1 - \epsilon.
        ! We split the cascaded systems decay operator:
        ! ja = \sqrt{\gamma} \sigma_{-} + \sqrt{\kappa} A,
        ! into separate components correponding to atomic decay, cavity decay
        ! and two cross terms. We also include a separate decay term to account
        ! for the detector placed outside the ring cavity.

        !--------------------!
        !    Cavity Decay    !
        !--------------------!
        ! Lindblad is (accounting for both cascaded and cavity decay)
        ! 2.0 * 0.5 \kappa (2 A \rho A^{\dagger} - A^{\dagger} A \rho
        !                                        - \rho A^{\dagger} A_),
        ! where A is a sum of all mode annihilation operators.

        ! \kappa * A * \rho * A^{\dagger} term
        ! For |nja>_{ja}|0>_{ka} <nja|_{ja'} <0|_{ka'}
        IF (nja > 0 .AND. njac > 0 .AND. ja == jac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k3(pn, pnc, :) = k3(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * sqrtl(njac) * &
                         & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))
        END IF
        IF (nja > 0 .AND. nkac > 0 .AND. ja == kac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k3(pn, pnc, :) = k3(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * &
                         & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))
        END IF
        IF (nka > 0 .AND. njac > 0 .AND. ka == jac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k3(pn, pnc, :) = k3(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(njac) * &
                         & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))
        END IF
        IF (nka > 0 .AND. nkac > 0 .AND. ka == kac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k3(pn, pnc, :) = k3(pn, pnc, :) + dt * 2.0d0 * kappa * &
                         & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))
        END IF

        ! -a^{\dagger}_{j} a_{j} \rho - \rho a^{\dagger}_{j} a_{j}
        k3(pa, pac, :) = k3(pa, pac, :) - dt * kappa * (DBLE(nja+nka) + DBLE(njac+nkac)) * &
                       & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))

        !--------------------!
        !   Cascade Decay    !
        !--------------------!
        ! \sqrt(2 epsilon gamma kappa) e^{i \phi_{j}} \sigma_{-} \rho a^{\dagger}_{j}
        ! Lower atom |\pm> and annihilate photon from mode jac/kac
        IF (njac > 0) THEN
          ! |1_{j},0_{k}> --> |0>, |2_{j},0_{k}> -> |1_{j},0_{k}>
          ! |1_{j},1_{k}> --> |0_{j},1_{k}>
          ! Annihilate a photon from mode jac
          pn = MAP_n2p(jac, kac, njac-1, nkac)
          ! |g><e| \rho
          k3(pa, pn, gg) = k3(pa, pn, gg) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg))
          k3(pa, pn, ge) = k3(pa, pn, ge) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
          k3(pa, pn, gf) = k3(pa, pn, gf) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
          ! \xi |e><f| \rho
          k3(pa, pn, eg) = k3(pa, pn, eg) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg))
          k3(pa, pn, ee) = k3(pa, pn, ee) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
          k3(pa, pn, ef) = k3(pa, pn, ef) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
          IF (nkac > 0) THEN
            ! |1_{j},1_{k}> --> |1_{j},0_{k}>
            ! Annihilate a photon from mode kac
            pn = MAP_n2p(jac, kac, njac, nkac-1)
            ! |g><e| \rho
            k3(pa, pn, gg) = k3(pa, pn, gg) + dt * gkal(kac) * &
                           & (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg))
            k3(pa, pn, ge) = k3(pa, pn, ge) + dt * gkal(kac) * &
                           & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
            k3(pa, pn, gf) = k3(pa, pn, gf) + dt * gkal(kac) * &
                           & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
            ! \xi |e><f| \rho
            k3(pa, pn, eg) = k3(pa, pn, eg) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg))
            k3(pa, pn, ee) = k3(pa, pn, ee) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
            k3(pa, pn, ef) = k3(pa, pn, ef) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
          END IF
        END IF

        ! \sqrt(2 epsilon gamma kappa) e^{-i \phi_{j}} a_{j} \rho \sigma_{+}
        ! Lower atom <\pm| and annihilate photon from mode ja/ka
        IF (nja > 0) THEN
          ! <1_{j}|<0_{k}| or <2_{j}|<0_{k}|
          ! Annihilate a photon from mode ja
          pn = MAP_n2p(ja, ka, nja-1, nka)
          ! \rho |e><g|
          k3(pn, pac, gg) = k3(pn, pac, gg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge))
          k3(pn, pac, eg) = k3(pn, pac, eg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
          k3(pn, pac, fg) = k3(pn, pac, fg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
          ! \xi \rho |f><e|
          k3(pn, pac, ge) = k3(pn, pac, ge) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf))
          k3(pn, pac, ee) = k3(pn, pac, ee) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
          k3(pn, pac, fe) = k3(pn, pac, fe) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
          IF (nka > 0) THEN
            ! |1_{j}>|1_{k}>
            ! Annihilate a photon from mode ka
            pn = MAP_n2p(ja, ka, nja, nka-1)
            ! \rho |e><g|
            k3(pn, pac, gg) = k3(pn, pac, gg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge))
            k3(pn, pac, eg) = k3(pn, pac, eg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
            k3(pn, pac, fg) = k3(pn, pac, fg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
            ! \xi \rho |f><e|
            k3(pn, pac, ge) = k3(pn, pac, ge) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf))
            k3(pn, pac, ee) = k3(pn, pac, ee) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
            k3(pn, pac, fe) = k3(pn, pac, fe) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
          END IF
        END IF

        ! Close p loop
      END DO
      ! Close pc loop
    END DO

    ! Calculate k4
    k4 = 0.0d0
    DO pac = 1, no_states
      DO pa = 1, no_states
        ! Get mode and photon number
        ja = MAP_p2n(pa, 1)
        jac = MAP_p2n(pac, 1)
        ka = MAP_p2n(pa, 2)
        kac = MAP_p2n(pac, 2)
        nja = MAP_p2n(pa, 3)
        njac = MAP_p2n(pac, 3)
        nka = MAP_p2n(pa, 4)
        nkac = MAP_p2n(pac, 4)

        ! Three-level equations of motion
        k4(pa, pac, gg) = k4(pa, pac, gg) + dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + k3(pa, pac, ge)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + k3(pa, pac, eg)) + &
                        & dt * gamma * (rho(pa, pac, ee) + k3(pa, pac, ee))
        k4(pa, pac, ge) = k4(pa, pac, ge) + dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + k3(pa, pac, gg)) - &
                        & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(pa, pac, ge) + k3(pa, pac, ge)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, gf) + k3(pa, pac, gf)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + k3(pa, pac, ee)) + &
                        & dt * gamma * xi * (rho(pa, pac, ef) + k3(pa, pac, ef))
        k4(pa, pac, gf) = k4(pa, pac, gf) + dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ge) + k3(pa, pac, ge)) - &
                        & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, gf) + k3(pa, pac, gf)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ef) + k3(pa, pac, ef))
        k4(pa, pac, eg) = k4(pa, pac, eg) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + k3(pa, pac, gg)) + &
                        & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(pa, pac, eg) + k3(pa, pac, eg)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + k3(pa, pac, ee)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fg) + k3(pa, pac, fg)) + &
                        & dt * gamma * xi * (rho(pa, pac, fe) + k3(pa, pac, fe))
        k4(pa, pac, ee) = k4(pa, pac, ee) - dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + k3(pa, pac, ge)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + k3(pa, pac, eg)) - &
                        & dt * gamma * (rho(pa, pac, ee) + k3(pa, pac, ee)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + k3(pa, pac, ef)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + k3(pa, pac, fe)) + &
                        & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + k3(pa, pac, ff))
        k4(pa, pac, ef) = k4(pa, pac, ef) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gf) + k3(pa, pac, gf)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + k3(pa, pac, ee)) + &
                        & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, ef) + k3(pa, pac, ef)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + k3(pa, pac, ff))
        k4(pa, pac, fg) = k4(pa, pac, fg) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, eg) + k3(pa, pac, eg)) + &
                        & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, fg) + k3(pa, pac, fg)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, fe) + k3(pa, pac, fe))
        k4(pa, pac, fe) = k4(pa, pac, fe) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + k3(pa, pac, ee)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, fg) + k3(pa, pac, fg)) - &
                        & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, fe) + k3(pa, pac, fe)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + k3(pa, pac, ff))
        k4(pa, pac, ff) = k4(pa, pac, ff) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + k3(pa, pac, ef)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + k3(pa, pac, fe)) - &
                        & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + k3(pa, pac, ff))

        ! Mode Hamiltonians
        k4(pa, pac, :) = k4(pa, pac, :) - i * dt * (&
                       & (wl(ja) * DBLE(nja)) - (wl(jac) * DBLE(njac)) + &
                       & (wl(ka) * DBLE(nka)) - (wl(kac) * DBLE(nkac))) * &
                       & (rho(pa, pac, :) + k3(pa, pac, :))

        !-----------------------!
        !  Commutator -i[H, p]  !
        !-----------------------!
        ! -i * H * rho
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{i \phi_{j}} a^{\dagger}_{j} \sigma_{-} \rho
        IF (nja < Fock .AND. nka == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == ja) THEN
              ! Add a photon to ja
              pn = MAP_n2p(ja, ka, nja+1, nka)
              ! |g><e| \rho
              k4(pn, pac, gg) = k4(pn, pac, gg) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, eg) + k3(pa, pac, eg))
              k4(pn, pac, ge) = k4(pn, pac, ge) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, ee) + k3(pa, pac, ee))
              k4(pn, pac, gf) = k4(pn, pac, gf) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, ef) + k3(pa, pac, ef))
              ! \xi |e><f| \rho
              k4(pn, pac, eg) = k4(pn, pac, eg) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, fg) + k3(pa, pac, fg))
              k4(pn, pac, ee) = k4(pn, pac, ee) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, fe) + k3(pa, pac, fe))
              k4(pn, pac, ef) = k4(pn, pac, ef) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, ff) + k3(pa, pac, ff))
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(ja, z, nja, nka+1)
              ! |g><e| \rho
              k4(pn, pac, gg) = k4(pn, pac, gg) - dt * gkal(z) * &
                              & (rho(pa, pac, eg) + k3(pa, pac, eg))
              k4(pn, pac, ge) = k4(pn, pac, ge) - dt * gkal(z) * &
                              & (rho(pa, pac, ee) + k3(pa, pac, ee))
              k4(pn, pac, gf) = k4(pn, pac, gf) - dt * gkal(z) * &
                              & (rho(pa, pac, ef) + k3(pa, pac, ef))
              ! \xi |e><f| \rho
              k4(pn, pac, eg) = k4(pn, pac, eg) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, fg) + k3(pa, pac, fg))
              k4(pn, pac, ee) = k4(pn, pac, ee) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, fe) + k3(pa, pac, fe))
              k4(pn, pac, ef) = k4(pn, pac, ef) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, ff) + k3(pa, pac, ff))
            END IF
          END DO
        END IF

        ! i * rho * H
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{-i \phi_{j}} \rho \sigma_{+} a_{j}
        IF (njac < Fock .AND. nkac == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == jac) THEN
              ! Add a photon to mode ja
              pn = MAP_n2p(jac, kac, njac+1, nkac)
              ! \rho |e><g|
              k4(pa, pn, gg) = k4(pa, pn, gg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, ge) + k3(pa, pac, ge))
              k4(pa, pn, eg) = k4(pa, pn, eg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, ee) + k3(pa, pac, ee))
              k4(pa, pn, fg) = k4(pa, pn, fg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, fe) + k3(pa, pac, fe))
              ! \xi \rho |f><e|
              k4(pa, pn, ge) = k4(pa, pn, ge) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, gf) + k3(pa, pac, gf))
              k4(pa, pn, ee) = k4(pa, pn, ee) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, ef) + k3(pa, pac, ef))
              k4(pa, pn, fe) = k4(pa, pn, fe) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, ff) + k3(pa, pac, ff))
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(jac, z, njac, nkac+1)
              ! \rho |e><g|
              k4(pa, pn, gg) = k4(pa, pn, gg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, ge) + k3(pa, pac, ge))
              k4(pa, pn, eg) = k4(pa, pn, eg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, ee) + k3(pa, pac, ee))
              k4(pa, pn, fg) = k4(pa, pn, fg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, fe) + k3(pa, pac, fe))
              ! \xi \rho |f><e|
              k4(pa, pn, ge) = k4(pa, pn, ge) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, gf) + k3(pa, pac, gf))
              k4(pa, pn, ee) = k4(pa, pn, ee) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, ef) + k3(pa, pac, ef))
              k4(pa, pn, fe) = k4(pa, pn, fe) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, ff) + k3(pa, pac, ff))
            END IF
          END DO
        END IF

        !-----------------------------!
        !      LINDBLAD DECAY         !
        !-----------------------------!
        ! The atom fluorecence is split by a beam splitter into two paths:
        ! one goes towards the cavity at a fraction \epsilon, and the other
        ! goes to a separate detector at a fraction 1 - \epsilon.
        ! We split the cascaded systems decay operator:
        ! ja = \sqrt{\gamma} \sigma_{-} + \sqrt{\kappa} A,
        ! into separate components correponding to atomic decay, cavity decay
        ! and two cross terms. We also include a separate decay term to account
        ! for the detector placed outside the ring cavity.

        !--------------------!
        !    Cavity Decay    !
        !--------------------!
        ! Lindblad is (accounting for both cascaded and cavity decay)
        ! 2.0 * 0.5 \kappa (2 A \rho A^{\dagger} - A^{\dagger} A \rho
        !                                        - \rho A^{\dagger} A_),
        ! where A is a sum of all mode annihilation operators.

        ! \kappa * A * \rho * A^{\dagger} term
        ! For |nja>_{ja}|0>_{ka} <nja|_{ja'} <0|_{ka'}
        IF (nja > 0 .AND. njac > 0 .AND. ja == jac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k4(pn, pnc, :) = k4(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * sqrtl(njac) * &
                         & (rho(pa, pac, :) + k3(pa, pac, :))
        END IF
        IF (nja > 0 .AND. nkac > 0 .AND. ja == kac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k4(pn, pnc, :) = k4(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * &
                         & (rho(pa, pac, :) + k3(pa, pac, :))
        END IF
        IF (nka > 0 .AND. njac > 0 .AND. ka == jac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k4(pn, pnc, :) = k4(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(njac) * &
                         & (rho(pa, pac, :) + k3(pa, pac, :))
        END IF
        IF (nka > 0 .AND. nkac > 0 .AND. ka == kac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k4(pn, pnc, :) = k4(pn, pnc, :) + dt * 2.0d0 * kappa * &
                         & (rho(pa, pac, :) + k3(pa, pac, :))
        END IF

        ! -a^{\dagger}_{j} a_{j} \rho - \rho a^{\dagger}_{j} a_{j}
        k4(pa, pac, :) = k4(pa, pac, :) - dt * kappa * (DBLE(nja+nka) + DBLE(njac+nkac)) * &
                       & (rho(pa, pac, :) + k3(pa, pac, :))

        !--------------------!
        !   Cascade Decay    !
        !--------------------!
        ! \sqrt(2 epsilon gamma kappa) e^{i \phi_{j}} \sigma_{-} \rho a^{\dagger}_{j}
        ! Lower atom |\pm> and annihilate photon from mode jac/kac
        IF (njac > 0) THEN
          ! |1_{j},0_{k}> --> |0>, |2_{j},0_{k}> -> |1_{j},0_{k}>
          ! |1_{j},1_{k}> --> |0_{j},1_{k}>
          ! Annihilate a photon from mode jac
          pn = MAP_n2p(jac, kac, njac-1, nkac)
          ! |g><e| \rho
          k4(pa, pn, gg) = k4(pa, pn, gg) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, eg) + k3(pa, pac, eg))
          k4(pa, pn, ge) = k4(pa, pn, ge) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, ee) + k3(pa, pac, ee))
          k4(pa, pn, gf) = k4(pa, pn, gf) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, ef) + k3(pa, pac, ef))
          ! \xi |e><f| \rho
          k4(pa, pn, eg) = k4(pa, pn, eg) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, fg) + k3(pa, pac, fg))
          k4(pa, pn, ee) = k4(pa, pn, ee) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, fe) + k3(pa, pac, fe))
          k4(pa, pn, ef) = k4(pa, pn, ef) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, ff) + k3(pa, pac, ff))
          IF (nkac > 0) THEN
            ! |1_{j},1_{k}> --> |1_{j},0_{k}>
            ! Annihilate a photon from mode kac
            pn = MAP_n2p(jac, kac, njac, nkac-1)
            ! |g><e| \rho
            k4(pa, pn, gg) = k4(pa, pn, gg) + dt * gkal(kac) * &
                           & (rho(pa, pac, eg) + k3(pa, pac, eg))
            k4(pa, pn, ge) = k4(pa, pn, ge) + dt * gkal(kac) * &
                           & (rho(pa, pac, ee) + k3(pa, pac, ee))
            k4(pa, pn, gf) = k4(pa, pn, gf) + dt * gkal(kac) * &
                           & (rho(pa, pac, ef) + k3(pa, pac, ef))
            ! \xi |e><f| \rho
            k4(pa, pn, eg) = k4(pa, pn, eg) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, fg) + k3(pa, pac, fg))
            k4(pa, pn, ee) = k4(pa, pn, ee) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, fe) + k3(pa, pac, fe))
            k4(pa, pn, ef) = k4(pa, pn, ef) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, ff) + k3(pa, pac, ff))
          END IF
        END IF

        ! \sqrt(2 epsilon gamma kappa) e^{-i \phi_{j}} a_{j} \rho \sigma_{+}
        ! Lower atom <\pm| and annihilate photon from mode ja/ka
        IF (nja > 0) THEN
          ! <1_{j}|<0_{k}| or <2_{j}|<0_{k}|
          ! Annihilate a photon from mode ja
          pn = MAP_n2p(ja, ka, nja-1, nka)
          ! \rho |e><g|
          k4(pn, pac, gg) = k4(pn, pac, gg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, ge) + k3(pa, pac, ge))
          k4(pn, pac, eg) = k4(pn, pac, eg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, ee) + k3(pa, pac, ee))
          k4(pn, pac, fg) = k4(pn, pac, fg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, fe) + k3(pa, pac, fe))
          ! \xi \rho |f><e|
          k4(pn, pac, ge) = k4(pn, pac, ge) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, gf) + k3(pa, pac, gf))
          k4(pn, pac, ee) = k4(pn, pac, ee) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, ef) + k3(pa, pac, ef))
          k4(pn, pac, fe) = k4(pn, pac, fe) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, ff) + k3(pa, pac, ff))
          IF (nka > 0) THEN
            ! |1_{j}>|1_{k}>
            ! Annihilate a photon from mode ka
            pn = MAP_n2p(ja, ka, nja, nka-1)
            ! \rho |e><g|
            k4(pn, pac, gg) = k4(pn, pac, gg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, ge) + k3(pa, pac, ge))
            k4(pn, pac, eg) = k4(pn, pac, eg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, ee) + k3(pa, pac, ee))
            k4(pn, pac, fg) = k4(pn, pac, fg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, fe) + k3(pa, pac, fe))
            ! \xi \rho |f><e|
            k4(pn, pac, ge) = k4(pn, pac, ge) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, gf) + k3(pa, pac, gf))
            k4(pn, pac, ee) = k4(pn, pac, ee) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, ef) + k3(pa, pac, ef))
            k4(pn, pac, fe) = k4(pn, pac, fe) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, ff) + k3(pa, pac, ff))
          END IF
        END IF

        ! Close p loop
      END DO
      ! Close pc loop
    END DO

    ! Update rho
    rho = rho + xis * (k1 + 2.0d0 * (k2 + k3) + k4)

    ! Check percentage
    IF (ten_percent /= 0) THEN
      IF (MOD(t, ten_percent) == 0 .AND. t /= 0) THEN
        CALL CPU_TIME(loop_check_time)
        percentage = NINT((100.0 * t) / (1.0 * tau_steps))
        loop_run_time = loop_check_time - loop_start_time
        loop_remaining_time = ((100.0 * (loop_check_time - loop_start_time)) / (1.0 * percentage)) - loop_run_time
        WRITE(*, FMT_corr) percentage, "%. Run time (g1 correlation): ", loop_run_time,&
                         & "s. Est. time left: ", loop_remaining_time, "s"
      END IF
    END IF

    ! Close time integration DO loop
  END DO

  ! Close file
  CLOSE(3)
END IF

!##############################################################################!
!                       SECOND-ORDER CORRELATION FUNCTION                      !
!##############################################################################!
IF (tau2_max > 0.0) THEN
  ! Propagate A \rho(ss) * A^{\dagger} (remove a photon from |ket> and <bra|)
  rho = 0.0
  DO pac = 1, no_states_two_photon
    DO pa = 1, no_states_two_photon
      ja = MAP_p2n(pa, 1)
      jac = MAP_p2n(pac, 1)
      ka = MAP_p2n(pa, 2)
      kac = MAP_p2n(pac, 2)
      nja = MAP_p2n(pa, 3)
      njac = MAP_p2n(pac, 3)
      nka = MAP_p2n(pa, 4)
      nkac = MAP_p2n(pac, 4)
      IF (nja > 0 .AND. njac > 0) THEN
        pn = MAP_n2p(ja, ka, nja-1, nka)
        pnc = MAP_n2p(jac, kac, njac-1, nkac)
        rho(pn, pnc, :) = rho(pn, pnc, :) + sqrtl(nja) * sqrtl(njac) * rho_ss(pa, pac, :)
      END IF
      IF (nja > 0 .AND. nkac > 0) THEN
        pn = MAP_n2p(ja, ka, nja-1, nka)
        pnc = MAP_n2p(jac, kac, njac, nkac-1)
        rho(pn, pnc, :) = rho(pn, pnc, :) + sqrtl(nja) * &
                        & rho_ss(pa, pac, :)
      END IF
      IF (nka > 0 .AND. njac > 0) THEN
        pn = MAP_n2p(ja, ka, nja, nka-1)
        pnc = MAP_n2p(jac, kac, njac-1, nkac)
        rho(pn, pnc, :) = rho(pn, pnc, :) + sqrtl(njac) * &
                        & rho_ss(pa, pac, :)
      END IF
      IF (nka > 0 .AND. nkac > 0) THEN
        pn = MAP_n2p(ja, ka, nja, nka-1)
        pnc = MAP_n2p(jac, kac, njac, nkac-1)
        rho(pn, pnc, :) = rho(pn, pnc, :) + rho_ss(pa, pac, :)
      END IF
    END DO
  END DO

  PRINT*, " "
  PRINT*, "Calculating second-order correlation ..."
  PRINT*, " "
  ! Change tau_steps for g2
  tau_steps = NINT(tau2_max / dt)

  ! Open file to write time and correlation to
  OPEN(UNIT=4, FILE=filename_g2, STATUS='REPLACE', ACTION='WRITE', RECL=4000)

  ! Ten percent of time steps
  ten_percent = NINT((1.0 * tau_steps / 10.0))
  ! Call CPU clock time
  CALL CPU_TIME(loop_start_time)

  ! Change number of states to single_photon
  no_states = no_states_one_photon
  Fock = 1

  ! Start time integration to find second-order correlation
  DO t = 0, tau_steps
    ! Calculate correlation
    corr = 0.0d0
    DO pac = 1, no_states_one_photon
      DO pa = 1, no_states_one_photon
        ja = MAP_p2n(pa, 1)
        jac = MAP_p2n(pac, 1)
        ka = MAP_p2n(pa, 2)
        kac = MAP_p2n(pac, 2)
        nja = MAP_p2n(pa, 3)
        njac = MAP_p2n(pac, 3)
        nka = MAP_p2n(pa, 4)
        nkac = MAP_p2n(pac, 4)

        ! To calculate < A^{\dagger} A > = <p' | A^{\dagger} A | p>, remove a
        ! photon from either side the, if the indice are the same, update the
        ! photona variable.

        IF (nja == 1 .AND. njac == 1 .AND. nka == 0 .AND. nkac == 0) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          IF (pn == pnc) THEN
            corr = corr + sqrtl(nja) * sqrtl(njac) * REAL(rho(pa, pac, gg) + rho(pa, pac, ee) + rho(pa, pac, ff))
          END IF
        END IF
      END DO
    END DO

    ! Normalise by steady-state mean-photon number
    IF (photona_ss /= 0.0) THEN
      corr = REAL(corr) / (photona_ss ** 2)
    END IF

    ! Write data
    WRITE(4,*) DBLE(t) * dt, REAL(corr)

    ! Calculate k1
    k1 = 0.0d0
    DO pac = 1, no_states
      DO pa = 1, no_states
        ! Get mode and photon number
        ja = MAP_p2n(pa, 1)
        jac = MAP_p2n(pac, 1)
        ka = MAP_p2n(pa, 2)
        kac = MAP_p2n(pac, 2)
        nja = MAP_p2n(pa, 3)
        njac = MAP_p2n(pac, 3)
        nka = MAP_p2n(pa, 4)
        nkac = MAP_p2n(pac, 4)

        ! Three-level equations of motion
        k1(pa, pac, gg) = k1(pa, pac, gg) + dt * i * 0.5d0 * Omega * rho(pa, pac, ge) - &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, eg) + &
                        & dt * gamma * rho(pa, pac, ee)
        k1(pa, pac, ge) = k1(pa, pac, ge) + dt * i * 0.5d0 * Omega * rho(pa, pac, gg) - &
                        & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * rho(pa, pac, ge) + &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, gf) - &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, ee) + &
                        & dt * gamma * xi * rho(pa, pac, ef)
        k1(pa, pac, gf) = k1(pa, pac, gf) + dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ge) - &
                        & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * rho(pa, pac, gf) - &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, ef)
        k1(pa, pac, eg) = k1(pa, pac, eg) - dt * i * 0.5d0 * Omega * rho(pa, pac, gg) + &
                        & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * rho(pa, pac, eg) + &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, ee) - &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, fg) + &
                        & dt * gamma * xi * rho(pa, pac, fe)
        k1(pa, pac, ee) = k1(pa, pac, ee) - dt * i * 0.5d0 * Omega * rho(pa, pac, ge) + &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, eg) - &
                        & dt * gamma * rho(pa, pac, ee) + &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ef) - &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, fe) + &
                        & dt * gamma * (xi ** 2) * rho(pa, pac, ff)
        k1(pa, pac, ef) = k1(pa, pac, ef) - dt * i * 0.5d0 * Omega * rho(pa, pac, gf) + &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ee) + &
                        & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * rho(pa, pac, ef) - &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ff)
        k1(pa, pac, fg) = k1(pa, pac, fg) - dt * i * 0.5d0 * xi * Omega * rho(pa, pac, eg) + &
                        & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * rho(pa, pac, fg) + &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, fe)
        k1(pa, pac, fe) = k1(pa, pac, fe) - dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ee) + &
                        & dt * i * 0.5d0 * Omega * rho(pa, pac, fg) - &
                        & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * rho(pa, pac, fe) + &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ff)
        k1(pa, pac, ff) = k1(pa, pac, ff) - dt * i * 0.5d0 * xi * Omega * rho(pa, pac, ef) + &
                        & dt * i * 0.5d0 * xi * Omega * rho(pa, pac, fe) - &
                        & dt * gamma * (xi ** 2) * rho(pa, pac, ff)

        ! Mode Hamiltonians
        k1(pa, pac, :) = k1(pa, pac, :) - i * dt * (&
                       & (wl(ja) * DBLE(nja)) - (wl(jac) * DBLE(njac)) + &
                       & (wl(ka) * DBLE(nka)) - (wl(kac) * DBLE(nkac))) * &
                       & rho(pa, pac, :)

        !-----------------------!
        !  Commutator -i[H, p]  !
        !-----------------------!
        ! -i * H * rho
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{i \phi_{j}} a^{\dagger}_{j} \sigma_{-} \rho
        IF (nja < Fock .AND. nka == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == ja) THEN
              ! Add a photon to ja
              pn = MAP_n2p(ja, ka, nja+1, nka)
              ! |g><e| \rho
              k1(pn, pac, gg) = k1(pn, pac, gg) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & rho(pa, pac, eg)
              k1(pn, pac, ge) = k1(pn, pac, ge) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & rho(pa, pac, ee)
              k1(pn, pac, gf) = k1(pn, pac, gf) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & rho(pa, pac, ef)
              ! \xi |e><f| \rho
              k1(pn, pac, eg) = k1(pn, pac, eg) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & rho(pa, pac, fg)
              k1(pn, pac, ee) = k1(pn, pac, ee) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & rho(pa, pac, fe)
              k1(pn, pac, ef) = k1(pn, pac, ef) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & rho(pa, pac, ff)
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(ja, z, nja, nka+1)
              ! |g><e| \rho
              k1(pn, pac, gg) = k1(pn, pac, gg) - dt * gkal(z) * &
                              & rho(pa, pac, eg)
              k1(pn, pac, ge) = k1(pn, pac, ge) - dt * gkal(z) * &
                              & rho(pa, pac, ee)
              k1(pn, pac, gf) = k1(pn, pac, gf) - dt * gkal(z) * &
                              & rho(pa, pac, ef)
              ! \xi |e><f| \rho
              k1(pn, pac, eg) = k1(pn, pac, eg) - dt * gkal(z) * xi * &
                              & rho(pa, pac, fg)
              k1(pn, pac, ee) = k1(pn, pac, ee) - dt * gkal(z) * xi * &
                              & rho(pa, pac, fe)
              k1(pn, pac, ef) = k1(pn, pac, ef) - dt * gkal(z) * xi * &
                              & rho(pa, pac, ff)
            END IF
          END DO
        END IF

        ! i * rho * H
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{-i \phi_{j}} \rho \sigma_{+} a_{j}
        IF (njac < Fock .AND. nkac == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == jac) THEN
              ! Add a photon to mode ja
              pn = MAP_n2p(jac, kac, njac+1, nkac)
              ! \rho |e><g|
              k1(pa, pn, gg) = k1(pa, pn, gg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & rho(pa, pac, ge)
              k1(pa, pn, eg) = k1(pa, pn, eg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & rho(pa, pac, ee)
              k1(pa, pn, fg) = k1(pa, pn, fg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & rho(pa, pac, fe)
              ! \xi \rho |f><e|
              k1(pa, pn, ge) = k1(pa, pn, ge) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & rho(pa, pac, gf)
              k1(pa, pn, ee) = k1(pa, pn, ee) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & rho(pa, pac, ef)
              k1(pa, pn, fe) = k1(pa, pn, fe) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & rho(pa, pac, ff)
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(jac, z, njac, nkac+1)
              ! \rho |e><g|
              k1(pa, pn, gg) = k1(pa, pn, gg) - dt * CONJG(gkal(z)) * &
                             & rho(pa, pac, ge)
              k1(pa, pn, eg) = k1(pa, pn, eg) - dt * CONJG(gkal(z)) * &
                             & rho(pa, pac, ee)
              k1(pa, pn, fg) = k1(pa, pn, fg) - dt * CONJG(gkal(z)) * &
                             & rho(pa, pac, fe)
              ! \xi \rho |f><e|
              k1(pa, pn, ge) = k1(pa, pn, ge) - dt * CONJG(gkal(z)) * xi * &
                             & rho(pa, pac, gf)
              k1(pa, pn, ee) = k1(pa, pn, ee) - dt * CONJG(gkal(z)) * xi * &
                             & rho(pa, pac, ef)
              k1(pa, pn, fe) = k1(pa, pn, fe) - dt * CONJG(gkal(z)) * xi * &
                             & rho(pa, pac, ff)
            END IF
          END DO
        END IF

        !-----------------------------!
        !      LINDBLAD DECAY         !
        !-----------------------------!
        ! The atom fluorecence is split by a beam splitter into two paths:
        ! one goes towards the cavity at a fraction \epsilon, and the other
        ! goes to a separate detector at a fraction 1 - \epsilon.
        ! We split the cascaded systems decay operator:
        ! ja = \sqrt{\gamma} \sigma_{-} + \sqrt{\kappa} A,
        ! into separate components correponding to atomic decay, cavity decay
        ! and two cross terms. We also include a separate decay term to account
        ! for the detector placed outside the ring cavity.

        !--------------------!
        !    Cavity Decay    !
        !--------------------!
        ! Lindblad is (accounting for both cascaded and cavity decay)
        ! 2.0 * 0.5 \kappa (2 A \rho A^{\dagger} - A^{\dagger} A \rho
        !                                        - \rho A^{\dagger} A_),
        ! where A is a sum of all mode annihilation operators.

        ! \kappa * A * \rho * A^{\dagger} term
        ! For |nja>_{ja}|0>_{ka} <nja|_{ja'} <0|_{ka'}
        IF (nja > 0 .AND. njac > 0 .AND. ja == jac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k1(pn, pnc, :) = k1(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * sqrtl(njac) * &
                         & rho(pa, pac, :)
        END IF
        IF (nja > 0 .AND. nkac > 0 .AND. ja == kac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k1(pn, pnc, :) = k1(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * &
                         & rho(pa, pac, :)
        END IF
        IF (nka > 0 .AND. njac > 0 .AND. ka == jac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k1(pn, pnc, :) = k1(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(njac) * &
                         & rho(pa, pac, :)
        END IF
        IF (nka > 0 .AND. nkac > 0 .AND. ka == kac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k1(pn, pnc, :) = k1(pn, pnc, :) + dt * 2.0d0 * kappa * &
                         & rho(pa, pac, :)
        END IF

        ! -a^{\dagger}_{j} a_{j} \rho - \rho a^{\dagger}_{j} a_{j}
        k1(pa, pac, :) = k1(pa, pac, :) - dt * kappa * (DBLE(nja+nka) + DBLE(njac+nkac)) * &
                       & rho(pa, pac, :)

        !--------------------!
        !   Cascade Decay    !
        !--------------------!
        ! \sqrt(2 epsilon gamma kappa) e^{i \phi_{j}} \sigma_{-} \rho a^{\dagger}_{j}
        ! Lower atom |\pm> and annihilate photon from mode jac/kac
        IF (njac > 0) THEN
          ! |1_{j},0_{k}> --> |0>, |2_{j},0_{k}> -> |1_{j},0_{k}>
          ! |1_{j},1_{k}> --> |0_{j},1_{k}>
          ! Annihilate a photon from mode jac
          pn = MAP_n2p(jac, kac, njac-1, nkac)
          ! |g><e| \rho
          k1(pa, pn, gg) = k1(pa, pn, gg) + dt * gkal(jac) * sqrtl(njac) * &
                         & rho(pa, pac, eg)
          k1(pa, pn, ge) = k1(pa, pn, ge) + dt * gkal(jac) * sqrtl(njac) * &
                         & rho(pa, pac, ee)
          k1(pa, pn, gf) = k1(pa, pn, gf) + dt * gkal(jac) * sqrtl(njac) * &
                         & rho(pa, pac, ef)
          ! \xi |e><f| \rho
          k1(pa, pn, eg) = k1(pa, pn, eg) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & rho(pa, pac, fg)
          k1(pa, pn, ee) = k1(pa, pn, ee) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & rho(pa, pac, fe)
          k1(pa, pn, ef) = k1(pa, pn, ef) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & rho(pa, pac, ff)
          IF (nkac > 0) THEN
            ! |1_{j},1_{k}> --> |1_{j},0_{k}>
            ! Annihilate a photon from mode kac
            pn = MAP_n2p(jac, kac, njac, nkac-1)
            ! |g><e| \rho
            k1(pa, pn, gg) = k1(pa, pn, gg) + dt * gkal(kac) * &
                           & rho(pa, pac, eg)
            k1(pa, pn, ge) = k1(pa, pn, ge) + dt * gkal(kac) * &
                           & rho(pa, pac, ee)
            k1(pa, pn, gf) = k1(pa, pn, gf) + dt * gkal(kac) * &
                           & rho(pa, pac, ef)
            ! \xi |e><f| \rho
            k1(pa, pn, eg) = k1(pa, pn, eg) + dt * gkal(kac) * xi * &
                           & rho(pa, pac, fg)
            k1(pa, pn, ee) = k1(pa, pn, ee) + dt * gkal(kac) * xi * &
                           & rho(pa, pac, fe)
            k1(pa, pn, ef) = k1(pa, pn, ef) + dt * gkal(kac) * xi * &
                           & rho(pa, pac, ff)
          END IF
        END IF

        ! \sqrt(2 epsilon gamma kappa) e^{-i \phi_{j}} a_{j} \rho \sigma_{+}
        ! Lower atom <\pm| and annihilate photon from mode ja/ka
        IF (nja > 0) THEN
          ! <1_{j}|<0_{k}| or <2_{j}|<0_{k}|
          ! Annihilate a photon from mode ja
          pn = MAP_n2p(ja, ka, nja-1, nka)
          ! \rho |e><g|
          k1(pn, pac, gg) = k1(pn, pac, gg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & rho(pa, pac, ge)
          k1(pn, pac, eg) = k1(pn, pac, eg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & rho(pa, pac, ee)
          k1(pn, pac, fg) = k1(pn, pac, fg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & rho(pa, pac, fe)
          ! \xi \rho |f><e|
          k1(pn, pac, ge) = k1(pn, pac, ge) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & rho(pa, pac, gf)
          k1(pn, pac, ee) = k1(pn, pac, ee) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & rho(pa, pac, ef)
          k1(pn, pac, fe) = k1(pn, pac, fe) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & rho(pa, pac, ff)
          IF (nka > 0) THEN
            ! |1_{j}>|1_{k}>
            ! Annihilate a photon from mode ka
            pn = MAP_n2p(ja, ka, nja, nka-1)
            ! \rho |e><g|
            k1(pn, pac, gg) = k1(pn, pac, gg) + dt * CONJG(gkal(ka)) * &
                            & rho(pa, pac, ge)
            k1(pn, pac, eg) = k1(pn, pac, eg) + dt * CONJG(gkal(ka)) * &
                            & rho(pa, pac, ee)
            k1(pn, pac, fg) = k1(pn, pac, fg) + dt * CONJG(gkal(ka)) * &
                            & rho(pa, pac, fe)
            ! \xi \rho |f><e|
            k1(pn, pac, ge) = k1(pn, pac, ge) + dt * CONJG(gkal(ka)) * xi * &
                            & rho(pa, pac, gf)
            k1(pn, pac, ee) = k1(pn, pac, ee) + dt * CONJG(gkal(ka)) * xi * &
                            & rho(pa, pac, ef)
            k1(pn, pac, fe) = k1(pn, pac, fe) + dt * CONJG(gkal(ka)) * xi * &
                            & rho(pa, pac, ff)
          END IF
        END IF

        ! Close p loop
      END DO
      ! Close pc loop
    END DO

    ! Calculate k2
    k2 = 0.0d0
    DO pac = 1, no_states
      DO pa = 1, no_states
        ! Get mode and photon number
        ja = MAP_p2n(pa, 1)
        jac = MAP_p2n(pac, 1)
        ka = MAP_p2n(pa, 2)
        kac = MAP_p2n(pac, 2)
        nja = MAP_p2n(pa, 3)
        njac = MAP_p2n(pac, 3)
        nka = MAP_p2n(pa, 4)
        nkac = MAP_p2n(pac, 4)

        ! Three-level equations of motion
        k2(pa, pac, gg) = k2(pa, pac, gg) + dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg)) + &
                        & dt * gamma * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
        k2(pa, pac, ge) = k2(pa, pac, ge) + dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + 0.5d0 * k1(pa, pac, gg)) - &
                        & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) + &
                        & dt * gamma * xi * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
        k2(pa, pac, gf) = k2(pa, pac, gf) + dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge)) - &
                        & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
        k2(pa, pac, eg) = k2(pa, pac, eg) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + 0.5d0 * k1(pa, pac, gg)) + &
                        & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg)) + &
                        & dt * gamma * xi * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
        k2(pa, pac, ee) = k2(pa, pac, ee) - dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg)) - &
                        & dt * gamma * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe)) + &
                        & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
        k2(pa, pac, ef) = k2(pa, pac, ef) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) + &
                        & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
        k2(pa, pac, fg) = k2(pa, pac, fg) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg)) + &
                        & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
        k2(pa, pac, fe) = k2(pa, pac, fe) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg)) - &
                        & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
        k2(pa, pac, ff) = k2(pa, pac, ff) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe)) - &
                        & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))

        ! Mode Hamiltonians
        k2(pa, pac, :) = k2(pa, pac, :) - i * dt * (&
                       & (wl(ja) * DBLE(nja)) - (wl(jac) * DBLE(njac)) + &
                       & (wl(ka) * DBLE(nka)) - (wl(kac) * DBLE(nkac))) * &
                       & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))

        !-----------------------!
        !  Commutator -i[H, p]  !
        !-----------------------!
        ! -i * H * rho
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{i \phi_{j}} a^{\dagger}_{j} \sigma_{-} \rho
        IF (nja < Fock .AND. nka == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == ja) THEN
              ! Add a photon to ja
              pn = MAP_n2p(ja, ka, nja+1, nka)
              ! |g><e| \rho
              k2(pn, pac, gg) = k2(pn, pac, gg) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg))
              k2(pn, pac, ge) = k2(pn, pac, ge) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
              k2(pn, pac, gf) = k2(pn, pac, gf) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
              ! \xi |e><f| \rho
              k2(pn, pac, eg) = k2(pn, pac, eg) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg))
              k2(pn, pac, ee) = k2(pn, pac, ee) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
              k2(pn, pac, ef) = k2(pn, pac, ef) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(ja, z, nja, nka+1)
              ! |g><e| \rho
              k2(pn, pac, gg) = k2(pn, pac, gg) - dt * gkal(z) * &
                              & (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg))
              k2(pn, pac, ge) = k2(pn, pac, ge) - dt * gkal(z) * &
                              & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
              k2(pn, pac, gf) = k2(pn, pac, gf) - dt * gkal(z) * &
                              & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
              ! \xi |e><f| \rho
              k2(pn, pac, eg) = k2(pn, pac, eg) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg))
              k2(pn, pac, ee) = k2(pn, pac, ee) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
              k2(pn, pac, ef) = k2(pn, pac, ef) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
            END IF
          END DO
        END IF

        ! i * rho * H
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{-i \phi_{j}} \rho \sigma_{+} a_{j}
        IF (njac < Fock .AND. nkac == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == jac) THEN
              ! Add a photon to mode ja
              pn = MAP_n2p(jac, kac, njac+1, nkac)
              ! \rho |e><g|
              k2(pa, pn, gg) = k2(pa, pn, gg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge))
              k2(pa, pn, eg) = k2(pa, pn, eg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
              k2(pa, pn, fg) = k2(pa, pn, fg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
              ! \xi \rho |f><e|
              k2(pa, pn, ge) = k2(pa, pn, ge) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf))
              k2(pa, pn, ee) = k2(pa, pn, ee) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
              k2(pa, pn, fe) = k2(pa, pn, fe) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(jac, z, njac, nkac+1)
              ! \rho |e><g|
              k2(pa, pn, gg) = k2(pa, pn, gg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge))
              k2(pa, pn, eg) = k2(pa, pn, eg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
              k2(pa, pn, fg) = k2(pa, pn, fg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
              ! \xi \rho |f><e|
              k2(pa, pn, ge) = k2(pa, pn, ge) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf))
              k2(pa, pn, ee) = k2(pa, pn, ee) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
              k2(pa, pn, fe) = k2(pa, pn, fe) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
            END IF
          END DO
        END IF

        !-----------------------------!
        !      LINDBLAD DECAY         !
        !-----------------------------!
        ! The atom fluorecence is split by a beam splitter into two paths:
        ! one goes towards the cavity at a fraction \epsilon, and the other
        ! goes to a separate detector at a fraction 1 - \epsilon.
        ! We split the cascaded systems decay operator:
        ! ja = \sqrt{\gamma} \sigma_{-} + \sqrt{\kappa} A,
        ! into separate components correponding to atomic decay, cavity decay
        ! and two cross terms. We also include a separate decay term to account
        ! for the detector placed outside the ring cavity.

        !--------------------!
        !    Cavity Decay    !
        !--------------------!
        ! Lindblad is (accounting for both cascaded and cavity decay)
        ! 2.0 * 0.5 \kappa (2 A \rho A^{\dagger} - A^{\dagger} A \rho
        !                                        - \rho A^{\dagger} A_),
        ! where A is a sum of all mode annihilation operators.

        ! \kappa * A * \rho * A^{\dagger} term
        ! For |nja>_{ja}|0>_{ka} <nja|_{ja'} <0|_{ka'}
        IF (nja > 0 .AND. njac > 0 .AND. ja == jac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k2(pn, pnc, :) = k2(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * sqrtl(njac) * &
                         & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))
        END IF
        IF (nja > 0 .AND. nkac > 0 .AND. ja == kac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k2(pn, pnc, :) = k2(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * &
                         & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))
        END IF
        IF (nka > 0 .AND. njac > 0 .AND. ka == jac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k2(pn, pnc, :) = k2(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(njac) * &
                         & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))
        END IF
        IF (nka > 0 .AND. nkac > 0 .AND. ka == kac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k2(pn, pnc, :) = k2(pn, pnc, :) + dt * 2.0d0 * kappa * &
                         & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))
        END IF

        ! -a^{\dagger}_{j} a_{j} \rho - \rho a^{\dagger}_{j} a_{j}
        k2(pa, pac, :) = k2(pa, pac, :) - dt * kappa * (DBLE(nja+nka) + DBLE(njac+nkac)) * &
                       & (rho(pa, pac, :) + 0.5d0 * k1(pa, pac, :))

        !--------------------!
        !   Cascade Decay    !
        !--------------------!
        ! \sqrt(2 epsilon gamma kappa) e^{i \phi_{j}} \sigma_{-} \rho a^{\dagger}_{j}
        ! Lower atom |\pm> and annihilate photon from mode jac/kac
        IF (njac > 0) THEN
          ! |1_{j},0_{k}> --> |0>, |2_{j},0_{k}> -> |1_{j},0_{k}>
          ! |1_{j},1_{k}> --> |0_{j},1_{k}>
          ! Annihilate a photon from mode jac
          pn = MAP_n2p(jac, kac, njac-1, nkac)
          ! |g><e| \rho
          k2(pa, pn, gg) = k2(pa, pn, gg) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg))
          k2(pa, pn, ge) = k2(pa, pn, ge) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
          k2(pa, pn, gf) = k2(pa, pn, gf) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
          ! \xi |e><f| \rho
          k2(pa, pn, eg) = k2(pa, pn, eg) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg))
          k2(pa, pn, ee) = k2(pa, pn, ee) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
          k2(pa, pn, ef) = k2(pa, pn, ef) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
          IF (nkac > 0) THEN
            ! |1_{j},1_{k}> --> |1_{j},0_{k}>
            ! Annihilate a photon from mode kac
            pn = MAP_n2p(jac, kac, njac, nkac-1)
            ! |g><e| \rho
            k2(pa, pn, gg) = k2(pa, pn, gg) + dt * gkal(kac) * &
                           & (rho(pa, pac, eg) + 0.5d0 * k1(pa, pac, eg))
            k2(pa, pn, ge) = k2(pa, pn, ge) + dt * gkal(kac) * &
                           & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
            k2(pa, pn, gf) = k2(pa, pn, gf) + dt * gkal(kac) * &
                           & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
            ! \xi |e><f| \rho
            k2(pa, pn, eg) = k2(pa, pn, eg) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, fg) + 0.5d0 * k1(pa, pac, fg))
            k2(pa, pn, ee) = k2(pa, pn, ee) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
            k2(pa, pn, ef) = k2(pa, pn, ef) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
          END IF
        END IF

        ! \sqrt(2 epsilon gamma kappa) e^{-i \phi_{j}} a_{j} \rho \sigma_{+}
        ! Lower atom <\pm| and annihilate photon from mode ja/ka
        IF (nja > 0) THEN
          ! <1_{j}|<0_{k}| or <2_{j}|<0_{k}|
          ! Annihilate a photon from mode ja
          pn = MAP_n2p(ja, ka, nja-1, nka)
          ! \rho |e><g|
          k2(pn, pac, gg) = k2(pn, pac, gg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge))
          k2(pn, pac, eg) = k2(pn, pac, eg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
          k2(pn, pac, fg) = k2(pn, pac, fg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
          ! \xi \rho |f><e|
          k2(pn, pac, ge) = k2(pn, pac, ge) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf))
          k2(pn, pac, ee) = k2(pn, pac, ee) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
          k2(pn, pac, fe) = k2(pn, pac, fe) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
          IF (nka > 0) THEN
            ! |1_{j}>|1_{k}>
            ! Annihilate a photon from mode ka
            pn = MAP_n2p(ja, ka, nja, nka-1)
            ! \rho |e><g|
            k2(pn, pac, gg) = k2(pn, pac, gg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, ge) + 0.5d0 * k1(pa, pac, ge))
            k2(pn, pac, eg) = k2(pn, pac, eg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, ee) + 0.5d0 * k1(pa, pac, ee))
            k2(pn, pac, fg) = k2(pn, pac, fg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, fe) + 0.5d0 * k1(pa, pac, fe))
            ! \xi \rho |f><e|
            k2(pn, pac, ge) = k2(pn, pac, ge) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, gf) + 0.5d0 * k1(pa, pac, gf))
            k2(pn, pac, ee) = k2(pn, pac, ee) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, ef) + 0.5d0 * k1(pa, pac, ef))
            k2(pn, pac, fe) = k2(pn, pac, fe) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, ff) + 0.5d0 * k1(pa, pac, ff))
          END IF
        END IF

        ! Close p loop
      END DO
      ! Close pc loop
    END DO

    ! Calculate k3
    k3 = 0.0d0
    DO pac = 1, no_states
      DO pa = 1, no_states
        ! Get mode and photon number
        ja = MAP_p2n(pa, 1)
        jac = MAP_p2n(pac, 1)
        ka = MAP_p2n(pa, 2)
        kac = MAP_p2n(pac, 2)
        nja = MAP_p2n(pa, 3)
        njac = MAP_p2n(pac, 3)
        nka = MAP_p2n(pa, 4)
        nkac = MAP_p2n(pac, 4)

        ! Three-level equations of motion
        k3(pa, pac, gg) = k3(pa, pac, gg) + dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg)) + &
                        & dt * gamma * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
        k3(pa, pac, ge) = k3(pa, pac, ge) + dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + 0.5d0 * k2(pa, pac, gg)) - &
                        & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) + &
                        & dt * gamma * xi * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
        k3(pa, pac, gf) = k3(pa, pac, gf) + dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge)) - &
                        & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
        k3(pa, pac, eg) = k3(pa, pac, eg) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + 0.5d0 * k2(pa, pac, gg)) + &
                        & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg)) + &
                        & dt * gamma * xi * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
        k3(pa, pac, ee) = k3(pa, pac, ee) - dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg)) - &
                        & dt * gamma * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe)) + &
                        & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
        k3(pa, pac, ef) = k3(pa, pac, ef) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) + &
                        & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
        k3(pa, pac, fg) = k3(pa, pac, fg) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg)) + &
                        & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
        k3(pa, pac, fe) = k3(pa, pac, fe) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg)) - &
                        & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
        k3(pa, pac, ff) = k3(pa, pac, ff) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe)) - &
                        & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))

        ! Mode Hamiltonians
        k3(pa, pac, :) = k3(pa, pac, :) - i * dt * (&
                       & (wl(ja) * DBLE(nja)) - (wl(jac) * DBLE(njac)) + &
                       & (wl(ka) * DBLE(nka)) - (wl(kac) * DBLE(nkac))) * &
                       & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))

        !-----------------------!
        !  Commutator -i[H, p]  !
        !-----------------------!
        ! -i * H * rho
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{i \phi_{j}} a^{\dagger}_{j} \sigma_{-} \rho
        IF (nja < Fock .AND. nka == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == ja) THEN
              ! Add a photon to ja
              pn = MAP_n2p(ja, ka, nja+1, nka)
              ! |g><e| \rho
              k3(pn, pac, gg) = k3(pn, pac, gg) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg))
              k3(pn, pac, ge) = k3(pn, pac, ge) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
              k3(pn, pac, gf) = k3(pn, pac, gf) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
              ! \xi |e><f| \rho
              k3(pn, pac, eg) = k3(pn, pac, eg) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg))
              k3(pn, pac, ee) = k3(pn, pac, ee) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
              k3(pn, pac, ef) = k3(pn, pac, ef) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(ja, z, nja, nka+1)
              ! |g><e| \rho
              k3(pn, pac, gg) = k3(pn, pac, gg) - dt * gkal(z) * &
                              & (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg))
              k3(pn, pac, ge) = k3(pn, pac, ge) - dt * gkal(z) * &
                              & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
              k3(pn, pac, gf) = k3(pn, pac, gf) - dt * gkal(z) * &
                              & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
              ! \xi |e><f| \rho
              k3(pn, pac, eg) = k3(pn, pac, eg) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg))
              k3(pn, pac, ee) = k3(pn, pac, ee) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
              k3(pn, pac, ef) = k3(pn, pac, ef) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
            END IF
          END DO
        END IF

        ! i * rho * H
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{-i \phi_{j}} \rho \sigma_{+} a_{j}
        IF (njac < Fock .AND. nkac == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == jac) THEN
              ! Add a photon to mode ja
              pn = MAP_n2p(jac, kac, njac+1, nkac)
              ! \rho |e><g|
              k3(pa, pn, gg) = k3(pa, pn, gg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge))
              k3(pa, pn, eg) = k3(pa, pn, eg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
              k3(pa, pn, fg) = k3(pa, pn, fg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
              ! \xi \rho |f><e|
              k3(pa, pn, ge) = k3(pa, pn, ge) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf))
              k3(pa, pn, ee) = k3(pa, pn, ee) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
              k3(pa, pn, fe) = k3(pa, pn, fe) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(jac, z, njac, nkac+1)
              ! \rho |e><g|
              k3(pa, pn, gg) = k3(pa, pn, gg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge))
              k3(pa, pn, eg) = k3(pa, pn, eg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
              k3(pa, pn, fg) = k3(pa, pn, fg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
              ! \xi \rho |f><e|
              k3(pa, pn, ge) = k3(pa, pn, ge) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf))
              k3(pa, pn, ee) = k3(pa, pn, ee) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
              k3(pa, pn, fe) = k3(pa, pn, fe) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
            END IF
          END DO
        END IF

        !-----------------------------!
        !      LINDBLAD DECAY         !
        !-----------------------------!
        ! The atom fluorecence is split by a beam splitter into two paths:
        ! one goes towards the cavity at a fraction \epsilon, and the other
        ! goes to a separate detector at a fraction 1 - \epsilon.
        ! We split the cascaded systems decay operator:
        ! ja = \sqrt{\gamma} \sigma_{-} + \sqrt{\kappa} A,
        ! into separate components correponding to atomic decay, cavity decay
        ! and two cross terms. We also include a separate decay term to account
        ! for the detector placed outside the ring cavity.

        !--------------------!
        !    Cavity Decay    !
        !--------------------!
        ! Lindblad is (accounting for both cascaded and cavity decay)
        ! 2.0 * 0.5 \kappa (2 A \rho A^{\dagger} - A^{\dagger} A \rho
        !                                        - \rho A^{\dagger} A_),
        ! where A is a sum of all mode annihilation operators.

        ! \kappa * A * \rho * A^{\dagger} term
        ! For |nja>_{ja}|0>_{ka} <nja|_{ja'} <0|_{ka'}
        IF (nja > 0 .AND. njac > 0 .AND. ja == jac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k3(pn, pnc, :) = k3(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * sqrtl(njac) * &
                         & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))
        END IF
        IF (nja > 0 .AND. nkac > 0 .AND. ja == kac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k3(pn, pnc, :) = k3(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * &
                         & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))
        END IF
        IF (nka > 0 .AND. njac > 0 .AND. ka == jac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k3(pn, pnc, :) = k3(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(njac) * &
                         & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))
        END IF
        IF (nka > 0 .AND. nkac > 0 .AND. ka == kac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k3(pn, pnc, :) = k3(pn, pnc, :) + dt * 2.0d0 * kappa * &
                         & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))
        END IF

        ! -a^{\dagger}_{j} a_{j} \rho - \rho a^{\dagger}_{j} a_{j}
        k3(pa, pac, :) = k3(pa, pac, :) - dt * kappa * (DBLE(nja+nka) + DBLE(njac+nkac)) * &
                       & (rho(pa, pac, :) + 0.5d0 * k2(pa, pac, :))

        !--------------------!
        !   Cascade Decay    !
        !--------------------!
        ! \sqrt(2 epsilon gamma kappa) e^{i \phi_{j}} \sigma_{-} \rho a^{\dagger}_{j}
        ! Lower atom |\pm> and annihilate photon from mode jac/kac
        IF (njac > 0) THEN
          ! |1_{j},0_{k}> --> |0>, |2_{j},0_{k}> -> |1_{j},0_{k}>
          ! |1_{j},1_{k}> --> |0_{j},1_{k}>
          ! Annihilate a photon from mode jac
          pn = MAP_n2p(jac, kac, njac-1, nkac)
          ! |g><e| \rho
          k3(pa, pn, gg) = k3(pa, pn, gg) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg))
          k3(pa, pn, ge) = k3(pa, pn, ge) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
          k3(pa, pn, gf) = k3(pa, pn, gf) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
          ! \xi |e><f| \rho
          k3(pa, pn, eg) = k3(pa, pn, eg) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg))
          k3(pa, pn, ee) = k3(pa, pn, ee) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
          k3(pa, pn, ef) = k3(pa, pn, ef) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
          IF (nkac > 0) THEN
            ! |1_{j},1_{k}> --> |1_{j},0_{k}>
            ! Annihilate a photon from mode kac
            pn = MAP_n2p(jac, kac, njac, nkac-1)
            ! |g><e| \rho
            k3(pa, pn, gg) = k3(pa, pn, gg) + dt * gkal(kac) * &
                           & (rho(pa, pac, eg) + 0.5d0 * k2(pa, pac, eg))
            k3(pa, pn, ge) = k3(pa, pn, ge) + dt * gkal(kac) * &
                           & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
            k3(pa, pn, gf) = k3(pa, pn, gf) + dt * gkal(kac) * &
                           & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
            ! \xi |e><f| \rho
            k3(pa, pn, eg) = k3(pa, pn, eg) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, fg) + 0.5d0 * k2(pa, pac, fg))
            k3(pa, pn, ee) = k3(pa, pn, ee) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
            k3(pa, pn, ef) = k3(pa, pn, ef) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
          END IF
        END IF

        ! \sqrt(2 epsilon gamma kappa) e^{-i \phi_{j}} a_{j} \rho \sigma_{+}
        ! Lower atom <\pm| and annihilate photon from mode ja/ka
        IF (nja > 0) THEN
          ! <1_{j}|<0_{k}| or <2_{j}|<0_{k}|
          ! Annihilate a photon from mode ja
          pn = MAP_n2p(ja, ka, nja-1, nka)
          ! \rho |e><g|
          k3(pn, pac, gg) = k3(pn, pac, gg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge))
          k3(pn, pac, eg) = k3(pn, pac, eg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
          k3(pn, pac, fg) = k3(pn, pac, fg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
          ! \xi \rho |f><e|
          k3(pn, pac, ge) = k3(pn, pac, ge) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf))
          k3(pn, pac, ee) = k3(pn, pac, ee) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
          k3(pn, pac, fe) = k3(pn, pac, fe) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
          IF (nka > 0) THEN
            ! |1_{j}>|1_{k}>
            ! Annihilate a photon from mode ka
            pn = MAP_n2p(ja, ka, nja, nka-1)
            ! \rho |e><g|
            k3(pn, pac, gg) = k3(pn, pac, gg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, ge) + 0.5d0 * k2(pa, pac, ge))
            k3(pn, pac, eg) = k3(pn, pac, eg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, ee) + 0.5d0 * k2(pa, pac, ee))
            k3(pn, pac, fg) = k3(pn, pac, fg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, fe) + 0.5d0 * k2(pa, pac, fe))
            ! \xi \rho |f><e|
            k3(pn, pac, ge) = k3(pn, pac, ge) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, gf) + 0.5d0 * k2(pa, pac, gf))
            k3(pn, pac, ee) = k3(pn, pac, ee) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, ef) + 0.5d0 * k2(pa, pac, ef))
            k3(pn, pac, fe) = k3(pn, pac, fe) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, ff) + 0.5d0 * k2(pa, pac, ff))
          END IF
        END IF

        ! Close p loop
      END DO
      ! Close pc loop
    END DO

    ! Calculate k4
    k4 = 0.0d0
    DO pac = 1, no_states
      DO pa = 1, no_states
        ! Get mode and photon number
        ja = MAP_p2n(pa, 1)
        jac = MAP_p2n(pac, 1)
        ka = MAP_p2n(pa, 2)
        kac = MAP_p2n(pac, 2)
        nja = MAP_p2n(pa, 3)
        njac = MAP_p2n(pac, 3)
        nka = MAP_p2n(pa, 4)
        nkac = MAP_p2n(pac, 4)

        ! Three-level equations of motion
        k4(pa, pac, gg) = k4(pa, pac, gg) + dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + k3(pa, pac, ge)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + k3(pa, pac, eg)) + &
                        & dt * gamma * (rho(pa, pac, ee) + k3(pa, pac, ee))
        k4(pa, pac, ge) = k4(pa, pac, ge) + dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + k3(pa, pac, gg)) - &
                        & dt * (i * ((0.5d0 * alpha) + delta) + 0.5d0 * gamma) * (rho(pa, pac, ge) + k3(pa, pac, ge)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, gf) + k3(pa, pac, gf)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + k3(pa, pac, ee)) + &
                        & dt * gamma * xi * (rho(pa, pac, ef) + k3(pa, pac, ef))
        k4(pa, pac, gf) = k4(pa, pac, gf) + dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ge) + k3(pa, pac, ge)) - &
                        & dt * ((2.0d0 * i * delta) + 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, gf) + k3(pa, pac, gf)) - &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ef) + k3(pa, pac, ef))
        k4(pa, pac, eg) = k4(pa, pac, eg) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gg) + k3(pa, pac, gg)) + &
                        & dt * (i * ((0.5d0 * alpha) + delta) - 0.5d0 * gamma) * (rho(pa, pac, eg) + k3(pa, pac, eg)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, ee) + k3(pa, pac, ee)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fg) + k3(pa, pac, fg)) + &
                        & dt * gamma * xi * (rho(pa, pac, fe) + k3(pa, pac, fe))
        k4(pa, pac, ee) = k4(pa, pac, ee) - dt * i * 0.5d0 * Omega * (rho(pa, pac, ge) + k3(pa, pac, ge)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, eg) + k3(pa, pac, eg)) - &
                        & dt * gamma * (rho(pa, pac, ee) + k3(pa, pac, ee)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + k3(pa, pac, ef)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + k3(pa, pac, fe)) + &
                        & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + k3(pa, pac, ff))
        k4(pa, pac, ef) = k4(pa, pac, ef) - dt * i * 0.5d0 * Omega * (rho(pa, pac, gf) + k3(pa, pac, gf)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + k3(pa, pac, ee)) + &
                        & dt * (i * ((0.5d0 * alpha) - delta) - 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, ef) + k3(pa, pac, ef)) - &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + k3(pa, pac, ff))
        k4(pa, pac, fg) = k4(pa, pac, fg) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, eg) + k3(pa, pac, eg)) + &
                        & dt * (2.0d0 * i * delta - 0.5d0 * gamma * (xi ** 2)) * (rho(pa, pac, fg) + k3(pa, pac, fg)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, fe) + k3(pa, pac, fe))
        k4(pa, pac, fe) = k4(pa, pac, fe) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ee) + k3(pa, pac, ee)) + &
                        & dt * i * 0.5d0 * Omega * (rho(pa, pac, fg) + k3(pa, pac, fg)) - &
                        & dt * (i * ((0.5d0 * alpha) - delta) + 0.5d0 * gamma * (1.0 + (xi ** 2))) * (rho(pa, pac, fe) + k3(pa, pac, fe)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ff) + k3(pa, pac, ff))
        k4(pa, pac, ff) = k4(pa, pac, ff) - dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, ef) + k3(pa, pac, ef)) + &
                        & dt * i * 0.5d0 * xi * Omega * (rho(pa, pac, fe) + k3(pa, pac, fe)) - &
                        & dt * gamma * (xi ** 2) * (rho(pa, pac, ff) + k3(pa, pac, ff))

        ! Mode Hamiltonians
        k4(pa, pac, :) = k4(pa, pac, :) - i * dt * (&
                       & (wl(ja) * DBLE(nja)) - (wl(jac) * DBLE(njac)) + &
                       & (wl(ka) * DBLE(nka)) - (wl(kac) * DBLE(nkac))) * &
                       & (rho(pa, pac, :) + k3(pa, pac, :))

        !-----------------------!
        !  Commutator -i[H, p]  !
        !-----------------------!
        ! -i * H * rho
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{i \phi_{j}} a^{\dagger}_{j} \sigma_{-} \rho
        IF (nja < Fock .AND. nka == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == ja) THEN
              ! Add a photon to ja
              pn = MAP_n2p(ja, ka, nja+1, nka)
              ! |g><e| \rho
              k4(pn, pac, gg) = k4(pn, pac, gg) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, eg) + k3(pa, pac, eg))
              k4(pn, pac, ge) = k4(pn, pac, ge) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, ee) + k3(pa, pac, ee))
              k4(pn, pac, gf) = k4(pn, pac, gf) - dt * gkal(ja) * sqrtl(nja+1) * &
                              & (rho(pa, pac, ef) + k3(pa, pac, ef))
              ! \xi |e><f| \rho
              k4(pn, pac, eg) = k4(pn, pac, eg) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, fg) + k3(pa, pac, fg))
              k4(pn, pac, ee) = k4(pn, pac, ee) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, fe) + k3(pa, pac, fe))
              k4(pn, pac, ef) = k4(pn, pac, ef) - dt * gkal(ja) * xi * sqrtl(nja+1) * &
                              & (rho(pa, pac, ff) + k3(pa, pac, ff))
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(ja, z, nja, nka+1)
              ! |g><e| \rho
              k4(pn, pac, gg) = k4(pn, pac, gg) - dt * gkal(z) * &
                              & (rho(pa, pac, eg) + k3(pa, pac, eg))
              k4(pn, pac, ge) = k4(pn, pac, ge) - dt * gkal(z) * &
                              & (rho(pa, pac, ee) + k3(pa, pac, ee))
              k4(pn, pac, gf) = k4(pn, pac, gf) - dt * gkal(z) * &
                              & (rho(pa, pac, ef) + k3(pa, pac, ef))
              ! \xi |e><f| \rho
              k4(pn, pac, eg) = k4(pn, pac, eg) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, fg) + k3(pa, pac, fg))
              k4(pn, pac, ee) = k4(pn, pac, ee) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, fe) + k3(pa, pac, fe))
              k4(pn, pac, ef) = k4(pn, pac, ef) - dt * gkal(z) * xi * &
                              & (rho(pa, pac, ff) + k3(pa, pac, ff))
            END IF
          END DO
        END IF

        ! i * rho * H
        ! Cascade Hamiltonian
        ! - SQRT(2 gamma kappa) \sum_{j=-N}^{N} e^{-i \phi_{j}} \rho \sigma_{+} a_{j}
        IF (njac < Fock .AND. nkac == 0) THEN
          ! Cycle through modes and create a photon
          DO z = -N, N
            IF (z == jac) THEN
              ! Add a photon to mode ja
              pn = MAP_n2p(jac, kac, njac+1, nkac)
              ! \rho |e><g|
              k4(pa, pn, gg) = k4(pa, pn, gg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, ge) + k3(pa, pac, ge))
              k4(pa, pn, eg) = k4(pa, pn, eg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, ee) + k3(pa, pac, ee))
              k4(pa, pn, fg) = k4(pa, pn, fg) - dt * CONJG(gkal(jac)) * sqrtl(njac+1) * &
                             & (rho(pa, pac, fe) + k3(pa, pac, fe))
              ! \xi \rho |f><e|
              k4(pa, pn, ge) = k4(pa, pn, ge) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, gf) + k3(pa, pac, gf))
              k4(pa, pn, ee) = k4(pa, pn, ee) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, ef) + k3(pa, pac, ef))
              k4(pa, pn, fe) = k4(pa, pn, fe) - dt * CONJG(gkal(jac)) * xi * sqrtl(njac+1) * &
                             & (rho(pa, pac, ff) + k3(pa, pac, ff))
            ELSE
              ! Add a photon to other modes
              pn = MAP_n2p(jac, z, njac, nkac+1)
              ! \rho |e><g|
              k4(pa, pn, gg) = k4(pa, pn, gg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, ge) + k3(pa, pac, ge))
              k4(pa, pn, eg) = k4(pa, pn, eg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, ee) + k3(pa, pac, ee))
              k4(pa, pn, fg) = k4(pa, pn, fg) - dt * CONJG(gkal(z)) * &
                             & (rho(pa, pac, fe) + k3(pa, pac, fe))
              ! \xi \rho |f><e|
              k4(pa, pn, ge) = k4(pa, pn, ge) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, gf) + k3(pa, pac, gf))
              k4(pa, pn, ee) = k4(pa, pn, ee) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, ef) + k3(pa, pac, ef))
              k4(pa, pn, fe) = k4(pa, pn, fe) - dt * CONJG(gkal(z)) * xi * &
                             & (rho(pa, pac, ff) + k3(pa, pac, ff))
            END IF
          END DO
        END IF

        !-----------------------------!
        !      LINDBLAD DECAY         !
        !-----------------------------!
        ! The atom fluorecence is split by a beam splitter into two paths:
        ! one goes towards the cavity at a fraction \epsilon, and the other
        ! goes to a separate detector at a fraction 1 - \epsilon.
        ! We split the cascaded systems decay operator:
        ! ja = \sqrt{\gamma} \sigma_{-} + \sqrt{\kappa} A,
        ! into separate components correponding to atomic decay, cavity decay
        ! and two cross terms. We also include a separate decay term to account
        ! for the detector placed outside the ring cavity.

        !--------------------!
        !    Cavity Decay    !
        !--------------------!
        ! Lindblad is (accounting for both cascaded and cavity decay)
        ! 2.0 * 0.5 \kappa (2 A \rho A^{\dagger} - A^{\dagger} A \rho
        !                                        - \rho A^{\dagger} A_),
        ! where A is a sum of all mode annihilation operators.

        ! \kappa * A * \rho * A^{\dagger} term
        ! For |nja>_{ja}|0>_{ka} <nja|_{ja'} <0|_{ka'}
        IF (nja > 0 .AND. njac > 0 .AND. ja == jac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k4(pn, pnc, :) = k4(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * sqrtl(njac) * &
                         & (rho(pa, pac, :) + k3(pa, pac, :))
        END IF
        IF (nja > 0 .AND. nkac > 0 .AND. ja == kac) THEN
          pn = MAP_n2p(ja, ka, nja-1, nka)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k4(pn, pnc, :) = k4(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(nja) * &
                         & (rho(pa, pac, :) + k3(pa, pac, :))
        END IF
        IF (nka > 0 .AND. njac > 0 .AND. ka == jac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac-1, nkac)
          k4(pn, pnc, :) = k4(pn, pnc, :) + dt * 2.0d0 * kappa * sqrtl(njac) * &
                         & (rho(pa, pac, :) + k3(pa, pac, :))
        END IF
        IF (nka > 0 .AND. nkac > 0 .AND. ka == kac) THEN
          pn = MAP_n2p(ja, ka, nja, nka-1)
          pnc = MAP_n2p(jac, kac, njac, nkac-1)
          k4(pn, pnc, :) = k4(pn, pnc, :) + dt * 2.0d0 * kappa * &
                         & (rho(pa, pac, :) + k3(pa, pac, :))
        END IF

        ! -a^{\dagger}_{j} a_{j} \rho - \rho a^{\dagger}_{j} a_{j}
        k4(pa, pac, :) = k4(pa, pac, :) - dt * kappa * (DBLE(nja+nka) + DBLE(njac+nkac)) * &
                       & (rho(pa, pac, :) + k3(pa, pac, :))

        !--------------------!
        !   Cascade Decay    !
        !--------------------!
        ! \sqrt(2 epsilon gamma kappa) e^{i \phi_{j}} \sigma_{-} \rho a^{\dagger}_{j}
        ! Lower atom |\pm> and annihilate photon from mode jac/kac
        IF (njac > 0) THEN
          ! |1_{j},0_{k}> --> |0>, |2_{j},0_{k}> -> |1_{j},0_{k}>
          ! |1_{j},1_{k}> --> |0_{j},1_{k}>
          ! Annihilate a photon from mode jac
          pn = MAP_n2p(jac, kac, njac-1, nkac)
          ! |g><e| \rho
          k4(pa, pn, gg) = k4(pa, pn, gg) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, eg) + k3(pa, pac, eg))
          k4(pa, pn, ge) = k4(pa, pn, ge) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, ee) + k3(pa, pac, ee))
          k4(pa, pn, gf) = k4(pa, pn, gf) + dt * gkal(jac) * sqrtl(njac) * &
                         & (rho(pa, pac, ef) + k3(pa, pac, ef))
          ! \xi |e><f| \rho
          k4(pa, pn, eg) = k4(pa, pn, eg) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, fg) + k3(pa, pac, fg))
          k4(pa, pn, ee) = k4(pa, pn, ee) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, fe) + k3(pa, pac, fe))
          k4(pa, pn, ef) = k4(pa, pn, ef) + dt * gkal(jac) * xi * sqrtl(njac) * &
                         & (rho(pa, pac, ff) + k3(pa, pac, ff))
          IF (nkac > 0) THEN
            ! |1_{j},1_{k}> --> |1_{j},0_{k}>
            ! Annihilate a photon from mode kac
            pn = MAP_n2p(jac, kac, njac, nkac-1)
            ! |g><e| \rho
            k4(pa, pn, gg) = k4(pa, pn, gg) + dt * gkal(kac) * &
                           & (rho(pa, pac, eg) + k3(pa, pac, eg))
            k4(pa, pn, ge) = k4(pa, pn, ge) + dt * gkal(kac) * &
                           & (rho(pa, pac, ee) + k3(pa, pac, ee))
            k4(pa, pn, gf) = k4(pa, pn, gf) + dt * gkal(kac) * &
                           & (rho(pa, pac, ef) + k3(pa, pac, ef))
            ! \xi |e><f| \rho
            k4(pa, pn, eg) = k4(pa, pn, eg) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, fg) + k3(pa, pac, fg))
            k4(pa, pn, ee) = k4(pa, pn, ee) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, fe) + k3(pa, pac, fe))
            k4(pa, pn, ef) = k4(pa, pn, ef) + dt * gkal(kac) * xi * &
                           & (rho(pa, pac, ff) + k3(pa, pac, ff))
          END IF
        END IF

        ! \sqrt(2 epsilon gamma kappa) e^{-i \phi_{j}} a_{j} \rho \sigma_{+}
        ! Lower atom <\pm| and annihilate photon from mode ja/ka
        IF (nja > 0) THEN
          ! <1_{j}|<0_{k}| or <2_{j}|<0_{k}|
          ! Annihilate a photon from mode ja
          pn = MAP_n2p(ja, ka, nja-1, nka)
          ! \rho |e><g|
          k4(pn, pac, gg) = k4(pn, pac, gg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, ge) + k3(pa, pac, ge))
          k4(pn, pac, eg) = k4(pn, pac, eg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, ee) + k3(pa, pac, ee))
          k4(pn, pac, fg) = k4(pn, pac, fg) + dt * CONJG(gkal(ja)) * sqrtl(nja) * &
                          & (rho(pa, pac, fe) + k3(pa, pac, fe))
          ! \xi \rho |f><e|
          k4(pn, pac, ge) = k4(pn, pac, ge) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, gf) + k3(pa, pac, gf))
          k4(pn, pac, ee) = k4(pn, pac, ee) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, ef) + k3(pa, pac, ef))
          k4(pn, pac, fe) = k4(pn, pac, fe) + dt * CONJG(gkal(ja)) * xi * sqrtl(nja) * &
                          & (rho(pa, pac, ff) + k3(pa, pac, ff))
          IF (nka > 0) THEN
            ! |1_{j}>|1_{k}>
            ! Annihilate a photon from mode ka
            pn = MAP_n2p(ja, ka, nja, nka-1)
            ! \rho |e><g|
            k4(pn, pac, gg) = k4(pn, pac, gg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, ge) + k3(pa, pac, ge))
            k4(pn, pac, eg) = k4(pn, pac, eg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, ee) + k3(pa, pac, ee))
            k4(pn, pac, fg) = k4(pn, pac, fg) + dt * CONJG(gkal(ka)) * &
                            & (rho(pa, pac, fe) + k3(pa, pac, fe))
            ! \xi \rho |f><e|
            k4(pn, pac, ge) = k4(pn, pac, ge) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, gf) + k3(pa, pac, gf))
            k4(pn, pac, ee) = k4(pn, pac, ee) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, ef) + k3(pa, pac, ef))
            k4(pn, pac, fe) = k4(pn, pac, fe) + dt * CONJG(gkal(ka)) * xi * &
                            & (rho(pa, pac, ff) + k3(pa, pac, ff))
          END IF
        END IF

        ! Close p loop
      END DO
      ! Close pc loop
    END DO

    ! Update rho
    rho = rho + xis * (k1 + 2.0d0 * (k2 + k3) + k4)

    ! Check percentage
    IF (ten_percent /= 0) THEN
      IF (MOD(t, ten_percent) == 0 .AND. t /= 0) THEN
        CALL CPU_TIME(loop_check_time)
        percentage = NINT((100.0 * t) / (1.0 * tau_steps))
        loop_run_time = loop_check_time - loop_start_time
        loop_remaining_time = ((100.0 * (loop_check_time - loop_start_time)) / (1.0 * percentage)) - loop_run_time
        WRITE(*, FMT_corr) percentage, "%. Run time (g2 correlation): ", loop_run_time, &
                         & "s. Est. time left: ", loop_remaining_time, "s"
      END IF
    END IF

    ! Close time integration DO loop
  END DO

  ! Close file
  CLOSE(4)
END IF

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, " "
PRINT*, "Runtime: ", end_time - start_time, "seconds"

END PROGRAM THREE_LEVEL_BANDPASS_MASTER_EQUATION_TWO_PHOTON
