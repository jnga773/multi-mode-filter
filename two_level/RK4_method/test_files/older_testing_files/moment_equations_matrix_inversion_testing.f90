! This program solves for the operator averages (moments) for a various range
! of operator combinations. The system modelled is a multi-mode cavity array
! driven by a spontaenous emitting, resonantly driven atom. The coupling is
! achieved via a cascaded systems approach. Using a non-Hermitian Hamiltonian
! (from trajectory theory), we can disregard some of the effects of the cascaded
! decay operator.

! The Hamiltonian is in the form H = H_a + H_c + H_I, where
! -  H_a = 0.5 * Omega * (\sigma_{+} + \sigma_{-}),
! -  H_c = \sum_{ja=-N}^{N} \Delta_{ja} a^{\dagger}_{N} a_{N},
! and
! -  H_I = -0.5 i \sqrt{\gamma \kappa} *
! \sum_{ja=-N}^{N} (e^{i \phi_{j}} a^{\dagger}_{ja} \sigma_{-} - e^{-i \phi_{j}} a_{ja} \sigma_{+})
! with
! -  \Delta_{ja} = \Delta_{0} + (ja \delta\omega)
! -  \phi_{ja} = \ja \pi / N

! The master equation is in the Lindblad form:
! - d/dt \rho = 1/(i\hbar) [H, \rho] +
!   (\gamma(1 - \epsilon(2N + 1))/2) (2 |g><e| \rho |e><g| - |e><e| \rho - \rho |e><e|) +
!    \frac{\kappa}{2} \sum_{j=-N}^{N} (2 a_{j} \rho a_{j}^{\dagger} -
!                             a_{j}^{\dagger} a_{j} \rho - \rho a_{j}^{\dagger} a_{j}) +
!    \sum_{j} 1/2 (2 C_{j} \rho C^{\dagger}_{j} - C^{\dagger}_{j} C_{j} \rho - \rho C^{\dagger}_{j} C_{j}),
! where
! - C_{j} = \sqrt{\epsilon \gamma} e^{i \phi_{j}} \sigma_{-} + \sqrt{2\kappa} a_{j}.

! The input parameters (gamma, Omega, w0, kappa, epsilon, N, dw, phase, dt,
! t_max, tau1_max, tau2_max) are taken from a NameList file (ParamList.nml).

! The program outputs the atomic population of the ground and excited states,
! the mean photon number in the cavity (<A^{\dagger} A>), and a parameter output
! file:
! - ./parameters.txt,
! - ./states.txt in columns (time, |g> population, |e> population, <A^{\dagger} A>)
! - ./g1_corr.txt in columns (tau, REAL(g1), IMAG(g1)),
! - ./g2_corr.txt in columns (tau, g2)

PROGRAM TWO_LEVEL_MULTI_MODE_FILTER_MOMENT_EQUATIONS

IMPLICIT NONE

! Parameters in terms of decay rate gamma
! Atom decay rate
REAL(KIND=8) :: gamma
! Drive strength (Rabi Frequency)
REAL(KIND=8) :: Omega

! Filter parameter stuff
! Central mode frequency of the filter cavity, with n mode frequencie either side
REAL(KIND=8) :: w0
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8) :: kappa
! Percentage of fluorecence aimed at cavity
REAL(KIND=8) :: epsilon
! Number of mode either side of w0, 2N + 1 total mode
INTEGER :: N
! Frequency spacing of modes
REAL(KIND=8) :: dw
! Phase modulation of mode coupling
INTEGER :: phase

! Quantum Object stuff

! First-order moments: Bloch equations (<\sigma_{-}>, <\sigma_{+}>, <\sigma_{z}>)
!       bloch(sm) = <\sigma_{-}>,
!       bloch(sp) = <\sigma_{+}>,
!       bloch(sz) = <\sigma_{z}>.
COMPLEX(KIND=8), DIMENSION(3) :: bloch, bloch_ss
COMPLEX(KIND=8), DIMENSION(3) :: k1_bloch, k2_bloch, k3_bloch, k4_bloch
! First-order moments: Cavity (<a>, <a^{\dagger})
!       cav1(j, a)  = <a_{j}>,
!       cav1(j, at) = <a^{\dagger}_{j}>.
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE :: cav1, cav1_ss
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE :: k1_cav1, k2_cav1, k3_cav1, k4_cav1

! Second-order moments: Cavity (<a^{\dagger} a>)
!       cav2(j, k, a) = <a_{j} a_{k}>
!       cav2(j, k, at) = <a^{\dagger}_{j} a^{\dagger}_{k}>
!       cav2(j, k, ata) = <a^{\dagger}_{j} a_{k}>
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE :: cav2, cav2_ss
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE :: k1_cav2, k2_cav2, k3_cav2, k4_cav2
! Second-order moments: Cavity and atom (<a \sigma_{-+z}>, <a^{\dagger} \sigma_{-+z}>)
!       cavsig2(j, a, sm) = <a_{j} \sigma_{-}>,
!       cavsig2(j, a, sp) = <a_{j} \sigma_{+}>,
!       cavsig2(j, a, sz) = <a_{j} \sigma_{z}>,
!       cavsig2(j, at, sm) = <a^{\dagger}_{j} \sigma_{-}>,
!       cavsig2(j, at, sm) = <a^{\dagger}_{j} \sigma_{+}>,
!       cavsig2(j, at, sz) = <a^{\dagger}_{j} \sigma_{z}>.
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE :: cavsig2, cavsig2_ss
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE :: k1_cavsig2, k2_cavsig2, k3_cavsig2, k4_cavsig2

! Third-order moments: Cavity (<a^{2} a^{\dagger}>, <a^{\dagger 2} a>)
!       cav3(j, k, l, a)   = <a^{\dagger}_{j} a_{k} a_{l}>,
!       cav3(j, k, l, at)  = <a^{\dagger}_{j} a^{\dagger}_{k} a_{l}>,
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE :: cav3, cav3_ss
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE :: k1_cav3, k2_cav3, k3_cav3, k4_cav3
! Third-order moments: Cavity and atom (<aa \sigma_{-+z}>, <a^{\dagger 2} \sigma_{-+z}>, <a^{\dagger} a \sigma_{-+z}>)
!       cavsig3(j, k, a, sm) = <a_{j} a_{k} \sigma_{-}>,
!       cavsig3(j, k, a, sp) = <a_{j} a_{k} \sigma_{+}>,
!       cavsig3(j, k, a, sz) = <a_{j} a_{k} \sigma_{z}>,
!       cavsig3(j, k, at, sm) = <a^{\dagger}_{j} a^{\dagger}_{l} \sigma_{-}>,
!       cavsig3(j, k, at, sp) = <a^{\dagger}_{j} a^{\dagger}_{l} \sigma_{+}>,
!       cavsig3(j, k, at, sz) = <a^{\dagger}_{j} a^{\dagger}_{l} \sigma_{z}>.
!       cavsig3(j, k, ata, sm) = <a^{\dagger}_{j} a_{l} \sigma_{-}>,
!       cavsig3(j, k, ata, sp) = <a^{\dagger}_{j} a_{l} \sigma_{+}>,
!       cavsig3(j, k, ata, sz) = <a^{\dagger}_{j} a_{l} \sigma_{z}>.
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE :: cavsig3, cavsig3_ss
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE :: k1_cavsig3, k2_cavsig3, k3_cavsig3, k4_cavsig3

! Fourth-order moments: Cavity (<a^{\dagger 2} a^{2}>)
!       cav4(j, k, l, m) = <a^{\dagger}_{j} a^{\dagger}_{k} a_{l} a_{m}>
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE :: cav4, cav4_ss
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE :: k1_cav4, k2_cav4, k3_cav4, k4_cav4
! Fourth-order moments: Cavity and atom ( <a^{\dagger} a^{2} \sigma_{-+z}>, <a^{\dagger 2} a \sigma_{-+z}>)
!       cavsig4(j, k, l, a, sm) = <a^{\dagger}_{j} a_{l} a_{m} \sigma_{-}>,
!       cavsig4(j, k, l, a, sp) = <a^{\dagger}_{j} a_{l} a_{m} \sigma_{+}>,
!       cavsig4(j, k, l, a, sz) = <a^{\dagger}_{j} a_{l} a_{m} \sigma_{z}>,
!       cavsig4(j, k, l, at, sm) = <a^{\dagger}_{j} a^{\dagger}_{l} a_{m} \sigma_{-}>,
!       cavsig4(j, k, l, at, sp) = <a^{\dagger}_{j} a^{\dagger}_{l} a_{m} \sigma_{+}>,
!       cavsig4(j, k, l, at, sz) = <a^{\dagger}_{j} a^{\dagger}_{l} a_{m} \sigma_{z}>,
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: cavsig4, cavsig4_ss
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: k1_cavsig4, k2_cavsig4, k3_cavsig4, k4_cavsig4

! Integer indices for: sigma_{-}, sigma_{+}, sigma_{z}
INTEGER, PARAMETER :: sm = 1, sp = 2, sz = 3
! Integer indices for: a, a^{\dagger}, a^{\dagger} a
INTEGER, PARAMETER :: a = 1, at = 2, ata = 3

! Time stuff
! Time step
REAL(KIND=8) :: dt
! Maximum time to integrate for
REAL(KIND=8) :: t_max
! Maximum number of steps to integrate for
INTEGER :: t_steps
! Maximum tau time to calculate correlation function for
REAL(KIND=8) :: tau1_max, tau2_max
! Max number of tau steps
INTEGER :: tau_steps
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
! Filter mode integers
INTEGER :: j, k, l, m
! Temporal complex values
COMPLEX(KIND=8) :: alpha, beta
! Inverse matrix stuff
COMPLEX(KIND=8), DIMENSION(3, 3) :: M_inverse
COMPLEX(KIND=8), DIMENSION(3) :: B_ss
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
! Print completion time check
! LOGICAL, PARAMETER :: progress_bar = .TRUE.
LOGICAL, PARAMETER :: progress_bar = .FALSE.

! Data stuff
! Sample rate for state populations
INTEGER :: sample_rate
! State population
REAL(KIND=8) :: rho_gg, rho_ee
! Mean photon number in cavity
REAL(KIND=8) :: photona, photona_ss
! First-order correlation
COMPLEX(KIND=8) :: corr

! Filename stuff
! I/O Format for percentage
CHARACTER(LEN=28), PARAMETER :: FMT_ss = "(T2,I3,A28,F9.2,A19,F9.2,A1)"
CHARACTER(LEN=28), PARAMETER :: FMT_corr = "(T2,I3,A27,F9.2,A19,F9.2,A1)"
! Filename of parameters
CHARACTER(LEN=27), PARAMETER :: filename_parameters = "./data_files/parameters.txt"
! Filename for state population
CHARACTER(LEN=23), PARAMETER :: filename_state = "./data_files/states.txt"
! Filename for first-order correlation of cavity a
CHARACTER(LEN=24), PARAMETER :: filename_g1 = "./data_files/g1_corr.txt"
! Filename for second-order correlation of cavity a
CHARACTER(LEN=24), PARAMETER :: filename_g2 = "./data_files/g2_corr.txt"
! Paramert Name List
CHARACTER(LEN=15), PARAMETER :: filename_ParamList = "./ParamList.nml"

!------------------------------------------------------------------------------!
!                    Calling parameters from ParamList_me.nml                  !
!------------------------------------------------------------------------------!
INTEGER :: IUNIT, ISTAT
CHARACTER(LEN=512) :: LINE
NAMELIST /ATOM/ gamma, Omega
NAMELIST /CAVITY_A/ w0, kappa, epsilon, N, dw, phase
NAMELIST /TIME/ dt, t_max, tau1_max, tau2_max

! Call start time from CPU_TIME
CALL CPU_TIME(start_time)

! read the PARAMS namelist
IUNIT = 420
OPEN(UNIT=IUNIT, FILE=filename_paramlist, STATUS="OLD", DELIM="QUOTE")
READ(UNIT=IUNIT, NML=ATOM, IOSTAT=istat)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') line
  CLOSE(IUNIT)
  PRINT *, "Invalid line in ATOM namelist: " // TRIM(LINE)
  CALL EXIT(1)
END IF
READ(UNIT=IUNIT, NML=CAVITY_A, IOSTAT=istat)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') line
  CLOSE(IUNIT)
  PRINT *, "Invalid line in CAVITY_A namelist: " // TRIM(LINE)
  CALL EXIT(1)
END IF
READ(UNIT=IUNIT, NML=TIME, IOSTAT=istat)
IF (ISTAT .NE. 0) THEN
  BACKSPACE(IUNIT)
  READ(IUNIT, FMT='(A)') line
  CLOSE(IUNIT)
  PRINT *, "Invalid line in TIME namelist: " // TRIM(LINE)
  CALL EXIT(1)
END IF
CLOSE(IUNIT)

! Max number of time steps
t_steps = NINT(t_max / dt)

!------------------------------------------------------------------------------!
!                       Allocate and initialise arrays                         !
!------------------------------------------------------------------------------!
! Allocate array of Delta and gka values
ALLOCATE(wl(-N:N))
wl = 0.0d0
ALLOCATE(gkal(-N:N))
gkal = 0.0d0
DO j = -N, N
  IF (N == 0) THEN
    wl(j) = w0
    gkal(j) = DSQRT(epsilon * gamma * kappa)
  ELSE
    wl(j) = w0 + DBLE(j) * dw
    ! Blackman window coefficient
    blackman = 1.0d0
    ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + ja) / (2.0d0 * DBLE(N))) + &
    !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + ja) / (2.0d0 * DBLE(N)))
    ! Mode dependent phase difference
    gkal(j) = DSQRT((epsilon / DBLE(2*N + 1)) * gamma * kappa) * blackman * EXP(i * DBLE(phase) * DBLE(j) * pi / DBLE(N))
  END IF
END DO

! Initialise moments
! First-order
bloch = 0.0d0; bloch_ss = 0.0d0
ALLOCATE(cav1(-N:N, 2)); cav1 = 0.0d0; ALLOCATE(cav1_ss(-N:N, 2)); cav1_ss = 0.0d0
! Second-order
ALLOCATE(cav2(-N:N, -N:N, 3)); cav2 = 0.0d0; ALLOCATE(cav2_ss(-N:N, -N:N, 3)); cav2_ss = 0.0d0
ALLOCATE(cavsig2(-N:N, 2, 3)); cavsig2 = 0.0d0; ALLOCATE(cavsig2_ss(-N:N, 2, 3)); cavsig2_ss = 0.0d0
! Third-order
ALLOCATE(cav3(-N:N, -N:N, -N:N, 2)); cav3 = 0.0d0; ALLOCATE(cav3_ss(-N:N, -N:N, -N:N, 2)); cavsig3 = 0.0d0
ALLOCATE(cavsig3(-N:N, -N:N, 3, 3)); cavsig3 = 0.0d0; ALLOCATE(cavsig3_ss(-N:N, -N:N, 3, 3)); cavsig3_ss = 0.0d0;
! Fourth-order
ALLOCATE(cav4(-N:N, -N:N, -N:N, -N:N)); cav4 = 0.0d0; ALLOCATE(cav4_ss(-N:N, -N:N, -N:N, -N:N)); cav4_ss = 0.0d0
ALLOCATE(cavsig4(-N:N, -N:N, -N:N, 2, 3)); cavsig4 = 0.0d0; ALLOCATE(cavsig4_ss(-N:N, -N:N, -N:N, 2, 3)); cavsig4_ss = 0.0d0

! Initialise Runge-Kutta vectors
! First-order
ALLOCATE(k1_cav1(-N:N, 2)); k1_bloch = 0.0d0
ALLOCATE(k2_cav1(-N:N, 2)); k2_bloch = 0.0d0
ALLOCATE(k3_cav1(-N:N, 2)); k3_bloch = 0.0d0
ALLOCATE(k4_cav1(-N:N, 2)); k4_bloch = 0.0d0
! Second-order
ALLOCATE(k1_cav2(-N:N, -N:N, 3)); k1_cav2 = 0.0d0
ALLOCATE(k2_cav2(-N:N, -N:N, 3)); k2_cav2 = 0.0d0
ALLOCATE(k3_cav2(-N:N, -N:N, 3)); k3_cav2 = 0.0d0
ALLOCATE(k4_cav2(-N:N, -N:N, 3)); k4_cav2 = 0.0d0
ALLOCATE(k1_cavsig2(-N:N, 2, 3)); k1_cavsig2 = 0.0d0
ALLOCATE(k2_cavsig2(-N:N, 2, 3)); k2_cavsig2 = 0.0d0
ALLOCATE(k3_cavsig2(-N:N, 2, 3)); k3_cavsig2 = 0.0d0
ALLOCATE(k4_cavsig2(-N:N, 2, 3)); k4_cavsig2 = 0.0d0
! Third-order
ALLOCATE(k1_cav3(-N:N, -N:N, -N:N, 2)); k1_cav3 = 0.0d0
ALLOCATE(k2_cav3(-N:N, -N:N, -N:N, 2)); k2_cav3 = 0.0d0
ALLOCATE(k3_cav3(-N:N, -N:N, -N:N, 2)); k3_cav3 = 0.0d0
ALLOCATE(k4_cav3(-N:N, -N:N, -N:N, 2)); k4_cav3 = 0.0d0
ALLOCATE(k1_cavsig3(-N:N, -N:N, 3, 3)); k1_cavsig3 = 0.0d0
ALLOCATE(k2_cavsig3(-N:N, -N:N, 3, 3)); k2_cavsig3 = 0.0d0
ALLOCATE(k3_cavsig3(-N:N, -N:N, 3, 3)); k3_cavsig3 = 0.0d0
ALLOCATE(k4_cavsig3(-N:N, -N:N, 3, 3)); k4_cavsig3 = 0.0d0
! Fourth-order
ALLOCATE(k1_cav4(-N:N, -N:N, -N:N, -N:N)); k1_cav4 = 0.0d0
ALLOCATE(k2_cav4(-N:N, -N:N, -N:N, -N:N)); k2_cav4 = 0.0d0
ALLOCATE(k3_cav4(-N:N, -N:N, -N:N, -N:N)); k3_cav4 = 0.0d0
ALLOCATE(k4_cav4(-N:N, -N:N, -N:N, -N:N)); k4_cav4 = 0.0d0
ALLOCATE(k1_cavsig4(-N:N, -N:N, -N:N, 2, 3)); k1_cavsig4 = 0.0d0
ALLOCATE(k2_cavsig4(-N:N, -N:N, -N:N, 2, 3)); k2_cavsig4 = 0.0d0
ALLOCATE(k3_cavsig4(-N:N, -N:N, -N:N, 2, 3)); k3_cavsig4 = 0.0d0
ALLOCATE(k4_cavsig4(-N:N, -N:N, -N:N, 2, 3)); k4_cavsig4 = 0.0d0

! Atom is initially in ground state with zero photons in cavity: |psi(0)> = |g, 0>
! Initial condition is then <\sigma_{z}> = -1.0, and <a^{\dagger} a> = 0.0
bloch(sz) = -1.0d0
! cav2(ata) = 0.0d0

! Initialise data variables
rho_gg = 0.0d0
rho_ee = 0.0d0
photona = 0.0d0
corr = 0.0d0

! Open file to write time to
OPEN(UNIT=1, FILE=filename_parameters, STATUS='REPLACE', ACTION='WRITE')
! Write parameter
WRITE(1,*) "Parameters are in the following order:"
WRITE(1,"(A10,F25.15)") "gamma =", gamma
WRITE(1,"(A10,F25.15)") "Omega =", Omega
WRITE(1,"(A10,F25.15)") "w0 =", w0
WRITE(1,"(A10,F25.15)") "kappa =", kappa
WRITE(1,"(A10,F25.15)") "epsilon =", epsilon
WRITE(1,"(A11,I9)") "N = ", N
WRITE(1,"(A11,F25.15)") "dw = ", dw
WRITE(1,"(A11,I9)") "phase = ", phase
WRITE(1,"(A10,F25.15)") "dt =", dt
IF (t_max <= 1000.0d0) THEN
  WRITE(1,"(A10,F25.15)") "Max time =", t_max
ELSE
  WRITE(1,"(A10,ES25.11E2)") "Max time =", t_max
END IF
WRITE(1,"(A10,F25.15)") "Max tau1 =", tau1_max
WRITE(1,"(A10,F25.15)") "Max tau2 =", tau2_max
! Close file
CLOSE(1)

!------------------------------------------------------------------------------!
!                                STATE-SOLVER                                  !
!------------------------------------------------------------------------------!

IF (t_max > 0.0) THEN
  ! Open file to write time and data to
  OPEN(UNIT=2, FILE=filename_state, STATUS='REPLACE', ACTION='WRITE', RECL=4000)

  ! Ten percent of time steps
  ten_percent = NINT((1.0 * t_steps / 10.0))
  ! Let max number of data points be 100,000 ~ 10mb file.
  IF (t_steps > 100000) THEN
    sample_rate = NINT(DBLE(t_steps) / 1d5)
  ELSE
    sample_rate = 1
  END IF
  ! Call CPU clock time
  CALL CPU_TIME(loop_start_time)

  ! Start time integration to find steady state
  DO t = 0, t_steps
    ! Calculate state probabilities and mean photon number
    rho_gg = 0.0d0
    rho_gg = REAL(0.5d0 * (1.0d0 - bloch(sz)))
    rho_ee = 0.0d0
    rho_ee = REAL(0.5d0 * (bloch(sz) + 1.0d0))
    photona = 0.0d0
    DO k = -N, N
      DO j = -N, N
        photona = photona + cav2(j, k, ata)
      END DO
    END DO

    ! If t_max is really big, only take a sample of results to write to file
    ! so file size isn't huge-mongous.
    IF (MOD(t, sample_rate) == 0) THEN
      ! Write data
      WRITE(2,*) DBLE(t) * dt, rho_gg, rho_ee, photona
    END IF

    ! Solve the operator moment differential equations using Runge-Kutta 4th Order

    !-------------------!
    !     CALCULATE     !
    !-------------------!
    ! Initialise Runge-Kutta Vectors
    k1_bloch = 0.0d0; k2_bloch = 0.0d0; k3_bloch = 0.0d0; k4_bloch = 0.0d0
    k1_cav1 = 0.0d0; k2_cav1 = 0.0d0; k3_cav1 = 0.0d0; k4_cav1 = 0.0d0
    k1_cav2 = 0.0d0; k2_cav2 = 0.0d0; k3_cav2 = 0.0d0; k4_cav2 = 0.0d0
    k1_cavsig2 = 0.0d0; k2_cavsig2 = 0.0d0; k3_cavsig2 = 0.0d0; k4_cavsig2 = 0.0d0
    k1_cav3 = 0.0d0; k2_cav3 = 0.0d0; k3_cav2 = 0.0d0; k4_cav3 = 0.0d0
    k1_cavsig3 = 0.0d0; k2_cavsig3 = 0.0d0; k3_cavsig3 = 0.0d0; k4_cavsig3 = 0.0d0
    k1_cav4 = 0.0d0; k2_cav4 = 0.0d0; k3_cav4 = 0.0d0; k4_cav4 = 0.0d0
    k1_cavsig4 = 0.0d0; k2_cavsig4 = 0.0d0; k3_cavsig4 = 0.0d0; k4_cavsig4 = 0.0d0

    !----------------------------!
    !     First-Order: Bloch     !
    !----------------------------!

    ! k1 - First-order: Bloch equations (<sm>, <sp>, <sz>)
    k1_bloch(sm) = -dt * 0.5d0 * gamma * bloch(sm) + &
                 & dt * i * 0.5d0 * Omega * bloch(sz)
    k1_bloch(sp) = -dt * 0.5d0 * gamma * bloch(sp) - &
                 & dt * i * 0.5d0 * Omega * bloch(sz)
    k1_bloch(sz) = -dt * gamma * bloch(sz) + &
                 & dt * i * Omega * bloch(sm) - &
                 & dt * i * Omega * bloch(sp) - &
                 & dt * gamma

    ! k2 - First-order: Bloch equations (<sm>, <sp>, <sz>)
    k2_bloch(sm) = -dt * 0.5d0 * gamma * (bloch(sm) + 0.5d0 * k1_bloch(sm)) + &
                 & dt * i * 0.5d0 * Omega * (bloch(sz) + 0.5d0 * k1_bloch(sz))
    k2_bloch(sp) = -dt * 0.5d0 * gamma * (bloch(sp) + 0.5d0 * k1_bloch(sp)) - &
                 & dt * i * 0.5d0 * Omega * (bloch(sz) + 0.5d0 * k1_bloch(sz))
    k2_bloch(sz) = -dt * gamma * (bloch(sz) + 0.5d0 * k1_bloch(sz)) + &
                 & dt * i * Omega * (bloch(sm) + 0.5d0 * k1_bloch(sm)) - &
                 & dt * i * Omega * (bloch(sp) + 0.5d0 * k1_bloch(sp)) - &
                 & dt * gamma

    ! k3 - First-order: Bloch equations (<sm>, <sp>, <sz>)
    k3_bloch(sm) = -dt * 0.5d0 * gamma * (bloch(sm) + 0.5d0 * k2_bloch(sm)) + &
                 & dt * i * 0.5d0 * Omega * (bloch(sz) + 0.5d0 * k2_bloch(sz))
    k3_bloch(sp) = -dt * 0.5d0 * gamma * (bloch(sp) + 0.5d0 * k2_bloch(sp)) - &
                 & dt * i * 0.5d0 * Omega * (bloch(sz) + 0.5d0 * k2_bloch(sz))
    k3_bloch(sz) = -dt * gamma * (bloch(sz) + 0.5d0 * k2_bloch(sz)) + &
                 & dt * i * Omega * (bloch(sm) + 0.5d0 * k2_bloch(sm)) - &
                 & dt * i * Omega * (bloch(sp) + 0.5d0 * k2_bloch(sp)) - &
                 & dt * gamma

    ! k4 - First-order: Bloch equations (<sm>, <sp>, <sz>)
    k4_bloch(sm) = -dt * 0.5d0 * gamma * (bloch(sm) + k3_bloch(sm)) + &
                 & dt * i * 0.5d0 * Omega * (bloch(sz) + k3_bloch(sz))
    k4_bloch(sp) = -dt * 0.5d0 * gamma * (bloch(sp) + k3_bloch(sp)) - &
                 & dt * i * 0.5d0 * Omega * (bloch(sz) + k3_bloch(sz))
    k4_bloch(sz) = -dt * gamma * (bloch(sz) + k3_bloch(sz)) + &
                 & dt * i * Omega * (bloch(sm) + k3_bloch(sm)) - &
                 & dt * i * Omega * (bloch(sp) + k3_bloch(sp)) - &
                 & dt * gamma

    ! Cycle through mode numbers
    DO j = -N, N
      !-----------------------------!
      !     First-Order: Cavity     !
      !-----------------------------!
      ! k1 - First-order: Cavity mode (<a>, <at>)
      k1_cav1(j, a) = -dt * (kappa + i * wl(j)) * cav1(j, a) - &
                    & dt * gkal(j) * bloch(sm)
      k1_cav1(j, at) = -dt * (kappa - i * wl(j)) * cav1(j, at) - &
                     & dt * CONJG(gkal(j)) * bloch(sp)

      ! k2 - First-order: Cavity mode (<a>, <at>)
      k2_cav1(j, a) = -dt * (kappa + i * wl(j)) * (cav1(j, a) + 0.5d0 * k1_cav1(j, a)) - &
                    & dt * gkal(j) * (bloch(sm) + 0.5d0 * k1_bloch(sm))
      k2_cav1(j, at) = -dt * (kappa - i * wl(j)) * (cav1(j, at) + 0.5d0 * k1_cav1(j, at)) - &
                     & dt * CONJG(gkal(j)) * (bloch(sp) + 0.5d0 * k1_bloch(sp))

      ! k3 - First-order: Cavity mode (<a>, <at>)
      k3_cav1(j, a) = -dt * (kappa + i * wl(j)) * (cav1(j, a) + 0.5d0 * k2_cav1(j, a)) - &
                    & dt * gkal(j) * (bloch(sm) + 0.5d0 * k2_bloch(sm))
      k3_cav1(j, at) = -dt * (kappa - i * wl(j)) * (cav1(j, at) + 0.5d0 * k2_cav1(j, at)) - &
                     & dt * CONJG(gkal(j)) * (bloch(sp) + 0.5d0 * k2_bloch(sp))

      ! k4 - First-order: Cavity mode (<a>, <at>)
      k4_cav1(j, a) = -dt * (kappa + i * wl(j)) * (cav1(j, a) + k3_cav1(j, a)) - &
                    & dt * gkal(j) * (bloch(sm) + k3_bloch(sm))
      k4_cav1(j, at) = -dt * (kappa - i * wl(j)) * (cav1(j, at) + k3_cav1(j, at)) - &
                     & dt * CONJG(gkal(j)) * (bloch(sp) + k3_bloch(sp))

      !---------------------------------------!
      !     Second-Order: Cavity and Atom     !
      !---------------------------------------!
      ! k1 - Second-order: a and \sigma (<a sm>, <a sp>, <a sz>)
      k1_cavsig2(j, a, sm) = -dt * (0.5d0 * gamma + kappa + i * wl(j)) * cavsig2(j, a, sm) + &
                           & dt * i * 0.5d0 * Omega * cavsig2(j, a, sz)
      k1_cavsig2(j, a, sp) = -dt * (0.5d0 * gamma + kappa + i * wl(j)) * cavsig2(j, a, sp) - &
                           & dt * i * 0.5d0 * Omega * cavsig2(j, a, sz) - &
                           & dt * 0.5d0 * gkal(j) * (bloch(sz) + 1.0d0)
      k1_cavsig2(j, a, sz) = -dt * (gamma + kappa + i * wl(j)) * cavsig2(j, a, sz) + &
                           & dt * i * Omega * cavsig2(j, a, sm) - &
                           & dt * i * Omega * cavsig2(j, a, sp) - &
                           & dt * gamma * cav1(j, a) + &
                           & dt * gkal(j) * bloch(sm)
      ! k1 - Second-order: at and \sigma (<at sm>, <at sp>, <at sz>)
      k1_cavsig2(j, at, sm) = -dt * (0.5d0 * gamma + kappa - i * wl(j)) * cavsig2(j, at, sm) + &
                            & dt * i * 0.5d0 * Omega * cavsig2(j, at, sz) - &
                            & dt * 0.5d0 * CONJG(gkal(j)) * (bloch(sz) + 1.0d0)
      k1_cavsig2(j, at, sp) = -dt * (0.5d0 * gamma + kappa - i * wl(j)) * cavsig2(j, at, sp) - &
                            & dt * i * 0.5d0 * Omega * cavsig2(j, at, sz)
      k1_cavsig2(j, at, sz) = -dt * (gamma + kappa - i * wl(j)) * cavsig2(j, at, sz) + &
                            & dt * i * Omega * cavsig2(j, at, sm) - &
                            & dt * i * Omega * cavsig2(j, at, sp) - &
                            & dt * gamma * cav1(j, at) + &
                            & dt * CONJG(gkal(j)) * bloch(sp)

      ! k2 - Second-order: a and \sigma (<a sm>, <a sp>, <a sz>)
      k2_cavsig2(j, a, sm) = -dt * (0.5d0 * gamma + kappa + i * wl(j)) * (cavsig2(j, a, sm) + 0.5d0 * k1_cavsig2(j, a, sm)) + &
                           & dt * i * 0.5d0 * Omega * (cavsig2(j, a, sz) + 0.5d0 * k1_cavsig2(j, a, sz))
      k2_cavsig2(j, a, sp) = -dt * (0.5d0 * gamma + kappa + i * wl(j)) * (cavsig2(j, a, sp) + 0.5d0 * k1_cavsig2(j, a, sp)) - &
                           & dt * i * 0.5d0 * Omega * (cavsig2(j, a, sz) + 0.5d0 * k1_cavsig2(j, a, sz)) - &
                           & dt * 0.5d0 * gkal(j) * ((bloch(sz) + 0.5d0 * k1_bloch(sz)) + 1.0d0)
      k2_cavsig2(j, a, sz) = -dt * (gamma + kappa + i * wl(j)) * (cavsig2(j, a, sz) + 0.5d0 * k1_cavsig2(j, a, sz)) + &
                           & dt * i * Omega * (cavsig2(j, a, sm) + 0.5d0 * k1_cavsig2(j, a, sm)) - &
                           & dt * i * Omega * (cavsig2(j, a, sp) + 0.5d0 * k1_cavsig2(j, a, sp)) - &
                           & dt * gamma * (cav1(j, a) + 0.5d0 * k1_cav1(j, a)) + &
                           & dt * gkal(j) * (bloch(sm) + 0.5d0 * k1_bloch(sm))
      ! k2 - Second-order: at and \sigma (<at sm>, <at sp>, <at sz>)
      k2_cavsig2(j, at, sm) = -dt * (0.5d0 * gamma + kappa - i * wl(j)) * (cavsig2(j, at, sm) + 0.5d0 * k1_cavsig2(j, at, sm)) + &
                            & dt * i * 0.5d0 * Omega * (cavsig2(j, at, sz) + 0.5d0 * k1_cavsig2(j, at, sz)) - &
                            & dt * 0.5d0 * CONJG(gkal(j)) * ((bloch(sz) + 0.5d0 * k1_bloch(sz)) + 1.0d0)
      k2_cavsig2(j, at, sp) = -dt * (0.5d0 * gamma + kappa - i * wl(j)) * (cavsig2(j, at, sp) + 0.5d0 * k1_cavsig2(j, at, sp)) - &
                            & dt * i * 0.5d0 * Omega * (cavsig2(j, at, sz) + 0.5d0 * k1_cavsig2(j, at, sz))
      k2_cavsig2(j, at, sz) = -dt * (gamma + kappa - i * wl(j)) * (cavsig2(j, at, sz) + 0.5d0 * k1_cavsig2(j, at, sz)) + &
                            & dt * i * Omega * (cavsig2(j, at, sm) + 0.5d0 * k1_cavsig2(j, at, sm)) - &
                            & dt * i * Omega * (cavsig2(j, at, sp) + 0.5d0 * k1_cavsig2(j, at, sp)) - &
                            & dt * gamma * (cav1(j, at) + 0.5d0 * k1_cav1(j, at)) + &
                            & dt * CONJG(gkal(j)) * (bloch(sp) + 0.5d0 * k1_bloch(sp))

      ! k3 - Second-order: a and \sigma (<a sm>, <a sp>, <a sz>)
      k3_cavsig2(j, a, sm) = -dt * (0.5d0 * gamma + kappa + i * wl(j)) * (cavsig2(j, a, sm) + 0.5d0 * k2_cavsig2(j, a, sm)) + &
                           & dt * i * 0.5d0 * Omega * (cavsig2(j, a, sz) + 0.5d0 * k2_cavsig2(j, a, sz))
      k3_cavsig2(j, a, sp) = -dt * (0.5d0 * gamma + kappa + i * wl(j)) * (cavsig2(j, a, sp) + 0.5d0 * k2_cavsig2(j, a, sp)) - &
                           & dt * i * 0.5d0 * Omega * (cavsig2(j, a, sz) + 0.5d0 * k2_cavsig2(j, a, sz)) - &
                           & dt * 0.5d0 * gkal(j) * ((bloch(sz) + 0.5d0 * k2_bloch(sz)) + 1.0d0)
      k3_cavsig2(j, a, sz) = -dt * (gamma + kappa + i * wl(j)) * (cavsig2(j, a, sz) + 0.5d0 * k2_cavsig2(j, a, sz)) + &
                           & dt * i * Omega * (cavsig2(j, a, sm) + 0.5d0 * k2_cavsig2(j, a, sm)) - &
                           & dt * i * Omega * (cavsig2(j, a, sp) + 0.5d0 * k2_cavsig2(j, a, sp)) - &
                           & dt * gamma * (cav1(j, a) + 0.5d0 * k2_cav1(j, a)) + &
                           & dt * gkal(j) * (bloch(sm) + 0.5d0 * k2_bloch(sm))
      ! k3 - Second-order: at and \sigma (<at sm>, <at sp>, <at sz>)
      k3_cavsig2(j, at, sm) = -dt * (0.5d0 * gamma + kappa - i * wl(j)) * (cavsig2(j, at, sm) + 0.5d0 * k2_cavsig2(j, at, sm)) + &
                            & dt * i * 0.5d0 * Omega * (cavsig2(j, at, sz) + 0.5d0 * k2_cavsig2(j, at, sz)) - &
                            & dt * 0.5d0 * CONJG(gkal(j)) * ((bloch(sz) + 0.5d0 * k2_bloch(sz)) + 1.0d0)
      k3_cavsig2(j, at, sp) = -dt * (0.5d0 * gamma + kappa - i * wl(j)) * (cavsig2(j, at, sp) + 0.5d0 * k2_cavsig2(j, at, sp)) - &
                            & dt * i * 0.5d0 * Omega * (cavsig2(j, at, sz) + 0.5d0 * k2_cavsig2(j, at, sz))
      k3_cavsig2(j, at, sz) = -dt * (gamma + kappa - i * wl(j)) * (cavsig2(j, at, sz) + 0.5d0 * k2_cavsig2(j, at, sz)) + &
                            & dt * i * Omega * (cavsig2(j, at, sm) + 0.5d0 * k2_cavsig2(j, at, sm)) - &
                            & dt * i * Omega * (cavsig2(j, at, sp) + 0.5d0 * k2_cavsig2(j, at, sp)) - &
                            & dt * gamma * (cav1(j, at) + 0.5d0 * k2_cav1(j, at)) + &
                            & dt * CONJG(gkal(j)) * (bloch(sp) + 0.5d0 * k2_bloch(sp))

      ! k4 - Second-order: a and \sigma (<a sm>, <a sp>, <a sz>)
      k4_cavsig2(j, a, sm) = -dt * (0.5d0 * gamma + kappa + i * wl(j)) * (cavsig2(j, a, sm) + k3_cavsig2(j, a, sm)) + &
                           & dt * i * 0.5d0 * Omega * (cavsig2(j, a, sz) + k3_cavsig2(j, a, sz))
      k4_cavsig2(j, a, sp) = -dt * (0.5d0 * gamma + kappa + i * wl(j)) * (cavsig2(j, a, sp) + k3_cavsig2(j, a, sp)) - &
                           & dt * i * 0.5d0 * Omega * (cavsig2(j, a, sz) + k3_cavsig2(j, a, sz)) - &
                           & dt * 0.5d0 * gkal(j) * ((bloch(sz) + k3_bloch(sz)) + 1.0d0)
      k4_cavsig2(j, a, sz) = -dt * (gamma + kappa + i * wl(j)) * (cavsig2(j, a, sz) + k3_cavsig2(j, a, sz)) + &
                           & dt * i * Omega * (cavsig2(j, a, sm) + k3_cavsig2(j, a, sm)) - &
                           & dt * i * Omega * (cavsig2(j, a, sp) + k3_cavsig2(j, a, sp)) - &
                           & dt * gamma * (cav1(j, a) + k3_cav1(j, a)) + &
                           & dt * gkal(j) * (bloch(sm) + k3_bloch(sm))
      ! k4 - Second-order: at and \sigma (<at sm>, <at sp>, <at sz>)
      k4_cavsig2(j, at, sm) = -dt * (0.5d0 * gamma + kappa - i * wl(j)) * (cavsig2(j, at, sm) + k3_cavsig2(j, at, sm)) + &
                            & dt * i * 0.5d0 * Omega * (cavsig2(j, at, sz) + k3_cavsig2(j, at, sz)) - &
                            & dt * 0.5d0 * CONJG(gkal(j)) * ((bloch(sz) + k3_bloch(sz)) + 1.0d0)
      k4_cavsig2(j, at, sp) = -dt * (0.5d0 * gamma + kappa - i * wl(j)) * (cavsig2(j, at, sp) + k3_cavsig2(j, at, sp)) - &
                            & dt * i * 0.5d0 * Omega * (cavsig2(j, at, sz) + k3_cavsig2(j, at, sz))
      k4_cavsig2(j, at, sz) = -dt * (gamma + kappa - i * wl(j)) * (cavsig2(j, at, sz) + k3_cavsig2(j, at, sz)) + &
                            & dt * i * Omega * (cavsig2(j, at, sm) + k3_cavsig2(j, at, sm)) - &
                            & dt * i * Omega * (cavsig2(j, at, sp) + k3_cavsig2(j, at, sp)) - &
                            & dt * gamma * (cav1(j, at) + k3_cav1(j, at)) + &
                            & dt * CONJG(gkal(j)) * (bloch(sp) + k3_bloch(sp))

      ! Cycle through mode numbers
      DO k = -N, N
        !------------------------------!
        !     Second-Order: Cavity     !
        !------------------------------!
        ! k1 - Second-order: Cavity (<a a>, <at at>, <at a>)
        k1_cav2(j, k, a) = -dt * (2.0d0 * kappa + i * (wl(j) + wl(k))) * cav2(j, k, a) - &
                         & dt * gkal(j) * cavsig2(k, a, sm) - &
                         & dt * gkal(k) * cavsig2(j, a, sm)
        k1_cav2(j, k, at) = -dt * (2.0d0 * kappa - i * (wl(j) + wl(k))) * cav2(j, k, at) - &
                          & dt * CONJG(gkal(j)) * cavsig2(k, at, sp) - &
                          & dt * CONJG(gkal(k)) * cavsig2(j, at, sp)
        k1_cav2(j, k, ata) = -dt * (2.0d0 * kappa - i * (wl(j) - wl(k))) * cav2(j, k, ata) - &
                           & dt * CONJG(gkal(j)) * cavsig2(k, a, sp) - &
                           & dt * gkal(k) * cavsig2(j, at, sm)

        ! k2 -  Second-order: Cavity (<a a>, <at at>, <at a>)
        k2_cav2(j, k, a) = -dt * (2.0d0 * kappa + i * (wl(j) + wl(k))) * (cav2(j, k, a) + 0.5d0 * k1_cav2(j, k, a)) - &
                         & dt * gkal(j) * (cavsig2(k, a, sm) + 0.5d0 * k1_cavsig2(k, a, sm)) - &
                         & dt * gkal(k) * (cavsig2(j, a, sm) + 0.5d0 * k1_cavsig2(j, a, sm))
        k2_cav2(j, k, at) = -dt * (2.0d0 * kappa - i * (wl(j) + wl(k))) * (cav2(j, k, at) + 0.5d0 * k1_cav2(j, k, at)) - &
                          & dt * CONJG(gkal(j)) * (cavsig2(k, at, sp) + 0.5d0 * k1_cavsig2(k, at, sp)) - &
                          & dt * CONJG(gkal(k)) * (cavsig2(j, at, sp) + 0.5d0 * k1_cavsig2(j, at, sp))
        k2_cav2(j, k, ata) = -dt * (2.0d0 * kappa - i * (wl(j) - wl(k))) * (cav2(j, k, ata) + 0.5d0 * k1_cav2(j, k, ata)) - &
                           & dt * CONJG(gkal(j)) * (cavsig2(k, a, sp) + 0.5d0 * k1_cavsig2(k, a, sp)) - &
                           & dt * gkal(k) * (cavsig2(j, at, sm) + 0.5d0 * k1_cavsig2(j, at, sm))

        ! k3 -  Second-order: Cavity (<a a>, <at at>, <at a>)
        k3_cav2(j, k, a) = -dt * (2.0d0 * kappa + i * (wl(j) + wl(k))) * (cav2(j, k, a) + 0.5d0 * k2_cav2(j, k, a)) - &
                         & dt * gkal(j) * (cavsig2(k, a, sm) + 0.5d0 * k2_cavsig2(k, a, sm)) - &
                         & dt * gkal(k) * (cavsig2(j, a, sm) + 0.5d0 * k2_cavsig2(j, a, sm))
        k3_cav2(j, k, at) = -dt * (2.0d0 * kappa - i * (wl(j) + wl(k))) * (cav2(j, k, at) + 0.5d0 * k2_cav2(j, k, at)) - &
                          & dt * CONJG(gkal(j)) * (cavsig2(k, at, sp) + 0.5d0 * k2_cavsig2(k, at, sp)) - &
                          & dt * CONJG(gkal(k)) * (cavsig2(j, at, sp) + 0.5d0 * k2_cavsig2(j, at, sp))
        k3_cav2(j, k, ata) = -dt * (2.0d0 * kappa - i * (wl(j) - wl(k))) * (cav2(j, k, ata) + 0.5d0 * k2_cav2(j, k, ata)) - &
                           & dt * CONJG(gkal(j)) * (cavsig2(k, a, sp) + 0.5d0 * k2_cavsig2(k, a, sp)) - &
                           & dt * gkal(k) * (cavsig2(j, at, sm) + 0.5d0 * k2_cavsig2(j, at, sm))

        ! k4 -  Second-order: Cavity (<a a>, <at at>, <at a>)
        k4_cav2(j, k, a) = -dt * (2.0d0 * kappa + i * (wl(j) + wl(k))) * (cav2(j, k, a) + k3_cav2(j, k, a)) - &
                         & dt * gkal(j) * (cavsig2(k, a, sm) + k3_cavsig2(k, a, sm)) - &
                         & dt * gkal(k) * (cavsig2(j, a, sm) + k3_cavsig2(j, a, sm))
        k4_cav2(j, k, at) = -dt * (2.0d0 * kappa - i * (wl(j) + wl(k))) * (cav2(j, k, at) + k3_cav2(j, k, at)) - &
                          & dt * CONJG(gkal(j)) * (cavsig2(k, at, sp) + k3_cavsig2(k, at, sp)) - &
                          & dt * CONJG(gkal(k)) * (cavsig2(j, at, sp) + k3_cavsig2(j, at, sp))
        k4_cav2(j, k, ata) = -dt * (2.0d0 * kappa - i * (wl(j) - wl(k))) * (cav2(j, k, ata) + k3_cav2(j, k, ata)) - &
                           & dt * CONJG(gkal(j)) * (cavsig2(k, a, sp) + k3_cavsig2(k, a, sp)) - &
                           & dt * gkal(k) * (cavsig2(j, at, sm) + k3_cavsig2(j, at, sm))

        !! Only calculate these for second-order correlation function
        ! IF (tau2_max > 0.0) THEN
          !--------------------------------------!
          !     Third-Order: Cavity and Atom     !
          !--------------------------------------!
          ! k1 - Third-order: a and \siga (<a^{2}, \sigma_{-}>, <a^{2} \sigma_{+}>, <a^{2} \sigma_{z}>)
          k1_cavsig3(j, k, a, sm) = -dt * (0.5d0 * gamma + 2.0d0 * kappa + i * (wl(j) + wl(k))) * cavsig3(j, k, a, sm) + &
                                  & dt * i * 0.5d0 * Omega * cavsig3(j, k, a, sz)
          k1_cavsig3(j, k, a, sp) = -dt * (0.5d0 * gamma + 2.0d0 * kappa + i * (wl(j) + wl(k))) * cavsig3(j, k, a, sp) - &
                                  & dt * i * 0.5d0 * Omega * cavsig3(j, k, a, sz) - &
                                  & dt * 0.5d0 * gkal(j) * cavsig2(k, a, sz) - &
                                  & dt * 0.5d0 * gkal(k) * cavsig2(j, a, sz) - &
                                  & dt * 0.5d0 * gkal(j) * cav1(k, a) - &
                                  & dt * 0.5d0 * gkal(k) * cav1(j, a)
          k1_cavsig3(j, k, a, sz) = -dt * (gamma + 2.0d0 * kappa + i * (wl(j) + wl(k))) * cavsig3(j, k, a, sz) + &
                                  & dt * i * Omega * cavsig3(j, k, a, sm) - &
                                  & dt * i * Omega * cavsig3(j, k, a, sp) - &
                                  & dt * gamma * cav2(j, k, a) + &
                                  & dt * gkal(j) * cavsig2(k, a, sm) + &
                                  & dt * gkal(k) * cavsig2(j, a, sm)
          ! k1 - Third-order: at and \siga (<at^{2}, \sigma_{-}>, <at^{2} \sigma_{+}>, <at^{2} \sigma_{z}>)
          k1_cavsig3(j, k, at, sm) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) + wl(k))) * cavsig3(j, k, at, sm) + &
                                   & dt * i * 0.5d0 * Omega * cavsig3(j, k, at, sz) - &
                                   & dt * 0.5d0 * CONJG(gkal(j)) * cavsig2(k, at, sz) - &
                                   & dt * 0.5d0 * CONJG(gkal(k)) * cavsig2(j, at, sz) - &
                                   & dt * 0.5d0 * CONJG(gkal(j)) * cav1(k, at) - &
                                   & dt * 0.5d0 * CONJG(gkal(k)) * cav1(j, at)
          k1_cavsig3(j, k, at, sp) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) + wl(k))) * cavsig3(j, k, at, sp) - &
                                   & dt * i * 0.5d0 * Omega * cavsig3(j, k, at, sz)
          k1_cavsig3(j, k, at, sz) =  -dt * (gamma + 2.0d0 * kappa - i * (wl(j) + wl(k))) * cavsig3(j, k, at, sz) + &
                                   & dt * i * Omega * cavsig3(j, k, at, sm) - &
                                   & dt * i * Omega * cavsig3(j, k, at, sp) - &
                                   & dt * gamma * cav2(j, k, at) + &
                                   & dt * CONJG(gkal(j)) * cavsig2(k, at, sp) + &
                                   & dt * CONJG(gkal(k)) * cavsig2(j, at, sp)
          ! k1 - Third-order: ata and \sig (<ata, \sigma_{-}>, <ata \sigma_{+}>, <ata \sigma_{z}>)
          k1_cavsig3(j, k, ata, sm) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) - wl(k))) * cavsig3(j, k, ata, sm) + &
                                    & dt * i * 0.5d0 * Omega * cavsig3(j, k, ata, sz) - &
                                    & dt * 0.5d0 * CONJG(gkal(j)) * cavsig2(k, a, sz) - &
                                    & dt * 0.5d0 * CONJG(gkal(j)) * cav1(k, a)
          k1_cavsig3(j, k, ata, sp) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) - wl(k))) * cavsig3(j, k, ata, sp) - &
                                    & dt * i * 0.5d0 * Omega * cavsig3(j, k, ata, sz) - &
                                    & dt * 0.5d0 * gkal(k) * cavsig2(j, at, sz) - &
                                    & dt * 0.5d0 * gkal(k) * cav1(j, at)
          k1_cavsig3(j, k, ata, sz) = -dt * (gamma + 2.0d0 * kappa - i * (wl(j) - wl(k))) * cavsig3(j, k, ata, sz) + &
                                    & dt * i * Omega * cavsig3(j, k, ata, sm) - &
                                    & dt * i * Omega * cavsig3(j, k, ata, sp) - &
                                    & dt * gamma * cav2(j, k, ata) + &
                                    & dt * CONJG(gkal(j)) * cavsig2(k, a, sp) + &
                                    & dt * gkal(k) * cavsig2(j, at, sm)

          ! k2 - Third-order: a and \siga (<a^{2}, \sigma_{-}>, <a^{2} \sigma_{+}>, <a^{2} \sigma_{z}>)
          k2_cavsig3(j, k, a, sm) = -dt * (0.5d0 * gamma + 2.0d0 * kappa + i * (wl(j) + wl(k))) * (cavsig3(j, k, a, sm) + 0.5d0 * k1_cavsig3(j, k, a, sm)) + &
                                  & dt * i * 0.5d0 * Omega * (cavsig3(j, k, a, sz) + 0.5d0 * k1_cavsig3(j, k, a, sz))
          k2_cavsig3(j, k, a, sp) = -dt * (0.5d0 * gamma + 2.0d0 * kappa + i * (wl(j) + wl(k))) * (cavsig3(j, k, a, sp) + 0.5d0 * k1_cavsig3(j, k, a, sp)) - &
                                  & dt * i * 0.5d0 * Omega * (cavsig3(j, k, a, sz) + 0.5d0 * k1_cavsig3(j, k, a, sz)) - &
                                  & dt * 0.5d0 * gkal(j) * (cavsig2(k, a, sz) + 0.5d0 * k1_cavsig2(k, a, sz)) - &
                                  & dt * 0.5d0 * gkal(k) * (cavsig2(j, a, sz) + 0.5d0 * k1_cavsig2(j, a, sz)) - &
                                  & dt * 0.5d0 * gkal(j) * (cav1(k, a) + 0.5d0 * k1_cav1(k, a)) - &
                                  & dt * 0.5d0 * gkal(k) * (cav1(j, a) + 0.5d0 * k1_cav1(j, a))
          k2_cavsig3(j, k, a, sz) = -dt * (gamma + 2.0d0 * kappa + i * (wl(j) + wl(k))) * (cavsig3(j, k, a, sz) + 0.5d0 * k1_cavsig3(j, k, a, sz)) + &
                                  & dt * i * Omega * (cavsig3(j, k, a, sm) + 0.5d0 * k1_cavsig3(j, k, a, sm)) - &
                                  & dt * i * Omega * (cavsig3(j, k, a, sp) + 0.5d0 * k1_cavsig3(j, k, a, sp)) - &
                                  & dt * gamma * (cav2(j, k, a) + 0.5d0 * k1_cav2(j, k, a)) + &
                                  & dt * gkal(j) * (cavsig2(k, a, sm) + 0.5d0 * k1_cavsig2(k, a, sm)) + &
                                  & dt * gkal(k) * (cavsig2(j, a, sm) + 0.5d0 * k1_cavsig2(j, a, sm))
          ! k2 - Third-order: at and \siga (<at^{2}, \sigma_{-}>, <at^{2} \sigma_{+}>, <at^{2} \sigma_{z}>)
          k2_cavsig3(j, k, at, sm) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) + wl(k))) * (cavsig3(j, k, at, sm) + 0.5d0 * k1_cavsig3(j, k, at, sm)) + &
                                   & dt * i * 0.5d0 * Omega * (cavsig3(j, k, at, sz) + 0.5d0 * k1_cavsig3(j, k, at, sz)) - &
                                   & dt * 0.5d0 * CONJG(gkal(j)) * (cavsig2(k, at, sz) + 0.5d0 * k1_cavsig2(k, at, sz)) - &
                                   & dt * 0.5d0 * CONJG(gkal(k)) * (cavsig2(j, at, sz) + 0.5d0 * k1_cavsig2(j, at, sz)) - &
                                   & dt * 0.5d0 * CONJG(gkal(j)) * (cav1(k, at) + 0.5d0 * k1_cav1(k, at)) - &
                                   & dt * 0.5d0 * CONJG(gkal(k)) * (cav1(j, at) + 0.5d0 * k1_cav1(j, at))
          k2_cavsig3(j, k, at, sp) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) + wl(k))) * (cavsig3(j, k, at, sp) + 0.5d0 * k1_cavsig3(j, k, at, sp)) - &
                                   & dt * i * 0.5d0 * Omega * (cavsig3(j, k, at, sz) + 0.5d0 * k1_cavsig3(j, k, at, sz))
          k2_cavsig3(j, k, at, sz) =  -dt * (gamma + 2.0d0 * kappa - i * (wl(j) + wl(k))) * (cavsig3(j, k, at, sz) + 0.5d0 * k1_cavsig3(j, k, at, sz)) + &
                                   & dt * i * Omega * (cavsig3(j, k, at, sm) + 0.5d0 * k1_cavsig3(j, k, at, sm)) - &
                                   & dt * i * Omega * (cavsig3(j, k, at, sp) + 0.5d0 * k1_cavsig3(j, k, at, sp)) - &
                                   & dt * gamma * (cav2(j, k, at) + 0.5d0 * k1_cav2(j, k, at)) + &
                                   & dt * CONJG(gkal(j)) * (cavsig2(k, at, sp) + 0.5d0 * k1_cavsig2(k, at, sp)) + &
                                   & dt * CONJG(gkal(k)) * (cavsig2(j, at, sp) + 0.5d0 * k1_cavsig2(j, at, sp))
          ! k2 - Third-order: ata and \sig (<ata, \sigma_{-}>, <ata \sigma_{+}>, <ata \sigma_{z}>)
          k2_cavsig3(j, k, ata, sm) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) - wl(k))) * (cavsig3(j, k, ata, sm) + 0.5d0 * k1_cavsig3(j, k, ata, sm)) + &
                                    & dt * i * 0.5d0 * Omega * (cavsig3(j, k, ata, sz) + 0.5d0 * k1_cavsig3(j, k, ata, sz)) - &
                                    & dt * 0.5d0 * CONJG(gkal(j)) * (cavsig2(k, a, sz) + 0.5d0 * k1_cavsig2(k, a, sz)) - &
                                    & dt * 0.5d0 * CONJG(gkal(j)) * (cav1(k, a) + 0.5d0 * k1_cav1(k, a))
          k2_cavsig3(j, k, ata, sp) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) - wl(k))) * (cavsig3(j, k, ata, sp) + 0.5d0 * k1_cavsig3(j, k, ata, sp)) - &
                                    & dt * i * 0.5d0 * Omega * (cavsig3(j, k, ata, sz) + 0.5d0 * k1_cavsig3(j, k, ata, sz)) - &
                                    & dt * 0.5d0 * gkal(k) * (cavsig2(j, at, sz) + 0.5d0 * k1_cavsig2(j, at, sz)) - &
                                    & dt * 0.5d0 * gkal(k) * (cav1(j, at) + 0.5d0 * k1_cav1(j, at))
          k2_cavsig3(j, k, ata, sz) = -dt * (gamma + 2.0d0 * kappa - i * (wl(j) - wl(k))) * (cavsig3(j, k, ata, sz) + 0.5d0 * k1_cavsig3(j, k, ata, sz)) + &
                                    & dt * i * Omega * (cavsig3(j, k, ata, sm) + 0.5d0 * k1_cavsig3(j, k, ata, sm)) - &
                                    & dt * i * Omega * (cavsig3(j, k, ata, sp) + 0.5d0 * k1_cavsig3(j, k, ata, sp)) - &
                                    & dt * gamma * (cav2(j, k, ata) + 0.5d0 * k1_cav2(j, k, ata)) + &
                                    & dt * CONJG(gkal(j)) * (cavsig2(k, a, sp) + 0.5d0 * k1_cavsig2(k, a, sp)) + &
                                    & dt * gkal(k) * (cavsig2(j, at, sm) + 0.5d0 * k1_cavsig2(j, at, sm))

          ! k3 - Third-order: a and \siga (<a^{2}, \sigma_{-}>, <a^{2} \sigma_{+}>, <a^{2} \sigma_{z}>)
          k3_cavsig3(j, k, a, sm) = -dt * (0.5d0 * gamma + 2.0d0 * kappa + i * (wl(j) + wl(k))) * (cavsig3(j, k, a, sm) + 0.5d0 * k2_cavsig3(j, k, a, sm)) + &
                                  & dt * i * 0.5d0 * Omega * (cavsig3(j, k, a, sz) + 0.5d0 * k2_cavsig3(j, k, a, sz))
          k3_cavsig3(j, k, a, sp) = -dt * (0.5d0 * gamma + 2.0d0 * kappa + i * (wl(j) + wl(k))) * (cavsig3(j, k, a, sp) + 0.5d0 * k2_cavsig3(j, k, a, sp)) - &
                                  & dt * i * 0.5d0 * Omega * (cavsig3(j, k, a, sz) + 0.5d0 * k2_cavsig3(j, k, a, sz)) - &
                                  & dt * 0.5d0 * gkal(j) * (cavsig2(k, a, sz) + 0.5d0 * k2_cavsig2(k, a, sz)) - &
                                  & dt * 0.5d0 * gkal(k) * (cavsig2(j, a, sz) + 0.5d0 * k2_cavsig2(j, a, sz)) - &
                                  & dt * 0.5d0 * gkal(j) * (cav1(k, a) + 0.5d0 * k2_cav1(k, a)) - &
                                  & dt * 0.5d0 * gkal(k) * (cav1(j, a) + 0.5d0 * k2_cav1(j, a))
          k3_cavsig3(j, k, a, sz) = -dt * (gamma + 2.0d0 * kappa + i * (wl(j) + wl(k))) * (cavsig3(j, k, a, sz) + 0.5d0 * k2_cavsig3(j, k, a, sz)) + &
                                  & dt * i * Omega * (cavsig3(j, k, a, sm) + 0.5d0 * k2_cavsig3(j, k, a, sm)) - &
                                  & dt * i * Omega * (cavsig3(j, k, a, sp) + 0.5d0 * k2_cavsig3(j, k, a, sp)) - &
                                  & dt * gamma * (cav2(j, k, a) + 0.5d0 * k2_cav2(j, k, a)) + &
                                  & dt * gkal(j) * (cavsig2(k, a, sm) + 0.5d0 * k2_cavsig2(k, a, sm)) + &
                                  & dt * gkal(k) * (cavsig2(j, a, sm) + 0.5d0 * k2_cavsig2(j, a, sm))
          ! k3 - Third-order: at and \siga (<at^{2}, \sigma_{-}>, <at^{2} \sigma_{+}>, <at^{2} \sigma_{z}>)
          k3_cavsig3(j, k, at, sm) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) + wl(k))) * (cavsig3(j, k, at, sm) + 0.5d0 * k2_cavsig3(j, k, at, sm)) + &
                                   & dt * i * 0.5d0 * Omega * (cavsig3(j, k, at, sz) + 0.5d0 * k2_cavsig3(j, k, at, sz)) - &
                                   & dt * 0.5d0 * CONJG(gkal(j)) * (cavsig2(k, at, sz) + 0.5d0 * k2_cavsig2(k, at, sz)) - &
                                   & dt * 0.5d0 * CONJG(gkal(k)) * (cavsig2(j, at, sz) + 0.5d0 * k2_cavsig2(j, at, sz)) - &
                                   & dt * 0.5d0 * CONJG(gkal(j)) * (cav1(k, at) + 0.5d0 * k2_cav1(k, at)) - &
                                   & dt * 0.5d0 * CONJG(gkal(k)) * (cav1(j, at) + 0.5d0 * k2_cav1(j, at))
          k3_cavsig3(j, k, at, sp) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) + wl(k))) * (cavsig3(j, k, at, sp) + 0.5d0 * k2_cavsig3(j, k, at, sp)) - &
                                   & dt * i * 0.5d0 * Omega * (cavsig3(j, k, at, sz) + 0.5d0 * k2_cavsig3(j, k, at, sz))
          k3_cavsig3(j, k, at, sz) =  -dt * (gamma + 2.0d0 * kappa - i * (wl(j) + wl(k))) * (cavsig3(j, k, at, sz) + 0.5d0 * k2_cavsig3(j, k, at, sz)) + &
                                   & dt * i * Omega * (cavsig3(j, k, at, sm) + 0.5d0 * k2_cavsig3(j, k, at, sm)) - &
                                   & dt * i * Omega * (cavsig3(j, k, at, sp) + 0.5d0 * k2_cavsig3(j, k, at, sp)) - &
                                   & dt * gamma * (cav2(j, k, at) + 0.5d0 * k2_cav2(j, k, at)) + &
                                   & dt * CONJG(gkal(j)) * (cavsig2(k, at, sp) + 0.5d0 * k2_cavsig2(k, at, sp)) + &
                                   & dt * CONJG(gkal(k)) * (cavsig2(j, at, sp) + 0.5d0 * k2_cavsig2(j, at, sp))
          ! k3 - Third-order: ata and \sig (<ata, \sigma_{-}>, <ata \sigma_{+}>, <ata \sigma_{z}>)
          k3_cavsig3(j, k, ata, sm) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) - wl(k))) * (cavsig3(j, k, ata, sm) + 0.5d0 * k2_cavsig3(j, k, ata, sm)) + &
                                    & dt * i * 0.5d0 * Omega * (cavsig3(j, k, ata, sz) + 0.5d0 * k2_cavsig3(j, k, ata, sz)) - &
                                    & dt * 0.5d0 * CONJG(gkal(j)) * (cavsig2(k, a, sz) + 0.5d0 * k2_cavsig2(k, a, sz)) - &
                                    & dt * 0.5d0 * CONJG(gkal(j)) * (cav1(k, a) + 0.5d0 * k2_cav1(k, a))
          k3_cavsig3(j, k, ata, sp) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) - wl(k))) * (cavsig3(j, k, ata, sp) + 0.5d0 * k2_cavsig3(j, k, ata, sp)) - &
                                    & dt * i * 0.5d0 * Omega * (cavsig3(j, k, ata, sz) + 0.5d0 * k2_cavsig3(j, k, ata, sz)) - &
                                    & dt * 0.5d0 * gkal(k) * (cavsig2(j, at, sz) + 0.5d0 * k2_cavsig2(j, at, sz)) - &
                                    & dt * 0.5d0 * gkal(k) * (cav1(j, at) + 0.5d0 * k2_cav1(j, at))
          k3_cavsig3(j, k, ata, sz) = -dt * (gamma + 2.0d0 * kappa - i * (wl(j) - wl(k))) * (cavsig3(j, k, ata, sz) + 0.5d0 * k2_cavsig3(j, k, ata, sz)) + &
                                    & dt * i * Omega * (cavsig3(j, k, ata, sm) + 0.5d0 * k2_cavsig3(j, k, ata, sm)) - &
                                    & dt * i * Omega * (cavsig3(j, k, ata, sp) + 0.5d0 * k2_cavsig3(j, k, ata, sp)) - &
                                    & dt * gamma * (cav2(j, k, ata) + 0.5d0 * k2_cav2(j, k, ata)) + &
                                    & dt * CONJG(gkal(j)) * (cavsig2(k, a, sp) + 0.5d0 * k2_cavsig2(k, a, sp)) + &
                                    & dt * gkal(k) * (cavsig2(j, at, sm) + 0.5d0 * k2_cavsig2(j, at, sm))

          ! k4 - Third-order: a and \siga (<a^{2}, \sigma_{-}>, <a^{2} \sigma_{+}>, <a^{2} \sigma_{z}>)
          k4_cavsig3(j, k, a, sm) = -dt * (0.5d0 * gamma + 2.0d0 * kappa + i * (wl(j) + wl(k))) * (cavsig3(j, k, a, sm) + k3_cavsig3(j, k, a, sm)) + &
                                  & dt * i * 0.5d0 * Omega * (cavsig3(j, k, a, sz) + k3_cavsig3(j, k, a, sz))
          k4_cavsig3(j, k, a, sp) = -dt * (0.5d0 * gamma + 2.0d0 * kappa + i * (wl(j) + wl(k))) * (cavsig3(j, k, a, sp) + k3_cavsig3(j, k, a, sp)) - &
                                  & dt * i * 0.5d0 * Omega * (cavsig3(j, k, a, sz) + k3_cavsig3(j, k, a, sz)) - &
                                  & dt * 0.5d0 * gkal(j) * (cavsig2(k, a, sz) + k3_cavsig2(k, a, sz)) - &
                                  & dt * 0.5d0 * gkal(k) * (cavsig2(j, a, sz) + k3_cavsig2(j, a, sz)) - &
                                  & dt * 0.5d0 * gkal(j) * (cav1(k, a) + k3_cav1(k, a)) - &
                                  & dt * 0.5d0 * gkal(k) * (cav1(j, a) + k3_cav1(j, a))
          k4_cavsig3(j, k, a, sz) = -dt * (gamma + 2.0d0 * kappa + i * (wl(j) + wl(k))) * (cavsig3(j, k, a, sz) + k3_cavsig3(j, k, a, sz)) + &
                                  & dt * i * Omega * (cavsig3(j, k, a, sm) + k3_cavsig3(j, k, a, sm)) - &
                                  & dt * i * Omega * (cavsig3(j, k, a, sp) + k3_cavsig3(j, k, a, sp)) - &
                                  & dt * gamma * (cav2(j, k, a) + k3_cav2(j, k, a)) + &
                                  & dt * gkal(j) * (cavsig2(k, a, sm) + k3_cavsig2(k, a, sm)) + &
                                  & dt * gkal(k) * (cavsig2(j, a, sm) + k3_cavsig2(j, a, sm))
          ! k4 - Third-order: at and \siga (<at^{2}, \sigma_{-}>, <at^{2} \sigma_{+}>, <at^{2} \sigma_{z}>)
          k4_cavsig3(j, k, at, sm) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) + wl(k))) * (cavsig3(j, k, at, sm) + k3_cavsig3(j, k, at, sm)) + &
                                   & dt * i * 0.5d0 * Omega * (cavsig3(j, k, at, sz) + k3_cavsig3(j, k, at, sz)) - &
                                   & dt * 0.5d0 * CONJG(gkal(j)) * (cavsig2(k, at, sz) + k3_cavsig2(k, at, sz)) - &
                                   & dt * 0.5d0 * CONJG(gkal(k)) * (cavsig2(j, at, sz) + k3_cavsig2(j, at, sz)) - &
                                   & dt * 0.5d0 * CONJG(gkal(j)) * (cav1(k, at) + k3_cav1(k, at)) - &
                                   & dt * 0.5d0 * CONJG(gkal(k)) * (cav1(j, at) + k3_cav1(j, at))
          k4_cavsig3(j, k, at, sp) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) + wl(k))) * (cavsig3(j, k, at, sp) + k3_cavsig3(j, k, at, sp)) - &
                                   & dt * i * 0.5d0 * Omega * (cavsig3(j, k, at, sz) + k3_cavsig3(j, k, at, sz))
          k4_cavsig3(j, k, at, sz) =  -dt * (gamma + 2.0d0 * kappa - i * (wl(j) + wl(k))) * (cavsig3(j, k, at, sz) + k3_cavsig3(j, k, at, sz)) + &
                                   & dt * i * Omega * (cavsig3(j, k, at, sm) + k3_cavsig3(j, k, at, sm)) - &
                                   & dt * i * Omega * (cavsig3(j, k, at, sp) + k3_cavsig3(j, k, at, sp)) - &
                                   & dt * gamma * (cav2(j, k, at) + k3_cav2(j, k, at)) + &
                                   & dt * CONJG(gkal(j)) * (cavsig2(k, at, sp) + k3_cavsig2(k, at, sp)) + &
                                   & dt * CONJG(gkal(k)) * (cavsig2(j, at, sp) + k3_cavsig2(j, at, sp))
          ! k4 - Third-order: ata and \sig (<ata, \sigma_{-}>, <ata \sigma_{+}>, <ata \sigma_{z}>)
          k4_cavsig3(j, k, ata, sm) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) - wl(k))) * (cavsig3(j, k, ata, sm) + k3_cavsig3(j, k, ata, sm)) + &
                                    & dt * i * 0.5d0 * Omega * (cavsig3(j, k, ata, sz) + k3_cavsig3(j, k, ata, sz)) - &
                                    & dt * 0.5d0 * CONJG(gkal(j)) * (cavsig2(k, a, sz) + k3_cavsig2(k, a, sz)) - &
                                    & dt * 0.5d0 * CONJG(gkal(j)) * (cav1(k, a) + k3_cav1(k, a))
          k4_cavsig3(j, k, ata, sp) = -dt * (0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) - wl(k))) * (cavsig3(j, k, ata, sp) + k3_cavsig3(j, k, ata, sp)) - &
                                    & dt * i * 0.5d0 * Omega * (cavsig3(j, k, ata, sz) + k3_cavsig3(j, k, ata, sz)) - &
                                    & dt * 0.5d0 * gkal(k) * (cavsig2(j, at, sz) + k3_cavsig2(j, at, sz)) - &
                                    & dt * 0.5d0 * gkal(k) * (cav1(j, at) + k3_cav1(j, at))
          k4_cavsig3(j, k, ata, sz) = -dt * (gamma + 2.0d0 * kappa - i * (wl(j) - wl(k))) * (cavsig3(j, k, ata, sz) + k3_cavsig3(j, k, ata, sz)) + &
                                    & dt * i * Omega * (cavsig3(j, k, ata, sm) + k3_cavsig3(j, k, ata, sm)) - &
                                    & dt * i * Omega * (cavsig3(j, k, ata, sp) + k3_cavsig3(j, k, ata, sp)) - &
                                    & dt * gamma * (cav2(j, k, ata) + k3_cav2(j, k, ata)) + &
                                    & dt * CONJG(gkal(j)) * (cavsig2(k, a, sp) + k3_cavsig2(k, a, sp)) + &
                                    & dt * gkal(k) * (cavsig2(j, at, sm) + k3_cavsig2(j, at, sm))
          ! Cycle through modes
          DO l = -N, N
            !-----------------------------!
            !     Third-Order: Cavity     !
            !-----------------------------!
            ! k1 - Third-order: Cavity (<a^{\dagger} a^{2}>, <a^{\dagger 2} a>)
            k1_cav3(j, k, l, a) = -dt * (3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * cav3(j, k, l, a) - &
                                & dt * CONJG(gkal(j)) * cavsig3(k, l, a, sp) - &
                                & dt * gkal(k) * cavsig3(j, l, ata, sm) - &
                                & dt * gkal(l) * cavsig3(j, k, ata, sm)
            k1_cav3(j, k, l, at) = -dt * (3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * cav3(j, k, l, at) - &
                                 & dt * CONJG(gkal(j)) * cavsig3(k, l, ata, sp) - &
                                 & dt * CONJG(gkal(k)) * cavsig3(j, l, ata, sp) - &
                                 & dt * gkal(l) * cavsig3(j, k, at, sm)
            ! k2 - Third-order: Cavity (<a^{\dagger} a^{2}>, <a^{\dagger 2} a>)
            k2_cav3(j, k, l, a) = -dt * (3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cav3(j, k, l, a) + 0.5d0 * k1_cav3(j, k, l, a)) - &
                                & dt * CONJG(gkal(j)) * (cavsig3(k, l, a, sp) + 0.5d0 * k1_cavsig3(k, l, a, sp)) - &
                                & dt * gkal(k) * (cavsig3(j, l, ata, sm) + 0.5d0 * k1_cavsig3(j, l, ata, sm)) - &
                                & dt * gkal(l) * (cavsig3(j, k, ata, sm) + 0.5d0 * k1_cavsig3(j, k, ata, sm))
            k2_cav3(j, k, l, at) = -dt * (3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cav3(j, k, l, at) + 0.5d0 * k1_cav3(j, k, l, at)) - &
                                 & dt * CONJG(gkal(j)) * (cavsig3(k, l, ata, sp) + 0.5d0 * k1_cavsig3(k, l, ata, sp)) - &
                                 & dt * CONJG(gkal(k)) * (cavsig3(j, l, ata, sp) + 0.5d0 * k1_cavsig3(j, l, ata, sp)) - &
                                 & dt * gkal(l) * (cavsig3(j, k, at, sm) + 0.5d0 * k1_cavsig3(j, k, at, sm))

            ! k3 - Third-order: Cavity (<a^{\dagger} a^{2}>, <a^{\dagger 2} a>)
            k3_cav3(j, k, l, a) = -dt * (3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cav3(j, k, l, a) + 0.5d0 * k2_cav3(j, k, l, a)) - &
                                & dt * CONJG(gkal(j)) * (cavsig3(k, l, a, sp) + 0.5d0 * k2_cavsig3(k, l, a, sp)) - &
                                & dt * gkal(k) * (cavsig3(j, l, ata, sm) + 0.5d0 * k2_cavsig3(j, l, ata, sm)) - &
                                & dt * gkal(l) * (cavsig3(j, k, ata, sm) + 0.5d0 * k2_cavsig3(j, k, ata, sm))
            k3_cav3(j, k, l, at) = -dt * (3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cav3(j, k, l, at) + 0.5d0 * k2_cav3(j, k, l, at)) - &
                                 & dt * CONJG(gkal(j)) * (cavsig3(k, l, ata, sp) + 0.5d0 * k2_cavsig3(k, l, ata, sp)) - &
                                 & dt * CONJG(gkal(k)) * (cavsig3(j, l, ata, sp) + 0.5d0 * k2_cavsig3(j, l, ata, sp)) - &
                                 & dt * gkal(l) * (cavsig3(j, k, at, sm) + 0.5d0 * k2_cavsig3(j, k, at, sm))

            ! k4 - Third-order: Cavity (<a^{\dagger} a^{2}>, <a^{\dagger 2} a>)
            k4_cav3(j, k, l, a) = -dt * (3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cav3(j, k, l, a) + k3_cav3(j, k, l, a)) - &
                                & dt * CONJG(gkal(j)) * (cavsig3(k, l, a, sp) + k3_cavsig3(k, l, a, sp)) - &
                                & dt * gkal(k) * (cavsig3(j, l, ata, sm) + k3_cavsig3(j, l, ata, sm)) - &
                                & dt * gkal(l) * (cavsig3(j, k, ata, sm) + k3_cavsig3(j, k, ata, sm))
            k4_cav3(j, k, l, at) = -dt * (3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cav3(j, k, l, at) + k3_cav3(j, k, l, at)) - &
                                 & dt * CONJG(gkal(j)) * (cavsig3(k, l, ata, sp) + k3_cavsig3(k, l, ata, sp)) - &
                                 & dt * CONJG(gkal(k)) * (cavsig3(j, l, ata, sp) + k3_cavsig3(j, l, ata, sp)) - &
                                 & dt * gkal(l) * (cavsig3(j, k, at, sm) + k3_cavsig3(j, k, at, sm))

            !---------------------------------------!
            !     Fourth-Order: Cavity and Atom     !
            !---------------------------------------!
            ! Fourth-order: a and sigma (<a^{\dagger} a a \sigma_{-+z}>)
            k1_cavsig4(j, k, l, a, sm) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * cavsig4(j, k, l, a, sm) + &
                                       & dt * i * 0.5d0 * Omega * cavsig4(j, k, l, a, sz) - &
                                       & dt * 0.5d0 * CONJG(gkal(j)) * cavsig3(k, l, a, sz) - &
                                       & dt * 0.5d0 * CONJG(gkal(j)) * cav2(k, l, a)
            k1_cavsig4(j, k, l, a, sp) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * cavsig4(j, k, l, a, sp) - &
                                       & dt * i * 0.5d0 * Omega * cavsig4(j, k, l, a, sz) - &
                                       & dt * 0.5d0 * gkal(k) * cavsig3(j, l, ata, sz) - &
                                       & dt * 0.5d0 * gkal(k) * cav2(j, l, ata) - &
                                       & dt * 0.5d0 * gkal(l) * cavsig3(j, k, ata, sz) - &
                                       & dt * 0.5d0 * gkal(l) * cav2(j, k, ata)
            k1_cavsig4(j, k, l, a, sz) = -dt * (gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * cavsig4(j, k, l, a, sz) + &
                                       & dt * i * Omega * cavsig4(j, k, l, a, sm) - &
                                       & dt * i * Omega * cavsig4(j, k, l, a, sp) - &
                                       & dt * gamma * cav3(j, k, l, a) + &
                                       & dt * CONJG(gkal(j)) * cavsig3(k, l, a, sp) + &
                                       & dt * gkal(k) * cavsig3(j, l, ata, sm) + &
                                       & dt * gkal(l) * cavsig3(j, k, ata, sm)
            ! k1 - Fourth-order: at and sigma (<a^{\dagger} a^{\dagger} a \sigma_{-+z}>)
            k1_cavsig4(j, k, l, at, sm) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * cavsig4(j, k, l, at, sm) + &
                                        & dt * i * 0.5d0 * Omega * cavsig4(j, k, l, at, sz) - &
                                        & dt * 0.5d0 * CONJG(gkal(j)) * cavsig3(k, l, ata, sz) - &
                                        & dt * 0.5d0 * CONJG(gkal(j)) * cav2(k, l, ata) - &
                                        & dt * 0.5d0 * CONJG(gkal(k)) * cavsig3(j, l, ata, sz) - &
                                        & dt * 0.5d0 * CONJG(gkal(k)) * cav2(j, l, ata)
            k1_cavsig4(j, k, l, at, sp) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * cavsig4(j, k, l, at, sp) - &
                                        & dt * i * 0.5d0 * Omega * cavsig4(j, k, l, at, sz) - &
                                        & dt * 0.5d0 * gkal(l) * cavsig3(j, k, at, sz) - &
                                        & dt * 0.5d0 * gkal(l) * cav2(j, k, at)
            k1_cavsig4(j, k, l, at, sz) = -dt * (gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * cavsig4(j, k, l, at, sz) + &
                                        & dt * i * Omega * cavsig4(j, k, l, at, sm) - &
                                        & dt * i * Omega * cavsig4(j, k, l, at, sp) - &
                                        & dt * gamma * cav3(j, k, l, at) + &
                                        & dt * CONJG(gkal(j)) * cavsig3(k, l, ata, sp) + &
                                        & dt * CONJG(gkal(k)) * cavsig3(j, l, ata, sp) + &
                                        & dt * gkal(l) * cavsig3(j, k, at, sm)

            ! k2 - Fourth-order: a and sigma (<a^{\dagger} a a \sigma_{-+z}>)
            k2_cavsig4(j, k, l, a, sm) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cavsig4(j, k, l, a, sm) + 0.5d0 * k1_cavsig4(j, k, l, a, sm)) + &
                                       & dt * i * 0.5d0 * Omega * (cavsig4(j, k, l, a, sz) + 0.5d0 * k1_cavsig4(j, k, l, a, sz)) - &
                                       & dt * 0.5d0 * CONJG(gkal(j)) * (cavsig3(k, l, a, sz) + 0.5d0 * k1_cavsig3(k, l, a, sz)) - &
                                       & dt * 0.5d0 * CONJG(gkal(j)) * (cav2(k, l, a) + 0.5d0 * k1_cav2(k, l, a))
            k2_cavsig4(j, k, l, a, sp) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cavsig4(j, k, l, a, sp) + 0.5d0 * k1_cavsig4(j, k, l, a, sp)) - &
                                       & dt * i * 0.5d0 * Omega * (cavsig4(j, k, l, a, sz) + 0.5d0 * k1_cavsig4(j, k, l, a, sz)) - &
                                       & dt * 0.5d0 * gkal(k) * (cavsig3(j, l, ata, sz) + 0.5d0 * k1_cavsig3(j, l, ata, sz)) - &
                                       & dt * 0.5d0 * gkal(k) * (cav2(j, l, ata) + 0.5d0 * k1_cav2(j, l, ata)) - &
                                       & dt * 0.5d0 * gkal(l) * (cavsig3(j, k, ata, sz) + 0.5d0 * k1_cavsig3(j, k, ata, sz)) - &
                                       & dt * 0.5d0 * gkal(l) * (cav2(j, k, ata) + 0.5d0 * k1_cav2(j, k, ata))
            k2_cavsig4(j, k, l, a, sz) = -dt * (gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cavsig4(j, k, l, a, sz) + 0.5d0 * k1_cavsig4(j, k, l, a, sz)) + &
                                       & dt * i * Omega * (cavsig4(j, k, l, a, sm) + 0.5d0 * k1_cavsig4(j, k, l, a, sm)) - &
                                       & dt * i * Omega * (cavsig4(j, k, l, a, sp) + 0.5d0 * k1_cavsig4(j, k, l, a, sp)) - &
                                       & dt * gamma * (cav3(j, k, l, a) + 0.5d0 * k1_cav3(j, k, l, a)) + &
                                       & dt * CONJG(gkal(j)) * (cavsig3(k, l, a, sp) + 0.5d0 * k1_cavsig3(k, l, a, sp)) + &
                                       & dt * gkal(k) * (cavsig3(j, l, ata, sm) + 0.5d0 * k1_cavsig3(j, l, ata, sm)) + &
                                       & dt * gkal(l) * (cavsig3(j, k, ata, sm) + 0.5d0 * k1_cavsig3(j, k, ata, sm))
            ! k2 - Fourth-order: at and sigma (<a^{\dagger} a^{\dagger} a \sigma_{-+z}>)
            k2_cavsig4(j, k, l, at, sm) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cavsig4(j, k, l, at, sm) + 0.5d0 * k1_cavsig4(j, k, l, at, sm)) + &
                                        & dt * i * 0.5d0 * Omega * (cavsig4(j, k, l, at, sz) + 0.5d0 * k1_cavsig4(j, k, l, at, sz)) - &
                                        & dt * 0.5d0 * CONJG(gkal(j)) * (cavsig3(k, l, ata, sz) + 0.5d0 * k1_cavsig3(k, l, ata, sz)) - &
                                        & dt * 0.5d0 * CONJG(gkal(j)) * (cav2(k, l, ata) + 0.5d0 * k1_cav2(k, l, ata)) - &
                                        & dt * 0.5d0 * CONJG(gkal(k)) * (cavsig3(j, l, ata, sz) + 0.5d0 * k1_cavsig3(j, l, ata, sz)) - &
                                        & dt * 0.5d0 * CONJG(gkal(k)) * (cav2(j, l, ata) + 0.5d0 * k1_cav2(j, l, ata))
            k2_cavsig4(j, k, l, at, sp) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cavsig4(j, k, l, at, sp) + 0.5d0 * k1_cavsig4(j, k, l, at, sp)) - &
                                        & dt * i * 0.5d0 * Omega * (cavsig4(j, k, l, at, sz) + 0.5d0 * k1_cavsig4(j, k, l, at, sz)) - &
                                        & dt * 0.5d0 * gkal(l) * (cavsig3(j, k, at, sz) + 0.5d0 * k1_cavsig3(j, k, at, sz)) - &
                                        & dt * 0.5d0 * gkal(l) * (cav2(j, k, at) + 0.5d0 * k1_cav2(j, k, at))
            k2_cavsig4(j, k, l, at, sz) = -dt * (gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cavsig4(j, k, l, at, sz) + 0.5d0 * k1_cavsig4(j, k, l, at, sz)) + &
                                        & dt * i * Omega * (cavsig4(j, k, l, at, sm) + 0.5d0 * k1_cavsig4(j, k, l, at, sm)) - &
                                        & dt * i * Omega * (cavsig4(j, k, l, at, sp) + 0.5d0 * k1_cavsig4(j, k, l, at, sp)) - &
                                        & dt * gamma * (cav3(j, k, l, at) + 0.5d0 * k1_cav3(j, k, l, at)) + &
                                        & dt * CONJG(gkal(j)) * (cavsig3(k, l, ata, sp) + 0.5d0 * k1_cavsig3(k, l, ata, sp)) + &
                                        & dt * CONJG(gkal(k)) * (cavsig3(j, l, ata, sp) + 0.5d0 * k1_cavsig3(j, l, ata, sp)) + &
                                        & dt * gkal(l) * (cavsig3(j, k, at, sm) + 0.5d0 * k1_cavsig3(j, k, at, sm))

            ! k3 - Fourth-order: a and sigma (<a^{\dagger} a a \sigma_{-+z}>)
            k3_cavsig4(j, k, l, a, sm) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cavsig4(j, k, l, a, sm) + 0.5d0 * k2_cavsig4(j, k, l, a, sm)) + &
                                       & dt * i * 0.5d0 * Omega * (cavsig4(j, k, l, a, sz) + 0.5d0 * k2_cavsig4(j, k, l, a, sz)) - &
                                       & dt * 0.5d0 * CONJG(gkal(j)) * (cavsig3(k, l, a, sz) + 0.5d0 * k2_cavsig3(k, l, a, sz)) - &
                                       & dt * 0.5d0 * CONJG(gkal(j)) * (cav2(k, l, a) + 0.5d0 * k2_cav2(k, l, a))
            k3_cavsig4(j, k, l, a, sp) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cavsig4(j, k, l, a, sp) + 0.5d0 * k2_cavsig4(j, k, l, a, sp)) - &
                                       & dt * i * 0.5d0 * Omega * (cavsig4(j, k, l, a, sz) + 0.5d0 * k2_cavsig4(j, k, l, a, sz)) - &
                                       & dt * 0.5d0 * gkal(k) * (cavsig3(j, l, ata, sz) + 0.5d0 * k2_cavsig3(j, l, ata, sz)) - &
                                       & dt * 0.5d0 * gkal(k) * (cav2(j, l, ata) + 0.5d0 * k2_cav2(j, l, ata)) - &
                                       & dt * 0.5d0 * gkal(l) * (cavsig3(j, k, ata, sz) + 0.5d0 * k2_cavsig3(j, k, ata, sz)) - &
                                       & dt * 0.5d0 * gkal(l) * (cav2(j, k, ata) + 0.5d0 * k2_cav2(j, k, ata))
            k3_cavsig4(j, k, l, a, sz) = -dt * (gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cavsig4(j, k, l, a, sz) + 0.5d0 * k2_cavsig4(j, k, l, a, sz)) + &
                                       & dt * i * Omega * (cavsig4(j, k, l, a, sm) + 0.5d0 * k2_cavsig4(j, k, l, a, sm)) - &
                                       & dt * i * Omega * (cavsig4(j, k, l, a, sp) + 0.5d0 * k2_cavsig4(j, k, l, a, sp)) - &
                                       & dt * gamma * (cav3(j, k, l, a) + 0.5d0 * k2_cav3(j, k, l, a)) + &
                                       & dt * CONJG(gkal(j)) * (cavsig3(k, l, a, sp) + 0.5d0 * k2_cavsig3(k, l, a, sp)) + &
                                       & dt * gkal(k) * (cavsig3(j, l, ata, sm) + 0.5d0 * k2_cavsig3(j, l, ata, sm)) + &
                                       & dt * gkal(l) * (cavsig3(j, k, ata, sm) + 0.5d0 * k2_cavsig3(j, k, ata, sm))
            ! k3 - Fourth-order: at and sigma (<a^{\dagger} a^{\dagger} a \sigma_{-+z}>)
            k3_cavsig4(j, k, l, at, sm) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cavsig4(j, k, l, at, sm) + 0.5d0 * k2_cavsig4(j, k, l, at, sm)) + &
                                        & dt * i * 0.5d0 * Omega * (cavsig4(j, k, l, at, sz) + 0.5d0 * k2_cavsig4(j, k, l, at, sz)) - &
                                        & dt * 0.5d0 * CONJG(gkal(j)) * (cavsig3(k, l, ata, sz) + 0.5d0 * k2_cavsig3(k, l, ata, sz)) - &
                                        & dt * 0.5d0 * CONJG(gkal(j)) * (cav2(k, l, ata) + 0.5d0 * k2_cav2(k, l, ata)) - &
                                        & dt * 0.5d0 * CONJG(gkal(k)) * (cavsig3(j, l, ata, sz) + 0.5d0 * k2_cavsig3(j, l, ata, sz)) - &
                                        & dt * 0.5d0 * CONJG(gkal(k)) * (cav2(j, l, ata) + 0.5d0 * k2_cav2(j, l, ata))
            k3_cavsig4(j, k, l, at, sp) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cavsig4(j, k, l, at, sp) + 0.5d0 * k2_cavsig4(j, k, l, at, sp)) - &
                                        & dt * i * 0.5d0 * Omega * (cavsig4(j, k, l, at, sz) + 0.5d0 * k2_cavsig4(j, k, l, at, sz)) - &
                                        & dt * 0.5d0 * gkal(l) * (cavsig3(j, k, at, sz) + 0.5d0 * k2_cavsig3(j, k, at, sz)) - &
                                        & dt * 0.5d0 * gkal(l) * (cav2(j, k, at) + 0.5d0 * k2_cav2(j, k, at))
            k3_cavsig4(j, k, l, at, sz) = -dt * (gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cavsig4(j, k, l, at, sz) + 0.5d0 * k2_cavsig4(j, k, l, at, sz)) + &
                                        & dt * i * Omega * (cavsig4(j, k, l, at, sm) + 0.5d0 * k2_cavsig4(j, k, l, at, sm)) - &
                                        & dt * i * Omega * (cavsig4(j, k, l, at, sp) + 0.5d0 * k2_cavsig4(j, k, l, at, sp)) - &
                                        & dt * gamma * (cav3(j, k, l, at) + 0.5d0 * k2_cav3(j, k, l, at)) + &
                                        & dt * CONJG(gkal(j)) * (cavsig3(k, l, ata, sp) + 0.5d0 * k2_cavsig3(k, l, ata, sp)) + &
                                        & dt * CONJG(gkal(k)) * (cavsig3(j, l, ata, sp) + 0.5d0 * k2_cavsig3(j, l, ata, sp)) + &
                                        & dt * gkal(l) * (cavsig3(j, k, at, sm) + 0.5d0 * k2_cavsig3(j, k, at, sm))

            ! k4 - Fourth-order: a and sigma (<a^{\dagger} a a \sigma_{-+z}>)
            k4_cavsig4(j, k, l, a, sm) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cavsig4(j, k, l, a, sm) + k3_cavsig4(j, k, l, a, sm)) + &
                                       & dt * i * 0.5d0 * Omega * (cavsig4(j, k, l, a, sz) + k3_cavsig4(j, k, l, a, sz)) - &
                                       & dt * 0.5d0 * CONJG(gkal(j)) * (cavsig3(k, l, a, sz) + k3_cavsig3(k, l, a, sz)) - &
                                       & dt * 0.5d0 * CONJG(gkal(j)) * (cav2(k, l, a) + k3_cav2(k, l, a))
            k4_cavsig4(j, k, l, a, sp) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cavsig4(j, k, l, a, sp) + k3_cavsig4(j, k, l, a, sp)) - &
                                       & dt * i * 0.5d0 * Omega * (cavsig4(j, k, l, a, sz) + k3_cavsig4(j, k, l, a, sz)) - &
                                       & dt * 0.5d0 * gkal(k) * (cavsig3(j, l, ata, sz) + k3_cavsig3(j, l, ata, sz)) - &
                                       & dt * 0.5d0 * gkal(k) * (cav2(j, l, ata) + k3_cav2(j, l, ata)) - &
                                       & dt * 0.5d0 * gkal(l) * (cavsig3(j, k, ata, sz) + k3_cavsig3(j, k, ata, sz)) - &
                                       & dt * 0.5d0 * gkal(l) * (cav2(j, k, ata) + k3_cav2(j, k, ata))
            k4_cavsig4(j, k, l, a, sz) = -dt * (gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l))) * (cavsig4(j, k, l, a, sz) + k3_cavsig4(j, k, l, a, sz)) + &
                                       & dt * i * Omega * (cavsig4(j, k, l, a, sm) + k3_cavsig4(j, k, l, a, sm)) - &
                                       & dt * i * Omega * (cavsig4(j, k, l, a, sp) + k3_cavsig4(j, k, l, a, sp)) - &
                                       & dt * gamma * (cav3(j, k, l, a) + k3_cav3(j, k, l, a)) + &
                                       & dt * CONJG(gkal(j)) * (cavsig3(k, l, a, sp) + k3_cavsig3(k, l, a, sp)) + &
                                       & dt * gkal(k) * (cavsig3(j, l, ata, sm) + k3_cavsig3(j, l, ata, sm)) + &
                                       & dt * gkal(l) * (cavsig3(j, k, ata, sm) + k3_cavsig3(j, k, ata, sm))
            ! k4 - Fourth-order: at and sigma (<a^{\dagger} a^{\dagger} a \sigma_{-+z}>)
            k4_cavsig4(j, k, l, at, sm) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cavsig4(j, k, l, at, sm) + k3_cavsig4(j, k, l, at, sm)) + &
                                        & dt * i * 0.5d0 * Omega * (cavsig4(j, k, l, at, sz) + k3_cavsig4(j, k, l, at, sz)) - &
                                        & dt * 0.5d0 * CONJG(gkal(j)) * (cavsig3(k, l, ata, sz) + k3_cavsig3(k, l, ata, sz)) - &
                                        & dt * 0.5d0 * CONJG(gkal(j)) * (cav2(k, l, ata) + k3_cav2(k, l, ata)) - &
                                        & dt * 0.5d0 * CONJG(gkal(k)) * (cavsig3(j, l, ata, sz) + k3_cavsig3(j, l, ata, sz)) - &
                                        & dt * 0.5d0 * CONJG(gkal(k)) * (cav2(j, l, ata) + k3_cav2(j, l, ata))
            k4_cavsig4(j, k, l, at, sp) = -dt * (0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cavsig4(j, k, l, at, sp) + k3_cavsig4(j, k, l, at, sp)) - &
                                        & dt * i * 0.5d0 * Omega * (cavsig4(j, k, l, at, sz) + k3_cavsig4(j, k, l, at, sz)) - &
                                        & dt * 0.5d0 * gkal(l) * (cavsig3(j, k, at, sz) + k3_cavsig3(j, k, at, sz)) - &
                                        & dt * 0.5d0 * gkal(l) * (cav2(j, k, at) + k3_cav2(j, k, at))
            k4_cavsig4(j, k, l, at, sz) = -dt * (gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l))) * (cavsig4(j, k, l, at, sz) + k3_cavsig4(j, k, l, at, sz)) + &
                                        & dt * i * Omega * (cavsig4(j, k, l, at, sm) + k3_cavsig4(j, k, l, at, sm)) - &
                                        & dt * i * Omega * (cavsig4(j, k, l, at, sp) + k3_cavsig4(j, k, l, at, sp)) - &
                                        & dt * gamma * (cav3(j, k, l, at) + k3_cav3(j, k, l, at)) + &
                                        & dt * CONJG(gkal(j)) * (cavsig3(k, l, ata, sp) + k3_cavsig3(k, l, ata, sp)) + &
                                        & dt * CONJG(gkal(k)) * (cavsig3(j, l, ata, sp) + k3_cavsig3(j, l, ata, sp)) + &
                                        & dt * gkal(l) * (cavsig3(j, k, at, sm) + k3_cavsig3(j, k, at, sm))
            ! Cycle through modes
            DO m = -N, N
              !------------------------------!
              !     Fourth-Order: Cavity     !
              !------------------------------!
              ! k1 - Fourth-Order: Cavity (<at at a a>)
              k1_cav4(j, k, l, m) = -dt * (4.0d0 * kappa - i * (wl(j) + wl(k) - wl(l) - wl(m))) * cav4(j, k, l, m) - &
                                  & dt * CONJG(gkal(j)) * cavsig4(k, l, m, a, sp) - &
                                  & dt * CONJG(gkal(k)) * cavsig4(j, l, m, a, sp) - &
                                  & dt * gkal(l) * cavsig4(j, k, m, at, sm) - &
                                  & dt * gkal(m) * cavsig4(j, k, l, at, sm)

              ! k2 - Fourth-Order: Cavity (<at at a a>)
              k2_cav4(j, k, l, m) = -dt * (4.0d0 * kappa - i * (wl(j) + wl(k) - wl(l) - wl(m))) * (cav4(j, k, l, m) + 0.5d0 * k1_cav4(j, k, l, m)) - &
                                  & dt * CONJG(gkal(j)) * (cavsig4(k, l, m, a, sp) + 0.5d0 * k1_cavsig4(k, l, m, a, sp)) - &
                                  & dt * CONJG(gkal(k)) * (cavsig4(j, l, m, a, sp) + 0.5d0 * k1_cavsig4(j, l, m, a, sp)) - &
                                  & dt * gkal(l) * (cavsig4(j, k, m, at, sm) + 0.5d0 * k1_cavsig4(j, k, m, at, sm)) - &
                                  & dt * gkal(m) * (cavsig4(j, k, l, at, sm) + 0.5d0 * k1_cavsig4(j, k, l, at, sm))

              ! k3 - Fourth-Order: Cavity (<at at a a>)
              k3_cav4(j, k, l, m) = -dt * (4.0d0 * kappa - i * (wl(j) + wl(k) - wl(l) - wl(m))) * (cav4(j, k, l, m) + 0.5d0 * k2_cav4(j, k, l, m)) - &
                                  & dt * CONJG(gkal(j)) * (cavsig4(k, l, m, a, sp) + 0.5d0 * k2_cavsig4(k, l, m, a, sp)) - &
                                  & dt * CONJG(gkal(k)) * (cavsig4(j, l, m, a, sp) + 0.5d0 * k2_cavsig4(j, l, m, a, sp)) - &
                                  & dt * gkal(l) * (cavsig4(j, k, m, at, sm) + 0.5d0 * k2_cavsig4(j, k, m, at, sm)) - &
                                  & dt * gkal(m) * (cavsig4(j, k, l, at, sm) + 0.5d0 * k2_cavsig4(j, k, l, at, sm))

              ! k4 - Fourth-Order: Cavity (<at at a a>)
              k4_cav4(j, k, l, m) = -dt * (4.0d0 * kappa - i * (wl(j) + wl(k) - wl(l) - wl(m))) * (cav4(j, k, l, m) + k3_cav4(j, k, l, m)) - &
                                  & dt * CONJG(gkal(j)) * (cavsig4(k, l, m, a, sp) + k3_cavsig4(k, l, m, a, sp)) - &
                                  & dt * CONJG(gkal(k)) * (cavsig4(j, l, m, a, sp) + k3_cavsig4(j, l, m, a, sp)) - &
                                  & dt * gkal(l) * (cavsig4(j, k, m, at, sm) + k3_cavsig4(j, k, m, at, sm)) - &
                                  & dt * gkal(m) * (cavsig4(j, k, l, at, sm) + k3_cavsig4(j, k, l, at, sm))
            END DO
          END DO
        ! END IF

      END DO
    END DO

    !---------------------------------!
    !     Update operator moments     !
    !---------------------------------!
    ! Update first-order moments
    bloch = bloch + xis * (k1_bloch + 2.0d0 * (k2_bloch + k3_bloch) + k4_bloch)
    cav1 = cav1 + xis * (k1_cav1 + 2.0d0 * (k2_cav1 + k3_cav1) + k4_cav1)
    ! Update second-order moments
    cav2 = cav2 + xis * (k1_cav2 + 2.0d0 * (k2_cav2 + k3_cav2) + k4_cav2)
    cavsig2 = cavsig2 + xis * (k1_cavsig2 + 2.0d0 * (k2_cavsig2 + k3_cavsig2) + k4_cavsig2)
    ! Update third-order moments
    cav3 = cav3 + xis * (k1_cav3 + 2.0d0 * (k2_cav3 + k3_cav3) + k4_cav3)
    cavsig3 = cavsig3 + xis * (k1_cavsig3 + 2.0d0 * (k2_cavsig3 + k3_cavsig3) + k4_cavsig3)
    ! ! Update fourt-order moments
    cav4 = cav4 + xis * (k1_cav4 + 2.0d0 * (k2_cav4 + k3_cav4) + k4_cav4)
    cavsig4 = cavsig4 + xis * (k1_cavsig4 + 2.0d0 * (k2_cavsig4 + k3_cavsig4) + k4_cavsig4)

    ! Check percentage
    IF (progress_bar .EQV. .TRUE.) THEN
      IF (MOD(t, ten_percent) == 0 .AND. t /= 0) THEN
        CALL CPU_TIME(loop_check_time)
        percentage = NINT((100.0 * t) / (1.0 * t_steps))
        loop_run_time = loop_check_time - loop_start_time
        loop_remaining_time = ((100.0 * (loop_check_time - loop_start_time)) / (1.0 * percentage)) - loop_run_time
        WRITE(*, FMT_ss) percentage, "%. Run time (steady state): ", loop_run_time, "s. Est. time left: ", &
                       & loop_remaining_time, "s"
      END IF
    END IF

    ! Close time integration DO loop
  END DO

  ! Close state file
  CLOSE(2)
END IF

!------------------------------------------------------------------------------!
!                               STEADY STATES                                  !
!------------------------------------------------------------------------------!

! Calculate steady state expressions
! First-order: Bloch equations
bloch_ss(sm) = -i * gamma * Omega / ((2.0d0 * Omega ** 2) + (gamma ** 2))
bloch_ss(sp) = i * gamma * Omega / ((2.0d0 * Omega ** 2) + (gamma ** 2))
bloch_ss(sz) = (gamma ** 2) / ((2.0d0 * Omega ** 2) + (gamma ** 2))
! Cycle through modes
DO j = -N, N
  !-----------------------------!
  !     First-order: Cavity     !
  !-----------------------------!
  cav1_ss(j, a) = - gkal(j) * bloch_ss(sm) / (kappa + i * wl(j))
  cav1_ss(j, at) = - CONJG(gkal(j)) * bloch_ss(sp) / (kappa - i * wl(j))
  !---------------------------------------!
  !     Second-order: Cavity and atom     !
  !---------------------------------------!
  ! <a \sigma>
  alpha = -(0.5d0 * gamma + kappa + i * wl(j))
  beta = -(gamma + kappa + i * wl(j))
  ! Set matrix elements
  M_inverse = INVERT_MATRIX(alpha, beta)
  B_ss = 0.0d0
  B_ss(1) = 0.0d0
  B_ss(2) = -0.5 * gkal(j) * (bloch_ss(sz) + 1.0d0)
  B_ss(3) = -gamma * cav1_ss(j, a) + &
          & gkal(j) * bloch_ss(sm)
  ! Matrix multiplication
  cavsig2_ss(j, a, sm) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
  cavsig2_ss(j, a, sp) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
  cavsig2_ss(j, a, sz) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))

  ! <at \sigma>
  alpha = -(0.5d0 * gamma + kappa - i * wl(j))
  beta = -(gamma + kappa - i * wl(j))
  ! Set matrix elements
  M_inverse = INVERT_MATRIX(alpha, beta)
  B_ss = 0.0d0
  B_ss(1) = -0.5 * CONJG(gkal(j)) * (bloch_ss(sz) + 1.0d0)
  B_ss(2) = 0.0d0
  B_ss(3) = -gamma * cav1_ss(j, at) + &
          & CONJG(gkal(j)) * bloch_ss(sp)
  ! Matrix multiplication
  cavsig2_ss(j, at, sm) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
  cavsig2_ss(j, at, sp) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
  cavsig2_ss(j, at, sz) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))
END DO

! Cycle through modes
DO j = -N, N
  DO k = -N, N
    !------------------------------!
    !     Second-order: Cavity     !
    !------------------------------!
    ! <a a>
    cav2_ss(j, k, a) = -(gkal(j) * cavsig2_ss(k, a, sm) + gkal(k) * cavsig2_ss(j, a, sm)) / &
                     & (2.0d0 * kappa + i * (wl(j) + wl(k)))
    ! <at at>
    cav2_ss(j, k, at) = -(CONJG(gkal(j)) * cavsig2_ss(k, at, sp) + CONJG(gkal(k)) * cavsig2_ss(j, at, sp)) / &
                      & (2.0d0 * kappa - i * (wl(j) + wl(k)))
    ! <at a>
    cav2_ss(j, k, ata) = -(CONJG(gkal(j)) * cavsig2_ss(k, a, sp) + gkal(k) * cavsig2_ss(j, at, sm)) / &
                       & (2.0d0 * kappa - i * (wl(j) - wl(k)))

    !--------------------------------------!
    !     Third-Order: Cavity and Atom     !
    !--------------------------------------!
    ! <a a \sigma>
    alpha = -(0.5d0 * gamma + 2.0d0 * kappa + i * (wl(j) + wl(k)))
    beta = -(gamma + 2.0d0 * kappa + i * (wl(j) + wl(k)))
    ! Set Matrix elements
    M_inverse = INVERT_MATRIX(alpha, beta)
    B_ss(1) = 0.0d0
    B_ss(2) = -0.5d0 * gkal(j) * (cavsig2_ss(k, a, sz) + cav1_ss(k, a)) + &
            & -0.5d0 * gkal(k) * (cavsig2_ss(j, a, sz) + cav1_ss(j, a))
    B_ss(3) = -gamma * cav2_ss(j, k, a) + &
            & gkal(j) * cavsig2_ss(k, a, sm) + &
            & gkal(k) * cavsig2_ss(j, a, sm)
    ! Matrix multiplication
    cavsig3_ss(j, k, a, sm) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
    cavsig3_ss(j, k, a, sp) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
    cavsig3_ss(j, k, a, sz) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))

    ! <at at \sigma>
    alpha = -(0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) + wl(k)))
    beta = -(gamma + 2.0d0 * kappa - i * (wl(j) + wl(k)))
    ! Set Matrix elements
    M_inverse = INVERT_MATRIX(alpha, beta)
    B_ss(1) = -0.5d0 * CONJG(gkal(j)) * (cavsig2_ss(k, at, sz) + cav1_ss(k, at)) + &
            & -0.5d0 * CONJG(gkal(k)) * (cavsig2_ss(j, at, sz) + cav1_ss(j, at))
    B_ss(2) = 0.0d0
    B_ss(3) = -gamma * cav2_ss(j, k, at) + &
            & CONJG(gkal(j)) * cavsig2_ss(k, at, sp) + &
            & CONJG(gkal(k)) * cavsig2_ss(j, at, sp)
    ! Matrix multiplication
    cavsig3_ss(j, k, at, sm) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
    cavsig3_ss(j, k, at, sp) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
    cavsig3_ss(j, k, at, sz) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))

    ! <at a \sigma>
    alpha = -(0.5d0 * gamma + 2.0d0 * kappa - i * (wl(j) - wl(k)))
    beta = -(gamma + 2.0d0 * kappa - i * (wl(j) - wl(k)))
    ! Set Matrix elements
    M_inverse = INVERT_MATRIX(alpha, beta)
    B_ss(1) = -0.5d0 * CONJG(gkal(j)) * (cavsig2_ss(k, a, sz) + cav1_ss(k, a))
    B_ss(2) = -0.5d0 * gkal(k) * (cavsig2_ss(j, at, sz) + cav1_ss(j, at))
    B_ss(3) = -gamma * cav2_ss(j, k, ata) + &
            & CONJG(gkal(j)) * cavsig2_ss(k, a, sp) + &
            & gkal(k) * cavsig2_ss(j, at, sm)
    ! Matrix multiplication
    cavsig3_ss(j, k, ata, sm) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
    cavsig3_ss(j, k, ata, sp) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
    cavsig3_ss(j, k, ata, sz) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))
  END DO
END DO

! Cycle through modes
DO j = -N, N
  DO k = -N, N
    DO l = -N, N
      !-----------------------------!
      !     Third-Order: Cavity     !
      !-----------------------------!
      ! <at a a>
      cav3_ss(j, k, l, a) = -(CONJG(gkal(j)) * cavsig3_ss(k, l, a, sp) + &
                          & gkal(k) * cavsig3_ss(j, l, ata, sm) + &
                          & gkal(l) * cavsig3_ss(j, k, ata, sm)) / &
                          & (3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l)))
      ! <at at a>
      cav3_ss(j, k, l, at) = -(CONJG(gkal(j)) * cavsig3_ss(k, l, ata, sp) + &
                           & CONJG(gkal(k)) * cavsig3_ss(j, l, ata, sp) + &
                           & gkal(l) * cavsig3_ss(j, k, at, sm)) / &
                           & (3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l)))

      !---------------------------------------!
      !     Fourth-Order: Cavity and Atom     !
      !---------------------------------------!
      ! <at a a \sigma>
      alpha = -(0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l)))
      beta = -(gamma + 3.0d0 * kappa - i * (wl(j) - wl(k) - wl(l)))
      ! Set Matrix elements
      M_inverse = INVERT_MATRIX(alpha, beta)
      B_ss(1) = -0.5d0 * CONJG(gkal(j)) * (cavsig3_ss(k, l, a, sz) + cav2_ss(k, l, a))
      B_ss(2) = -0.5d0 * gkal(k) * (cavsig3_ss(j, l, ata, sz) + cav2_ss(j, l, ata)) + &
              & -0.5d0 * gkal(l) * (cavsig3_ss(j, k, ata, sz) + cav2_ss(j, k, ata))
      B_ss(3) = -gamma * cav3_ss(j, k, l, a) + &
              & CONJG(gkal(j)) * cavsig3_ss(k, l, a, sp) + &
              & gkal(k) * cavsig3_ss(j, l, ata, sm) + &
              & gkal(l) * cavsig3_ss(j, k, ata, sm)
      ! Matrix multiplication
      cavsig4_ss(j, k, l, a, sm) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
      cavsig4_ss(j, k, l, a, sp) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
      cavsig4_ss(j, k, l, a, sz) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))

      ! <at at a \sigma>
      alpha = -(0.5d0 * gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l)))
      beta = -(gamma + 3.0d0 * kappa - i * (wl(j) + wl(k) - wl(l)))
      ! Set Matrix elements
      M_inverse = INVERT_MATRIX(alpha, beta)
      B_ss(1) = -0.5d0 * CONJG(gkal(j)) * (cavsig3_ss(k, l, ata, sz) + cav2_ss(k, l, ata)) + &
              & -0.5d0 * CONJG(gkal(k)) * (cavsig3_ss(j, l, ata, sz) + cav2_ss(j, l, ata))
      B_ss(2) = -0.5d0 * gkal(l) * (cavsig3_ss(j, k, at, sz) + cav2_ss(j, k, at))
      B_ss(3) = -gamma * cav3_ss(j, k, l, at) + &
              & CONJG(gkal(j)) * cavsig3_ss(k, l, ata, sp) + &
              & CONJG(gkal(k)) * cavsig3_ss(j, l, ata, sp) + &
              & gkal(l) * cavsig3_ss(j, k, at, sm)
      ! Matrix multiplication
      cavsig4_ss(j, k, l, at, sm) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
      cavsig4_ss(j, k, l, at, sp) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
      cavsig4_ss(j, k, l, at, sz) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))

    END DO
  END DO
END DO

! Cycle through modes
DO j = -N, N
  DO k = -N, N
    DO l = -N, N
      DO m = -N, N
        !------------------------------!
        !     Fourth-Order: Cavity     !
        !------------------------------!
        ! <at at a a>
        cav4_ss(j, k, l, m) = -(CONJG(gkal(j)) * cavsig4_ss(k, l, m, a, sp) + &
                            & CONJG(gkal(k)) * cavsig4_ss(j, l, m, a, sp) + &
                            & gkal(l) * cavsig4_ss(j, k, m, at, sm) + &
                            & gkal(m) * cavsig4_ss(j, k, l, at, sm)) / &
                            & (4.0d0 * kappa - i * (wl(j) + wl(k) - wl(l) - wl(m)))
      END DO
    END DO
  END DO
END DO

IF (N == 0) THEN
  ! PRINT*, "<a sm>_{ss} = ", cavsig2(0, a, sm)
  ! PRINT*, "<a sm>_{ss} = ", cavsig2_ss(0, a, sm)
  ! PRINT*, " "
  ! PRINT*, "<a sp>_{ss} = ", cavsig2(0, a, sp)
  ! PRINT*, "<a sp>_{ss} = ", cavsig2_ss(0, a, sp)
  ! PRINT*, " "
  ! PRINT*, "<a a>_ss = ", cav2(0, 0, a)
  ! PRINT*, "<a a>_ss = ", cav2_ss(0, 0, a)
  ! PRINT*, " "
  ! PRINT*, "<at at>_ss = ", cav2(0, 0, at)
  ! PRINT*, "<at at>_ss = ", cav2_ss(0, 0, at)
  ! PRINT*, " "
  ! PRINT*, "<at a>_ss = ", cav2(0, 0, ata)
  ! PRINT*, "<at a>_ss = ", cav2_ss(0, 0, ata)
ELSE
  ! First-Order: Bloch
  PRINT*, "===================================================="
  PRINT*, "First-Order: Bloch"
  PRINT*, "===================================================="
  PRINT*, "<sm>_{ss} = ", bloch(sm)
  PRINT*, "<sm>_{ss} = ", bloch_ss(sm)
  PRINT*, " "
  PRINT*, "<sp>_{ss} = ", bloch(sp)
  PRINT*, "<sp>_{ss} = ", bloch_ss(sp)
  PRINT*, " "
  PRINT*, "<sm>_{sz} = ", bloch(sz)
  PRINT*, "<sm>_{sz} = ", bloch_ss(sz)
  PRINT*, "===================================================="

  PRINT*, " "

  PRINT*, "===================================================="
  PRINT*, "First-Order: Cavity"
  PRINT*, "===================================================="
  PRINT*, "<a_{-1}>_{ss} = ", cav1(-1, a)
  PRINT*, "<a_{-1}>_{ss} = ", cav1_ss(-1, a)
  PRINT*, " "
  PRINT*, "<at_{-1}>_{ss} = ", cav1(-1, at)
  PRINT*, "<at_{-1}>_{ss} = ", cav1_ss(-1, at)
  PRINT*, "===================================================="
  PRINT*, " "

  ! Second-Order: Cavity and Atom
  PRINT*, "===================================================="
  PRINT*, "Second-Order: Cavity and Atom"
  PRINT*, "===================================================="
  PRINT*, "<a_{-1} sm>_{ss} = ", cavsig2(-1, a, sm)
  PRINT*, "<a_{-1} sm>_{ss} = ", cavsig2_ss(-1, a, sm)
  PRINT*, " "
  PRINT*, "<a_{-1} sp>_{ss} = ", cavsig2(-1, a, sp)
  PRINT*, "<a_{-1} sp>_{ss} = ", cavsig2_ss(-1, a, sp)
  PRINT*, " "
  PRINT*, "<a_{-1} sz>_{ss} = ", cavsig2(-1, a, sz)
  PRINT*, "<a_{-1} sz>_{ss} = ", cavsig2_ss(-1, a, sz)
  PRINT*, "----------------------------------------------------"
  PRINT*, "<at_{-1} sm>_{ss} = ", cavsig2(-1, at, sm)
  PRINT*, "<at_{-1} sm>_{ss} = ", cavsig2_ss(-1, at, sm)
  PRINT*, " "
  PRINT*, "<at_{-1} sp>_{ss} = ", cavsig2(-1, at, sp)
  PRINT*, "<at_{-1} sp>_{ss} = ", cavsig2_ss(-1, at, sp)
  PRINT*, " "
  PRINT*, "<at_{-1} sz>_{ss} = ", cavsig2(-1, at, sz)
  PRINT*, "<at_{-1} sz>_{ss} = ", cavsig2_ss(-1, at, sz)
  PRINT*, "===================================================="

  PRINT*, " "

  ! Second-Order: Cavity
  PRINT*, "===================================================="
  PRINT*, "Second-Order: Cavity"
  PRINT*, "===================================================="
  PRINT*, "<a_{-1} a_{0}>_ss = ", cav2(-1, 0, a)
  PRINT*, "<a_{-1} a_{0}>_ss = ", cav2_ss(-1, 0, a)
  PRINT*, " "
  PRINT*, "<at_{-1} at_{0}>_ss = ", cav2(-1, 0, at)
  PRINT*, "<at_{-1} at_{0}>_ss = ", cav2_ss(-1, 0, at)
  PRINT*, " "
  PRINT*, "<at_{-1} a_{0}>_ss = ", cav2(-1, 0, ata)
  PRINT*, "<at_{-1} a_{0}>_ss = ", cav2_ss(-1, 0, ata)
  PRINT*, "===================================================="

  PRINT*, " "

  ! Third-Order: Cavity and Atom
  PRINT*, "===================================================="
  PRINT*, "Third-Order: Cavity and Atom"
  PRINT*, "===================================================="
  PRINT*, "<a_{-1} a_{0} sm>_{ss}> = ", cavsig3(-1, 0, a, sm)
  PRINT*, "<a_{-1} a_{0} sm>_{ss}> = ", cavsig3_ss(-1, 0, a, sm)
  PRINT*, " "
  PRINT*, "<a_{-1} a_{0} sp>_{ss}> = ", cavsig3(-1, 0, a, sp)
  PRINT*, "<a_{-1} a_{0} sp>_{ss}> = ", cavsig3_ss(-1, 0, a, sp)
  PRINT*, " "
  PRINT*, "<a_{-1} a_{0} sz>_{ss}> = ", cavsig3(-1, 0, a, sz)
  PRINT*, "<a_{-1} a_{0} sz>_{ss}> = ", cavsig3_ss(-1, 0, a, sz)
  PRINT*, "----------------------------------------------------"
  PRINT*, "<at_{-1} at_{0} sm>_{ss}> = ", cavsig3(-1, 0, at, sm)
  PRINT*, "<at_{-1} at_{0} sm>_{ss}> = ", cavsig3_ss(-1, 0, at, sm)
  PRINT*, " "
  PRINT*, "<at_{-1} at_{0} sp>_{ss}> = ", cavsig3(-1, 0, at, sp)
  PRINT*, "<at_{-1} at_{0} sp>_{ss}> = ", cavsig3_ss(-1, 0, at, sp)
  PRINT*, " "
  PRINT*, "<at_{-1} at_{0} sz>_{ss}> = ", cavsig3(-1, 0, at, sz)
  PRINT*, "<at_{-1} at_{0} sz>_{ss}> = ", cavsig3_ss(-1, 0, at, sz)
  PRINT*, "----------------------------------------------------"
  PRINT*, "<at_{-1} a_{0} sm>_{ss} = ", cavsig3(-1, 0, ata, sm)
  PRINT*, "<at_{-1} a_{0} sm>_{ss} = ", cavsig3_ss(-1, 0, ata, sm)
  PRINT*, " "
  PRINT*, "<at_{-1} a_{0} sp>_{ss} = ", cavsig3(-1, 0, ata, sp)
  PRINT*, "<at_{-1} a_{0} sp>_{ss} = ", cavsig3_ss(-1, 0, ata, sp)
  PRINT*, " "
  PRINT*, "<at_{-1} a_{0} sz>_{ss} = ", cavsig3(-1, 0, ata, sz)
  PRINT*, "<at_{-1} a_{0} sz>_{ss} = ", cavsig3_ss(-1, 0, ata, sz)
  PRINT*, "===================================================="

  PRINT*, " "

  ! Third-Order: Cavity
  PRINT*, "===================================================="
  PRINT*, "Third-Order: Cavity"
  PRINT*, "===================================================="
  PRINT*, "<at_{-1} a_{0} a_{1}>_{ss} = ", cav3(-1, 0, 1, a)
  PRINT*, "<at_{-1} a_{0} a_{1}>_{ss} = ", cav3_ss(-1, 0, 1, a)
  PRINT*, " "
  PRINT*, "<at_{-1} at_{0} a_{1}>_{ss} = ", cav3(-1, 0, 1, at)
  PRINT*, "<at_{-1} at_{0} a_{1}>_{ss} = ", cav3_ss(-1, 0, 1, at)
  PRINT*, "===================================================="

  PRINT*, " "

  ! Fourth-Order: Cavity and Atom
  PRINT*, "===================================================="
  PRINT*, "Fourth-Order: Cavity and Atom"
  PRINT*, "===================================================="
  PRINT*, "<at_{-1} a_{0} a_{1} sm>_{ss} = ", cavsig4(-1, 0, 1, a, sm)
  PRINT*, "<at_{-1} a_{0} a_{1} sm>_{ss} = ", cavsig4_ss(-1, 0, 1, a, sm)
  PRINT*, " "
  PRINT*, "<at_{-1} a_{0} a_{1} sp>_{ss} = ", cavsig4(-1, 0, 1, a, sp)
  PRINT*, "<at_{-1} a_{0} a_{1} sp>_{ss} = ", cavsig4_ss(-1, 0, 1, a, sp)
  PRINT*, " "
  PRINT*, "<at_{-1} a_{0} a_{1} sz>_{ss} = ", cavsig4(-1, 0, 1, a, sz)
  PRINT*, "<at_{-1} a_{0} a_{1} sz>_{ss} = ", cavsig4_ss(-1, 0, 1, a, sz)
  PRINT*, "----------------------------------------------------"
  PRINT*, "<at_{-1} at_{0} a_{1} sm>_{ss} = ", cavsig4(-1, 0, 1, at, sm)
  PRINT*, "<at_{-1} at_{0} a_{1} sm>_{ss} = ", cavsig4_ss(-1, 0, 1, at, sm)
  PRINT*, " "
  PRINT*, "<at_{-1} at_{0} a_{1} sp>_{ss} = ", cavsig4(-1, 0, 1, at, sp)
  PRINT*, "<at_{-1} at_{0} a_{1} sp>_{ss} = ", cavsig4_ss(-1, 0, 1, at, sp)
  PRINT*, " "
  PRINT*, "<at_{-1} at_{0} a_{1} sz>_{ss} = ", cavsig4(-1, 0, 1, at, sz)
  PRINT*, "<at_{-1} at_{0} a_{1} sz>_{ss} = ", cavsig4_ss(-1, 0, 1, at, sz)
  PRINT*, "===================================================="

  PRINT*, " "

  ! Fourth-Order: Cavity
  PRINT*, "===================================================="
  PRINT*, "Fourth-Order: Cavity"
  PRINT*, "===================================================="
  PRINT*, "<at_{-1} at_{0} a_{1} a_{0}>_{ss} = ", cav4(-1, 0, 1, 0)
  PRINT*, "<at_{-1} at_{0} a_{1} a_{0}>_{ss} = ", cav4_ss(-1, 0, 1, 0)
  PRINT*, "===================================================="
END IF

! Steady state mean photon number
photona_ss = 0.0d0
DO j = -N, N
  DO k = -N, N
    photona_ss = photona_ss + cav2_ss(j, k, ata)
  END DO
END DO

PRINT*, " "
PRINT*, "Steady state found!"
! PRINT*, "Mean photon number in cavity A =", photona
PRINT*, "Mean photon number in cavity A =", photona_ss
PRINT*, " "

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
PRINT*, "Runtime: ", end_time - start_time, "seconds"

CONTAINS

FUNCTION INVERT_MATRIX(alpha_in, beta_in)
  IMPLICIT NONE
  ! Input arguments
  COMPLEX(KIND=8), INTENT(IN) :: alpha_in, beta_in
  COMPLEX(KIND=8), DIMENSION(3, 3) :: INVERT_MATRIX
  ! Temporal denominator
  COMPLEX(KIND=8) :: denominator
  denominator = 0.5d0 / ((alpha * beta) + (Omega ** 2))
  ! Fill matrix
  !----------------------------------------------------------------------------!
  INVERT_MATRIX(1, 1) = denominator * (2.0d0 * alpha * beta + (Omega ** 2)) / alpha
  INVERT_MATRIX(1, 2) = denominator * ((Omega ** 2) / alpha)
  INVERT_MATRIX(1, 3) = denominator * -i * Omega
  !----------------------------------------------------------------------------!
  INVERT_MATRIX(2, 1) = denominator * ((Omega ** 2) / alpha)
  INVERT_MATRIX(2, 2) = denominator * (2.0d0 * alpha * beta + (Omega ** 2)) / alpha
  INVERT_MATRIX(2, 3) = denominator * i * Omega
  !----------------------------------------------------------------------------!
  INVERT_MATRIX(3, 1) = denominator * -2.0d0 * i * Omega
  INVERT_MATRIX(3, 2) = denominator * 2.0d0 * i * Omega
  INVERT_MATRIX(3, 3) = denominator * 2.0d0 * alpha
  !----------------------------------------------------------------------------!
END FUNCTION INVERT_MATRIX

END PROGRAM TWO_LEVEL_MULTI_MODE_FILTER_MOMENT_EQUATIONS