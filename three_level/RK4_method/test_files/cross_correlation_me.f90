PROGRAM Two_Filer_Master_Equation_Cross_Correlation

  IMPLICIT NONE

  ! Parameters in terms of decay rate gamma
  ! Atom decay rate
  REAL(KIND=8) :: gamma
  ! Drive strength (Rabi Frequency)
  REAL(KIND=8) :: omega
  ! Drive detuning from two-photon resonance, \omega_{d} - \omega_{gf}
  REAL(KIND=8) :: delta
  ! Drive strength ratio of the two levels, |g> <-> |e> and |e> <-> |f>
  REAL(KIND=8) :: xi
  ! Difference between two energy levels, \omega_{ef} - \omega_{ge}
  REAL(KIND=8) :: alpha

  ! Filter parameter stuff
  ! Detuning of cavity "a" resonance frequency with drive frequency
  ! \Delta_{f} = \omega_{0} - \omega_{d}. \Delta_{f} = 0 is resonant with
  ! \omega_{gf} / 2 if \delta = 0.
  REAL(KIND=8) :: D_a
  ! Detuning of cavity "b" resonance frequency with drive frequency
  REAL(KIND=8) :: D_b
  ! Cavity linewidth/transmission of cavity a
  REAL(KIND=8) :: kappa_a
  ! Cavity linewidth/transmission of cavity b
  REAL(KIND=8) :: kappa_b

  ! Eigenfrequencies for position of spectrum peaks for delta = 0
  ! eigen drive strength
  REAL(KIND=8) :: Omega_t
  ! Positive eigenfrequency
  REAL(KIND=8) :: wp
  ! Negative eigenfrequency
  REAL(KIND=8) :: wm

  ! Quantum object stuff
  ! Hilbert Space - max number of photons in cavity 0 -> N
  INTEGER, PARAMETER :: N = 2
  ! Hilbert truncated dimension. 3 atomic states x (N+1)^2 Fock states for
  ! two cavities.
  INTEGER, PARAMETER :: d_Hilb = INT(((N + 1) ** 2) * 3)
  ! Density operator matrix
  COMPLEX(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: rho
  ! Hamiltonian in extended basis
  COMPLEX(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: H
  ! Atom decay operators in extended basis
  REAL(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: sigmam, sigmap, sigmapm
  ! sigma_gg = |g><g|, sigma_ee = |e><e|, sigma_ff = |f><f|
  REAL(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: gg, ee, ff
  ! Lowering and raising operators for cavity a
  REAL(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: a, a_dag
  ! Lowering and raising operators for cavity b
  REAL(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: b, b_dag
  ! Lindblad operator J_a = \sqrt(gamma)\Sigma + \sqrt{\kappa_a}a for cavity a
  REAL(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: J_a, J_a_dag, JJ_a
  ! Lindblad operator J_b for cavity b coupling
  REAL(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: J_b, J_b_dag, JJ_b
  ! Number operator for cavity a N_s = a^{\dagger} a
  REAL(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: N_a
  ! Number operator for cavity b N_b = b^{\dagger} b
  REAL(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: N_b

  ! Correlation stuff
  ! Steady state density rho
  COMPLEX(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: rho_ss
  ! a^{\dagger}a * rho for correlation calculation
  COMPLEX(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: Nrho
  ! Nrho copy for calculating TRACE(MATMUL(a_dag, Nrho))
  COMPLEX(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: corr_mat
  ! Steady state mean photon number in cavity a and b
  REAL(KIND=8) :: mean_photon_a, mean_photon_b
  ! Correlation value for cavity a
  REAL(KIND=8) :: corr_a
  ! Cross correlation value for cavity b
  REAL(KIND=8) :: corr_b

  ! Time step stuff
  ! Time step
  REAL(KIND=8), PARAMETER :: dt = 0.001
  ! Max time
  REAL(KIND=8), PARAMETER :: steady_time = 50
  ! Max number of time steps for tau calculations
  INTEGER, PARAMETER :: steady_steps = INT(steady_time / dt)
  ! Maximum tau time to calculate correlation function for
  REAL(KIND=8) :: tau_max
  ! Max number of tau steps
  INTEGER :: tau_steps

  ! Useful stuff
  ! Time step integer for loop
  INTEGER(KIND=8) :: k
  ! General integer counters
  INTEGER :: m
  ! Integer counters for a photon number and b photon number
  INTEGER :: na, nb
  ! Integer place for photon number in psi vector (na + 1) * (nb + 1)
  INTEGER :: nplace, nplace_pm1
  ! Complex i = SQRT(-1)
  COMPLEX(KIND=8), PARAMETER :: i = CMPLX(0,1)
  ! Temporal value
  REAL(KIND=8) :: temp
  ! Runge-Kutta 4th Order Vectors/Matrices
  COMPLEX(KIND=8), DIMENSION(d_Hilb, d_Hilb) :: k1, k2, k3, k4
  ! 1.0 / 6.0 to maybe save on calculation time
  REAL(KIND=8), PARAMETER :: xis = 1.0 / 6.0

  ! Data stuff
  ! State population values
  REAL(KIND=8) :: popgg, popee, popff
  ! Mean photon number for cavity a and cavity b
  REAL(KIND=8) :: photon_a, photon_b
  ! Filename for saving data
  CHARACTER(LEN=44) :: filename_correlation = "./data_files/me_correlation.txt"
  ! CHARACTER(LEN=36) :: filename_time = "./data_files/master_equation/tau.txt"

  !----------------------------------------------------------------------------!
  !                   Calling parameters from ParamList_me.nml                 !
  !----------------------------------------------------------------------------!
  ! namelists
  INTEGER :: iunit
  INTEGER :: istat
  CHARACTER(LEN=512) :: line
  NAMELIST /PARAMS/ gamma, omega, delta, xi, alpha
  ! NAMELIST /EIGENVALUES/ Omega_t, wm, wp
  NAMELIST /CAVITY/ D_a, kappa_a, D_b, kappa_b
  NAMELIST /TIME/ tau_max

  ! Parameters in terms of decay rate gamma
  ! Atom decay rate
  gamma = 1.0
  ! Drive strength (Rabi Frequency)
  omega = 40.0
  ! Drive detuning from two-photon resonance, \omega_{d} - \omega_{gf}
  delta = 0.0
  ! Drive strength ratio of the two levels, |g> <-> |e> and |e> <-> |f>
  xi = 1.0
  ! Difference between two energy levels, \omega_{ef} - \omega_{ge}
  alpha = -120.0

  ! read the PARAMS namelist
  iunit = 31
  OPEN(iunit, FILE="test_files/cross_correlation_me.nml", STATUS="OLD", DELIM="QUOTE")
  READ(iunit, NML=PARAMS, IOSTAT=istat)
  IF (istat .NE. 0) THEN
    BACKSPACE(iunit)
    READ(iunit, FMT='(A)') line
    CLOSE(iunit)
    PRINT *, "Invalid line in PARAMS namelist: " // TRIM(line)
    CALL EXIT(1)
  END IF
  CLOSE(IUNIT)
  print *, "DEBUG: gamma is: ", gamma

  ! ! Set eigenvalues of system for \delta = 0
  ! ! eigen drive strength
  ! Omega_t = SQRT(((0.25*alpha)**2) + ((0.5*omega) &
  !                                    & ** 2) * (1 + (xi ** 2)))
  ! ! Positive eigenfrequency
  ! wp = -(0.25 * alpha) + Omega_t
  ! ! Negative eigenfrequency
  ! wm = -(0.25 * alpha) - Omega_t
  !
  ! iunit = 31
  ! OPEN(iunit, FILE="ParamList.nml", STATUS="OLD", DELIM="QUOTE")
  ! READ(iunit, NML=EIGENVALUES, IOSTAT=istat)
  ! IF (istat .NE. 0) THEN
  !   BACKSPACE(iunit)
  !   READ(iunit, FMT='(A)') line
  !   CLOSE(iunit)
  !   PRINT *, "Invalid line in PARAMS namelist: " // TRIM(line)
  !   CALL EXIT(1)
  ! END IF
  ! CLOSE(IUNIT)
  ! print *, "DEBUG: gamma is: ", gamma

  ! Cavity param list
  ! Default value for cavities
  D_a = 0.0
  kappa_a = 10.0
  D_b = 0.0
  kappa_b = 10.0

  ! read the CAVITY namelist
  iunit = 31
  OPEN(iunit, FILE="test_files/cross_correlation_me.nml", STATUS="OLD", DELIM="QUOTE")
  READ(iunit, NML=CAVITY, IOSTAT=istat)
  IF (istat .NE. 0) THEN
    BACKSPACE(iunit)
    READ(iunit, FMT='(A)') line
    CLOSE(iunit)
    PRINT *, "Invalid line in TIME namelist: " // TRIM(line)
    CALL EXIT(1)
  END IF
  CLOSE(iunit)

  ! Time param list
  ! Default value for max time, in units of \gamma \tau
  tau_max = 20.0

  ! read the TIME namelist
  iunit = 31
  OPEN(iunit, FILE="test_files/cross_correlation_me.nml", STATUS="OLD", DELIM="QUOTE")
  READ(iunit, NML=TIME, IOSTAT=istat)
  IF (istat .NE. 0) THEN
    BACKSPACE(iunit)
    READ(iunit, FMT='(A)') line
    CLOSE(iunit)
    PRINT *, "Invalid line in TIME namelist: " // TRIM(line)
    CALL EXIT(1)
  END IF
  CLOSE(iunit)

  ! Set number of steps based on tau_max values
  tau_steps = INT(tau_max / dt) + 1

  !----------------------------------------------------------------------------!
  !                    END parameters from ParamList_me.nml                    !
  !----------------------------------------------------------------------------!

  ! Initialising matrices
  ! Operators in extended basis {|1,g>, |1,e>, |1,f>, |2,g>, ..., |N,f>}

  ! In the atom basis the kets are: |g> = (1, 0, 0)^T
  !                                 |e> = (0, 1, 0)^T
  !                                 \f> = (0, 0, 1)^T

  ! Atom Raising and Lowering Operators
  sigmam = 0
  sigmap = 0

  DO na=0,N
    DO nb=0,N
      ! photon number index
      nplace = (3 * (N + 1) * na) + (3 * nb)
      ! \Sigma^{\dagger} = |g><e| + \xi|e><f
      sigmam(nplace + 1, nplace + 2) = 1.0
      sigmam(nplace + 2, nplace + 3) = xi
      ! \Sigma^{\dagger} = |e><g| + \xi|f><e|
      sigmap(nplace + 2, nplace + 1) = 1.0
      sigmap(nplace + 3, nplace + 2) = xi
    END DO
  END DO

  ! sigmapm = \Sigma^{\dagger}\Sigma
  sigmapm = 0
  sigmapm = MATMUL(sigmap, sigmam)

  ! State operators gg = |g><g|, ee = |e><e|, ff = |f><f|
  gg = 0
  ee = 0
  ff = 0

  DO na=0,N
    DO nb=0,N
      ! photon number index
      nplace = (3 * (N + 1) * na) + (3 * nb)
      gg(nplace + 1, nplace + 1) = 1.0
      ee(nplace + 2, nplace + 2) = 1.0
      ff(nplace + 3, nplace + 3) = 1.0
    END DO
  END DO

  ! Cavity annihilation operator a: a|n> = SQRT(n)|n-1>
  a = 0
  DO na=0,N
    DO nb=0,N
      ! photon number index
      nplace = (3 * (N + 1) * na) + (3 * nb)
      ! photon number index minus one na - 1
      nplace_pm1 = (3 * (N + 1) * (na - 1)) + (3 * nb)
      DO m=1,3
        a(nplace_pm1 + m, nplace + m) = SQRT(1.0 * na)
      END DO
    END DO
  END DO
  ! Cavity creation operator a^{\dagger}: a^{\dagger}|n> = SQRT(n+1)|n+1>
  a_dag = 0
  FORALL (na=1:d_Hilb, nb=1:d_Hilb)
    a_dag(na,nb) = a(nb,na)
  END FORALL

  ! Number operator N = a^{\dagger} a: N|n> = n|n>
  N_a = MATMUL(a_dag, a)

  ! Cavity lowering operator for cavity b
  b = 0
  DO na=0,N
    DO nb=0,N
      ! photon number index
      nplace = (3 * (N + 1) * na) + (3 * nb)
      ! photon number index minus one na - 1
      nplace_pm1 = (3 * (N + 1) * na) + (3 * (nb - 1))
      DO m=1,3
        b(nplace_pm1 + m, nplace + m) = SQRT(1.0 * nb)
      END DO
    END DO
  END DO

  b_dag = 0
  FORALL (na=1:d_Hilb, nb=1:d_Hilb)
    b_dag(na,nb) = b(nb,na)
  END FORALL

  N_b = MATMUL(b_dag, b)

  ! Hamiltonian for the atom
  H = 0
  DO na=0,N
    DO nb=0,N
      ! photon number index
      nplace = (3 * (N + 1) * na) + (3 * nb)
      H(nplace + 1, nplace + 2) = 0.5 * omega
      H(nplace + 2, nplace + 1) = 0.5 * omega
      H(nplace + 2, nplace + 2) = -(0.5 * alpha + delta)
      H(nplace + 2, nplace + 3) = 0.5 * omega * xi
      H(nplace + 3, nplace + 2) = 0.5 * omega * xi
      H(nplace + 3, nplace + 3) = -2.0 * delta
    END DO
  END DO

  ! H = H_atom + H_filter + i\hbar\sqrt{\gamma\kappa_a}(a\Sigma^{\dagger} +
  ! \Sigma a^{\dagger})
  H = H + (D_a * N_a) + (D_b * N_b)
  H = H + 0.5 * i * SQRT(0.5 * gamma * kappa_a) * (MATMUL(a, sigmap) - &
    & MATMUL(sigmam, a_dag))
  H = H + 0.5 * i * SQRT(0.5 * gamma * kappa_b) * (MATMUL(b, sigmap) - &
    & MATMUL(sigmam, b_dag))

  ! Linblad Decay operator J = \sqrt(gamma)\Sigma + \sqrt{\kappa_a}a
  J_a = (SQRT(0.5 * gamma) * sigmam) + (SQRT(kappa_a) * a)
  J_a_dag = (SQRT(0.5 * gamma) * sigmap) + (SQRT(kappa_a) * a_dag)

  JJ_a = MATMUL(J_a_dag, J_a)

  J_b = (SQRT(0.5 * gamma) * sigmam) + (SQRT(kappa_b) * b)
  J_b_dag = (SQRT(0.5 * gamma) * sigmap) + (SQRT(kappa_b) * b_dag)

  JJ_b = MATMUL(J_b_dag, J_b)

  ! Initialise density operator in the ground state <g|rho|g> = 1
  rho = 0
  rho(1,1) = 1.0

  ! Initialise Runge-Kutta 4-th Order Vectors
  k1 = 0
  k2 = 0
  k3 = 0
  k4 = 0

  ! Solve the master equation for a long time to find the system's steady state
  ! density matrix
  DO k=0,steady_steps
    k1 = -i * dt * ( MATMUL(H, rho) - MATMUL(rho, H) ) + &
       ! Cavity a cascaded evolution
       & dt * MATMUL(MATMUL(J_a, rho), J_a_dag) - &
       & 0.5 * dt * (MATMUL(JJ_a, rho) + MATMUL(rho, JJ_a)) + &
       ! Cavity a emission
       & kappa_a * dt * MATMUL(MATMUL(a, rho), a_dag) - &
       & 0.5 * kappa_a * dt * (MATMUL(N_a, rho) + MATMUL(rho, N_a)) + &
       ! Cavity b cascaded evolution
       & dt * MATMUL(MATMUL(J_b, rho), J_b_dag) - &
       & 0.5 * dt * (MATMUL(JJ_b, rho) + MATMUL(rho, JJ_b)) + &
       ! Cavity b emission
       & kappa_b * dt * MATMUL(MATMUL(b, rho), b_dag) - &
       & 0.5 * kappa_b * dt * (MATMUL(N_b, rho) + MATMUL(rho, N_b))

    k2 = -i * dt * ( MATMUL(H, (rho + 0.5*k1)) - MATMUL((rho + 0.5*k1), H) ) + &
       ! Cavity a cascaded evolution
       & dt * MATMUL(MATMUL(J_a, (rho + 0.5*k1)), J_a_dag) - &
       & 0.5 * dt * (MATMUL(JJ_a, (rho + 0.5*k1)) + MATMUL((rho + 0.5*k1), JJ_a)) + &
       ! Cavity a emission
       & kappa_a * dt * MATMUL(MATMUL(a, (rho + 0.5*k1)), a_dag) - &
       & 0.5 * kappa_a * dt * (MATMUL(N_a, (rho + 0.5*k1)) + MATMUL((rho + 0.5*k1), N_a)) + &
       ! Cavity b cascaded evolution
       & dt * MATMUL(MATMUL(J_b, (rho + 0.5*k1)), J_b_dag) - &
       & 0.5 * dt * (MATMUL(JJ_b, (rho + 0.5*k1)) + MATMUL((rho + 0.5*k1), JJ_b)) + &
       ! Cavity b emission
       & kappa_b * dt * MATMUL(MATMUL(b, (rho + 0.5*k1)), b_dag) - &
       & 0.5 * kappa_b * dt * (MATMUL(N_b, (rho + 0.5*k1)) + MATMUL((rho + 0.5*k1), N_b))

    k3 = -i * dt * ( MATMUL(H, (rho + 0.5*k2)) - MATMUL((rho + 0.5*k2), H) ) + &
       ! Cavity a cascaded evolution
       & dt * MATMUL(MATMUL(J_a, (rho + 0.5*k2)), J_a_dag) - &
       & 0.5 * dt * (MATMUL(JJ_a, (rho + 0.5*k2)) + MATMUL((rho + 0.5*k2), JJ_a)) + &
       ! Cavity a emission
       & kappa_a * dt * MATMUL(MATMUL(a, (rho + 0.5*k2)), a_dag) - &
       & 0.5 * kappa_a * dt * (MATMUL(N_a, (rho + 0.5*k2)) + MATMUL((rho + 0.5*k2), N_a)) + &
       ! Cavity b cascaded evolution
       & dt * MATMUL(MATMUL(J_b, (rho + 0.5*k2)), J_b_dag) - &
       & 0.5 * dt * (MATMUL(JJ_b, (rho + 0.5*k2)) + MATMUL((rho + 0.5*k2), JJ_b)) + &
       ! Cavity b emission
       & kappa_b * dt * MATMUL(MATMUL(b, (rho + 0.5*k2)), b_dag) - &
       & 0.5 * kappa_b * dt * (MATMUL(N_b, (rho + 0.5*k2)) + MATMUL((rho + 0.5*k2), N_b))

    k4 = -i * dt * ( MATMUL(H, (rho + k3)) - MATMUL((rho + k3), H) ) + &
       ! Cavity a cascaded evolution
       & dt * MATMUL(MATMUL(J_a, (rho + k3)), J_a_dag) - &
       & 0.5 * dt * (MATMUL(JJ_a, (rho + k3)) + MATMUL((rho + k3), JJ_a)) + &
       ! Cavity a emission
       & kappa_a * dt * MATMUL(MATMUL(a, (rho + k3)), a_dag) - &
       & 0.5 * kappa_a * dt * (MATMUL(N_a, (rho + k3)) + MATMUL((rho + k3), N_a)) + &
       ! Cavity b cascaded evolution
       & dt * MATMUL(MATMUL(J_b, (rho + k3)), J_b_dag) - &
       & 0.5 * dt * (MATMUL(JJ_b, (rho + k3)) + MATMUL((rho + k3), JJ_b)) + &
       ! Cavity b emission
       & kappa_b * dt * MATMUL(MATMUL(b, (rho + k3)), b_dag) - &
       & 0.5 * kappa_b * dt * (MATMUL(N_b, (rho + k3)) + MATMUL((rho + k3), N_b))

    rho = rho + (1.0 / 6.0) * (k1 + 2.0*(k2 + k3) + k4)

    PRINT*, (100.0 * k) / (1.0 * steady_steps), "% complete"
  END DO

  PRINT*, "steady state found"
  rho_ss = rho

  ! Calculate mean photon number in cavity for steady state
  mean_photon_a = EXPECTATION(N_a, rho_ss)
  mean_photon_b = EXPECTATION(N_b, rho_ss)

  ! Now propagate a * rho_ss * a^{\dagger}
  Nrho = 0
  Nrho = MATMUL(MATMUL(a, rho_ss), a_dag)

  ! ! Open the file used to write data
  ! OPEN(UNIT = 1, file = filename_time, STATUS = 'replace', ACTION = 'write')!
  ! ! Write the values of the parameters to the first 9 lines. The time values
  ! ! will start on the 10th line
  ! WRITE(1,*) omega
  ! WRITE(1,*) delta
  ! WRITE(1,*) xi
  ! WRITE(1,*) alpha
  ! WRITE(1,*) D_a
  ! WRITE(1,*) kappa_a
  ! WRITE(1,*) D_b
  ! WRITE(1,*) kappa_b
  ! WRITE(1,*) ' '

  ! Open files for state probabilities to be written to
  OPEN(UNIT = 2, file = filename_correlation, STATUS = 'replace', ACTION = 'write')

  ! Calculate correlation for tau time
  DO k=0,tau_steps
    ! Calculate correlation value of cavity a
    corr_mat = MATMUL(N_a, Nrho)
    corr_a = EXPECTATION(N_a, Nrho)
    corr_a = corr_a / (mean_photon_a ** 2)

    ! Calculate cross-correlation value of cavity b
    corr_mat = MATMUL(N_b, Nrho)
    corr_b = EXPECTATION(N_b, Nrho)
    ! corr_b = corr_b / (mean_photon_b ** 2)
    corr_b = corr_b / (mean_photon_a * mean_photon_b)

    ! Write values to files
    WRITE(2,*) k * dt, REAL(corr_a), REAL(corr_b)

    ! Propagate Nrho
    k1 = -i * dt * ( MATMUL(H, Nrho) - MATMUL(Nrho, H) ) + &
       ! Cavity a cascaded evolution
       & dt * MATMUL(MATMUL(J_a, Nrho), J_a_dag) - &
       & 0.5 * dt * (MATMUL(JJ_a, Nrho) + MATMUL(Nrho, JJ_a)) + &
       ! Cavity a emission
       & kappa_a * dt * MATMUL(MATMUL(a, Nrho), a_dag) - &
       & 0.5 * kappa_a * dt * (MATMUL(N_a, Nrho) + MATMUL(Nrho, N_a)) + &
       ! Cavity b cascaded evolution
       & dt * MATMUL(MATMUL(J_b, Nrho), J_b_dag) - &
       & 0.5 * dt * (MATMUL(JJ_b, Nrho) + MATMUL(Nrho, JJ_b)) + &
       ! Cavity b emission
       & kappa_b * dt * MATMUL(MATMUL(b, Nrho), b_dag) - &
       & 0.5 * kappa_b * dt * (MATMUL(N_b, Nrho) + MATMUL(Nrho, N_b))

    k2 = -i * dt * ( MATMUL(H, (Nrho + 0.5*k1)) - MATMUL((Nrho + 0.5*k1), H) ) + &
       ! Cavity a cascaded evolution
       & dt * MATMUL(MATMUL(J_a, (Nrho + 0.5*k1)), J_a_dag) - &
       & 0.5 * dt * (MATMUL(JJ_a, (Nrho + 0.5*k1)) + MATMUL((Nrho + 0.5*k1), JJ_a)) + &
       ! Cavity a emission
       & kappa_a * dt * MATMUL(MATMUL(a, (Nrho + 0.5*k1)), a_dag) - &
       & 0.5 * kappa_a * dt * (MATMUL(N_a, (Nrho + 0.5*k1)) + MATMUL((Nrho + 0.5*k1), N_a)) + &
       ! Cavity b cascaded evolution
       & dt * MATMUL(MATMUL(J_b, (Nrho + 0.5*k1)), J_b_dag) - &
       & 0.5 * dt * (MATMUL(JJ_b, (Nrho + 0.5*k1)) + MATMUL((Nrho + 0.5*k1), JJ_b)) + &
       ! Cavity b emission
       & kappa_b * dt * MATMUL(MATMUL(b, (Nrho + 0.5*k1)), b_dag) - &
       & 0.5 * kappa_b * dt * (MATMUL(N_b, (Nrho + 0.5*k1)) + MATMUL((Nrho + 0.5*k1), N_b))

    k3 = -i * dt * ( MATMUL(H, (Nrho + 0.5*k2)) - MATMUL((Nrho + 0.5*k2), H) ) + &
       ! Cavity a cascaded evolution
       & dt * MATMUL(MATMUL(J_a, (Nrho + 0.5*k2)), J_a_dag) - &
       & 0.5 * dt * (MATMUL(JJ_a, (Nrho + 0.5*k2)) + MATMUL((Nrho + 0.5*k2), JJ_a)) + &
       ! Cavity a emission
       & kappa_a * dt * MATMUL(MATMUL(a, (Nrho + 0.5*k2)), a_dag) - &
       & 0.5 * kappa_a * dt * (MATMUL(N_a, (Nrho + 0.5*k2)) + MATMUL((Nrho + 0.5*k2), N_a)) + &
       ! Cavity b cascaded evolution
       & dt * MATMUL(MATMUL(J_b, (Nrho + 0.5*k2)), J_b_dag) - &
       & 0.5 * dt * (MATMUL(JJ_b, (Nrho + 0.5*k2)) + MATMUL((Nrho + 0.5*k2), JJ_b)) + &
       ! Cavity b emission
       & kappa_b * dt * MATMUL(MATMUL(b, (Nrho + 0.5*k2)), b_dag) - &
       & 0.5 * kappa_b * dt * (MATMUL(N_b, (Nrho + 0.5*k2)) + MATMUL((Nrho + 0.5*k2), N_b))

    k4 = -i * dt * ( MATMUL(H, (Nrho + k3)) - MATMUL((Nrho + k3), H) ) + &
       ! Cavity a cascaded evolution
       & dt * MATMUL(MATMUL(J_a, (Nrho + k3)), J_a_dag) - &
       & 0.5 * dt * (MATMUL(JJ_a, (Nrho + k3)) + MATMUL((Nrho + k3), JJ_a)) + &
       ! Cavity a emission
       & kappa_a * dt * MATMUL(MATMUL(a, (Nrho + k3)), a_dag) - &
       & 0.5 * kappa_a * dt * (MATMUL(N_a, (Nrho + k3)) + MATMUL((Nrho + k3), N_a)) + &
       ! Cavity b cascaded evolution
       & dt * MATMUL(MATMUL(J_b, (Nrho + k3)), J_b_dag) - &
       & 0.5 * dt * (MATMUL(JJ_b, (Nrho + k3)) + MATMUL((Nrho + k3), JJ_b)) + &
       ! Cavity b emission
       & kappa_b * dt * MATMUL(MATMUL(b, (Nrho + k3)), b_dag) - &
       & 0.5 * kappa_b * dt * (MATMUL(N_b, (Nrho + k3)) + MATMUL((Nrho + k3), N_b))

    Nrho = Nrho + (1.0 / 6.0) * (k1 + 2.0*(k2 + k3) + k4)

    PRINT*, (100.0 * k) / (1.0 * tau_steps), "% complete"
  END DO

CONTAINS

FUNCTION EXPECTATION(AtA, rho)
  ! Expectation of an operator <A> = Tr{A\rho}. This can be used to find
  ! the mean photon number by using a^{\dagger}a.
  IMPLICIT NONE
  ! Input: operators to find expectation
  REAL(KIND=8), DIMENSION(:,:) :: AtA
  ! Input: density matrix
  COMPLEX(KIND=8), DIMENSION(:,:) :: rho
  ! Multiplication of input matrices to calculate the trace
  COMPLEX(KIND=8), DIMENSION(:,:), ALLOCATABLE :: trace_mat
  ! Dimension of matrices and counter integer
  INTEGER :: dimen, i
  ! Output expectation value
  REAL(KIND=8) :: EXPECTATION

  dimen = INT(SQRT(SIZE(rho) * 1.0))
  ALLOCATE(trace_mat(dimen, dimen))
  trace_mat = 0
  trace_mat = MATMUL(AtA, rho)
  EXPECTATION = 0.0
  DO i = 1,dimen
    EXPECTATION = EXPECTATION + REAL(trace_mat(i,i))
  END DO
END FUNCTION EXPECTATION

END PROGRAM Two_Filer_Master_Equation_Cross_Correlation
