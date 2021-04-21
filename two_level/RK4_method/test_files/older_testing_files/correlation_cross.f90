! This program solves for the population inversion <\sigma_{z}(t)> and the mean
! photon number <A^{\dagger} A(t)> for the multi-mode filter array. This version
! of the program uses an analytic expression, derived in some notes somewhere in
! my folders.

! The Hamiltonian is in the form H = H_a + H_c + H_I, where
! -  H_a = 0.5 * Omega * (\sigma_{+} + \sigma_{-}),
! -  H_c = \sum_{ja=-N}^{N} \omega_{ja} a^{\dagger}_{N} a_{N},
! and
! -  H_I = -0.5 i \sqrt{\gamma \kappa} *
! \sum_{ja=-N}^{N} (e^{i \phi_{j}} a^{\dagger}_{ja} \sigma_{-} - e^{-i \phi_{j}} a_{ja} \sigma_{+})
! with
! -  \omega_{ja} = \omega_{0} + (ja \delta\omega)
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

PROGRAM TWO_LEVEL_MULTI_MODE_CROSS_CORRELATIONS

IMPLICIT NONE

! Parameters in terms of decay rate gamma
! Atom decay rate
REAL(KIND=8) :: gamma
! Drive strength (Rabi Frequency)
REAL(KIND=8) :: Omega

! Filter parameter stuff
! Central mode frequency of the filter cavity, with n mode frequencie either side
REAL(KIND=8) :: w0a, w0b
! Cavity linewidth/transmission of cavity mode
REAL(KIND=8) :: kappaa, kappab
! Percentage of fluorecence aimed at cavity
REAL(KIND=8) :: epsilon
! Number of mode either side of w0, 2N + 1 total mode
INTEGER :: N
! Frequency spacing of modes
REAL(KIND=8) :: dwa, dwb
! Phase modulation of mode coupling
INTEGER :: phase
! Kappa values for both cavities
REAL(KIND=8), DIMENSION(2) :: kappa
! List of omega values
REAL(KIND=8), DIMENSION(:, :), ALLOCATABLE :: wl
! List of mode dependent cascade coupling values
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE :: gkl

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
! Time step integer
INTEGER :: t
! Time unit
REAL(KIND=8) :: dtime

! Quantum Object stuff
! First-order moments: Bloch equations (<\sigma_{-}>, <\sigma_{+}>, <\sigma_{z}>)
COMPLEX(KIND=8), DIMENSION(3) :: bloch_ss
! First-order moments: Cavity (<a>, <a^{\dagger})
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE :: cav1_ss
! Second-order moments: Cavity (<a^{\dagger} a>)
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :), ALLOCATABLE :: cav2_ss
! Second-order moments: Cavity and atom (<a \sigma_{-+z}>, <a^{\dagger}, \sigma_{-+z}>)
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE :: cavsig2_ss
! Third-order moments: Cavity (<a^{2} a^{\dagger}>, <a^{\dagger 2} a>)
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :, :), ALLOCATABLE :: cav3_ss
! Third-order moments: Cavity and atom
!   (<aa \sigma_{-+z}>, <a^{\dagger 2} \sigma_{-+z}>, <a^{\dagger} a \sigma_{-+z}>)
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: cavsig3_ss
! Fourth-order moments: Cavity (<a^{\dagger 2} a^{2}>)
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :, :, :), ALLOCATABLE :: cav4_ss
! Fourth-order moments: Cavity and atom ( <a^{\dagger} a^{2} \sigma_{-+z}>, <a^{\dagger 2} a \sigma_{-+z}>)
COMPLEX(KIND=8), DIMENSION(:, :, :, :, :, :, :, :), ALLOCATABLE :: cavsig4_ss
! Integer indices for: sigma_{-}, sigma_{+}, sigma_{z}
INTEGER, PARAMETER :: sm = 1, sp = 2, sz = 3
! Integer indices for: a, a^{\dagger}, a^{\dagger} a
INTEGER, PARAMETER :: a = 1, at = 2, ata = 3, B = 2
! Integer placeholders for loops
INTEGER :: Cj, Ck, Cl, Cm

! Analytic Things
! omega shortcut
COMPlEX(KIND=8) :: delta
! Temporal complex values
COMPLEX(KIND=8) :: alpha
! Inverse matrix stuff
COMPLEX(KIND=8), DIMENSION(3, 3) :: M_inverse
COMPLEX(KIND=8), DIMENSION(3) :: B_ss
! Initial moments and steady state moments
COMPLEX(KIND=8), DIMENSION(3) :: sm0, sp0, sz0, smss, spss, szss
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE :: aj0, ajss
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE :: ajsm0, ajsp0, ajsz0, ajsmss, ajspss, ajszss
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE :: atjak0, atjakss
! Differential Equation Coefficients
COMPLEX(KIND=8), DIMENSION(3) :: C1, C2, C3
COMPLEX(KIND=8), DIMENSION(:, :), ALLOCATABLE :: C4, C6, C7, C8
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE :: C9
! Cheaty variables
COMPLEX(KIND=8), DIMENSION(3) :: v1, v2, v3
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE :: va, vb, vc
! Cheaty function
COMPLEX(KIND=8), DIMENSION(:, :, :), ALLOCATABLE :: f_jk0
COMPLEX(KIND=8), DIMENSION(:, :, :, :), ALLOCATABLE :: f_jk
! Second-order correlation
COMPLEX(KIND=8), DIMENSION(3) :: g2

! Useful things
! Filter mode integers
INTEGER :: j, k, l, m
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
! Steady state mean photon number in cavity
REAL(KIND=8) :: photonA_ss, photonB_ss

! Filename stuff
! I/O Format for percentage
CHARACTER(LEN=28), PARAMETER :: FMT_ss = "(T2,I3,A28,F9.2,A19,F9.2,A1)"
CHARACTER(LEN=28), PARAMETER :: FMT_corr = "(T2,I3,A27,F9.2,A19,F9.2,A1)"
! Filename of parameters
CHARACTER(LEN=27), PARAMETER :: filename_parameters = "./data_files/parameters.txt"
! Filename for state population
! CHARACTER(LEN=23), PARAMETER :: filename_state = "./data_files/states.txt"
! Filename for first-order correlation of cavity a
! CHARACTER(LEN=24), PARAMETER :: filename_g1 = "./data_files/g1_corr.txt"
! Filename for second-order correlation of cavity a
CHARACTER(LEN=34), PARAMETER :: filename_g2 = "./data_files/g2_corr_and_cross.txt"
! Paramert Name List
CHARACTER(LEN=15), PARAMETER :: filename_ParamList = "./ParamList.nml"

!------------------------------------------------------------------------------!
!                    Calling parameters from ParamList_me.nml                  !
!------------------------------------------------------------------------------!
! NameList things
! Status and unit integers
INTEGER :: ISTAT, IUNIT
! Line to be read from file
CHARACTER(LEN=512) :: LINE
! Namelist parameters
NAMELIST /ATOM/ Gamma, Omega
NAMELIST /FILTER/ epsilon, N, phase
NAMELIST /CAVITYA/ kappaa, w0a, dwa
NAMELIST /CAVITYB/ kappab, w0b, dwb
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

! Max number of time steps
t_steps = NINT(t_max / dt)

! omega shortcut
IF (0.25 * gamma > Omega) THEN
  delta = SQRT(((0.25d0 * gamma) ** 2) - (Omega ** 2))
ELSE
  delta = i * SQRT((Omega ** 2) - ((0.25d0 * gamma) ** 2))
END IF

!------------------------------------------------------------------------------!
!                       Allocate and initialise arrays                         !
!------------------------------------------------------------------------------!
! Kappa values
kappa = 0.0d0
kappa(1) = kappaa; kappa(2) = kappab
! Allocate array of omega and gka values
ALLOCATE(wl(-N:N, 2)); wl = 0.0d0
ALLOCATE(gkl(-N:N, 2)); gkl = 0.0d0
DO j = -N, N
  IF (N == 0) THEN
    ! Cavity a
    wl(j, a) = w0a
    gkl(j, a) = DSQRT(0.5d0 * epsilon * gamma * kappa(a))
    ! Cavity b
    wl(j, b) = w0b
    gkl(j, b) = DSQRT(0.5d0 * epsilon * gamma * kappa(b))
  ELSE
    ! Blackman window coefficient
    blackman = 1.0d0
    ! blackman = 0.42d0 - 0.5d0 * COS(2.0d0 * pi * DBLE(N + ja) / (2.0d0 * DBLE(N))) + &
    !          & 0.08d0 * COS(4.0d0 * pi * DBLE(N + ja) / (2.0d0 * DBLE(N)))
    ! Cavity a
    wl(j, a) = w0a + DBLE(j) * dwa
    gkl(j, a) = DSQRT(0.5d0 * (epsilon / DBLE(2*N + 1)) * gamma * kappa(a)) * blackman * EXP(i * DBLE(phase) * DBLE(j) * pi / DBLE(N))
    ! Cavity a
    wl(j, b) = w0b + DBLE(j) * dwb
    gkl(j, b) = DSQRT(0.5d0 * (epsilon / DBLE(2*N + 1)) * gamma * kappa(b)) * blackman * EXP(i * DBLE(phase) * DBLE(j) * pi / DBLE(N))
  END IF
END DO

! Initialise moments
! First-order
bloch_ss = 0.0d0
ALLOCATE(cav1_ss(-N:N, 4, 2)); cav1_ss = 0.0d0
! Second-order
ALLOCATE(cav2_ss(-N:N, -N:N, 6, 2, 2)); cav2_ss = 0.0d0
ALLOCATE(cavsig2_ss(-N:N, 4, 3, 2)); cavsig2_ss = 0.0d0
! Third-order
ALLOCATE(cav3_ss(-N:N, -N:N, -N:N, 2, 2, 2, 2)); cav3_ss = 0.0d0
ALLOCATE(cavsig3_ss(-N:N, -N:N, 3, 3, 2, 2)); cavsig3_ss = 0.0d0;
! Fourth-order
ALLOCATE(cav4_ss(-N:N, -N:N, -N:N, -N:N, 2, 2, 2, 2)); cav4_ss = 0.0d0
ALLOCATE(cavsig4_ss(-N:N, -N:N, -N:N, 2, 3, 2, 2, 2)); cavsig4_ss = 0.0d0

! Open file to write time to
OPEN(UNIT=1, FILE=filename_parameters, STATUS='REPLACE', ACTION='WRITE')
! Write parameter
! WRITE(1,*) "Parameters are in the following order:"
! WRITE(1,"(A10,F25.15)") "gamma =", gamma
! WRITE(1,"(A10,F25.15)") "Omega =", Omega
! WRITE(1,"(A10,F25.15)") "w0_a =", w0a
! WRITE(1,"(A10,F25.15)") "w0_b =", w0b
! WRITE(1,"(A10,F25.15)") "kappa_a =", kappaa
! WRITE(1,"(A10,F25.15)") "kappa_b =", kappab
! WRITE(1,"(A10,F25.15)") "epsilon =", epsilon
! WRITE(1,"(A11,I9)") "N = ", N
! WRITE(1,"(A11,F25.15)") "dw_a = ", dwa
! WRITE(1,"(A11,F25.15)") "dw_b = ", dwb
! WRITE(1,"(A11,I9)") "phase = ", phase
! WRITE(1,"(A10,F25.15)") "dt =", dt
! IF (t_max <= 1000.0d0) THEN
!   WRITE(1,"(A10,F25.15)") "Max time =", t_max
! ELSE
!   WRITE(1,"(A10,ES25.11E2)") "Max time =", t_max
! END IF
! WRITE(1,"(A10,F25.15)") "Max tau1 =", tau1_max
! WRITE(1,"(A10,F25.15)") "Max tau2 =", tau2_max
! ! Close file
! CLOSE(1)

WRITE(1,*) "Parameters are in the following order:"
WRITE(1,"(A10,F25.15)") "Gamma =", Gamma
WRITE(1,"(A10,F25.15)") "Omega =", Omega
WRITE(1,"(A10,F25.15)") "w0a =", w0a
WRITE(1,"(A10,F25.15)") "kappaa =", kappaa
WRITE(1,"(A11,F25.15)") "dwa = ", dwa
WRITE(1,"(A10,F25.15)") "w0b =", w0b
WRITE(1,"(A10,F25.15)") "kappab =", kappab
WRITE(1,"(A11,F25.15)") "dwb = ", dwb
WRITE(1,"(A10,F25.15)") "epsilon =", epsilon
WRITE(1,"(A11,I9)") "N = ", N
WRITE(1,"(A11,I9)") "phase = ", phase
WRITE(1,"(A10,F25.15)") "dt =", dt
WRITE(1,"(A10,F25.15)") "Max tau1 =", tau1_max
! Close file
CLOSE(1)

!------------------------------------------------------------------------------!
!                               STEADY STATES                                  !
!------------------------------------------------------------------------------!

! Calculate steady state expressions
! First-order: Bloch equations
bloch_ss(sm) = -i * gamma * Omega / ((2.0d0 * Omega ** 2) + (gamma ** 2))
bloch_ss(sp) = i * gamma * Omega / ((2.0d0 * Omega ** 2) + (gamma ** 2))
bloch_ss(sz) = -(gamma ** 2) / ((2.0d0 * Omega ** 2) + (gamma ** 2))
! Cycle through modes
DO Cj = 1, 2
  DO j = -N, N
    !-----------------------------!
    !     First-order: Cavity     !
    !-----------------------------!
    cav1_ss(j, a, Cj) = - gkl(j, Cj) * bloch_ss(sm) / (kappa(Cj) + i * wl(j, Cj))
    cav1_ss(j, at, Cj) = - CONJG(gkl(j, Cj)) * bloch_ss(sp) / (kappa(Cj) - i * wl(j, Cj))
    !---------------------------------------!
    !     Second-order: Cavity and atom     !
    !---------------------------------------!
    ! <a \sigma>
    alpha = -(0.5d0 * gamma + kappa(Cj) + i * wl(j, Cj))
    ! Set matrix elements
    M_inverse = INVERT_MATRIX(gamma, Omega, alpha)
    B_ss = 0.0d0
    B_ss(1) = 0.0d0
    B_ss(2) = -0.5 * gkl(j, Cj) * (bloch_ss(sz) + 1.0d0)
    B_ss(3) = -gamma * cav1_ss(j, a, Cj) + &
            & gkl(j, Cj) * bloch_ss(sm)
    ! Matrix multiplication
    cavsig2_ss(j, a, sm, Cj) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
    cavsig2_ss(j, a, sp, Cj) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
    cavsig2_ss(j, a, sz, Cj) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))

    ! <at \sigma>
    alpha = -(0.5d0 * gamma + kappa(Cj) - i * wl(j, Cj))
    ! Set matrix elements
    M_inverse = INVERT_MATRIX(gamma, Omega, alpha)
    B_ss = 0.0d0
    B_ss(1) = -0.5 * CONJG(gkl(j, Cj)) * (bloch_ss(sz) + 1.0d0)
    B_ss(2) = 0.0d0
    B_ss(3) = -gamma * cav1_ss(j, at, Cj) + &
            & CONJG(gkl(j, Cj)) * bloch_ss(sp)
    ! Matrix multiplication
    cavsig2_ss(j, at, sm, Cj) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
    cavsig2_ss(j, at, sp, Cj) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
    cavsig2_ss(j, at, sz, Cj) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))
  END DO
END DO

! Cycle through modes
DO Ck = 1, 2
  DO Cj = 1, 2
    ! Cycle through modes
    DO j = -N, N
      DO k = -N, N
        !------------------------------!
        !     Second-order: Cavity     !
        !------------------------------!
        ! <a a>
        cav2_ss(j, k, a, Cj, Ck) = -(gkl(j, Cj) * cavsig2_ss(k, a, sm, Ck) + gkl(k, Ck) * cavsig2_ss(j, a, sm, Cj)) / &
                         & ((kappa(Cj) + kappa(Ck)) + i * (wl(j, Cj) + wl(k, Ck)))
        ! <at at>
        cav2_ss(j, k, at, Cj, Ck) = -(CONJG(gkl(j, Cj)) * cavsig2_ss(k, at, sp, Ck) + CONJG(gkl(k, Ck)) * cavsig2_ss(j, at, sp, Cj)) / &
                          & ((kappa(Cj) + kappa(Ck)) - i * (wl(j, Cj) + wl(k, Ck)))
        ! <at a>
        cav2_ss(j, k, ata, Cj, Ck) = -(CONJG(gkl(j, Cj)) * cavsig2_ss(k, a, sp, Ck) + gkl(k, Ck) * cavsig2_ss(j, at, sm, Cj)) / &
                           & ((kappa(Cj) + kappa(Ck)) - i * (wl(j, Cj) - wl(k, Ck)))

        !--------------------------------------!
        !     Third-Order: Cavity and Atom     !
        !--------------------------------------!
        ! <a a \sigma>
        alpha = -(0.5d0 * gamma + (kappa(Cj) + kappa(Ck)) + i * (wl(j, Cj) + wl(k, Ck)))
        ! Set Matrix elements
        M_inverse = INVERT_MATRIX(gamma, Omega, alpha)
        B_ss(1) = 0.0d0
        B_ss(2) = -0.5d0 * gkl(j, Cj) * (cavsig2_ss(k, a, sz, Ck) + cav1_ss(k, a, Ck)) + &
                & -0.5d0 * gkl(k, Ck) * (cavsig2_ss(j, a, sz, Cj) + cav1_ss(j, a, Cj))
        B_ss(3) = -gamma * cav2_ss(j, k, a, Cj, Ck) + &
                & gkl(j, Cj) * cavsig2_ss(k, a, sm, Ck) + &
                & gkl(k, Ck) * cavsig2_ss(j, a, sm, Cj)
        ! Matrix multiplication
        cavsig3_ss(j, k, a, sm, Cj, Ck) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
        cavsig3_ss(j, k, a, sp, Cj, Ck) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
        cavsig3_ss(j, k, a, sz, Cj, Ck) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))

        ! <at at \sigma>
        alpha = -(0.5d0 * gamma + (kappa(Cj) + kappa(Ck)) - i * (wl(j, Cj) + wl(k, Ck)))
        ! Set Matrix elements
        M_inverse = INVERT_MATRIX(gamma, Omega, alpha)
        B_ss(1) = -0.5d0 * CONJG(gkl(j, Cj)) * (cavsig2_ss(k, at, sz, Ck) + cav1_ss(k, at, Ck)) + &
                & -0.5d0 * CONJG(gkl(k, Ck)) * (cavsig2_ss(j, at, sz, Cj) + cav1_ss(j, at, Cj))
        B_ss(2) = 0.0d0
        B_ss(3) = -gamma * cav2_ss(j, k, at, Cj, Ck) + &
                & CONJG(gkl(j, Cj)) * cavsig2_ss(k, at, sp, Ck) + &
                & CONJG(gkl(k, Ck)) * cavsig2_ss(j, at, sp, Cj)
        ! Matrix multiplication
        cavsig3_ss(j, k, at, sm, Cj, Ck) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
        cavsig3_ss(j, k, at, sp, Cj, Ck) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
        cavsig3_ss(j, k, at, sz, Cj, Ck) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))

        ! <at a \sigma>
        alpha = -(0.5d0 * gamma + (kappa(Cj) + kappa(Ck)) - i * (wl(j, Cj) - wl(k, Ck)))
        ! Set Matrix elements
        M_inverse = INVERT_MATRIX(gamma, Omega, alpha)
        B_ss(1) = -0.5d0 * CONJG(gkl(j, Cj)) * (cavsig2_ss(k, a, sz, Ck) + cav1_ss(k, a, Ck))
        B_ss(2) = -0.5d0 * gkl(k, Ck) * (cavsig2_ss(j, at, sz, Cj) + cav1_ss(j, at, Cj))
        B_ss(3) = -gamma * cav2_ss(j, k, ata, Cj, Ck) + &
                & CONJG(gkl(j, Cj)) * cavsig2_ss(k, a, sp, Ck) + &
                & gkl(k, Ck) * cavsig2_ss(j, at, sm, Cj)
        ! Matrix multiplication
        cavsig3_ss(j, k, ata, sm, Cj, Ck) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
        cavsig3_ss(j, k, ata, sp, Cj, Ck) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
        cavsig3_ss(j, k, ata, sz, Cj, Ck) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))
      END DO
    END DO
  END DO
END DO

DO Cl = 1, 2
  DO Ck = 1, 2
    DO Cj = 1, 2
      ! Cycle through modes
      DO j = -N, N
        DO k = -N, N
          DO l = -N, N
            !-----------------------------!
            !     Third-Order: Cavity     !
            !-----------------------------!
            ! <at a a>
            cav3_ss(j, k, l, a, Cj, Ck, Cl) = -(CONJG(gkl(j, Cj)) * cavsig3_ss(k, l, a, sp, Ck, Cl) + &
                                & gkl(k, Ck) * cavsig3_ss(j, l, ata, sm, Cj, Cl) + &
                                & gkl(l, Cl) * cavsig3_ss(j, k, ata, sm, Cj, Ck)) / &
                                & ((kappa(Cj) + kappa(Ck) + kappa(Cl)) - i * (wl(j, Cj) - wl(k, Ck) - wl(l, Cl)))
            ! <at at a>
            cav3_ss(j, k, l, at, Cj, Ck, Cl) = -(CONJG(gkl(j, Cj)) * cavsig3_ss(k, l, ata, sp, Ck, Cl) + &
                                 & CONJG(gkl(k, Ck)) * cavsig3_ss(j, l, ata, sp, Cj, Cl) + &
                                 & gkl(l, Cl) * cavsig3_ss(j, k, at, sm, Cj, Ck)) / &
                                 & ((kappa(Cj) + kappa(Ck) + kappa(Cl)) - i * (wl(j, Cj) + wl(k, Ck) - wl(l, Cl)))

            !---------------------------------------!
            !     Fourth-Order: Cavity and Atom     !
            !---------------------------------------!
            ! <at a a \sigma>
            alpha = -(0.5d0 * gamma + (kappa(Cj) + kappa(Ck) + kappa(Cl)) - i * (wl(j, Cj) - wl(k, Ck) - wl(l, Cl)))
            ! Set Matrix elements
            M_inverse = INVERT_MATRIX(gamma, Omega, alpha)
            B_ss(1) = -0.5d0 * CONJG(gkl(j, Cj)) * (cavsig3_ss(k, l, a, sz, Ck, Cl) + cav2_ss(k, l, a, Ck, Cl))
            B_ss(2) = -0.5d0 * gkl(k, Ck) * (cavsig3_ss(j, l, ata, sz, Cj, Cl) + cav2_ss(j, l, ata, Cj, Cl)) + &
                    & -0.5d0 * gkl(l, Cl) * (cavsig3_ss(j, k, ata, sz, Cj, Cl) + cav2_ss(j, k, ata, Cj, Ck))
            B_ss(3) = -gamma * cav3_ss(j, k, l, a, Cj, Ck, Cl) + &
                    & CONJG(gkl(j, Cj)) * cavsig3_ss(k, l, a, sp, Ck, Cl) + &
                    & gkl(k, Ck) * cavsig3_ss(j, l, ata, sm, Cj, Cl) + &
                    & gkl(l, Cl) * cavsig3_ss(j, k, ata, sm, Cj, Ck)
            ! Matrix multiplication
            cavsig4_ss(j, k, l, a, sm, Cj, Ck, Cl) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
            cavsig4_ss(j, k, l, a, sp, Cj, Ck, Cl) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
            cavsig4_ss(j, k, l, a, sz, Cj, Ck, Cl) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))

            ! <at at a \sigma>
            alpha = -(0.5d0 * gamma + (kappa(Cj) + kappa(Ck) + kappa(Cl)) - i * (wl(j, Cj) + wl(k, Ck) - wl(l, Cl)))
            ! Set Matrix elements
            M_inverse = INVERT_MATRIX(gamma, Omega, alpha)
            B_ss(1) = -0.5d0 * CONJG(gkl(j, Cj)) * (cavsig3_ss(k, l, ata, sz, Ck, Cl) + cav2_ss(k, l, ata, Ck, Cl)) + &
                    & -0.5d0 * CONJG(gkl(k, Ck)) * (cavsig3_ss(j, l, ata, sz, Cj, Cl) + cav2_ss(j, l, ata, Cj, Cl))
            B_ss(2) = -0.5d0 * gkl(l, Cl) * (cavsig3_ss(j, k, at, sz, Cj, Ck) + cav2_ss(j, k, at, Cj, Ck))
            B_ss(3) = -gamma * cav3_ss(j, k, l, at, Cj, Ck, Cl) + &
                    & CONJG(gkl(j, Cj)) * cavsig3_ss(k, l, ata, sp, Ck, Cl) + &
                    & CONJG(gkl(k, Ck)) * cavsig3_ss(j, l, ata, sp, Cj, Cl) + &
                    & gkl(l, Cl) * cavsig3_ss(j, k, at, sm, Cj, Ck)
            ! Matrix multiplication
            cavsig4_ss(j, k, l, at, sm, Cj, Ck, Cl) = -(M_inverse(1, 1) * B_ss(1) + M_inverse(1, 2) * B_ss(2) + M_inverse(1, 3) * B_ss(3))
            cavsig4_ss(j, k, l, at, sp, Cj, Ck, Cl) = -(M_inverse(2, 1) * B_ss(1) + M_inverse(2, 2) * B_ss(2) + M_inverse(2, 3) * B_ss(3))
            cavsig4_ss(j, k, l, at, sz, Cj, Ck, Cl) = -(M_inverse(3, 1) * B_ss(1) + M_inverse(3, 2) * B_ss(2) + M_inverse(3, 3) * B_ss(3))

          END DO
        END DO
      END DO
    END DO
  END DO
END DO

DO Cm = 1, 2
  DO Cl = 1, 2
    DO Ck = 1, 2
      DO Cj = 1, 2
        ! Cycle through modes
        DO j = -N, N
          DO k = -N, N
            DO l = -N, N
              DO m = -N, N
                !------------------------------!
                !     Fourth-Order: Cavity     !
                !------------------------------!
                ! <at at a a>
                cav4_ss(j, k, l, m, Cj, Ck, Cl, Cm) = -(CONJG(gkl(j, Cj)) * cavsig4_ss(k, l, m, a, sp, Ck, Cl, Cm) + &
                                    & CONJG(gkl(k, Ck)) * cavsig4_ss(j, l, m, a, sp, Cj, Cl, Cm) + &
                                    & gkl(l, Cl) * cavsig4_ss(j, k, m, at, sm, Cj, Ck, Cm) + &
                                    & gkl(m, Cm) * cavsig4_ss(j, k, l, at, sm, Cj, Ck, Cl)) / &
                                    & ((kappa(Cj) + kappa(Ck) + kappa(Cl) + kappa(Cm)) - i * (wl(j, Cj) + wl(k, Ck) - wl(l, Cl) - wl(m, Cm)))
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
END DO

! Steady state mean photon number
photonA_ss = 0.0d0
photonB_ss = 0.0d0
DO j = -N, N
  DO k = -N, N
    ! Steady state mean photon number
    photonA_ss = photonA_ss + cav2_ss(j, k, ata, a, a)
    photonB_ss = photonB_ss + cav2_ss(j, k, ata, b, b)
  END DO
END DO

PRINT*, "<A^{\dagger} A>_{ss} = ", photonA_ss
PRINT*, "<B^{\dagger} B>_{ss} = ", photonB_ss

!------------------------------------------------------------------------------!
!                                STATE-SOLVER                                  !
!------------------------------------------------------------------------------!
! Solve for the ground and excited state population and the mean photon number

! Calculate second-order correlation initial and steady states
! First-order: Atom
sm0 = 0.0d0; sp0 = 0.0d0; sz0 = 0.0d0
spss = 0.0d0; szss = 0.0d0; szss = 0.0d0
! First-order: Cavity
ALLOCATE(aj0(-N:N, 3)); aj0 = 0.0d0; ALLOCATE(ajss(-N:N, 3)); ajss = 0.0d0
! Second-order: Cavity and Atom
ALLOCATE(ajsm0(-N:N, 3)); ajsm0 = 0.0d0; ALLOCATE(ajsmss(-N:N, 3)); ajsmss = 0.0d0
ALLOCATE(ajsp0(-N:N, 3)); ajsp0 = 0.0d0; ALLOCATE(ajspss(-N:N, 3)); ajspss = 0.0d0
ALLOCATE(ajsz0(-N:N, 3)); ajsz0 = 0.0d0; ALLOCATE(ajszss(-N:N, 3)); ajszss = 0.0d0
! Second-order: Cavity
ALLOCATE(atjak0(-N:N, -N:N, 3)); atjak0 = 0.0d0; ALLOCATE(atjakss(-N:N, -N:N, 3)); atjakss = 0.0d0

! Steady States - Auto-correlation of Cavity A
smss(1) = bloch_ss(sm) * photonA_ss
spss(1) = bloch_ss(sp) * photonA_ss
szss(1) = bloch_ss(sz) * photonA_ss
! Steady States - Auto-correlation of Cavity B
smss(2) = bloch_ss(sm) * photonB_ss
spss(2) = bloch_ss(sp) * photonB_ss
szss(2) = bloch_ss(sz) * photonB_ss
! Steady States - Cross-correlation
smss(3) = smss(1); spss(3) = spss(1); szss(3) = szss(1)
! Cycle through modes
DO j = -N, N
  ! Steady States - Auto-correlation of Cavity A
  ajss(j, 1) = cav1_ss(j, a, A) * photonA_ss
  ajsmss(j, 1) = cavsig2_ss(j, a, sm, A) * photonA_ss
  ajspss(j, 1) = cavsig2_ss(j, a, sp, A) * photonA_ss
  ajszss(j, 1) = cavsig2_ss(j, a, sz, A) * photonA_ss
  ! Steady States - Auto-correlation of Cavity B
  ajss(j, 2) = cav1_ss(j, a, B) * photonB_ss
  ajsmss(j, 2) = cavsig2_ss(j, a, sm, B) * photonB_ss
  ajspss(j, 2) = cavsig2_ss(j, a, sp, B) * photonB_ss
  ajszss(j, 2) = cavsig2_ss(j, a, sz, B) * photonB_ss
  ! Steady States - Cross-correlation
  ajss(j, 3) = cav1_ss(j, a, B) * photonA_ss
  ajsmss(j, 3) = cavsig2_ss(j, a, sm, B) * photonA_ss
  ajspss(j, 3) = cavsig2_ss(j, a, sp, B) * photonA_ss
  ajszss(j, 3) = cavsig2_ss(j, a, sz, B) * photonA_ss
  DO m = -N, N
    ! Steady states - Auto-correlation of Cavity A
    atjakss(j, m, 1) = cav2_ss(j, m, ata, A, A) * photonA_ss
    ! Steady states - Auto-correlation of Cavity B
    atjakss(j, m, 2) = cav2_ss(j, m, ata, B, B) * photonB_ss
    ! Steady states - Cross-correlation
    atjakss(j, m, 3) = cav2_ss(j, m, ata, B, B) * photonA_ss
    ! Initial states - Auto-correlation of Cavity A
    sm0(1) = sm0(1) + cavsig3_ss(j, m, ata, sm, A, A)
    sp0(1) = sp0(1) + cavsig3_ss(j, m, ata, sp, A, A)
    sz0(1) = sz0(1) + cavsig3_ss(j, m, ata, sz, A, A)
    ! Initial states - Auto-correlation of Cavity B
    sm0(2) = sm0(2) + cavsig3_ss(j, m, ata, sm, B, B)
    sp0(2) = sp0(2) + cavsig3_ss(j, m, ata, sp, B, B)
    sz0(2) = sz0(2) + cavsig3_ss(j, m, ata, sz, B, B)
    ! Initial states - Cross-correlation
    sm0(3) = sm0(3) + cavsig3_ss(j, m, ata, sm, A, A)
    sp0(3) = sp0(3) + cavsig3_ss(j, m, ata, sp, A, A)
    sz0(3) = sz0(3) + cavsig3_ss(j, m, ata, sz, A, A)
    ! Cycle through modes
    DO l = -N, N
      ! Initial states - Auto-correlation of Cavity A
      aj0(l, 1) = aj0(l, 1) + cav3_ss(j, l, m, a, A, A, A)
      ajsm0(l, 1) = ajsm0(l, 1) + cavsig4_ss(j, l, m, a, sm, A, A, A)
      ajsp0(l, 1) = ajsp0(l, 1) + cavsig4_ss(j, l, m, a, sp, A, A, A)
      ajsz0(l, 1) = ajsz0(l, 1) + cavsig4_ss(j, l, m, a, sz, A, A, A)
      ! Initial states - Auto-correlation of Cavity B
      aj0(l, 2) = aj0(l, 2) + cav3_ss(j, l, m, a, B, B, B)
      ajsm0(l, 2) = ajsm0(l, 2) + cavsig4_ss(j, l, m, a, sm, B, B, B)
      ajsp0(l, 2) = ajsp0(l, 2) + cavsig4_ss(j, l, m, a, sp, B, B, B)
      ajsz0(l, 2) = ajsz0(l, 2) + cavsig4_ss(j, l, m, a, sz, B, B, B)
      ! Initial states - Cross-correlation
      aj0(l, 3) = aj0(l, 3) + cav3_ss(j, l, m, a, A, B, A)
      ajsm0(l, 3) = ajsm0(l, 3) + cavsig4_ss(j, l, m, a, sm, A, B, A)
      ajsp0(l, 3) = ajsp0(l, 3) + cavsig4_ss(j, l, m, a, sp, A, B, A)
      ajsz0(l, 3) = ajsz0(l, 3) + cavsig4_ss(j, l, m, a, sz, A, B, A)
      ! Cycle through modes
      DO k = -N, N
        ! Initial states - Auto-correlation of Cavity A
        atjak0(k, l, 1) = atjak0(k, l, 1) + cav4_ss(j, k, l, m, A, A, A, A)
        ! Initial states - Auto-correlation of Cavity B
        atjak0(k, l, 2) = atjak0(k, l, 2) + cav4_ss(j, k, l, m, B, B, B, B)
        ! Initial states - Cross-correlation
        atjak0(k, l, 3) = atjak0(k, l, 3) + cav4_ss(j, k, l, m, A, B, B, A)
      END DO
    END DO
  END DO
END DO

! Calculate Differential Coefficients
! Bloch equations
C1 = 0.0d0; C2 = 0.0d0; C3 = 0.0d0
! Auto-correlation of Cavity A
C1(1) = 0.5d0 * (sm0(1) - smss(1)) + 0.5d0 * (sp0(1) - spss(1))
C2(1) = (-i * (gamma + 4.0d0 * delta) * (gamma - 4.0d0 * delta) / (32.0d0 * delta * Omega)) * &
      & ((sm0(1) - smss(1)) - (sp0(1) - spss(1))) + &
      & ((gamma + 4.0d0 * delta) / (8.0d0 * delta)) * (sz0(1) - szss(1))
C3(1) = (i * (gamma + 4.0d0 * delta) * (gamma - 4.0d0 * delta) / (32.0d0 * delta * Omega)) * &
      & ((sm0(1) - smss(1)) - (sp0(1) - spss(1))) - &
      & ((gamma - 4.0d0 * delta) / (8.0d0 * delta)) * (sz0(1) - szss(1))
! Auto-correlation of Cavity B
C1(2) = 0.5d0 * (sm0(2) - smss(2)) + 0.5d0 * (sp0(2) - spss(2))
C2(2) = (-i * (gamma + 4.0d0 * delta) * (gamma - 4.0d0 * delta) / (32.0d0 * delta * Omega)) * &
      & ((sm0(2) - smss(2)) - (sp0(2) - spss(2))) + &
      & ((gamma + 4.0d0 * delta) / (8.0d0 * delta)) * (sz0(2) - szss(2))
C3(2) = (i * (gamma + 4.0d0 * delta) * (gamma - 4.0d0 * delta) / (32.0d0 * delta * Omega)) * &
      & ((sm0(2) - smss(2)) - (sp0(2) - spss(2))) - &
      & ((gamma - 4.0d0 * delta) / (8.0d0 * delta)) * (sz0(2) - szss(2))
! Cross-correlation
C1(3) = C1(1); C2(3) = C2(1); C3(3) = C3(1)

! Cheaty variables
ALLOCATE(va(-N:N, 0:3, 3)); va = 0.0d0
ALLOCATE(vb(-N:N, 0:4, 3)); vb = 0.0d0
ALLOCATE(vc(-N:N, 0:4, 3)); vc = 0.0d0
! First-order cavity
ALLOCATE(C4(-N:N, 3)); C4 = 0.0d0
! Second-order cavity and atom
ALLOCATE(C6(-N:N, 3)); C6 = 0.0d0
ALLOCATE(C7(-N:N, 3)); C7 = 0.0d0
ALLOCATE(C8(-N:N, 3)); C8 = 0.0d0
DO j = -N, N
  !----------------------!
  !     Coefficients     !
  !----------------------!
  ! C4 - Auto-correlation of Cavity A
  C4(j, 1) = aj0(j, 1) - ajss(j, 1) - &
           & (gkl(j, A) * C1(1) / (0.5d0 * gamma - (kappa(A) + i * wl(j, A)))) + &
           & (2.0d0 * i * Omega * gkl(j, A) * C2(1) / ((gamma + 4.0d0 * delta) * (0.75 * gamma + delta - (kappa(A) + i * wl(j, A))))) + &
           & (2.0d0 * i * Omega * gkl(j, A) * C3(1) / ((gamma - 4.0d0 * delta) * (0.75 * gamma - delta - (kappa(A) + i * wl(j, A)))))
  ! C4 - Auto-correlation of Cavity B
  C4(j, 2) = aj0(j, 2) - ajss(j, 2) - &
           & (gkl(j, B) * C1(2) / (0.5d0 * gamma - (kappa(B) + i * wl(j, B)))) + &
           & (2.0d0 * i * Omega * gkl(j, B) * C2(2) / ((gamma + 4.0d0 * delta) * (0.75 * gamma + delta - (kappa(B) + i * wl(j, B))))) + &
           & (2.0d0 * i * Omega * gkl(j, B) * C3(2) / ((gamma - 4.0d0 * delta) * (0.75 * gamma - delta - (kappa(B) + i * wl(j, B)))))
  ! C4 - Cross-correlation
  C4(j, 3) = aj0(j, 3) - ajss(j, 3) - &
           & (gkl(j, B) * C1(3) / (0.5d0 * gamma - (kappa(B) + i * wl(j, B)))) + &
           & (2.0d0 * i * Omega * gkl(j, B) * C2(3) / ((gamma + 4.0d0 * delta) * (0.75 * gamma + delta - (kappa(B) + i * wl(j, B))))) + &
           & (2.0d0 * i * Omega * gkl(j, B) * C3(3) / ((gamma - 4.0d0 * delta) * (0.75 * gamma - delta - (kappa(B) + i * wl(j, B)))))

  !-------------------------------------------------------------------!
  !     Cheaty Variables (a, b, c) - Auto-correlation of Cavity A     !
  !-------------------------------------------------------------------!
  ! Cheaty alpha
  alpha = -(0.5d0 * gamma + kappa(A) + i * wl(j, A))

  ! a variable
  va(j, 0, 1) = (szss(1) + 1.0d0) / (alpha)
  va(j, 2, 1) = 1.0d0 / (0.75d0 * gamma + delta + alpha)
  va(j, 3, 1) = 1.0d0 / (0.75d0 * gamma - delta + alpha)
  ! b variable
  vb(j, 0, 1) = (Omega * gkl(j, A) * smss(1) - &
              & 0.125d0 * i * gkl(j, A) * (gamma - 4.0d0 * delta) * (szss(1) + 1.0d0) - &
              & gamma * Omega * ajss(j, 1)) / &
              & (0.25d0 * gamma + delta - alpha)
  vb(j, 1, 1) = -Omega * gkl(j, A) * (1.0d0 - (gamma / (0.5d0 * gamma - (kappa(A) + i * wl(j, A))))) / &
              & (0.25d0 * gamma - delta + alpha)
  vb(j, 2, 1) = -i * gkl(j, A) * (0.125d0 * (gamma - 4.0d0 * delta) + &
              & (2.0d0 * (Omega ** 2) / (gamma + 4.0d0 * delta)) * &
              & (1.0d0 - (gamma / (0.75d0 * gamma + delta - (kappa(A) + i * wl(j, A)))))) / &
              & (kappa(A) + i * wl(j, A))
  vb(j, 3, 1) = -i * gkl(j, A) * (0.125d0 * (gamma - 4.0d0 * delta) + &
              & (2.0d0 * (Omega ** 2) / (gamma - 4.0d0 * delta)) * &
              & (1.0d0 - (gamma / (0.75d0 * gamma - delta - (kappa(A) + i * wl(j, A)))))) / &
              & (2.0d0 * delta + kappa(A) + i * wl(j, A))
  vb(j, 4, 1) = -gamma * Omega / (0.75d0 * gamma + delta)
  ! c variable
  vc(j, 0, 1) = (Omega * gkl(j, A) * smss(1) - &
              & 0.125d0 * i * gkl(j, A) * (gamma + 4.0d0 * delta) * (szss(1) + 1.0d0) - &
              & gamma * Omega * ajss(j, 1)) / &
              & (0.25d0 * gamma - delta - alpha)
  vc(j, 1, 1) = -Omega * gkl(j, A) * (1.0d0 - (gamma / (0.5d0 * gamma - (kappa(A) + i * wl(j, A))))) / &
              & (0.25d0 * gamma + delta + alpha)
  vc(j, 2, 1) = i * gkl(j, A) * (0.125d0 * (gamma + 4.0d0 * delta) + &
              & (2.0d0 * (Omega ** 2) / (gamma + 4.0d0 * delta)) * &
              & (1.0d0 - (gamma / (0.75d0 * gamma + delta - (kappa(A) + i * wl(j, A)))))) / &
              & (2.0d0 * delta - (kappa(A) + i * wl(j, A)))
  vc(j, 3, 1) = -i * gkl(j, A) * (0.125d0 * (gamma + 4.0d0 * delta) + &
              & (2.0d0 * (Omega ** 2) / (gamma - 4.0d0 * delta)) * &
              & (1.0d0 - (gamma / (0.75d0 * gamma - delta - (kappa(A) + i * wl(j, A)))))) / &
              & (kappa(A) + i * wl(j, A))
  vc(j, 4, 1) = -gamma * Omega / (0.75d0 * gamma - delta)
  ! Initial values for v1(t), v2(t), and v3(t)
  v1(1) = 0.0d0
  v1(1) = (-i * (vb(j, 1, 1) - vc(j, 1, 1)) * C1(1) + &
        & (gkl(j, A) * delta * va(j, 2, 1) - i * (vb(j, 2, 1) - vc(j, 2, 1))) * C2(1) + &
        & (gkl(j, A) * delta * va(j, 3, 1) - i * (vb(j, 3, 1) - vc(j, 3, 1))) * C3(1) + &
        & -i * (vb(j, 4, 1) - vc(j, 4, 1)) * C4(j, 1)) / &
        & (4.0d0 * delta)

  v2(1) = 0.0d0
  v2(1) = (i * (vb(j, 1, 1) - vc(j, 1, 1)) * C1(1) + &
        & (gkl(j, A) * delta * va(j, 2, 1) + i * (vb(j, 2, 1) - vc(j, 2, 1))) * C2(1) + &
        & (gkl(j, A) * delta * va(j, 3, 1) + i * (vb(j, 3, 1) - vc(j, 3, 1))) * C3(1) + &
        & i * (vb(j, 4, 1) - vc(j, 4, 1)) * C4(j, 1)) / &
        & (4.0d0 * delta)

  v3(1) = 0.0d0
  v3(1) = (vb(j, 1, 1) * (gamma + 4.0d0 * delta) - vc(j, 1, 1) * (gamma - 4.0d0 * delta)) * C1(1) + &
        & (vb(j, 2, 1) * (gamma + 4.0d0 * delta) - vc(j, 2, 1) * (gamma - 4.0d0 * delta)) * C2(1) + &
        & (vb(j, 3, 1) * (gamma + 4.0d0 * delta) - vc(j, 3, 1) * (gamma - 4.0d0 * delta)) * C3(1) + &
        & (vb(j, 4, 1) * (gamma + 4.0d0 * delta) - vc(j, 4, 1) * (gamma - 4.0d0 * delta)) * C4(j, 1)
  v3(1) = 0.125d0 * v3(1) / (Omega * delta)

  !-------------------------------------------------------------------!
  !     Cheaty Variables (a, b, c) - Auto-correlation of Cavity B     !
  !-------------------------------------------------------------------!
  ! Cheaty alpha
  alpha = -(0.5d0 * gamma + kappa(B) + i * wl(j, B))

  ! a variable
  va(j, 0, 2) = (szss(2) + 1.0d0) / (alpha)
  va(j, 2, 2) = 1.0d0 / (0.75d0 * gamma + delta + alpha)
  va(j, 3, 2) = 1.0d0 / (0.75d0 * gamma - delta + alpha)
  ! b variable
  vb(j, 0, 2) = (Omega * gkl(j, B) * smss(2) - &
              & 0.125d0 * i * gkl(j, B) * (gamma - 4.0d0 * delta) * (szss(2) + 1.0d0) - &
              & gamma * Omega * ajss(j, 2)) / &
              & (0.25d0 * gamma + delta - alpha)
  vb(j, 1, 2) = -Omega * gkl(j, B) * (1.0d0 - (gamma / (0.5d0 * gamma - (kappa(B) + i * wl(j, B))))) / &
              & (0.25d0 * gamma - delta + alpha)
  vb(j, 2, 2) = -i * gkl(j, B) * (0.125d0 * (gamma - 4.0d0 * delta) + &
              & (2.0d0 * (Omega ** 2) / (gamma + 4.0d0 * delta)) * &
              & (1.0d0 - (gamma / (0.75d0 * gamma + delta - (kappa(B) + i * wl(j, B)))))) / &
              & (kappa(B) + i * wl(j, B))
  vb(j, 3, 2) = -i * gkl(j, B) * (0.125d0 * (gamma - 4.0d0 * delta) + &
              & (2.0d0 * (Omega ** 2) / (gamma - 4.0d0 * delta)) * &
              & (1.0d0 - (gamma / (0.75d0 * gamma - delta - (kappa(B) + i * wl(j, B)))))) / &
              & (2.0d0 * delta + kappa(B) + i * wl(j, B))
  vb(j, 4, 2) = -gamma * Omega / (0.75d0 * gamma + delta)
  ! c variable
  vc(j, 0, 2) = (Omega * gkl(j, B) * smss(2) - &
              & 0.125d0 * i * gkl(j, B) * (gamma + 4.0d0 * delta) * (szss(2) + 1.0d0) - &
              & gamma * Omega * ajss(j, 2)) / &
              & (0.25d0 * gamma - delta - alpha)
  vc(j, 1, 2) = -Omega * gkl(j, B) * (1.0d0 - (gamma / (0.5d0 * gamma - (kappa(B) + i * wl(j, B))))) / &
              & (0.25d0 * gamma + delta + alpha)
  vc(j, 2, 2) = i * gkl(j, B) * (0.125d0 * (gamma + 4.0d0 * delta) + &
              & (2.0d0 * (Omega ** 2) / (gamma + 4.0d0 * delta)) * &
              & (1.0d0 - (gamma / (0.75d0 * gamma + delta - (kappa(B) + i * wl(j, B)))))) / &
              & (2.0d0 * delta - (kappa(B) + i * wl(j, B)))
  vc(j, 3, 2) = -i * gkl(j, B) * (0.125d0 * (gamma + 4.0d0 * delta) + &
              & (2.0d0 * (Omega ** 2) / (gamma - 4.0d0 * delta)) * &
              & (1.0d0 - (gamma / (0.75d0 * gamma - delta - (kappa(B) + i * wl(j, B)))))) / &
              & (kappa(B) + i * wl(j, B))
  vc(j, 4, 2) = -gamma * Omega / (0.75d0 * gamma - delta)
  ! Initial values for v1(t), v2(t), and v3(t)
  v1(2) = 0.0d0
  v1(2) = (-i * (vb(j, 1, 2) - vc(j, 1, 2)) * C1(2) + &
        & (gkl(j, B) * delta * va(j, 2, 2) - i * (vb(j, 2, 2) - vc(j, 2, 2))) * C2(2) + &
        & (gkl(j, B) * delta * va(j, 3, 2) - i * (vb(j, 3, 2) - vc(j, 3, 2))) * C3(2) + &
        & -i * (vb(j, 4, 2) - vc(j, 4, 2)) * C4(j, 2)) / &
        & (4.0d0 * delta)

  v2(2) = 0.0d0
  v2(2) = (i * (vb(j, 1, 2) - vc(j, 1, 2)) * C1(2) + &
        & (gkl(j, B) * delta * va(j, 2, 2) + i * (vb(j, 2, 2) - vc(j, 2, 2))) * C2(2) + &
        & (gkl(j, B) * delta * va(j, 3, 2) + i * (vb(j, 3, 2) - vc(j, 3, 2))) * C3(2) + &
        & i * (vb(j, 4, 2) - vc(j, 4, 2)) * C4(j, 2)) / &
        & (4.0d0 * delta)

  v3(2) = 0.0d0
  v3(2) = (vb(j, 1, 2) * (gamma + 4.0d0 * delta) - vc(j, 1, 2) * (gamma - 4.0d0 * delta)) * C1(2) + &
        & (vb(j, 2, 2) * (gamma + 4.0d0 * delta) - vc(j, 2, 2) * (gamma - 4.0d0 * delta)) * C2(2) + &
        & (vb(j, 3, 2) * (gamma + 4.0d0 * delta) - vc(j, 3, 2) * (gamma - 4.0d0 * delta)) * C3(2) + &
        & (vb(j, 4, 2) * (gamma + 4.0d0 * delta) - vc(j, 4, 2) * (gamma - 4.0d0 * delta)) * C4(j, 2)
  v3(2) = 0.125d0 * v3(2) / (Omega * delta)
  !--------------------------------------------------------!
  !     Cheaty Variables (a, b, c) - Cross-correlation     !
  !--------------------------------------------------------!
  ! Cheaty alpha
  alpha = -(0.5d0 * gamma + kappa(B) + i * wl(j, B))

  ! a variable
  va(j, 0, 3) = (szss(3) + 1.0d0) / (alpha)
  va(j, 2, 3) = 1.0d0 / (0.75d0 * gamma + delta + alpha)
  va(j, 3, 3) = 1.0d0 / (0.75d0 * gamma - delta + alpha)
  ! b variable
  vb(j, 0, 3) = (Omega * gkl(j, B) * smss(3) - &
              & 0.125d0 * i * gkl(j, B) * (gamma - 4.0d0 * delta) * (szss(3) + 1.0d0) - &
              & gamma * Omega * ajss(j, 3)) / &
              & (0.25d0 * gamma + delta - alpha)
  vb(j, 1, 3) = -Omega * gkl(j, B) * (1.0d0 - (gamma / (0.5d0 * gamma - (kappa(B) + i * wl(j, B))))) / &
              & (0.25d0 * gamma - delta + alpha)
  vb(j, 2, 3) = -i * gkl(j, B) * (0.125d0 * (gamma - 4.0d0 * delta) + &
              & (2.0d0 * (Omega ** 2) / (gamma + 4.0d0 * delta)) * &
              & (1.0d0 - (gamma / (0.75d0 * gamma + delta - (kappa(B) + i * wl(j, B)))))) / &
              & (kappa(B) + i * wl(j, B))
  vb(j, 3, 3) = -i * gkl(j, B) * (0.125d0 * (gamma - 4.0d0 * delta) + &
              & (2.0d0 * (Omega ** 2) / (gamma - 4.0d0 * delta)) * &
              & (1.0d0 - (gamma / (0.75d0 * gamma - delta - (kappa(B) + i * wl(j, B)))))) / &
              & (2.0d0 * delta + kappa(B) + i * wl(j, B))
  vb(j, 4, 3) = -gamma * Omega / (0.75d0 * gamma + delta)
  ! c variable
  vc(j, 0, 3) = (Omega * gkl(j, B) * smss(3) - &
              & 0.125d0 * i * gkl(j, B) * (gamma + 4.0d0 * delta) * (szss(3) + 1.0d0) - &
              & gamma * Omega * ajss(j, 3)) / &
              & (0.25d0 * gamma - delta - alpha)
  vc(j, 1, 3) = -Omega * gkl(j, B) * (1.0d0 - (gamma / (0.5d0 * gamma - (kappa(B) + i * wl(j, B))))) / &
              & (0.25d0 * gamma + delta + alpha)
  vc(j, 2, 3) = i * gkl(j, B) * (0.125d0 * (gamma + 4.0d0 * delta) + &
              & (2.0d0 * (Omega ** 2) / (gamma + 4.0d0 * delta)) * &
              & (1.0d0 - (gamma / (0.75d0 * gamma + delta - (kappa(B) + i * wl(j, B)))))) / &
              & (2.0d0 * delta - (kappa(B) + i * wl(j, B)))
  vc(j, 3, 3) = -i * gkl(j, B) * (0.125d0 * (gamma + 4.0d0 * delta) + &
              & (2.0d0 * (Omega ** 2) / (gamma - 4.0d0 * delta)) * &
              & (1.0d0 - (gamma / (0.75d0 * gamma - delta - (kappa(B) + i * wl(j, B)))))) / &
              & (kappa(B) + i * wl(j, B))
  vc(j, 4, 3) = -gamma * Omega / (0.75d0 * gamma - delta)
  ! Initial values for v1(t), v2(t), and v3(t)
  v1(3) = 0.0d0
  v1(3) = (-i * (vb(j, 1, 3) - vc(j, 1, 3)) * C1(3) + &
        & (gkl(j, B) * delta * va(j, 2, 3) - i * (vb(j, 2, 3) - vc(j, 2, 3))) * C2(3) + &
        & (gkl(j, B) * delta * va(j, 3, 3) - i * (vb(j, 3, 3) - vc(j, 3, 3))) * C3(3) + &
        & -i * (vb(j, 4, 3) - vc(j, 4, 3)) * C4(j, 3)) / &
        & (4.0d0 * delta)

  v2(3) = 0.0d0
  v2(3) = (i * (vb(j, 1, 3) - vc(j, 1, 3)) * C1(3) + &
        & (gkl(j, B) * delta * va(j, 2, 3) + i * (vb(j, 2, 3) - vc(j, 2, 3))) * C2(3) + &
        & (gkl(j, B) * delta * va(j, 3, 3) + i * (vb(j, 3, 3) - vc(j, 3, 3))) * C3(3) + &
        & i * (vb(j, 4, 3) - vc(j, 4, 3)) * C4(j, 3)) / &
        & (4.0d0 * delta)

  v3(3) = 0.0d0
  v3(3) = (vb(j, 1, 3) * (gamma + 4.0d0 * delta) - vc(j, 1, 3) * (gamma - 4.0d0 * delta)) * C1(3) + &
        & (vb(j, 2, 3) * (gamma + 4.0d0 * delta) - vc(j, 2, 3) * (gamma - 4.0d0 * delta)) * C2(3) + &
        & (vb(j, 3, 3) * (gamma + 4.0d0 * delta) - vc(j, 3, 3) * (gamma - 4.0d0 * delta)) * C3(3) + &
        & (vb(j, 4, 3) * (gamma + 4.0d0 * delta) - vc(j, 4, 3) * (gamma - 4.0d0 * delta)) * C4(j, 3)
  v3(3) = 0.125d0 * v3(3) / (Omega * delta)

  !-----------------------------------------------------!
  !     Coefficients - Auto-correlation of Cavity A     !
  !-----------------------------------------------------!
  ! C6
  C6(j, 1) = 0.5d0 * (ajsm0(j, 1) - ajsmss(j, 1) - v1(1)) + 0.5d0 * (ajsp0(j, 1) - ajspss(j, 1) - v2(1))
  ! C7
  C7(j, 1) = (-i * (gamma + 4.0d0 * delta) * (gamma - 4.0d0 * delta) / (32.0d0 * delta * Omega)) * &
           & ((ajsm0(j, 1) - ajsmss(j, 1) - v1(1)) - (ajsp0(j, 1) - ajspss(j, 1) - v2(1))) + &
           & (gamma + 4.0d0 * delta) * (ajsz0(j, 1) - ajszss(j, 1) - v3(1)) / (8.0d0 * delta)
  ! C8
  C8(j, 1) = (i * (gamma + 4.0d0 * delta) * (gamma - 4.0d0 * delta) / (32.0d0 * delta * Omega)) * &
           & ((ajsm0(j, 1) - ajsmss(j, 1) - v1(1)) - (ajsp0(j, 1) - ajspss(j, 1) - v2(1))) - &
           & (gamma - 4.0d0 * delta) * (ajsz0(j, 1) - ajszss(j, 1) - v3(1)) / (8.0d0 * delta)

  !-----------------------------------------------------!
  !     Coefficients - Auto-correlation of Cavity B     !
  !-----------------------------------------------------!
  ! C6
  C6(j, 2) = 0.5d0 * (ajsm0(j, 2) - ajsmss(j, 2) - v1(2)) + 0.5d0 * (ajsp0(j, 2) - ajspss(j, 2) - v2(2))
  ! C7
  C7(j, 2) = (-i * (gamma + 4.0d0 * delta) * (gamma - 4.0d0 * delta) / (32.0d0 * delta * Omega)) * &
           & ((ajsm0(j, 2) - ajsmss(j, 2) - v1(2)) - (ajsp0(j, 2) - ajspss(j, 2) - v2(2))) + &
           & (gamma + 4.0d0 * delta) * (ajsz0(j, 2) - ajszss(j, 2) - v3(2)) / (8.0d0 * delta)
  ! C8
  C8(j, 2) = (i * (gamma + 4.0d0 * delta) * (gamma - 4.0d0 * delta) / (32.0d0 * delta * Omega)) * &
           & ((ajsm0(j, 2) - ajsmss(j, 2) - v1(2)) - (ajsp0(j, 2) - ajspss(j, 2) - v2(2))) - &
           & (gamma - 4.0d0 * delta) * (ajsz0(j, 2) - ajszss(j, 2) - v3(2)) / (8.0d0 * delta)

  !------------------------------------------!
  !     Coefficients - Cross-correlation     !
  !------------------------------------------!
  ! C6
  C6(j, 3) = 0.5d0 * (ajsm0(j, 3) - ajsmss(j, 3) - v1(3)) + 0.5d0 * (ajsp0(j, 3) - ajspss(j, 3) - v2(3))
  ! C7
  C7(j, 3) = (-i * (gamma + 4.0d0 * delta) * (gamma - 4.0d0 * delta) / (32.0d0 * delta * Omega)) * &
           & ((ajsm0(j, 3) - ajsmss(j, 3) - v1(3)) - (ajsp0(j, 3) - ajspss(j, 3) - v2(3))) + &
           & (gamma + 4.0d0 * delta) * (ajsz0(j, 3) - ajszss(j, 3) - v3(3)) / (8.0d0 * delta)
  ! C8
  C8(j, 3) = (i * (gamma + 4.0d0 * delta) * (gamma - 4.0d0 * delta) / (32.0d0 * delta * Omega)) * &
           & ((ajsm0(j, 3) - ajsmss(j, 3) - v1(3)) - (ajsp0(j, 3) - ajspss(j, 3) - v2(3))) - &
           & (gamma - 4.0d0 * delta) * (ajsz0(j, 3) - ajszss(j, 3) - v3(3)) / (8.0d0 * delta)
END DO

! Cheaty function
ALLOCATE(f_jk(-N:N, -N:N, 1:8, 3)); f_jk = 0.0d0
ALLOCATE(f_jk0(-N:N, -N:N, 3)); f_jk0 = 0.0d0
DO j = -N, N
  DO k = -N, N
    ! Auto-correlation of Cavity A
    f_jk(j, k, 1, 1) = -0.25d0 * C1(1) * i * (vb(k, 1, 1) - vc(k, 1, 1)) / &
                     & (delta * (0.5d0 * gamma - 2.0d0 * kappa(A) + i * (wl(j, A) - wl(k, A))))
    f_jk(j, k, 2, 1) = -0.25d0 * C2(1) * (gkl(k, A) * delta * va(k, 2, 1) + i * (vb(k, 2, 1) - vc(k, 2, 1))) / &
                     & (delta * (0.75d0 * gamma + delta - 2.0d0 * kappa(A) + i * (wl(j, A) - wl(k, A))))
    f_jk(j, k, 3, 1) = -0.25d0 * C3(1) * (gkl(k, A) * delta * va(k, 3, 1) + i * (vb(k, 3, 1) - vc(k, 3, 1))) / &
                     & (delta * (0.75d0 * gamma - delta - 2.0d0 * kappa(A) + i * (wl(j, A) - wl(k, A))))
    f_jk(j, k, 4, 1) = 0.25d0 * C4(k, 1) * i * (vb(k, 4, 1) - vc(k, 4, 1)) / (delta * (kappa(A) - i * wl(j, A)))
    f_jk(j, k, 6, 1) = -C6(k, 1) / (0.5d0 * gamma - (kappa(A) - i * wl(j, A)))
    f_jk(j, k, 7, 1) = -C7(k, 1) * 2.0d0 * i * Omega / &
                     & ((gamma + 4.0d0 * delta) * (0.75d0 * gamma + delta - (kappa(A) - i * wl(j, A))))
    f_jk(j, k, 8, 1) = -C8(k, 1) * 2.0d0 * i * Omega / &
                     & ((gamma - 4.0d0 * delta) * (0.75d0 * gamma - delta - (kappa(A) - i * wl(j, A))))
    ! Initial value
    f_jk0(j, k, 1) = f_jk(j, k, 1, 1) + f_jk(j, k, 2, 1) + f_jk(j, k, 3, 1) + f_jk(j, k, 4, 1) + &
                   & f_jk(j, k, 6, 1) + f_jk(j, k, 7, 1) + f_jk(j, k, 8, 1)

    ! Auto-correlation of Cavity B
    f_jk(j, k, 1, 2) = -0.25d0 * C1(2) * i * (vb(k, 1, 2) - vc(k, 1, 2)) / &
                     & (delta * (0.5d0 * gamma - 2.0d0 * kappa(B) + i * (wl(j, B) - wl(k, B))))
    f_jk(j, k, 2, 2) = -0.25d0 * C2(2) * (gkl(k, B) * delta * va(k, 2, 2) + i * (vb(k, 2, 2) - vc(k, 2, 2))) / &
                     & (delta * (0.75d0 * gamma + delta - 2.0d0 * kappa(B) + i * (wl(j, B) - wl(k, B))))
    f_jk(j, k, 3, 2) = -0.25d0 * C3(2) * (gkl(k, B) * delta * va(k, 3, 2) + i * (vb(k, 3, 2) - vc(k, 3, 2))) / &
                     & (delta * (0.75d0 * gamma - delta - 2.0d0 * kappa(B) + i * (wl(j, B) - wl(k, B))))
    f_jk(j, k, 4, 2) = 0.25d0 * C4(k, 2) * i * (vb(k, 4, 2) - vc(k, 4, 2)) / (delta * (kappa(B) - i * wl(j, B)))
    f_jk(j, k, 6, 2) = -C6(k, 2) / (0.5d0 * gamma - (kappa(B) - i * wl(j, B)))
    f_jk(j, k, 7, 2) = -C7(k, 2) * 2.0d0 * i * Omega / &
                     & ((gamma + 4.0d0 * delta) * (0.75d0 * gamma + delta - (kappa(B) - i * wl(j, B))))
    f_jk(j, k, 8, 2) = -C8(k, 2) * 2.0d0 * i * Omega / &
                     & ((gamma - 4.0d0 * delta) * (0.75d0 * gamma - delta - (kappa(B) - i * wl(j, B))))
    ! Initial value
    f_jk0(j, k, 2) = f_jk(j, k, 1, 2) + f_jk(j, k, 2, 2) + f_jk(j, k, 3, 2) + f_jk(j, k, 4, 2) + &
                   & f_jk(j, k, 6, 2) + f_jk(j, k, 7, 2) + f_jk(j, k, 8, 2)

    ! Cross-correlation
    f_jk(j, k, 1, 3) = -0.25d0 * C1(3) * i * (vb(k, 1, 3) - vc(k, 1, 3)) / &
                     & (delta * (0.5d0 * gamma - 2.0d0 * kappa(B) + i * (wl(j, B) - wl(k, B))))
    f_jk(j, k, 2, 3) = -0.25d0 * C2(3) * (gkl(k, B) * delta * va(k, 2, 3) + i * (vb(k, 2, 3) - vc(k, 2, 3))) / &
                     & (delta * (0.75d0 * gamma + delta - 2.0d0 * kappa(B) + i * (wl(j, B) - wl(k, B))))
    f_jk(j, k, 3, 3) = -0.25d0 * C3(3) * (gkl(k, B) * delta * va(k, 3, 3) + i * (vb(k, 3, 3) - vc(k, 3, 3))) / &
                     & (delta * (0.75d0 * gamma - delta - 2.0d0 * kappa(B) + i * (wl(j, B) - wl(k, B))))
    f_jk(j, k, 4, 3) = 0.25d0 * C4(k, 3) * i * (vb(k, 4, 3) - vc(k, 4, 3)) / (delta * (kappa(B) - i * wl(j, B)))
    f_jk(j, k, 6, 3) = -C6(k, 3) / (0.5d0 * gamma - (kappa(B) - i * wl(j, B)))
    f_jk(j, k, 7, 3) = -C7(k, 3) * 2.0d0 * i * Omega / &
                     & ((gamma + 4.0d0 * delta) * (0.75d0 * gamma + delta - (kappa(B) - i * wl(j, B))))
    f_jk(j, k, 8, 3) = -C8(k, 3) * 2.0d0 * i * Omega / &
                     & ((gamma - 4.0d0 * delta) * (0.75d0 * gamma - delta - (kappa(B) - i * wl(j, B))))
    ! Initial value
    f_jk0(j, k, 3) = f_jk(j, k, 1, 3) + f_jk(j, k, 2, 3) + f_jk(j, k, 3, 3) + f_jk(j, k, 4, 3) + &
                   & f_jk(j, k, 6, 3) + f_jk(j, k, 7, 3) + f_jk(j, k, 8, 3)
  END DO
END DO

! Calculate Differential Coefficients
ALLOCATE(C9(-N:N, -N:N, 3)); C9 = 0.0d0
DO j = -N, N
  DO k = -N, N
    ! Auto-correlation of Cavity A
    C9(j, k, 1) = atjak0(j, k, 1) - atjakss(j, k, 1) + &
                & (CONJG(gkl(j, A)) * f_jk0(j, k, 1)) + (gkl(k, A) * CONJG(f_jk0(k, j, 1)))
    ! Auto-correlation of Cavity B
    C9(j, k, 2) = atjak0(j, k, 2) - atjakss(j, k, 2) + &
                & (CONJG(gkl(j, B)) * f_jk0(j, k, 2)) + (gkl(k, B) * CONJG(f_jk0(k, j, 2)))
    ! Cross-correlation
    C9(j, k, 3) = atjak0(j, k, 3) - atjakss(j, k, 3) + &
                & (CONJG(gkl(j, B)) * f_jk0(j, k, 3)) + (gkl(k, B) * CONJG(f_jk0(k, j, 3)))
  END DO
END DO

! Open file to write time and correlation to
OPEN(UNIT=4, FILE=filename_g2, STATUS='REPLACE', ACTION='WRITE', RECL=4000)

! Number of time steps to integrate over
tau_steps = NINT(tau2_max / dt)
! Ten percent of time steps
ten_percent = NINT((1.0 * tau_steps / 10.0))

! Sample rate
IF (tau_steps > 100000) THEN
 sample_rate = NINT(DBLE(tau_steps) / 1d5)
ELSE
 sample_rate = 1
END IF
! Call CPU clock time
CALL CPU_TIME(loop_start_time)

! Start time integration to find steady state
DO t = 0, tau_steps
  ! Update time
  dtime = DBLE(t) * dt

  !-------------------!
  !     CALCULATE     !
  !-------------------!
  ! Auto-correlation of Cavity A
  g2(1) = photonA_ss * photonA_ss
  ! Auto-correlation of Cavity B
  g2(2) = photonB_ss * photonB_ss
  ! Cross-correlation of Cavity A with Cavity B
  g2(3) = photonA_ss * photonB_ss

  ! Cycle through auto-A, auto-B, and cross-AB
  DO m = 1, 3
    IF (m .EQ. 1) THEN
      ! Auto-correlation of Cavity A
      Cj = A
    ELSE
      ! Auto-correlation of Cavity B or cross-correlation
      Cj = B
    END IF
    ! Cycle through modes
    DO j = -N, N
      DO k = -N, N
        g2(m) = g2(m) + C9(j, k, m) * exp(-(2.0d0 * kappa(Cj) - i * (wl(j, Cj) - wl(k, Cj))) * dtime)
        g2(m) = g2(m) - CONJG(gkl(j, Cj)) * ( &
              & f_jk(j, k, 1, m) * EXP(-0.5d0 * gamma * dtime) + &
              & f_jk(j, k, 2, m) * EXP(-(0.75d0 * gamma + delta) * dtime) + &
              & f_jk(j, k, 3, m) * EXP(-(0.75d0 * gamma - delta) * dtime) + &
              & f_jk(j, k, 4, m) * EXP(-(kappa(Cj) + i * wl(k, Cj)) * dtime) + &
              & f_jk(j, k, 6, m) * EXP(-(0.5d0 * gamma + (kappa(Cj) + i * wl(k, Cj))) * dtime) + &
              & f_jk(j, k, 7, m) * EXP(-(0.75d0 * gamma + delta + (kappa(Cj) + i * wl(k, Cj))) * dtime) + &
              & f_jk(j, k, 8, m) * EXP(-(0.75d0 * gamma - delta + (kappa(Cj) + i * wl(k, Cj))) * dtime))
        g2(m) = g2(m) - gkl(k, Cj) * CONJG( &
              & f_jk(k, j, 1, m) * EXP(-0.5d0 * gamma * dtime) + &
              & f_jk(k, j, 2, m) * EXP(-(0.75d0 * gamma + delta) * dtime) + &
              & f_jk(k, j, 3, m) * EXP(-(0.75d0 * gamma - delta) * dtime) + &
              & f_jk(k, j, 4, m) * EXP(-(kappa(Cj) + i * wl(j, Cj)) * dtime) + &
              & f_jk(k, j, 6, m) * EXP(-(0.5d0 * gamma + (kappa(Cj) + i * wl(j, Cj))) * dtime) + &
              & f_jk(k, j, 7, m) * EXP(-(0.75d0 * gamma + delta + (kappa(Cj) + i * wl(j, Cj))) * dtime) + &
              & f_jk(k, j, 8, m) * EXP(-(0.75d0 * gamma - delta + (kappa(Cj) + i * wl(j, Cj))) * dtime))
      END DO
    END DO
  END DO

  !-------------------!
  !     WRITE DATA    !
  !-------------------!
  ! Calculate state probabilities and mean photon number
  IF (photonA_ss .NE. 0.0d0) THEN
    g2(1) = g2(1) / (photonA_ss * photonA_ss)
  END IF
  IF (photonB_ss .NE. 0.0d0) THEN
    g2(2) = g2(2) / (photonB_ss * photonB_ss)
  END IF
  IF (photonA_ss .NE. 0.0d0 .AND. photonB_ss .NE. 0.0d0) THEN
    g2(3) = g2(3) / (photonA_ss * photonB_ss)
  END IF

  ! If t_max is really big, only take a sample of results to write to file
  ! so file size isn't huge-mongous.
  IF (MOD(t, sample_rate) == 0) THEN
    ! Write data
    WRITE(4,*) DBLE(t) * dt, REAL(g2(1)), REAL(g2(2)), REAL(g2(3))
  END IF


  ! Check percentage
  IF (progress_bar .EQV. .TRUE.) THEN
    IF (MOD(t, ten_percent) == 0 .AND. t /= 0) THEN
      CALL CPU_TIME(loop_check_time)
      percentage = NINT((100.0 * t) / (1.0 * tau_steps))
      loop_run_time = loop_check_time - loop_start_time
      loop_remaining_time = ((100.0 * (loop_check_time - loop_start_time)) / (1.0 * percentage)) - loop_run_time
      WRITE(*, FMT_ss) percentage, "%. Run time (steady state): ", loop_run_time, "s. Est. time left: ", &
                     & loop_remaining_time, "s"
    END IF
  END IF

  ! Close time integration DO loop
END DO

! Close state file
CLOSE(4)

! Call end time from CPU_TIME
CALL CPU_TIME(end_time)
! PRINT*, "Runtime: ", end_time - start_time, "seconds"

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

END PROGRAM TWO_LEVEL_MULTI_MODE_CROSS_CORRELATIONS
