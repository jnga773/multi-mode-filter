#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 13:44:32 2021

@author: jnga773
"""

#-----------------------------------------------------------------------------#
#                        STATE-STRING TO INDEX CONVERTER                      #
#-----------------------------------------------------------------------------#

def _state_str_to_index(state_in):
    """
    Converts the input state string {'0', '+', '-'} into an index {0, 1, 2}.

    Parameters
    ----------
    state_in : string
        Input state to convert

    Returns
    -------
    index_out : integer
        Index for array
    """
    # Dictionary for string to index
    state_dic = {'0': 0, '+': 1, '-': 2}
    # Check if state_in /= '0', '+', or '-'
    if state_in not in ['0', '+', '-']:
        from sys import exit
        exit("pee pee poo poo: state_in != '0', '+', or '-'!")
    # Find index
    index_out = state_dic[state_in]

    return index_out

def _which_transition(w0_in, Omega_in, alpha_in, delta_in, xi_in):
    """
    Checks where w0 is an compares to dressed-state transition
    frequency to see which dressed-state transition has occured.

    Parameters
    ----------
    w0_in : float
        Central resonance frequency of filter.
    Omega_in : float
        Driving amplitude.
    alpha_in : float
        Anharmonicity of atom.
    delta_in : float
        Driven detuning from two-photon resonance.
    xi_in : float
        Dipole moment ratio.

    Returns
    -------
    transition_out = list, string
        Two element list where first is the initial state and
        the second element is the final state ['0', '+', '-']
    """
    # Calculate eigenvalues
    w0, wp, wm = three_level_eig(Omega_in, alpha_in, delta_in, xi_in)

    # For each transition, numerically order transition frequency
    w_trans = [-(wp - wm), -(wp - w0), -(w0 - wm), 0.0,
               w0 - wm, wp - w0, wp - wm]

    # Transition labels
    labels_trans = [['-', '+'], ['0', '+'], ['-', '0'], ['+', '+'],
                    ['0', '-'], ['+', '0'], ['+', '-']]

    # Cycle through w_trans and compare with w0_in to find matching transition
    w_check = False
    for index, w in enumerate(w_trans):
        if round(w0_in, 2) == round(w, 2):
            # Set check to True so no fail
            w_check = True
            # Grab the transition labels
            transition_out = labels_trans[index]

    # If w0_in does not match any transition frequency, exit the function
    # with an error.
    if w_check is False:
        from sys import exit
        exit("pee pee poo poo: w0_in does not match any eigen-transition!")
    else:
        return transition_out

#-----------------------------------------------------------------------------#
#                             PARAMETER CALCULATOR                            #
#-----------------------------------------------------------------------------#

def three_level_eig(Omega_in, alpha_in, delta_in, xi_in, vals_or_vecs='vals'):
    """
    Calculates the dressed/eigenvalues of the Hamiltonian in matrix form.

    Parameters
    ----------
    Omega_in : float
        Driving amplitude.
    alpha_in : float
        Anharmonicity of atom.
    delta_in : float
        Driven detuning from two-photon resonance.
    xi_in : float
        Dipole moment ratio.
    vals_or_vecs : string, optional
        Choose output to be eigenvalues ('vals'), eigenvectors ('vecs') or
        both ('both'). The default is 'vals'.

    Returns
    -------
    eigvals_out : array
        Array of eigenvalues in order [w0, wp, wm].
    eigvecs_out : matrix
        S matrix from diagonalisation of H; containing eigenvectors in order
        [v0, vp, vm].
    """
    from numpy.linalg import eig
    from numpy import matrix
    # Set up matrix
    H = matrix([[0, 0.5 * Omega_in, 0],
                [0.5 * Omega_in, -(0.5 * alpha_in) - delta_in, 0.5 * xi_in * Omega_in],
                [0, 0.5 * xi_in * Omega_in, -2 * delta_in]])
    # Calculate eigenvalues
    # eigvals_out = eigvals(H)
    eigvals_out, eigvecs_out = eig(H)

    # Get the indicies that would sort them from big to small
    sorted_indices = eigvals_out.argsort()[::-1]
    # Return the sorted eigenvalues and eigenvectors
    eigvals_out = eigvals_out[sorted_indices]
    eigvecs_out = eigvecs_out[:, sorted_indices]

    # This sorting will be in the order [\omega_{+}, \omega_{0}, \omega_{-}].
    # To match it with the analytic calculations I've done, let's change it to
    # the the order of [\omega_{0}, \omega_{+}, \omega_{-}]
    sorted_indices = [1, 0, 2]
    eigvals_out = eigvals_out[sorted_indices]
    eigvecs_out = eigvecs_out[:, sorted_indices]

    # Return the output depending on vals_or_vecs
    if vals_or_vecs == 'vals':
        return eigvals_out
    if vals_or_vecs == 'vecs':
        return eigvecs_out
    if vals_or_vecs == 'both':
        return eigvals_out, eigvecs_out
    
def print_transition_frequencies(Omega_in, alpha_in, delta_in, xi_in,
                                 output=False):
    """
    Prints or outputs the frequencies for the dressed-state transitions.
    
    Parameters
    ----------
    Omega_in : float
        Driving amplitude.
    alpha_in : float
        Anharmonicity of atom.
    delta_in : float
        Driven detuning from two-photon resonance.
    xi_in : float
        Dipole moment ratio.
    vals_or_vecs : string, optional
        Choose output to be eigenvalues ('vals'), eigenvectors ('vecs') or
        both ('both'). The default is 'vals'.
    output : boolean (default = False)
        If True, outputs the frequencies as a list
        
    Returns
    -------
    transition_frequencies : array like
        List of transition frequencies in numerical order
    """
    # Calculate eigen-frequencies
    w0, wp, wm = three_level_eig(Omega_in, alpha_in, delta_in, xi_in,
                                 vals_or_vecs='vals')
    # Print as a pretty thingy
    print("|\omega_+ - \omega_-| = {}".format(abs(wp - wm)))
    print("|\omega_0 - \omega_+| = {}".format(abs(w0 - wp)))
    print("|\omega_0 - \omega_-| = {}".format(abs(w0 - wm)))
    
    if output:
        from numpy import array
        transition_frequencies = [wm - wp,
                                  w0 - wp,
                                  wm - w0,
                                  0.0,
                                  w0 - wm,
                                  wp - w0,
                                  wp - wm]
        transition_frequencies = array(transition_frequencies)
        
        return transition_frequencies

def _Sigma_matrix_elements(Gamma_in, Omega_in, alpha_in, delta_in, xi_in):
    """
    Calculates the matrix elements (a1 - a9) of the \Sigma_{-} operator in the
    dressed state basis.

    Parameters
    ----------
    Gamma_in : float
        Atomic decay rate.
    Omega_in : float
        Driving amplitude.
    alpha_in : float
        Anharmonicity of atom.
    delta_in : float
        Driven detuning from two-photon resonance.
    xi_in : float
        Dipole moment ratio.
    return_matrix : boolean, optional
        If True, returns the 3x3 matrix for calculating the diagonal dressed
        state moments.

    Returns
    -------
    a_out : array_like
        Array of the matrix elements
    """
    from numpy.linalg import inv
    from numpy import array

    # Calculate eigenvectors
    eigvecs = three_level_eig(Omega_in, alpha_in, delta_in, xi_in, "vecs")

    # Invert eigvec matrix
    Sinv = inv(eigvecs)

    # Grab each element from this matrix
    g0 = Sinv[0, 0]
    gp = Sinv[1, 0]
    gm = Sinv[2, 0]
    e0 = Sinv[0, 1]
    ep = Sinv[1, 1]
    em = Sinv[2, 1]
    f0 = Sinv[0, 2]
    fp = Sinv[1, 2]
    fm = Sinv[2, 2]

    # Nine matrix elements
    a1 = (g0 * e0) + xi_in * (e0 * f0)
    a2 = (g0 * ep) + xi_in * (e0 * fp)
    a3 = (g0 * em) + xi_in * (e0 * fm)
    a4 = (gp * e0) + xi_in * (ep * f0)
    a5 = (gp * ep) + xi_in * (ep * fp)
    a6 = (gp * em) + xi_in * (ep * fm)
    a7 = (gm * e0) + xi_in * (em * f0)
    a8 = (gm * ep) + xi_in * (em * fp)
    a9 = (gm * em) + xi_in * (em * fm)

    # Create array
    a_out = array([a1, a2, a3, a4, a5, a6, a7, a8, a9], dtype='complex')

    # Return
    return a_out

def Gamma_values(Gamma_in, Omega_in, alpha_in, delta_in, xi_in,
                 return_matrix=False):
    """
    Calculates the dressed state Gamma values for each transition

    Parameters
    ----------
    Gamma_in : float
        Atomic decay rate.
    Omega_in : float
        Driving amplitude.
    alpha_in : float
        Anharmonicity of atom.
    delta_in : float
        Driven detuning from two-photon resonance.
    xi_in : float
        Dipole moment ratio.
    return_matrix : boolean, optional
        If True, returns the 3x3 matrix for calculating the diagonal dressed
        state moments.

    Returns
    -------
    Gpz_out, Gmz_out, Gpm_out : float
        Gamma values for the off-diagonal moments.
    M_diag_out : matrix
        3x3 matrix for coupled diagonal moment equations
    """

    # Calculate matrix elements
    a_elements = _Sigma_matrix_elements(Gamma_in, Omega_in, alpha_in, delta_in, xi_in)
    a1, a2, a3, a4, a5, a6, a7, a8, a9 = a_elements

    # Different diagonal Gamma rates
    Gpz_out = Gamma_in * ((a1 ** 2) + (a2 ** 2) + (a4 ** 2) + (a5 ** 2) +
                      (a7 ** 2) + (a8 ** 2) - (a1 * a5))
    Gmz_out = Gamma_in * ((a1 ** 2) + (a3 ** 2) + (a4 ** 2) + (a6 ** 2) +
                      (a7 ** 2) + (a9 ** 2) - (a1 * a9))
    Gpm_out = Gamma_in * ((a2 ** 2) + (a3 ** 2) + (a5 ** 2) + (a6 ** 2) +
                      (a8 ** 2) + (a9 ** 2) - (a5 * a9))

    return Gpz_out, Gmz_out, Gpm_out

def _get_Gammaij(Gamma_in, Omega_in, alpha_in, delta_in, xi_in,
                 state_init, state_final):
    """
    Calculate the correct Gamma_ij value from the initial and final states.

    Parameters
    ----------
    Gamma_in : float
        Atomic decay rate.
    Omega_in : float
        Driving amplitude.
    alpha_in : float
        Anharmonicity of atom.
    delta_in : float
        Driven detuning from two-photon resonance.
    xi_in : float
        Dipole moment ratio.
    state_init : string [0', '+', '-']
        Initial state.
    state_final : string [0', '+', '-']
        Final state.
    
    Returns
    -------
    Gamma_out : float
        Corret Gamma rate for off-diagonal moment evolution.
    """
    from numpy import matrix

    # Get indices
    i_init = _state_str_to_index(state_init)
    i_final = _state_str_to_index(state_final)

    # Calculate dressed state Gamma rates (no matrix)
    Gpz, Gmz, Gpm = Gamma_values(Gamma_in, Omega_in, alpha_in, delta_in, xi_in)

    # Sort into matrix
    mat = matrix([[0.0, Gpz, Gmz],
                  [Gpz, 0.0, Gpm],
                  [Gmz, Gpm, 0.0]])
    
    # Grab gamma_value
    Gamma_out = mat[i_init, i_final]

    # Print warning if Gamma is 0.0
    if Gamma_out == 0.0:
        print("Gamma_ij = 0.0!")
    
    return Gamma_out

def _diagonal_matrix(Gamma_in, Omega_in, alpha_in, delta_in, xi_in,
                     matrix_dim=2):
    """
    Calculates the evolution matrix for the diagonal moments,
        < |0><0| >, < |+><+| >, and < |-><-| >,
    with either the 3x3 matrix or the 2x2 matrix and non-homogeneous vector.

    Parameters
    ----------
    Gamma_in : float
        Atomic decay rate.
    Omega_in : float
        Driving amplitude.
    alpha_in : float
        Anharmonicity of atom.
    delta_in : float
        Driven detuning from two-photon resonance.
    xi_in : float
        Dipole moment ratio.
    matrix_dim : integer (default = 2)
        3x3 matrix or 2x2 matrix and B Vector

    Returns
    -------
    matrix_out : matrix
        Evolution matrix of the operators
    B_out : matrix
        If 3_or_2 == 2, outputs the non-homogeneous B vector too.
    """
    from numpy import zeros
    
    # Sigma matrix elements
    a_elements = _Sigma_matrix_elements(
        Gamma_in, Omega_in, alpha_in, delta_in, xi_in)
    a1, a2, a3, a4, a5, a6, a7, a8, a9 = a_elements
    
    if matrix_dim == 3:
        # Matrix for evolution
        M = zeros(shape=(3, 3), dtype='complex')
        # Matrix elements
        # d/dt |0><0|
        M[0, 0] = -Gamma_in * ((a4 ** 2) + (a7 ** 2))
        M[0, 1] = Gamma_in * (a2 ** 2)
        M[0, 2] = Gamma_in * (a3 ** 2)
        # d/dt |+><+|
        M[1, 0] = Gamma_in * (a4 ** 2)
        M[1, 1] = -Gamma_in * ((a2 ** 2) + (a8 ** 2))
        M[1, 2] = Gamma_in * (a6 ** 2)
        # d/dt |-><-|
        M[2, 0] = Gamma_in * (a7 ** 2)
        M[2, 1] = Gamma_in * (a8 ** 2)
        M[2, 2] = -Gamma_in * ((a3 ** 2) + (a6 ** 2))
        
        # Return a zero vector
        B = zeros(shape=(3, 3), dtype='complex')
        
    elif matrix_dim == 2:
        # Matrix for evolution
        M = zeros(shape=(2, 2), dtype='complex')
        # Matrix elements
        # d/dt |0><0|
        M[0, 0] = -Gamma_in * ((a4 ** 2) + (a7 ** 2) + (a3 ** 2))
        M[0, 1] = Gamma_in * ((a2 ** 2) - (a3 ** 2))
        # d/dt |+><+|
        M[1, 0] = Gamma_in * ((a4 ** 2) - (a6 ** 2))
        M[1, 1] = -Gamma_in * ((a2 ** 2) + (a6 ** 2) + (a8 ** 2))
    
        # B vector thing
        B = zeros(shape=(2, 1), dtype='complex')
        B[0, 0] = Gamma_in * (a3 ** 2)
        B[1, 0] = Gamma_in * (a6 ** 2)
        
    # Return
    return M, B
    
#-----------------------------------------------------------------------------#
#                        DRESSED-STATE MOMENT EQUATIONS                       #
#-----------------------------------------------------------------------------#

def diagonal_ss(Gamma_in, Omega_in, alpha_in, delta_in, xi_in):
    """
    Calculates the steady states of the diagonal moments
    """
    from numpy import matmul
    from numpy.linalg import inv
    
    # Matrix for evolution
    M, B = _diagonal_matrix(Gamma_in, Omega_in, alpha_in, delta_in, xi_in,
                            matrix_dim=2)

    # Invert the matrix
    M_inv = inv(M)
    
    # Multiply by the non-homogeneous vector
    ss_out = -1.0 * matmul(M_inv, B)
    
    # Grab individual solutions
    zz_out = ss_out[0, 0].real
    pp_out = ss_out[1, 0].real
    mm_out = 1 - (zz_out + pp_out)
    
    return zz_out, pp_out, mm_out

def off_diagonal_moments(t_in, Gamma_in, Omega_in, alpha_in, delta_in, xi_in,
                         state_init, state_final, initial_in=1.0):
    """
    Calculates the time evolution of the off-diagonal moments:
        \sigma_{ij} = |i><j|, i /= j.
    
    Parameters
    ----------
    t_in : float, array
        Array of t times for evolution.
    Gamma_in : float
        Atomic decay rate.
    Omega_in : float
        Driving amplitude.
    alpha_in : float
        Anharmonicity of atom.
    delta_in : float
        Driven detuning from two-photon resonance.
    xi_in : float
        Dipole moment ratio.
    state_init : string [0', '+', '-']
        Initial state.
    state_final : string [0', '+', '-']
        Final state.
    initial_in : float (default = 1.0)
        Initial condition of moment.
    
    Returns
    -------
    moment_out : complex array
        Time evolution of operator < |i><j|(t) >
    """
    from numpy import exp

    # Calculate the dressed-state frequencies
    w0, wp, wm = three_level_eig(Omega_in, alpha_in, delta_in, xi_in)
    w_dressed = [w0, wp, wm]

    # Grab indices
    i_init = _state_str_to_index(state_init)
    i_final = _state_str_to_index(state_final)

    # Grab necessary dressed-state frequencies
    w_init = w_dressed[i_init]
    w_final = w_dressed[i_final]

    # Calculate dressed state Gamma rates
    Gammaij = _get_Gammaij(Gamma_in, Omega_in, alpha_in, delta_in, xi_in,
                           state_init, state_final)

    # Set initial condition
    moment_out = initial_in
    # Calculate moment evolution
    moment_out *= exp((-0.5 * Gammaij + 1j * (w_init - w_final)) * t_in)

    return moment_out

def diagonal_moments(t_in, Gamma_in, Omega_in, alpha_in, delta_in, xi_in,
                     initial_state, output_state='all', matrix_dim=2):
    """
    Calculates the time evolution of the off-diagonal moments:
        \sigma_{ij} = |i > < i | .

    Parameters
    ----------
    t_in: float, array
        Array of t times for evolution.
    Gamma_in: float
        Atomic decay rate.
    Omega_in: float
        Driving amplitude.
    alpha_in: float
        Anharmonicity of atom.
    delta_in: float
        Driven detuning from two-photon resonance.
    xi_in: float
        Dipole moment ratio.
    initial_state : string
        Initial state ['0', '+', '-']
    output_state : string
        Which operator we want to output ['all', '0', '+', '-']
    matrix_dim : integer (default = 2)
        3x3 matrix or 2x2 matrix and B Vector
    
    Returns
    -------
    moment_out : complex, array
        Time evolution of operator < |i><i| (t) >
    """
    from numpy import zeros, matmul

    # Time step
    dt = t_in[1] - t_in[0]
    
    # Calculate Matrix and non-homogeneous vector
    M, B = _diagonal_matrix(Gamma_in, Omega_in, alpha_in, delta_in, xi_in,
                            matrix_dim)

    # State vector for [s0, sp, sm]
    X = zeros(shape=(matrix_dim, 1), dtype='complex')
    
    if initial_state != '-':
        # Get index of initial state and output state
        i_init = _state_str_to_index(initial_state)
        # Set initial condition
        X[i_init] = 1.0

    # Runge-Kutta vectors
    k1 = zeros(shape=(matrix_dim, 1), dtype='complex')
    k2 = zeros(shape=(matrix_dim, 1), dtype='complex')
    k3 = zeros(shape=(matrix_dim, 1), dtype='complex')
    k4 = zeros(shape=(matrix_dim, 1), dtype='complex')

    # data arrays
    s0 = zeros(shape=len(t_in), dtype='complex')
    sp = zeros(shape=len(t_in), dtype='complex')
    sm = zeros(shape=len(t_in), dtype='complex')

    # Calculate X with RK4
    for step in range(len(t_in)):
        # for step in range(2):
        # Update data
        s0[step] = X[0, 0]
        sp[step] = X[1, 0]
        if matrix_dim == 3:
            sm[step] = X[2]
        elif matrix_dim == 2:
            sm[step] = 1 - (X[0, 0] + X[1, 0])

        # Calculate Runge-Kutta Vectors
        k1 = dt * (matmul(M, X) + B)
        k2 = dt * (matmul(M, X + 0.5 * k1) + B)
        k3 = dt * (matmul(M, X + 0.5 * k2) + B)
        k4 = dt * (matmul(M, X + k3) + B)

        # Update X vector
        X += (1/6) * (k1 + 2 * (k2 + k3) + k4)

    # Create list of all moments
    all_moments = [s0, sp, sm]
    
    if output_state == "all":
        return all_moments
    else:
        i_out = _state_str_to_index(output_state)
        return all_moments[i_out]
    
#-----------------------------------------------------------------------------#
#                           FIRST-ORDER CORRELATION                           #
#-----------------------------------------------------------------------------#

def dressed_g1_calc(tau_in, w0_in, Gamma_in, Omega_in, alpha_in, delta_in,
                    xi_in):
    """
    Calculates the approximate dressed state first-order correlation function
    based on input parameters.

    Parameters
    ----------
    tau_in : float, array
        Array of tau times for correlation function.
    w0_in : float
        Central frequency of filter.
    Gamma_in : float
        Atomic decay rate.
    Omega_in : float
        Driving amplitude.
    alpha_in : float
        Anharmonicity of atom.
    delta_in : float
        Driven detuning from two-photon resonance.
    xi_in : float
        Dipole moment ratio.

    Returns
    -------
    corr_out : complex array
    """
    # From the central resonance frequency w0_in, check which transition
    # has occured.
    transition = _which_transition(w0_in, Omega_in, alpha_in, delta_in,
                                   xi_in)
    
    # Calculate steady states
    steady_state = diagonal_ss(Gamma_in, Omega_in, alpha_in, delta_in, xi_in)
    # Transition is from |i> -> |j> = |j><i|, so steady state value is
    # |j><j| operator
    ss_norm = steady_state[_state_str_to_index(transition[1])]
    
    if transition[0] != transition[1]:
        # Off-diagonal transition
        corr_out = ss_norm * off_diagonal_moments(tau_in, Gamma_in, Omega_in,
                                                  alpha_in, delta_in, xi_in,
                                                  transition[0], transition[1])
    else:
        # Centre peak
        corr_out = ss_norm * diagonal_moments(tau_in, Gamma_in, Omega_in,
                                              alpha_in, delta_in, xi_in,
                                              transition[0], transition[1])

    # Normalise by steady state value
    if ss_norm != 0.0:
        corr_out *= 1 / ss_norm
        
    return corr_out
    
#-----------------------------------------------------------------------------#
#                           SECOND-ORDER CORRELATION                          #
#-----------------------------------------------------------------------------#

def dressed_g2_auto_calc(tau_in, w0_in, Gamma_in, Omega_in, alpha_in, delta_in,
                         xi_in):
    """
    Calculates the exponential decay g1 for any transition that isn't the
    central peak.
    
    Parameters
    ----------
    tau_in : float, array
        Array of tau times for correlation function.
    w0_in : float
        Central frequency of filter.
    Gamma_in : float
        Atomic decay rate.
    Omega_in : float
        Driving amplitude.
    alpha_in : float
        Anharmonicity of atom.
    delta_in : float
        Driven detuning from two-photon resonance.
    xi_in : float
        Dipole moment ratio.
    
    Returns
    -------
    corr_out : complex, array
    """
    # From the central resonance frequency w0_in, check which transition
    # has occured.
    transition = _which_transition(w0_in, Omega_in, alpha_in, delta_in,
                                   xi_in)
    
    # Calculate steady states
    steady_state = diagonal_ss(Gamma_in, Omega_in, alpha_in, delta_in, xi_in)
    # Transition is from |i> -> |j> = |j><i|, so steady state value is
    # |j><j| operator
    ss_norm = steady_state[_state_str_to_index(transition[1])]
    
    # Calculate the thing
    corr_out = diagonal_moments(tau_in, Gamma_in, Omega_in, alpha_in,
                                delta_in, xi_in,
                                transition[0], transition[1])
    # Normalise
    corr_out = corr_out.real / (ss_norm ** 2)
    return corr_out

def dressed_g2_calc(tau_in, w0_in, Gamma_in, Omega_in, alpha_in, delta_in,
                    xi_in, w0b_in=None):
    """
    Calculates the exponential decay g1 for any transition that isn't the
    central peak.
    
    Parameters
    ----------
    tau_in : float, array
        Array of tau times for correlation function.
    w0_in : float
        Central frequency of filter.
    Gamma_in : float
        Atomic decay rate.
    Omega_in : float
        Driving amplitude.
    alpha_in : float
        Anharmonicity of atom.
    delta_in : float
        Driven detuning from two-photon resonance.
    xi_in : float
        Dipole moment ratio.
    w0b_in : float (Default None)
        Used if calculating cross-correlations
    
    Returns
    -------
    corr_out : complex, array
    """
    # If no second-filter given, set second frequency to be the same as first
    if w0b_in is None:
        w0b_in = w0_in
    
    # From the central resonance frequency w0_in, check which transition
    # has occured first.
    transition1 = _which_transition(w0_in, Omega_in, alpha_in, delta_in,
                                    xi_in)
    # Now check what the second transition was from w0b
    transition2 = _which_transition(w0b_in, Omega_in, alpha_in, delta_in,
                                    xi_in)
    
    # Calculate steady states
    steady_state = diagonal_ss(Gamma_in, Omega_in, alpha_in, delta_in, xi_in)
    # Transition is from |i> -> |j> = |j><i|, so steady state value is
    # |j><j| operator
    ss_norm1 = steady_state[_state_str_to_index(transition1[0])]
    ss_norm2 = steady_state[_state_str_to_index(transition2[0])]
    
    # Calculate the thing
    corr_out = diagonal_moments(tau_in, Gamma_in, Omega_in, alpha_in, delta_in,
                                xi_in, initial_state = transition1[1],
                                output_state = transition2[0])
    # Normalise
    corr_out = (ss_norm1 * corr_out.real) / (ss_norm1 * ss_norm2)
    
    return corr_out
