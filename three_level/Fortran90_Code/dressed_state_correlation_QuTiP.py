# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 19:33:44 2023

@author: Jacob
"""

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
                                 print_frequencies=True, output=False):
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
    if print_frequencies:
        print("|\omega_+ - \omega_-| = {}".format(abs(wp - wm)))
        print("|\omega_0 - \omega_+| = {}".format(abs(w0 - wp)))
        print("|\omega_0 - \omega_-| = {}".format(abs(w0 - wm)))
    
    if output:
        from numpy import array
        transition_frequencies_out = [wm - wp,
                                     w0 - wp,
                                     wm - w0,
                                     0.0,
                                     w0 - wm,
                                     wp - w0,
                                     wp - wm]
        transition_frequencies_out = array(transition_frequencies_out)
        
        return transition_frequencies_out
    
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
    # from numpy import matrix, matmul
    # Calculate eigenvalues
    w0, wp, wm = three_level_eig(Omega_in, alpha_in, delta_in, xi_in)

    # For each transition, numerically order transition frequency
    w_trans = [-(wp - wm), -(wp - w0), -(w0 - wm), 0.0,
               w0 - wm, wp - w0, wp - wm]

    # # Transition labels
    # labels_trans = [['d', 'u'], ['m', 'u'],
    #                 ['d', 'm'], ['+', '+'], ['m', 'd'],
    #                 ['u', 'm'], ['u', 'd']]

    # # Turn states into vectors
    # m = matrix([[1], [0], [0]])
    # u = matrix([[0], [1], [0]])
    # l = matrix([[0], [0], [1]])
    
    operator_trans = ['ul_p', 'um_p',
                      'ml_p', 'sz', 'ml_m',
                      'um_m', 'ul_m']

    # Cycle through w_trans and compare with w0_in to find matching transition
    w_check = False
    for index, w in enumerate(w_trans):
        if round(w0_in, 2) == round(w, 2):
            # Set check to True so no fail
            w_check = True
            # Grab the transition labels
            # transition_out = labels_trans[index]
            transition_out = operator_trans[index]

    # If w0_in does not match any transition frequency, exit the function
    # with an error.
    if w_check is False:
        from sys import exit
        exit("pee pee poo poo: w0_in does not match any eigen-transition!")
    else:
        return transition_out
    
def _Sigma_matrix_elements(Omega_in, alpha_in, delta_in, xi_in):
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
    a_out = array([a1, a2, a3,
                   a4, a5, a6,
                   a7, a8, a9], dtype='complex')

    # Return
    return a_out

def _diagonal_matrix(Gamma_in, Omega_in, alpha_in, delta_in, xi_in,
                     matrix_dim=2):
    """
    Calculates the evolution matrix for the diagonal moments,
        < |m><m| >, < |u><u| >, and < |l><l| >,
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
    a1, a2, a3, a4, a5, a6, a7, a8, a9 = \
        _Sigma_matrix_elements(Omega_in, alpha_in, delta_in, xi_in)
    
    # Matrix for evolution
    M_out = zeros(shape=(3, 3), dtype='complex')
    # Matrix elements
    # d/dt |m><m|
    M_out[0, 0] = -Gamma_in * ((a4 ** 2) + (a7 ** 2))
    M_out[0, 1] = Gamma_in * (a2 ** 2)
    M_out[0, 2] = Gamma_in * (a3 ** 2)
    # d/dt |u><u|
    M_out[1, 0] = Gamma_in * (a4 ** 2)
    M_out[1, 1] = -Gamma_in * ((a2 ** 2) + (a8 ** 2))
    M_out[1, 2] = Gamma_in * (a6 ** 2)
    # d/dt |l><l|
    M_out[2, 0] = Gamma_in * (a7 ** 2)
    M_out[2, 1] = Gamma_in * (a8 ** 2)
    M_out[2, 2] = -Gamma_in * ((a3 ** 2) + (a6 ** 2))
        
    # Return
    return M_out

def steady_state_diagonal(Gamma_in, Omega_in, alpha_in, delta_in, xi_in):
    """
    Calculates the steady states of the diagonal moments
    """
    from numpy import where, abs
    # from numpy.linalg import inv
    from numpy.linalg import eig
    # from numpy import matrix
    
    # Matrix for evolution
    M = _diagonal_matrix(Gamma_in, Omega_in, alpha_in, delta_in, xi_in)

    # Calculate eigenvalues of matrix
    eigvals, eigvecs = eig(M)
    
    # Find the index in the eigvals array that corresponds to the zero
    # eigenvalue
    zero_index = where(abs(eigvals) < 1e-8)[0][0]
    
    # Usnig this index value, grab the corresponding eigenvector
    eigvec_ss = eigvecs[:, zero_index]
    
    # Normalise by sum
    steady_state = eigvec_ss.real / (sum(eigvec_ss.real))
    
    # Output them babes
    m_out = steady_state[0]
    u_out = steady_state[1]
    l_out = steady_state[2]
    
    return m_out, u_out, l_out

#-----------------------------------------------------------------------------#
#                               QuTiP FUNCTION                                #
#-----------------------------------------------------------------------------#
def dressed_g2_calc(tau_in, w0a_in, Gamma_in, Omega_in, alpha_in, delta_in,
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
        w0b_in = w0a_in
        
    # Import qutip functions
    from qutip import basis, correlation_3op_1t, lindblad_dissipator
    from numpy import sqrt
    
    #-------------------------------------------------------------------------#
    #                              PARAMETERS                                 #
    #-------------------------------------------------------------------------#
    # Calculate eigenfrequencies
    w_m, w_u, w_l = three_level_eig(Omega_in, alpha_in, delta_in, xi_in)
    
    # Sigma matrix elements
    a1, a2, a3, a4, a5, a6, a7, a8, a9 = \
        _Sigma_matrix_elements(Omega_in, alpha_in, delta_in, xi_in)
    
    #-------------------------------------------------------------------------#
    #                              QuTiP THINGS                               #
    #-------------------------------------------------------------------------#
    #-------------------#
    #     Operators     #
    #-------------------#
    # Dressed state: |m>, |u>, |l>
    m, u, l = (basis(3, 0), basis(3, 1), basis(3, 2))

    # |u><m| operators
    um_m = m * u.dag(); um_p = u * m.dag()
    ml_m = l * m.dag(); ml_p = m * l.dag()
    ul_m = l * u.dag(); ul_p = u * l.dag()

    mm = m * m.dag(); uu = u * u.dag(); ll = l * l.dag()
    
    # Try sz
    # sz = uu - ll
    # sz = uu - mm
    # sz = mm - ll
    # sz = mm + ll + uu
    sz = uu + ll
    
    # Sort into dictionary for quick and easy access
    operator_dict = {'um_m': um_m, 'um_p': um_p,
                     'ml_m': ml_m, 'ml_p': ml_p,
                     'ul_m': ul_m, 'ul_p': ul_p,
                     'sz': sz}

    #---------------------#
    #     Hamiltonian     #
    #---------------------#
    # Hamiltonian
    H_A = (w_m * mm) + (w_u * uu) + (w_l * ll)
    
    # Collapse operators
    c_ops = [sqrt(Gamma_in * (a1 ** 2)) * lindblad_dissipator(mm),
             sqrt(Gamma_in * (a5 ** 2)) * lindblad_dissipator(uu),
             sqrt(Gamma_in * (a9 ** 2)) * lindblad_dissipator(ll),
             #-----------------------#
             sqrt(Gamma_in * a1 * a5) * (lindblad_dissipator(mm, uu) + lindblad_dissipator(uu, mm)),
             sqrt(Gamma_in * a1 * a9) * (lindblad_dissipator(mm, ll) + lindblad_dissipator(ll, mm)),
             sqrt(Gamma_in * a5 * a9) * (lindblad_dissipator(uu, ll) + lindblad_dissipator(ll, uu)),
             #-----------------------#
             sqrt(Gamma_in * (a2 ** 2)) * lindblad_dissipator(um_m),
             sqrt(Gamma_in * (a4 ** 2)) * lindblad_dissipator(um_p),
             #-----------------------#
             sqrt(Gamma_in * (a3 ** 2)) * lindblad_dissipator(ml_p),
             sqrt(Gamma_in * (a7 ** 2)) * lindblad_dissipator(ml_m),
             #-----------------------#
             sqrt(Gamma_in * (a6 ** 2)) * lindblad_dissipator(ul_p),
             sqrt(Gamma_in * (a8 ** 2)) * lindblad_dissipator(ul_m)]
    
    #-------------------------------------------------------------------------#
    #                 SOLVE SECOND-ORDER CORRELATION FUNCTION                 #
    #-------------------------------------------------------------------------#
    # Turn frequency into operator index
    operator1_index = \
        _which_transition(w0a_in, Omega_in, alpha_in, delta_in, xi_in)
    
    operator2_index = \
        _which_transition(w0b_in, Omega_in, alpha_in, delta_in, xi_in)
        
    # First operator
    operator1 = operator_dict[operator1_index]
    operator2 = operator_dict[operator2_index]
    
    # Set operators
    A = operator1.dag()
    B = operator2.dag() * operator2
    C = operator1
    
    # Calculate G2 function
    G2 = correlation_3op_1t(H_A, None, tau_in, c_ops,
                            a_op=A, b_op=B, c_op=C)
    
    # Normalise
    g2_out = G2 / G2[-1]
    g2_out = g2_out.real
    
    return g2_out

if __name__ == '__main__':
    import numpy as np
    
    #-------------------------------------------------------------------------#
    #                              PARAMETERS                                 #
    #-------------------------------------------------------------------------#
    # Atomic decay rate
    Gamma = 1.0
    # Driving amplitude
    # Omega = 0.01
    Omega = 40.0
    # Drive detuning from two-photon resonance
    delta = 0.0
    # Atomic anharmonicity
    alpha = -120.0
    # Dipole ratio
    xi = 1.5

    # Time step
    dt = 0.001
    # Max tau
    tau_max = 10.0
    # Tau list
    tau = np.arange(0, tau_max + dt, dt)
    
    # Index for transition frequency
    trans_freq_index = 3
    
    # Transition frequencies
    trans_freq = print_transition_frequencies(Omega, alpha, delta, xi, False, True)
    
    # Set frequency of filter
    w0 = trans_freq[trans_freq_index]
    
    #-------------------------------------------------------------------------#
    #                            TEST FUNCTIONS                               #
    #-------------------------------------------------------------------------#
    # # Test steady states of diagonal moments
    # m_out, u_out, l_out = steady_state_diagonal(Gamma, Omega, alpha, delta, xi)
    
    # Test transition label
    target_transition = _which_transition(w0, Omega, alpha, delta, xi)
    print("Target transition operator (w0 = {:.3f}) : sigma^{}".format(w0, target_transition))
    
    # Calculate dressed g2
    g2_dressed = dressed_g2_calc(tau, w0, Gamma, Omega, alpha, delta, xi,
                                 w0b_in=None)
    
    #-------------------------------------------------------------------------#
    #                                  PLOT                                   #
    #-------------------------------------------------------------------------#
    import matplotlib.pyplot as plt
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = 'serif'
    
    plt.close('all')
    
    fig = plt.figure(num='Dressed G2', figsize=[8, 6])
    ax = plt.gca()

    ax.plot(tau, g2_dressed)

    ax.set_xlabel(r'$\Gamma \tau$', fontsize=11)
    ax.set_ylabel(r'$g^{(2)}_{\mathrm{dressed}}(\tau)$', fontsize=11)

    fig.tight_layout()
    fig.show()
    