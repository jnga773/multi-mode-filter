#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 10:40:03 2022

@author: jnga773
"""

def steady_state_moments(kappa_in, Delta_in, K_in):
    """
    Calculates the steady states of the operator moment equations of the
    degenerate parametric oscillator, up to fourth-order.
    
    Parameters
    ----------
    kappa_in : float
        Cavity decay rate
    Delta_in : float
        Drive detuning from resonance (\omega_{C} - 2 \omega_{P})
    K_in : float
        Coupling rate (\lambda in TeX)
    
    Returns
    -------
    ss_out : dictionary
         Dictionary of all the steady state moment values
    
    """
    from numpy import matrix, matmul
    from numpy.linalg import inv
        
    #=========================================================================#
    #              DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS             #
    #=========================================================================#
    #-----------------------#
    #     Moment Labels     #
    #-----------------------#
    # First-order
    # a = 0; at = 1
    # Second-order
    a2 = 0; ata = 1; at2 = 2
    # Third-order
    a3 = 0; ata2 = 1; at2a = 2; at3 = 3
    # Fourth-order
    a4 = 0; ata3 = 1; at2a2 = 2; at3a = 3; at4 = 4
    
    #=========================================================================#
    #                      CALCULATE STEADY-STATE MOMENTS                     #
    #=========================================================================#
    #---------------------#
    #     First-Order     #
    #---------------------#
    # # Evolution matrix
    Mat = matrix([[-(kappa_in + 1j * Delta_in), K_in],
                  [K_in, -(kappa_in - 1j * Delta_in)]
                  ], dtype='complex')
    # Steady state moments are actually just zero so
    ss_out = {'a': 0.0j, 'at': 0.0j}
    
    #----------------------#
    #     Second-Order     #
    #----------------------#
    # Evolution matrix
    Mat = matrix([[-2 * (kappa_in + 1j * Delta_in), 2 * K_in, 0],
                  [K_in, -2 * kappa_in, K_in],
                  [0, 2 * K_in, -2 * (kappa_in - 1j * Delta_in)]
                  ], dtype='complex')
    # Non homogeneous vector
    B = matrix([[K_in],
                [0],
                [K_in]
                ], dtype='complex')
    # Invert matrix
    Inv_temp = inv(Mat)
    # Multiply
    ss_temp = -matmul(Inv_temp, B)
    
    # Output array
    ss_out['a2'] = ss_temp[a2, 0]
    ss_out['ata'] = ss_temp[ata, 0]
    ss_out['at2'] = ss_temp[at2, 0]
    
    #---------------------#
    #     Third-Order     #
    #---------------------#
    # Evolution matrix
    Mat = matrix([[-3 * (kappa_in + 1j * Delta_in), 3 * K_in, 0, 0],
                  [K_in, -(3 * kappa_in + 1j * Delta_in), 2 * K_in, 0],
                  [0, 2 * K_in, -(3 * kappa_in - 1j * Delta_in), K_in],
                  [0, 0, 3 * K_in, -3 * (kappa_in - 1j * Delta_in)]
                  ], dtype='complex')
    # Non homogeneous vector
    B = matrix([[3 * K_in * ss_out['a']],
                [K_in * ss_out['at']],
                [K_in * ss_out['a']],
                [3 * K_in * ss_out['at']]
                ], dtype='complex')
    # Invert matrix
    Inv_temp = inv(Mat)
    # Multiply
    ss_temp = -matmul(Inv_temp, B)
    
    # Output array
    ss_out['a3'] = ss_temp[a3, 0]
    ss_out['ata2'] = ss_temp[ata2, 0]
    ss_out['at2a'] = ss_temp[at2a, 0]
    ss_out['at3'] = ss_temp[at3, 0]
    
    #----------------------#
    #     Fourth-Order     #
    #----------------------#
    # Evolution matrix
    Mat = matrix([[-4 * (kappa_in + 1j * Delta_in), 4 * K_in, 0, 0, 0],
                  [K_in, -(4 * kappa_in + 2j * Delta_in), 3 * K_in, 0, 0],
                  [0, 2 * K_in, -4 * kappa_in, 2 * K_in, 0],
                  [0, 0, 3 * K_in, -(4 * kappa_in - 2j * Delta_in), K_in],
                  [0, 0, 0, 4 * K_in, -4 * (kappa_in - 1j * Delta_in)]
                  ], dtype='complex')
    # Non homogeneous vector
    B = matrix([[6 * K_in * ss_out['a2']],
                [3 * K_in * ss_out['ata']],
                [K_in * (ss_out['a2'] + ss_out['at2'])],
                [3 * K_in * ss_out['ata']],
                [6 * K_in * ss_out['at2']]
                ], dtype='complex')
    # Invert matrix
    Inv_temp = inv(Mat)
    # Multiply
    ss_temp = -matmul(Inv_temp, B)
    
    # Output array
    ss_out['a4'] = ss_temp[a4, 0]
    ss_out['ata3'] = ss_temp[ata3, 0]
    ss_out['at2a2'] = ss_temp[at2a2, 0]
    ss_out['at3a'] = ss_temp[at3a, 0]
    ss_out['at4'] = ss_temp[at4, 0]
    
    #=========================================================================#
    #                                  OUTPUT                                 #
    #=========================================================================#
    return ss_out

def time_evolution_moments_RK4(kappa_in, Delta_in, K_in, t_in, output='all',
                               initial_conditions=None):
    """
    Calculates the time evolution of the operator moment equations of the
    degenerate parametric oscillator, up to fourth-order, using Runge-Kutta 
    4th order.
    
    Parameters
    ----------
    kappa_in : float
        Cavity decay rate
    Delta_in : float
        Drive detuning from resonance (\omega_{C} - 2 \omega_{P})
    K_in : float
        Coupling rate (\lambda in TeX)
    t_in : array-like
        List of times to integrate over
    output : string
        Either 'all' or 'g1', or 'g2', depending on what I want to calculate.
        'all' will time integrate every moment, 'g1' will only integrate the
        first-order moments, 'g2' will integrate up to second-order.
    initial_conditions : dictionary
        Initial conditions for moment equations. If None, will default to zero
        photons in the cavity <ata(0)> = 0. Otherwise, if supplied with a
        dictionary, will use those as initial conditions.
    
    Returns
    -------
    moments_out : dictionary
         Dictionary of all the steady state moment values
    
    """
    from numpy import zeros, matrix, matmul

    # Check if output is 'all', 'g1', or 'g2'
    if output not in ['all', 'g1', 'g2']:
        from sys import exit
        exit("69! output should be {'all', 'g1', 'g2'}!")
    
    # Check if initial_conditions is None or a dictionary
    if initial_conditions is not None:
        # Initial condition is not blank, so check if dictionary
        if type(initial_conditions) is not dict:
            exit("70! initial_conditions should be None or a dictionary!")
        else:
            initial_condition_check = True
    else:
        initial_condition_check = False

    # Time step
    dt = t_in[1] - t_in[0]
    
    # If calculating g1 or g2, the non-homogenous part is multipled by the
    # steady state value of <a>_ss or <ata>_ss
    if initial_condition_check is True and output == 'g1' or output == 'g2':
        multiplier = initial_conditions['multiplier']
    elif initial_condition_check is False:
        multiplier = 1.0

    #=========================================================================#
    #              DEFINING ANALYTIC MATRICES/EIGENVALUES/VECTORS             #
    #=========================================================================#
    #-----------------------#
    #     Moment Labels     #
    #-----------------------#
    # First-order
    a = 0; at = 1

    # Second-order
    a2 = 0; ata = 1; at2 = 2

    # Third-order
    a3 = 0; ata2 = 1; at2a = 2; at3 = 3

    # Fourth-order
    a4 = 0; ata3 = 1; at2a2 = 2; at3a = 3; at4 = 4

    #-----------------------#
    #     Moment Arrays     #
    #-----------------------#
    # First-order moments
    m1 = zeros(shape=(2, 1), dtype='complex')

    if output == 'all' or output == 'g2':
        # Second-order moments
        m2 = zeros(shape=(3, 1), dtype='complex')

    if output == 'all':
        # Third-order moments
        m3 = zeros(shape=(4, 1), dtype='complex')
        # Fourth-order moments
        m4 = zeros(shape=(5, 1), dtype='complex')   

    #---------------------#
    #     Data Arrays     #
    #---------------------#
    # First-order moments
    a_t = zeros(shape=len(t_in), dtype='complex')
    at_t = zeros(shape=len(t_in), dtype='complex')

    if output == 'all' or output == 'g2':
        # Second-order moments
        a2_t = zeros(shape=len(t_in), dtype='complex')
        ata_t = zeros(shape=len(t_in), dtype='complex')
        at2_t = zeros(shape=len(t_in), dtype='complex')
    
    if output == 'all':
        # Third-order moments
        a3_t = zeros(shape=len(t_in), dtype='complex')
        ata2_t = zeros(shape=len(t_in), dtype='complex')
        at2a_t = zeros(shape=len(t_in), dtype='complex')
        at3_t = zeros(shape=len(t_in), dtype='complex')
        
        # Fourth-order moments
        a4_t = zeros(shape=len(t_in), dtype='complex')
        ata3_t = zeros(shape=len(t_in), dtype='complex')
        at2a2_t = zeros(shape=len(t_in), dtype='complex')
        at3a_t = zeros(shape=len(t_in), dtype='complex')
        at4_t = zeros(shape=len(t_in), dtype='complex')

    #----------------------------#
    #     Evolution Matrices     #
    #----------------------------#
    # First-order moments
    Mat1 = matrix([[-(kappa_in + 1j * Delta_in), K_in],
                   [K_in, -(kappa_in - 1j * Delta_in)]
                   ], dtype='complex')

    if output == 'g2' or output =='all':
        # Second-order moments
        Mat2 = matrix([[-2 * (kappa_in + 1j * Delta_in), 2 * K_in, 0],
                       [K_in, -2 * kappa_in, K_in],
                       [0, 2 * K_in, -2 * (kappa_in - 1j * Delta_in)]
                       ], dtype='complex')
        # Non homogeneous vector
        B2 = matrix([[K_in * multiplier],
                     [0],
                     [K_in * multiplier]
                     ], dtype='complex')

    if output == 'all':
        # Third-order moments
        Mat3 = matrix([[-3 * (kappa_in + 1j * Delta_in), 3 * K_in, 0, 0],
                       [K_in, -(3 * kappa_in + 1j * Delta_in), 2 * K_in, 0],
                       [0, 2 * K_in, -(3 * kappa_in - 1j * Delta_in), K_in],
                       [0, 0, 3 * K_in, -3 * (kappa_in - 1j * Delta_in)]
                       ], dtype='complex')

        # Fourth-order moments
        Mat4 = matrix([[-4 * (kappa_in + 1j * Delta_in), 4 * K_in, 0, 0, 0],
                       [K_in, -(4 * kappa_in + 2j * Delta_in), 3 * K_in, 0, 0],
                       [0, 2 * K_in, -4 * kappa_in, 2 * K_in, 0],
                       [0, 0, 3 * K_in, -(4 * kappa_in - 2j * Delta_in), K_in],
                       [0, 0, 0, 4 * K_in, -4 * (kappa_in - 1j * Delta_in)]
                       ], dtype='complex') 

    #----------------------------#
    #     Initial Conditions     #
    #----------------------------#
    # If initial_conditions is None, all moments are 0.0 at t=0. Otherwise,
    # feed initial conditions from input dictionary
    if initial_conditions is not None:
        # First-order
        m1[a, 0] = initial_conditions['a']
        m1[at, 0] = initial_conditions['at']
        
        if output == 'all' or output == 'g2':
            # Second-order
            m2[a2, 0] = initial_conditions['a2']
            m2[ata, 0] = initial_conditions['ata']
            m2[at2, 0] = initial_conditions['at2']
        
        if output == 'all':
            # Third-order
            m3[a3, 0] = initial_conditions['a3']
            m3[ata2, 0] = initial_conditions['ata2']
            m3[at2a, 0] = initial_conditions['at2a']
            m3[at3, 0] = initial_conditions['at3']
            
            # Fourth-order
            m4[a4, 0] = initial_conditions['a4']
            m4[ata3, 0] = initial_conditions['ata3']
            m4[at2a2, 0] = initial_conditions['at2a2']
            m4[at3a, 0] = initial_conditions['at3a']
            m4[at4, 0] = initial_conditions['at4']
        
    #=========================================================================#
    #                  CALCULATE TIME INTEGRATE OF MOMENTS                    #
    #=========================================================================#
    # Loop through all time steps
    for index, time in enumerate(t_in):
        #=====================================================================#
        #                      CALCULATE AND WRITE DATA                       #
        #=====================================================================#
        # Output data to arrays
        # First-order
        a_t[index] = m1[a, 0]
        at_t[index] = m1[at, 0]
        
        if output == 'all' or output == 'g2':
            # Second-order
            a2_t[index] = m2[a2, 0]
            ata_t[index] = m2[ata, 0]
            at2_t[index] = m2[at2, 0]
        
        if output == 'all':
            # Third-order
            a3_t[index] = m3[a3, 0]
            ata2_t[index] = m3[ata2, 0]
            at2a_t[index] = m3[at2a, 0]
            at3_t[index] = m3[at3, 0]
            
            # Fourth-order
            a4_t[index] = m4[a4, 0]
            ata3_t[index] = m4[ata3, 0]
            at2a2_t[index] = m4[at2a2, 0]
            at3a_t[index] = m4[at3a, 0]
            at4_t[index] = m4[at4, 0]
        
        #=====================================================================#
        #              CALCULATE USING FOURTH-ORDER RUNGE-KUTTA               #
        #=====================================================================#
        #---------------------#
        #     First-order     #
        #---------------------#
        k1_m1 = dt * matmul(Mat1, m1)
        k2_m1 = dt * matmul(Mat1, (m1 + 0.5 * k1_m1))
        k3_m1 = dt * matmul(Mat1, (m1 + 0.5 * k2_m1))
        k4_m1 = dt * matmul(Mat1, (m1 + k3_m1))
        
        #----------------------#
        #     Second-order     #
        #----------------------#
        if output == 'all' or output == 'g2':
            k1_m2 = dt * matmul(Mat2, m2) + dt * B2
            k2_m2 = dt * matmul(Mat2, (m2 + 0.5 * k1_m2)) + dt * B2
            k3_m2 = dt * matmul(Mat2, (m2 + 0.5 * k2_m2)) + dt * B2
            k4_m2 = dt * matmul(Mat2, (m2 + k3_m2)) + dt * B2
        
        #---------------------#
        #     Third-order     #
        #---------------------#
        if output == 'all':
            B3 = matrix([[3 * K_in * m1[a, 0]],
                         [K_in * m1[at, 0]],
                         [K_in * m1[a, 0]],
                         [3 * K_in * m1[at, 0]]
                         ], dtype='complex')
            k1_m3 = dt * matmul(Mat3, m3) + dt * B3
            
            B3 = matrix([[3 * K_in * (m1[a, 0] + 0.5 * k1_m1[a, 0])],
                         [K_in * (m1[at, 0] + 0.5 * k1_m1[at, 0])],
                         [K_in * (m1[a, 0] + 0.5 * k1_m1[a, 0])],
                         [3 * K_in * (m1[at, 0] + 0.5 * k1_m1[at, 0])]
                         ], dtype='complex')
            k2_m3 = dt * matmul(Mat3, (m3 + 0.5 * k1_m3)) + dt * B3
            
            B3 = matrix([[3 * K_in * (m1[a, 0] + 0.5 * k2_m1[a, 0])],
                         [K_in * (m1[at, 0] + 0.5 * k2_m1[at, 0])],
                         [K_in * (m1[a, 0] + 0.5 * k2_m1[a, 0])],
                         [3 * K_in * (m1[at, 0] + 0.5 * k2_m1[at, 0])]
                         ],dtype='complex')
            k3_m3 = dt * matmul(Mat3, (m3 + 0.5 * k2_m3)) + dt * B3
            
            B3 = matrix([[3 * K_in * (m1[a, 0] + k3_m1[a, 0])],
                         [K_in * (m1[at, 0] + k3_m1[at, 0])],
                         [K_in * (m1[a, 0] + k3_m1[a, 0])],
                         [3 * K_in * (m1[at, 0] + k3_m1[at, 0])]
                         ], dtype='complex')
            k4_m3 = dt * matmul(Mat3, (m3 + k3_m3)) + dt * B3
        
        #----------------------#
        #     Fourth-order     #
        #----------------------#
        if output == 'all':
            B4 = matrix([[6 * K_in * m2[a2, 0]],
                         [3 * K_in * m2[ata, 0]],
                         [K_in * (m2[a2, 0] + m2[at2, 0])],
                         [3 * K_in * m2[ata, 0]],
                         [6 * K_in * m2[at2, 0]]
                         ], dtype='complex')
            k1_m4 = dt * matmul(Mat4, m4) + dt * B4
            
            B4 = matrix([[6 * K_in * (m2[a2, 0] + 0.5 * k1_m2[a2, 0])],
                         [3 * K_in * (m2[ata, 0] + 0.5 * k1_m2[ata, 0])],
                         [K_in * ((m2[a2, 0] + 0.5 * k1_m2[a2, 0]) + (m2[at2, 0] + 0.5 * k1_m2[at2, 0]))],
                         [3 * K_in * (m2[ata, 0] + 0.5 * k1_m2[ata, 0])],
                         [6 * K_in * (m2[at2, 0] + 0.5 * k1_m2[at2, 0])]
                         ], dtype='complex')
            k2_m4 = dt * matmul(Mat4, (m4 + 0.5 * k1_m4)) + dt * B4
            
            B4 = matrix([[6 * K_in * (m2[a2, 0] + 0.5 * k2_m2[a2, 0])],
                         [3 * K_in * (m2[ata, 0] + 0.5 * k2_m2[ata, 0])],
                         [K_in * ((m2[a2, 0] + 0.5 * k2_m2[a2, 0]) + (m2[at2, 0] + 0.5 * k2_m2[at2, 0]))],
                         [3 * K_in * (m2[ata, 0] + 0.5 * k2_m2[ata, 0])],
                         [6 * K_in * (m2[at2, 0] + 0.5 * k2_m2[at2, 0])]
                         ], dtype='complex')
            k3_m4 = dt * matmul(Mat4, (m4 + 0.5 * k2_m4)) + dt * B4
            
            B4 = matrix([[6 * K_in * (m2[a2, 0] + k3_m2[a2, 0])],
                         [3 * K_in * (m2[ata, 0] + k3_m2[ata, 0])],
                         [K_in * ((m2[a2, 0] + k3_m2[a2, 0]) + (m2[at2, 0] + k3_m2[at2, 0]))],
                         [3 * K_in * (m2[ata, 0] + k3_m2[ata, 0])],
                         [6 * K_in * (m2[at2, 0] + k3_m2[at2, 0])]
                         ], dtype='complex')
            k4_m4 = dt * matmul(Mat4, (m4 + k3_m4)) + dt * B4
        
        #=====================================================================#
        #               UPDATE ARRAYS FROM RUNGE-KUTTA ARRAYS                 #
        #=====================================================================#
        # First-order
        m1 += (1/6) * (k1_m1 + 2.0 * (k2_m1 + k3_m1) + k4_m1)
        if output == 'all' or output == 'g2':
            # Second-order
            m2 += (1/6) * (k1_m2 + 2.0 * (k2_m2 + k3_m2) + k4_m2)
        if output == 'all':
            # Third-order
            m3 += (1/6) * (k1_m3 + 2.0 * (k2_m3 + k3_m3) + k4_m3) 
            # Fourth-order
            m4 += (1/6) * (k1_m4 + 2.0 * (k2_m4 + k3_m4) + k4_m4)  
        
        
    #=========================================================================#
    #                                  OUTPUT                                 #
    #=========================================================================#
    if output == 'all':
        # output all moments
        moments_out = {'a': a_t, 'at': at_t,
                       'a2': a2_t, 'ata': ata_t, 'at2': at2_t,
                       'a3': a3_t, 'ata2': ata2_t, 'at2a': at2a_t, 'at3': at3_t,
                       'a4': a4_t, 'ata3': ata3_t, 'at2a2': at2a2_t, 'at3a': at3a_t, 'at4': at4_t}
    
    elif output =='g1':
        # Output first-order moment < a^{\dagger}(t) >
        moments_out = at_t
    
    elif output == 'g2':
        # Output second-order moment < a^{\dagger} a(t) >
        moments_out = ata_t
        
    return moments_out

def g1_RK4(kappa_in, Delta_in, K_in, tau_in):
    """
    Calculates the time evolution of the first-order correlation function,
        g^{(1)}(\tau) = < a^{\dagger}(\tau) a(0) >
    of the
    degenerate parametric oscillator, up to fourth-order, using Runge-Kutta 
    4th order.
    
    Parameters
    ----------
    kappa_in : float
        Cavity decay rate
    Delta_in : float
        Drive detuning from resonance (\omega_{C} - 2 \omega_{P})
    K_in : float
        Coupling rate (\lambda in TeX)
    tau_in : array-like
        List of times to integrate over
    
    Returns
    -------
    moments_out : dictionary
         Dictionary of all the steady state moment values
    
    """
    # Calculate steady states 
    steady_states = steady_state_moments(kappa_in, Delta_in, K_in)
    # Use these as the initial conditions for the time integration
    init_cond = {'a': steady_states['a2'],
                 'at': steady_states['ata'],
                 'multiplier': steady_states['a']}
    
    # Calculate G1
    G1 = time_evolution_moments_RK4(kappa_in, Delta_in, K_in, tau_in, output='g1',
                                    initial_conditions=init_cond)
    # Normalise by the steady state value of < a^{\dagger} a >
    g1 = G1 / steady_states['ata']
    
    return g1

def g2_RK4(kappa_in, Delta_in, K_in, tau_in):
    """
    Calculates the time evolution of the second-order correlation function,
        g^{(2)}(\tau) = < a^{\dagger}(0) a^{\dagger} a(\tau) a(0) >
    of the
    degenerate parametric oscillator, up to fourth-order, using Runge-Kutta 
    4th order.
    
    Parameters
    ----------
    kappa_in : float
        Cavity decay rate
    Delta_in : float
        Drive detuning from resonance (\omega_{C} - 2 \omega_{P})
    K_in : float
        Coupling rate (\lambda in TeX)
    tau_in : array-like
        List of times to integrate over
    
    Returns
    -------
    moments_out : dictionary
         Dictionary of all the steady state moment values
    
    """
    # Calculate steady states 
    steady_states = steady_state_moments(kappa_in, Delta_in, K_in)
    # Use these as the initial conditions for the time integration
    init_cond = {'a': steady_states['ata2'],
                 'at': steady_states['at2a'],
                 'a2': steady_states['ata3'],
                 'ata': steady_states['at2a2'],
                 'at2': steady_states['at3a'],
                 'multiplier': steady_states['ata']}
    # Calculate G1
    G2 = time_evolution_moments_RK4(kappa_in, Delta_in, K_in, tau_in, output='g2',
                                    initial_conditions=init_cond)
    # Normalise by the steady state value of < a^{\dagger} a >
    g2 = G2.real / (steady_states['ata'].real ** 2)
    
    print("Steady state photon number = {:.3f}".format(steady_states['ata'].real))
    
    return g2

if __name__ == '__main__':
    import numpy as np
    
    # Test parameters
    # Decay rate
    kappa = 1.0
    # dw = wc - wp
    dw = 0.5
    # Coupling rate
    K = 0.5
    # Max time
    t_max = 10.0
    # Time step
    dt = 0.001
    # Time list
    t = np.arange(0, t_max + dt, dt)

    # Calculate steady states
    steady_states = steady_state_moments(kappa, dw, K)
    
    # Calculate time integration
    test = time_evolution_moments_RK4(kappa, dw, K, t)