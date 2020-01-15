# -*- coding: utf-8 -*-
'''
Created on Mon Dec 9 12:45:36 2019

@author: sam
'''


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def RV2COE(mu, state):
    '''
    Converts a state vector to the classical orbital elements
    * does not include special cases *
    Args:
        mu - gravitational parameter of primary body
        state - position and velocity as a 6 element array
    Returns:
        h_mag - specific angular momentum
        E - specific mechanical energy 
        n_mag - magnitude of the node vector 
        e_mag - ecentricity 
        p - semiparameter 
        a - semimajor axis 
        i_deg - inclination in radians 
        Omega_deg - longitude of the ascending node in radians 
        omega_deg - argument of perigee in radians 
        true_deg - true anomally in radians
    '''
    
    tol = 1*10**-6
    K = [0, 0, 1]
    
    # position vector
    r = state[:3]
    r_mag = np.linalg.norm(r)
    
    # velocity vector
    v = state[3:]
    v_mag = np.linalg.norm(v)
    
    # specific angular momentum
    h = np.cross(r, v)
    h_mag = np.linalg.norm(h)
    
    # node vector
    n = np.cross(K, h)
    n_mag = np.linalg.norm(n)
    
    # eccentricity
    e = ((v_mag**2 - (mu/r_mag))*r - np.dot(r, v)*v)/mu
    e_mag = np.linalg.norm(e)
    
    # specific energy
    E = (v_mag**2/2) - (mu/r_mag)
    
    # semiparameter and semimajor axis depending on orbit type
    if 1 - tol < e_mag < 1 + tol:
        p = (h_mag**2/mu)
        a = np.inf
    else:
        a = -mu/(2*E)
        p = a*(1 - e_mag**2)
        
    # inclination
    i = np.arccos(h[2]/h_mag)
    
    # longitude of the ascending node
    Omega = np.arccos(n[0]/n_mag)
    if n[1] < tol:
        Omega = 2*np.pi - Omega
    
    # argument of perigee
    omega = np.arccos((np.dot(n, e))/(n_mag*e_mag))
    if e[2] < 0:
        omega = 2*np.pi - omega
       
    # true anomaly
    true = np.arccos(np.dot(e, r)/(e_mag*r_mag))
    if np.dot(r, v) < 0:
        true = 2*np.pi - true
    
    return h_mag, E, n_mag, e_mag, p, a, i, Omega, omega, true


def COE2RV(mu, p, e, i, Omega, omega, true):
    '''
    Converts the classical orbital elements to a state vector
    * does not include special cases *
    Args:
        mu - gravitational parameter of primary body
        p - semiparameter 
        e - ecentricity  
        i - inclination in radians 
        Omega - longitude of the ascending node in radians 
        omega - argument of perigee in radians 
        true - true anomally in radians
    Returns:
        state - state vector
    '''
    
    # position and velocity PQW vectors
    r_pqw = np.matrix([[(p*np.cos(true))/(1 + e*np.cos(true))],
                       [(p*np.sin(true))/(1 + e*np.cos(true))], 
                       [0]])
    v_pqw = np.matrix([[-np.sqrt(mu/p)*np.sin(true)],
                       [np.sqrt(mu/p)*(e + np.cos(true))], 
                       [0]])
    
    # PQW to IJK transformation matrix
    pqw2ijk = np.matrix([[ np.cos(Omega)*np.cos(omega) - np.sin(Omega)*np.sin(omega)*np.cos(i),
                          -np.cos(Omega)*np.sin(omega) - np.sin(Omega)*np.cos(omega)*np.cos(i),
                           np.sin(Omega)*np.sin(i)],
                         [ np.sin(Omega)*np.cos(omega) + np.cos(Omega)*np.sin(omega)*np.cos(i),
                          -np.sin(Omega)*np.sin(omega) + np.cos(Omega)*np.cos(omega)*np.cos(i),
                          -np.cos(Omega)*np.sin(i)],
                         [ np.sin(omega)*np.sin(i),
                           np.cos(omega)*np.sin(i),
                           np.cos(i)]])
    
    # position and velocity IJK vectors
    r_ijk = pqw2ijk*r_pqw
    v_ijk = pqw2ijk*v_pqw
    
    # convert to a single state vector
    r = r_ijk.getA1()          
    v = v_ijk.getA1()
    state = np.concatenate((r, v))
    
    return state


def state(r, v, gamma, inc):
    '''
    Converts position and velocity magnitudes into a state array
    Args:
        r - position vector
        v - velocity magnitude
        gamma - flight path angle in degrees
        inc - inclination angle in degrees
    Returns:
        state - state array
    '''
    g = gamma*(np.pi/180)
    i = inc*(np.pi/180)
    state = np.concatenate((r, [v*np.cos(i)*np.sin(g), 
                                v*np.cos(i)*np.cos(g), 
                                v*np.sin(i)]))
    return state


if __name__ == '__main__':
    
    # define some constants
    mu = 37.931*10**6 # saturns gravitational parameter
    saturn_equatorial = 60268 # km
    saturn_polar = 54364 # km
    r_titan = 1.2*10**6 # titans orbit radius 
    r_encel = 238000 # km
    v_titan = 5.57 # km/s
    
    gamma = np.linspace(1, 90, 20) # flight path angle array
    inc = np.linspace(1, 90, 20) # inclination angle array
    data = pd.DataFrame([]) # initialize an empty dataframe
    
    for g in gamma:
        for i in inc:
            v1 = 7.94 # max v1 to prevent escape from Saturn system in km/s
            r = r_encel + 10
            
            while r > r_encel:
                
                v1 = v1 - 0.01 # max v1 to prevent escape from Saturn system in km/s
                pos = np.array([r_titan, 0, 0]) # initial satellite position
                statevec = state(pos, v1, g, i)
                # calculate classical orbital elements
                h, E, n, e, p, a, incl, Omega, omega, true = RV2COE(mu, statevec)
                # True Anomally at Descending Node
                nu = np.pi - omega
                # Radius at Descending Node
                r = p/(1 + e*np.cos(nu))
                # print some stuff to see progress
                print('%-13s %-20s % -20s %-20s %-20s'  
                  %('v1', 'gamma', 'inclination', 'r', 'w'))
                print('%5.1f %20.10f %20.10f %20.10f %20.10f' 
                  %(v1, g, i, r, nu))
            
            # Add to the Dataframe
            data = data.append(pd.DataFrame({'Saturn fpa (deg)': g, 
                                             'Saturn inclination (deg)': i,
                                             'Saturn v1 max': abs(v1),
                                             'Radius at Descending Node': r,
                                             'Angular Momentum': h,
                                             'Mechanical Energy': E,
                                             'Eccentricity': e,
                                             'Semimajor Axis': a,
                                             'Inclination': incl*(180/np.pi),
                                             'Longitude of Ascending Node': Omega*(180/np.pi),
                                             'Argument of Perigee': omega*(180/np.pi),
                                             'True Anomally': nu*(180/np.pi)},
                                             index = [0]), ignore_index = True)
#    v1 = 3.24
#    g = 20
#    i = 20        
#    pos = np.array([r_titan, 0, 0]) # initial satellite position
#    statevec = state(pos, v1, g, i)
#    # calculate classical orbital elements
#    h, E, n, e, p, a, incl, Omega, omega, true = RV2COE(mu, statevec)
#    nu = np.pi - omega
#    r = p/(1 + e*np.cos(nu)) 
#    print(' ')
#    print('Position: {} km \
#          \nVelocity: {} km/s \
#          '.format(statevec[:3], statevec[3:]))
#    print(' ')
#    print('Specific angular momentum: {} km^2/s \
#          \nSpecific mechanical energy: {} J \
#          \nNode magnitude: {} km^2/s \
#          \nEccentricity magnitude: {} \
#          \nSemiparameter: {} km \
#          \nSemimajor axis: {} km \
#          \nInclination: {} deg \
#          \nLongitude of the ascending node: {} deg \
#          \nArgument of perigee: {} deg \
#          \nTrue anomally: {} deg \
#          '.format(h, E, n, e, p, a, incl*(180/np.pi), 
#                   Omega*(180/np.pi), omega*(180/np.pi), true*(180/np.pi))) 
#    nu = 180 - omega   
#    r = p/(1 + e*np.cos(nu))    
#    print(' ')
#    print('Radius at Descending Node: {} km'.format(r))
        
