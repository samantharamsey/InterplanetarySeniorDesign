# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 15:59:54 2021

@author: saman
"""


import numpy as np


def newton_method(x0, step_size, tolerance):
    '''
    Determines the roots of a non-linear single variable function using 
    derivative estimation and Taylor Series Expansion
    Args:
        x0 - initial condition estimate
        tolerance - the tolerance for convergence
        step_size - determines the size of delta x
    '''
    f0 = function(x0)
    residual = abs(f0) 
    iteration = 1
    fx = (function(x0 + step_size) - f0) / step_size
    print('%-13s %-20s %-20s %-20s' 
          %('Iterations', 'Function Value', 'Root Value', 'Residual'))
    while residual > tolerance:
        x1 = x0 + ((0 - f0)/fx)
        f1 = function(x1)
        residual = abs(f1)
        f0 = f1
        x0 = x1
        fx = (function(x0 + step_size) - f0) / step_size
        print('%5.1f %20.10f %20.10f %20.10f' 
              %(iteration, f1, x1, residual))
        iteration = iteration + 1
    return x1
        
def inclination():
    ''' calculates inclination using Eqn 5.28 from Tewari '''
    one = np.tan(dec)
    two = np.sin(RA - LAN)
    return np.arctan2(one, two)

def energy():
    ''' energy equation to find velocity '''
    energy = (1/2)*v_inf**2 - mu_saturn/r_titan
    return energy

def velocity():
    ''' probe velocity at Titan intercept wrt Saturn '''
    E = energy()
    vi = np.sqrt(2*abs(E))
    return vi
    
def semi_ax():
    ''' semi-major axis '''
    en = energy()
    return -mu_saturn/(2*en)

def function(e):
    a = semi_ax()
    w = np.arccos(-1/e) - np.arccos(np.cos(dec)*np.cos(LAN - RA))
    sol = (a*(1 - e**2))/(1 - e*np.cos(w)) + r_titan
    return sol

def omega():
    e = eccentricity
    w = w = np.arccos(-1/e) - np.arccos(np.cos(dec)*np.cos(LAN - RA))
    return w

def FPA():
    ''' computes flight path angle '''
    # eccentricity and argument of periapse
    w = omega()
    e = eccentricity
    FPA = np.arctan(e*np.sin(w/(1 + e*np.cos(w))))
    return FPA
    
def vel_components():
    ''' breaks up vi into its vector components '''
    vi = velocity()
    # break vi into the z and xy plane components
    i = inclination()
    vz = np.sin(i)*vi
    vxy = np.cos(i)*vi
    # flight path angle
    gamma = FPA()
    theta = np.pi/2 - LAN + gamma
    phi = np.pi/2 - theta
    # use phi to break into x and y components
    vx = vxy*np.sin(phi)
    vy = vxy*np.cos(phi)
    return vx, vy, vz
    
def reference_trans():
    ''' converts velocity in saturn reference frame to titan '''
    # components wrt Saturn
    vx, vy, vz = vel_components()
    # Titan orbital velocity components
    vxt = v_titan*np.sin(LAN)
    vyt = v_titan*np.cos(LAN)
    vzt = 0
    # relative velocity to Titan
    v_infx = vx - vxt
    v_infy = vy - vyt
    v_infz = vz - vzt
    return v_infx, v_infy, v_infz

def state():
    ''' creates state vector '''
    rx = r_titan*np.sin(LAN)
    ry = r_titan*np.cos(LAN)
    rz = 0
    vx, vy, vz = reference_trans()
    return np.array([rx, ry, rz, vx, vy, vz])

def RV2COE(mu, state):
    '''
    Converts a state vector to the classical orbital elements
    *** does not include special cases ***
    Args:
        mu - gravitational parameter of primary body
        state - position and velocity as a 6 element array
    Returns:
        h_mag - specific angular momentum
        E - specific mechanical energy
        n_mag - magnitude of the node vector
        e_mag - eccentricity
        p - semiparameter
        a - semimajor axis
        i_deg - inclination in radians
        Omega_deg - longitude of the ascending node in radians
        omega_deg - argument of perigee in radians
        true_deg - true anomaly in radians
    '''

    tol = 1*10**-6
    K = [0, 0, 1]

    # position vector
    r = state[:3]
    r_mag = np.linalg.norm(r)

    # velocity vector
    v = state[3:]
    v_mag = np.linalg.norm(v)
    print(v_mag)

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
        p = (h_mag**2 / mu)
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

if __name__ == '__main__':
    
    # define some constants
    v_titan = 5.57 # Titan's orbital velocity in km/s
    r_titan = 1221865 # Titan's orbital radius in km
    mu_saturn = 37931187.9
    mu_titan = 0.0225*(3.986*10**5)
    
    # MAnE arrival conditions from CASESMRY output file
    v_inf = 5.0 # km/s
    RA    =  72.66*(np.pi/180) # deg converted to rad
    dec   = -11.69*(np.pi/180) # deg converted to rad
    
    # node occurs at Titan intercept - specifies LAN
    LAN = 60*(np.pi/180) # deg converted to rad
    
    # newton method to calculate eccentricity
    eccentricity = newton_method(1.5, 10**-5, 10**-6)
    
    # do stuff
    titan_state = state()
    h, E, n, e, p, a, i, LANt, omegat, nu = RV2COE(mu_titan, titan_state)
    