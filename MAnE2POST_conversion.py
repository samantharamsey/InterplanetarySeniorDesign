# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 15:59:54 2021

@author: saman
"""


import numpy as np
from D_potato import RV2COE


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
    return np.array([vx, vy, vz, rx, ry, rz])

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
    
    # do stuff
    # titan_state = state()
    # h, E, n, e, p, a, i, LANt, omegat, nu = RV2COE(mu_titan, titan_state)
    
    eccentricity = newton_method(1.5, 10**-5, 10**-6)