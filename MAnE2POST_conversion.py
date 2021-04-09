# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 15:59:54 2021

@author: saman
"""


import numpy as np
import mpmath as mp


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
    one = np.tan(sat_dec)
    two = np.sin(sat_RA - sat_LAN)
    return np.arctan2(one, two)

def energy():
    ''' energy equation to find velocity '''
    energy = (1/2)*v_inf**2 + mu_saturn/r_titan
    return energy

def velocity():
    ''' probe velocity at Titan intercept wrt Saturn '''
    E = energy()
    vi = np.sqrt(2*E)
    return vi
    
def semi_ax():
    ''' semi-major axis '''
    en = (1/2)*v_inf**2
    return mu_saturn/(2*en)

def function(e):
    a = semi_ax()
    w = np.arccos(-1/e) - np.arccos(np.cos(sat_dec)*np.cos(sat_LAN - sat_RA))
    sol = (a*(1 - e**2))/(1 - e*np.cos(w)) + r_titan
    return sol

def omega():
    e = sat_e
    w = w = np.arccos(-1/e) - np.arccos(np.cos(sat_dec)*np.cos(sat_LAN - sat_RA))
    return w

def FPA():
    ''' computes flight path angle '''
    # eccentricity and argument of periapse
    w = omega()
    e = sat_e
    FPA = np.arctan(e*np.sin(w)/(1 + e*np.cos(w)))
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
    theta = np.pi/2 - sat_LAN + gamma
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
    vxt = v_titan*np.cos(sat_LAN)
    vyt = v_titan*np.sin(sat_LAN)
    vzt = 0
    # relative velocity to Titan
    v_infx = vx - vxt
    v_infy = vy - vyt
    v_infz = vz - vzt
    return v_infx, v_infy, v_infz

def COE(rp, az):
    ''' calculate COE wrt Titan '''
    vx, vy, vz = reference_trans()
    vmag = np.sqrt((vx**2) + (vy**2) + (vz**2))
    E = (vmag**2/2)
    a = -mu_titan/(2*E)
    e = (-rp/a) + 1
    vxy = np.sqrt((vx**2) + (vy**2))
    dec = abs(np.arctan(vz/vxy))
    i = np.arccos(np.cos(dec)*np.sin(az))
    RA = np.arctan(vy/vx)
    LAN = RA - np.arcsin(np.tan(dec)*(np.cos(i)/np.sin(i)))
    w = np.arccos(-1/e) - np.arccos(np.cos(dec)*np.cos(LAN - RA))
    return E, a, e, dec, i, RA, LAN, w


if __name__ == '__main__':
    
    
    ################################ CONSTANTS ################################

    # define some constants
    v_titan   = 5.57    # Titan's orbital velocity in km/s
    r_titan   = 1221865 # Titan's orbital radius in km
    mu_saturn = 37931187.9
    mu_titan  = 0.0225*(3.986*10**5)
    
    
    ################################# INPUTS #################################
    
    alt       = 300         # km
    azimuth   = np.pi/2     # radians
    alt_titan = 2575 + alt  # km
    
    # MAnE arrival conditions from CASESMRY output file
    v_inf   =   5.0              # km/s
    sat_RA  =  72.66*(np.pi/180) # deg converted to rad
    sat_dec = -11.69*(np.pi/180) # deg converted to rad
    
    # node occurs at Titan intercept - specifies LAN
    sat_LAN = 60*(np.pi/180) # deg converted to rad
    
    
    ############################## CALCULATIONS ###############################
    
    # newton method to calculate eccentricity wrt Saturn
    sat_e = newton_method(1.5, 10**-5, 10**-6)
    
    # calculate the orbital elements wrt Titan
    E, a, e, dec, i, RA, LAN, w = COE(alt_titan, azimuth)
    
    
    
    
    
    
    