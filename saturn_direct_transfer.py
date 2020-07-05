'''
Created on Sat May 30 13:00:24 2020
@author: sam

Attempting a patched conic trajectory to Saturn.

Required:
    SPICE Kernels
'''

import pandas as pd
import numpy as np
import spiceypy as spice
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D


def furnsh_kernels():
    '''
    Furnshes the kernels needed for a Saturn mission
    '''
    
    spice.tkvrsn('TOOLKIT')
    path = r'C:/Spice_Kernels/'
    spice.furnsh(path + r'de430.bsp')
    spice.furnsh(path + r'naif0009.tls')
    spice.furnsh(path + r'sat425.bsp')
    spice.furnsh(path + r'cpck05Mar2004.tpc')
    spice.furnsh(path + r'020514_SE_SAT105.bsp')
    spice.furnsh(path + r'981005_PLTEPH-DE405S.bsp')
    spice.furnsh(path + r'030201AP_SK_SM546_T45.bsp')
    spice.furnsh(path + r'04135_04171pc_psiv2.bc')
    spice.furnsh(path + r'sat425.inp')
    spice.furnsh(path + r'sat425.cmt')


def distance(r):
    '''
    Standard euclidean distance fromula returns the scalar
    '''
    return np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)


def cos_law(s1, s2, angle):
    '''
    Cosine law returns the final side
    '''
    return np.sqrt(s1**2.0 + s2**2.0 - 2.0*s1*s2*np.cos(angle))
    

def Lambert_Battin(mu, r1, r2, dt, tm = 1, orbit_type = 1):
    '''
    Lambert solver using Battin's method from Vallado's Fundamentals of 
    Astrodynamics and Applications
    Args:
        mu - primary body gravitational parameter
        r1 - initial position
        r2 - final position
        dt - change in time
        tm - transfer method term
             +/- 1 depending on short/long term transfer
    Returns:
        v1 - delta velocity 1
        v2 - delta velocity 2
    '''
    
    r1m = distance(r1)
    r2m = distance(r2)
    # angles and convenience parameters
    c_nu = np.dot(r1, r2)/(r1m*r2m)
    s_nu = tm*np.sqrt(1 - c_nu**2.0)
    nu = np.arctan2(s_nu, c_nu)
    if nu < 0:
        nu = nu + (2*np.pi)
    s_nu2 = np.sin(nu/4)**2.0
    c_nu2 = np.cos(nu/4)**2.0
    c = cos_law(r1m, r2m, nu)
    s = (r1m + r2m + c)/2
    eps = (r2m - r1m)/r1m
    t2w = (eps**2.0/4.0)/(np.sqrt(r2m/r1m) + r2m/r1m*(2.0 + np.sqrt(r2m/r1m)))
    rop = np.sqrt(r1m*r2m)*(c_nu2 + t2w)

    if nu < np.pi:
        tp = s_nu2 + t2w
        btm = s_nu2 + t2w + np.cos(nu/2.0)
        l = tp/btm
    else:
        tp = c_nu2 + t2w - np.cos(nu/2.0)
        btm = c_nu2 + t2w
        l = tp/btm

    m = mu*dt**2.0/(8*rop**3)

    if orbit_type == 1:
        x = l
    else:
        x = 0

    err = 10  # error
    tol = 1e-6  # tolerance
    iters = 0
    iters_max = 100
    # coefficients
    c1 = 8.0
    c2 = 1.0
    c3 = 9.0/7.0
    c4 = 16.0/63.0
    c5 = 25.0/99.0
    
    while err > tol:
        eta = x/(np.sqrt(1.0 + x) + 1.0)**2.0
        xi = c1*(np.sqrt(1.0 + x) + 1.0)/(3.0 + c2/(5.0 + eta + c3*eta/
                                         (1.0 + c4*eta/(1.0 + c5*eta))))
        bt1 = (1.0 + 2*x + l)
        bt2 = (4.0*x + xi*(3.0 + x))
        h1 = (l + x)**2.0*(1.0 + 3.0*x + xi)/(bt1*bt2)
        h2 = (m*(x - l + xi))/(bt1*bt2)
        B = 27.0*h2/(4.0*(1.0 + h1)**3.0)
        U = B/(2.0*(np.sqrt(1.0 + B) + 1.0))
        K = (1/3)/(1 + (4/27)*U/(1 + (8/27)*U/(1 + (700/2907)*U)))
        y = (1+h1)/3*(2 + np.sqrt(1 + B)/(1 + 2*U*K**2.0))
        xn = np.sqrt(((1-l)/2)**2.0 + m/y**2) - (1 + l)/2
        err = np.abs(xn - x)
        x = xn
        iters += 1
        if iters > iters_max:
            break
    # semi major axis
    a = mu*(dt**2.0)/(16*rop**2*x*y**2.0)
    if a > 0:
        s_beta = np.sqrt((s - c)/(2*a))
        beta = 2.0*np.arcsin(s_beta)
        if nu > np.pi:
            beta = -beta
        a_min = s/2
        tmin = np.sqrt(a_min**3.0/mu)*(np.pi - beta + np.sin(beta))
        alpha = 2.0*np.arcsin(np.sqrt(s/(2*a)))
        if dt > tmin:
            alpha = 2*np.pi - alpha
        dE = alpha - beta
        f = 1 - a/r1m*(1 - np.cos(dE))
        g = dt - np.sqrt(a**3.0/mu)*(dE - np.sin(dE))
        gdot = 1 - a/r2m*(1 - np.cos(dE))
    else:
        alpha = 2.0*np.arcsinh(np.sqrt(s/(-2*a)))
        beta = 2.0*np.arcsinh(np.sqrt((s-c)/(-2*a)))
        dH = alpha - beta
        f = 1 - a/r1m*(1 - np.cosh(dH))
        g = dt - np.sqrt(-a**3.0/mu)*(np.sinh(dH) - dH)
        gdot = 1 - a/r2m*(1 - np.cosh(dH))

    v1 = (r2 - f*r1)/g
    v2 = (gdot*r2 - r1)/g
    return v1, v2


def epoch(utc):
    '''
    Converts two UTC dates to julian time and calculate the change in epoch
    Args:
        utc - mission start and end dates as two strings in an array
    Returns:
        etOne - start date in julian time
        etTwo - end date in julian time
        dt - change in epoch in seconds
    '''
    # convert to julian time
    etOne = spice.str2et(utc[0])
    etTwo = spice.str2et(utc[1])
    # difference in epoch
    dt = etTwo - etOne
    return etOne, etTwo, dt


def find_state(body, et, observer):
    '''
    Determines the position and velocity vectors
    Args:
        body - desired body
        et - epoch
        observer - observing body
    Returns
        r - position vector
        v - velocity
    '''
    state, lightTimes = spice.spkezr(body, et, 'J2000', 'NONE', observer)
    r = state[:3]
    v = state[3:]
    return r, v        


def unit_vector(r):
    '''
    Converts a position vector to a unit vector
    '''
    return r/np.linalg.norm(r)


def intercept(r1, r2):
    '''
    Determines the spacecraft intercept angle
    Args:
        r1 - unit vector from the primary body to the Sun
        r2 - unit vector of the secondary body relative to the primary body
    Returns
        theta - intercept location with respect to the primary body
    '''
    return np.arccos(np.dot(r1, r2))


def doStuff(t1, dt):
    
    # initialize an empty dataframe
    data = pd.DataFrame([])
    
    # convert to julian seconds
    etOne = spice.str2et(t1)
    etTwo = etOne + dt
    
    # one day increments
    inc = 1*24*3600 
    
    for i in range(int(dt/inc)):
        # position and velocity relative to the Sun
        r_earth,  v_earth  = find_state('Earth',  etOne, 'Sun') # at departure
        r_saturn, v_saturn = find_state('Saturn', etTwo, 'Sun') # at arrival
        r_titan,  v_titan  = find_state('Titan',  etTwo, 'Sun') # at arrival
        
        # position and velocity relative to Saturn
        r_es,  v_es  = find_state('Earth', etOne, 'Saturn') # at departure
        r_ts,  v_ts  = find_state('Titan', etTwo, 'Saturn') # at arrival
        r_sun, v_sun = find_state('Sun',   etTwo, 'Saturn') # at arrival
        
        # unit vectors relative to Saturn
        u_sun   = unit_vector(r_sun)
        u_titan = unit_vector(r_ts)
        
        # find the intercept location at Saturn
        theta = intercept(u_sun, u_titan)
        
        # lambert solver
        v1, v2 = Lambert_Battin(mu_sun, r_earth, r_saturn, dt)        
        # hyperbolic excess velocity relative to Earth at the SOI
        v_inf_e = v1 - v_earth        
        # hyperbolic excess velocity relative to Saturn at the SOI
        v_inf_s = v2 - v_saturn        
        # departure c3
        c3_depart  = v_inf_e**2.0
        # arrival c3
        c3_arrival = v_inf_s**2.0
        etOne += inc

        # Add to the Dataframe
        data = data.append(pd.DataFrame({'Intercept (deg)': theta*(180/np.pi), 
                                         'Earth Departure v1 (km/s)': np.linalg.norm(v1),
                                         'Saturn Arrival v2 (km/s)': np.linalg.norm(v2),
                                         'Departure v inf (km/s)': np.linalg.norm(v_inf_e),
                                         'Arrival v inf (km/s)': np.linalg.norm(v_inf_s),
                                         'Departure c3': np.linalg.norm(c3_depart),
                                         'Arrival c3': np.linalg.norm(c3_arrival)},
                                         index = [0]), ignore_index = True)
    return data
        
if __name__ == '__main__':
    
    # define some constants
    mu_sun = 132712440018.0 # km^3s^2
    mu_earth = 398600.4415 # km^3s^2
    r_earth = 6378.137 # km
    
    # furnsh the SPICE kernels
    furnsh_kernels()
    
    # mission start and end dates
    utc = ['May 20, 2022', 'Dec 1, 2030']
    etOne, etTwo, dt = epoch(utc)
    dt = 6.5*365*24*3600
    
    doStuff(utc[0], dt)
#     # position and velocity relative to the Sun
#     r_earth,  v_earth  = find_state('Earth',  etOne, 'Sun') # at departure
#     r_saturn, v_saturn = find_state('Saturn', etTwo, 'Sun') # at arrival
#     r_titan,  v_titan  = find_state('Titan',  etTwo, 'Sun') # at arrival
    
#     # position and velocity relative to Saturn
#     r_es,  v_es  = find_state('Earth', etOne, 'Saturn') # at departure
#     r_ts,  v_ts  = find_state('Titan', etTwo, 'Saturn') # at arrival
#     r_sun, v_sun = find_state('Sun',   etTwo, 'Saturn') # at arrival
    
#     # unit vectors relative to Saturn
#     u_sun   = unit_vector(r_sun)
#     u_titan = unit_vector(r_ts)
    
#     # find the intercept location at Saturn
#     theta = intercept(u_sun, u_titan)
    
# ########################### INTERPLANETARY TRANSFER ###########################

#     # lambert solver
#     v1, v2 = Lambert_Battin(mu_sun, r_earth, r_saturn, dt)
#     # hyperbolic excess velocity relative to Earth at the SOI
#     v_inf_e = v1 - v_earth
#     # hyperbolic excess velocity relative to Saturn at the SOI
#     v_inf_s = v2 - v_saturn
#     # departure c3
#     c3_depart = v_inf_e**2.0
    
############################ DEPARTURE FROM EARTH #############################
    
    # determine perigee altitude
    # alt = 185 # km
    # # energy equation for hyperbolic trajectory from earth
    # epsilon = (np.linalg.norm(v_inf_e)**2)/2 
    # # required velocity at perigee
    # v_p = np.sqrt(2*(mu_earth/(r_earth + alt) + epsilon))
    # # circular parking orbit prior to interplanetary injection
    # v_circular = np.sqrt(mu_earth/(r_earth + alt))
    # # delta v to inject onto interplanetary orbit
    # delta_v1 = v_p - v_circular
    
    
    
    
    
    
        # fig = plt.figure()
    # ax  = fig.add_subplot(111, projection = '3d')
    # x  = [0, sun_unitV[0]]
    # x1 = [sun_unitV[0], titan_unitV[0]]
    # y  = [0, sun_unitV[1]]
    # y1 = [sun_unitV[1], titan_unitV[1]]
    # z  = [0, sun_unitV[2]]
    # z1 = [sun_unitV[2], titan_unitV[2]]
    # ax.plot(x,y,z)
    # ax.plot(x1,y1,z1)