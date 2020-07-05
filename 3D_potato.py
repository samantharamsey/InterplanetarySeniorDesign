# -*- coding: utf-8 -*-
'''
Created on Mon Dec 9 12:45:36 2019
@author: Samantha Ramsey

Calculates the family of acceptable velocity vectors post Titan aerogravity 
assist that will result in a science orbit about Saturn that will allow for 
frequent fly-bys of Enceladus.

Requires:
    C:\Senior_Design\TitanAGAMission\Data\PotatoData\3D_potato.hdf
'''

# imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



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


def COE2RV(mu, p, e, i, Omega, omega, true):
    '''
    Converts the classical orbital elements to a state vector
    *** does not include special cases ***
    Args:
        mu - gravitational parameter of primary body
        p - semiparameter
        e - eccentricity
        i - inclination in radians
        Omega - longitude of the ascending node in radians
        omega - argument of perigee in radians
        true - true anomaly in radians
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
    pqw2ijk = np.matrix([[np.cos(Omega)*np.cos(omega) - np.sin(Omega)*np.sin(omega)*np.cos(i),
                         -np.cos(Omega)*np.sin(omega) - np.sin(Omega)*np.cos(omega)*np.cos(i),
                          np.sin(Omega)*np.sin(i)],
                         [np.sin(Omega)*np.cos(omega) + np.cos(Omega)*np.sin(omega)*np.cos(i),
                         -np.sin(Omega)*np.sin(omega) + np.cos(Omega)*np.cos(omega)*np.cos(i),
                         -np.cos(Omega)*np.sin(i)],
                         [np.sin(omega)*np.sin(i),
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
    i =   inc*(np.pi/180)
    state = np.concatenate((r, [v*np.cos(i)*np.sin(g),
                                v*np.cos(i)*np.cos(g),
                                v*np.sin(i)]))
    return state


def calc_new_v(gamma, inclination):
    '''
    Calculates necessary v1
    Args:
        gamma - flight path angle array
        inclination - inclination angle array
    Returns:
        v1 - post aerocapture velocity that results in satisfactory orbit
        orbital elements
    '''
    
    # constants
    escape_v = 7.94 # km/s
    increment = 0.01 # km/s
    
    # initialize an empty dataframe
    data = pd.DataFrame([])
    
    for g in gamma:
        for i in inclination:
            v1 = escape_v # max v1 to prevent escape from Saturn system in km/s
            r = r_encel + 10
            
            while r > r_encel:
                v1 = v1 - increment # max v1 to prevent escape from Saturn system in km/s
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


def calculation_loop(inc, gamma):
    '''
    Evaluates based on quadrant
    *** BROKEN ***
    Args:
        inclination - array of inclinations
        gamma - array of flight path angles
    '''

    for g in gamma:
        for i in inc:

            if 0 < g < 90 or 270 < g < 360:
                if 0 < i < 180:
                    v1, r, h, E, n, e, p, a, incl, Omega, omega, nu = calc_new_v(7.94, -0.01, gamma, inc)
                else:
                    v1, r, h, E, n, e, p, a, incl, Omega, omega, nu = calc_new_v(-7.94, 0.01, gamma, inc)
            elif 90 < g < 270:
                if 0 < i < 180:
                    v1, r, h, E, n, e, p, a, incl, Omega, omega, nu = calc_new_v(7.94, -0.01, gamma, inc)
                else:
                    v1, r, h, E, n, e, p, a, incl, Omega, omega, nu = calc_new_v(-7.94, 0.01, gamma, inc)
            else:
                v1, r, h, E, n, e, p, a, incl, Omega, omega, nu = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
            
        return v1, r, h, E, n, e, p, a, incl, Omega, omega, nu
    
    
def plot_single_inclination(inclination, gamma):
    '''
    Plots a single degree of inclination for a range of flight path angles
    *** Use this to see a single inclination curve ***
    Args:
        inclination - the inclination you want to plot in degrees
                      allowable range: 0 < inclination < 360
        gamma - the final flight path angle to plot in degrees
                this function will plot from fpa = 0 - fpa = gamma
                allowable range: 0 < gamma < 360
    '''

    fpa   = data['Saturn fpa (deg)'] < gamma
    lower = data['Saturn inclination (deg)'] < inclination + 1
    upper = data['Saturn inclination (deg)'] > inclination - 1
    new_data = data[fpa]
    new_data = new_data[lower]
    new_data = new_data[upper]
    
    # determine the plot color based on the perigee radius
    green  = new_data['Radius at Descending Node'] > 230000
    orange = new_data['Radius at Descending Node'].between(61000, 230000)
    red    = new_data['Radius at Descending Node'] < 61000
    
    # create a mini dataframe of each color
    data3 = new_data[green]
    data4 = new_data[orange]
    data5 = new_data[red]
    
    # convert to Cartesian coordinates
    X_green  = data3['Saturn v1 max']*np.cos(
               data3['Saturn inclination (deg)']*(np.pi/180))*np.sin(
               data3['Saturn fpa (deg)']*(np.pi/180))
    Y_green  = data3['Saturn v1 max']*np.cos(
               data3['Saturn inclination (deg)']*(np.pi/180))*np.cos(
               data3['Saturn fpa (deg)']*(np.pi/180))
    Z_green  = data3['Saturn v1 max']*np.sin(
               data3['Saturn inclination (deg)']*(np.pi/180))
    X_orange = data4['Saturn v1 max']*np.cos(
               data4['Saturn inclination (deg)']*(np.pi/180))*np.sin(
               data4['Saturn fpa (deg)']*(np.pi / 180))
    Y_orange = data4['Saturn v1 max']*np.cos(
               data4['Saturn inclination (deg)']*(np.pi/180))*np.cos(
               data4['Saturn fpa (deg)']*(np.pi/180))
    Z_orange = data4['Saturn v1 max']*np.sin(
               data4['Saturn inclination (deg)']*(np.pi/180))
    X_red    = data5['Saturn v1 max']*np.cos(
               data5['Saturn inclination (deg)']*(np.pi/180))*np.sin(
               data5['Saturn fpa (deg)']*(np.pi/180))
    Y_red    = data5['Saturn v1 max']*np.cos(
               data5['Saturn inclination (deg)']*(np.pi/180))*np.cos(
               data5['Saturn fpa (deg)']*(np.pi/180))
    Z_red    = data5['Saturn v1 max']*np.sin(
               data5['Saturn inclination (deg)']*(np.pi/180))
    
    # do the plotting
    plt.rcParams.update({'font.size': 10})
    plt.rcParams['font.family'] = 'times new roman'
    fig = plt.figure()
    ax  = fig.add_subplot(111, projection = '3d')
    ax.scatter(X_green, Y_green, Z_green, c = 'green')
    ax.scatter(X_orange, Y_orange, Z_orange, c = 'orange')
    ax.scatter(X_red, Y_red, Z_red, c = 'red')
    ax.set_xlabel('Radial Velocity Component (km/s)')
    ax.set_ylabel('Tangential Velocity Component (km/s)')
    ax.set_zlabel('Vertical Velocity Component (km/s)')
    plt.legend(['Trajectories with $r_p$ about equal to $r_{enceladus}$',
                'Escape Velocity Constraint',
                'Trajectories with Perigee Radius Smaller than Saturn\'s Radius',], loc = 3)
    plt.show()


def plot_multiple_inclinations(inclination, gamma):
    '''
    Plots a selection of inclinations for a range of flight path angles
    *** Use this to see a small selection of inclination curves ***
    Args:
        inclination - an [array] of inclinations you want to plot in degrees
                      0 < inclination < 360
        gamma - the final flight path angle to plot in degrees
                as a single number; 0 < gamma < 360
    '''
    
    plt.rcParams.update({'font.size': 12})
    plt.rcParams['font.family'] = 'times new roman'

    fig = plt.figure()
    ax  = fig.add_subplot(111, projection='3d')
    
    for i in inclination:
        fpa   = data['Saturn fpa (deg)'] < gamma
        lower = data['Saturn inclination (deg)'] < i + 1
        upper = data['Saturn inclination (deg)'] > i - 1
        new_data = data[fpa]
        new_data = new_data[lower]
        new_data = new_data[upper]
        
        # determine the plot color based on the perigee radius
        green  = new_data['Radius at Descending Node'] > 230000
        orange = new_data['Radius at Descending Node'].between(61000, 230000)
        red    = new_data['Radius at Descending Node'] < 61000
        
        # create a mini dataframe of each color
        data3 = new_data[green]
        data4 = new_data[orange]
        data5 = new_data[red]
        
        # convert to Cartesian coordinates
        X_green  = data3['Saturn v1 max']*np.cos(
                   data3['Saturn inclination (deg)']*(np.pi/180))*np.sin(
                   data3['Saturn fpa (deg)']*(np.pi/180))
        Y_green  = data3['Saturn v1 max']*np.cos(
                   data3['Saturn inclination (deg)']*(np.pi/180))*np.cos(
                   data3['Saturn fpa (deg)']*(np.pi/180))
        Z_green  = data3['Saturn v1 max']*np.sin(
                   data3['Saturn inclination (deg)']*(np.pi/180))
        X_orange = data4['Saturn v1 max']*np.cos(
                   data4['Saturn inclination (deg)']*(np.pi/180))*np.sin(
                   data4['Saturn fpa (deg)']*(np.pi / 180))
        Y_orange = data4['Saturn v1 max']*np.cos(
                   data4['Saturn inclination (deg)']*(np.pi/180))*np.cos(
                   data4['Saturn fpa (deg)']*(np.pi/180))
        Z_orange = data4['Saturn v1 max']*np.sin(
                   data4['Saturn inclination (deg)']*(np.pi/180))
        X_red    = data5['Saturn v1 max']*np.cos(
                   data5['Saturn inclination (deg)']*(np.pi/180))*np.sin(
                   data5['Saturn fpa (deg)']*(np.pi/180))
        Y_red    = data5['Saturn v1 max']*np.cos(
                   data5['Saturn inclination (deg)']*(np.pi/180))*np.cos(
                   data5['Saturn fpa (deg)']*(np.pi/180))
        Z_red    = data5['Saturn v1 max']*np.sin(
                   data5['Saturn inclination (deg)']*(np.pi/180))
        
        # do the plotting 
        ax.scatter(X_green, Y_green, Z_green, c = 'green')
        ax.scatter(X_orange, Y_orange, Z_orange, c = 'orange')
        ax.scatter(X_red, Y_red, Z_red, c = 'red')
        
    ax.set_xlabel('Radial Velocity Component (km/s)')
    ax.set_ylabel('Tangential Velocity Component (km/s)')
    ax.set_zlabel('Vertical Velocity Component (km/s)')
    plt.title('Family of Acceptable Post-AGA Velocity Vectors wrt Saturn')
    plt.legend(['Trajectories with $r_p$ about equal to $r_{enceladus}$',
                'Escape Velocity Constraint',
                'Trajectories with Perigee Radius Smaller than Saturn\'s Radius',], loc = 3)
    plt.show()


def plot_inclination_range(inclination, gamma):
    '''
    Plots a range of inclinations for a range of flight path angles
    *** Use this to go from 0 up to the desired values ***
    Args:
        inclination - the final inclination you want to plot in degrees
                      as a single number; 0 < inclination < 360
        gamma - the final flight path angle to plot in degrees
                as a single number; 0 < gamma < 360   
    '''
    
    fpa = data['Saturn fpa (deg)'] < gamma
    inc = data['Saturn inclination (deg)'] < inclination
    new_data = data[fpa]
    new_data = new_data[inc]
    
    # determine the plot color based on the perigee radius
    green  = new_data['Radius at Descending Node'] > 230000
    orange = new_data['Radius at Descending Node'].between(61000, 230000)
    red    = new_data['Radius at Descending Node'] < 61000
    
    # create a mini dataframe of each color
    data3 = new_data[green]
    data4 = new_data[orange]
    data5 = new_data[red]
    
    # convert to Cartesian coordinates
    X_green  = data3['Saturn v1 max']*np.cos(
               data3['Saturn inclination (deg)']*(np.pi/180))*np.sin(
               data3['Saturn fpa (deg)']*(np.pi/180))
    Y_green  = data3['Saturn v1 max']*np.cos(
               data3['Saturn inclination (deg)']*(np.pi/180))*np.cos(
               data3['Saturn fpa (deg)']*(np.pi/180))
    Z_green  = data3['Saturn v1 max']*np.sin(
               data3['Saturn inclination (deg)']*(np.pi/180))
    X_orange = data4['Saturn v1 max']*np.cos(
               data4['Saturn inclination (deg)']*(np.pi/180))*np.sin(
               data4['Saturn fpa (deg)']*(np.pi / 180))
    Y_orange = data4['Saturn v1 max']*np.cos(
               data4['Saturn inclination (deg)']*(np.pi/180))*np.cos(
               data4['Saturn fpa (deg)']*(np.pi/180))
    Z_orange = data4['Saturn v1 max']*np.sin(
               data4['Saturn inclination (deg)']*(np.pi/180))
    X_red    = data5['Saturn v1 max']*np.cos(
               data5['Saturn inclination (deg)']*(np.pi/180))*np.sin(
               data5['Saturn fpa (deg)']*(np.pi/180))
    Y_red    = data5['Saturn v1 max']*np.cos(
               data5['Saturn inclination (deg)']*(np.pi/180))*np.cos(
               data5['Saturn fpa (deg)']*(np.pi/180))
    Z_red    = data5['Saturn v1 max']*np.sin(
               data5['Saturn inclination (deg)']*(np.pi/180))
    
    # do the plotting
    plt.rcParams.update({'font.size': 10})
    plt.rcParams['font.family'] = 'times new roman'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.scatter(X_green, Y_green, Z_green, c = 'green')
    ax.scatter(X_orange, Y_orange, Z_orange, c = 'orange')
    ax.scatter(X_red, Y_red, Z_red, c = 'red')
    ax.set_xlabel('Radial Velocity Component (km/s)')
    ax.set_ylabel('Tangential Velocity Component (km/s)')
    ax.set_zlabel('Vertical Velocity Component (km/s)')
    plt.title('Family of Acceptable Post-AGA Velocity Vectors wrt Saturn')
    plt.legend(['Trajectories with $r_p$ about equal to $r_{enceladus}$',
                'Escape Velocity Constraint',
                'Trajectories with Perigee Radius Smaller than Saturn\'s Radius'], loc = 8)
    plt.show()


def referenceframe_transformation(data):
    '''
    Converts from Saturn centered reference frame to Titan centered
    *** Appends to DataFrame - will not overwrite columns ***
    '''

    # wrt Saturn
    vx_max = data['Saturn v1 max']*np.cos(
             data['Saturn inclination (deg)']*(np.pi/180))*np.sin(
             data['Saturn fpa (deg)']*(np.pi/180))
    vy_max = data['Saturn v1 max']*np.cos(
             data['Saturn inclination (deg)']*(np.pi/180))*np.cos(
             data['Saturn fpa (deg)']*(np.pi/180))
    vz_max = data['Saturn v1 max']*np.sin(
             data['Saturn inclination (deg)']*(np.pi/180))

    # wrt Titan
    new_vymax = vy_max - v_titan

    vmax = []
    for i in range(len(data)):
        vmax_titan = np.linalg.norm([vx_max[i], new_vymax[i], vz_max[i]])
        vmax.append(vmax_titan)
    titan_fpa_max = np.arctan2(vx_max, new_vymax)

    data.insert(12, 'Titan fpa (deg)', titan_fpa_max*(180/np.pi), True)
    data.insert(13, 'Titan v1', vmax, True)
    data.to_hdf(filepath + filename + r'.hdf', key = 'df')


if __name__ == '__main__':
    
    ############################## CALCULATIONS ##############################

    # load pre-existing file into the dataframe
    # UPDATE FILEPATH BEFORE RUNNING
    filepath = r'C:\Senior_Design\TitanAGAMission\Data\PotatoData'
    filename = r'\3D_potato'
    data = pd.read_hdf(filepath + filename + r'.hdf')    
    
    # define some constants
    mu = 37.931*10**6  # saturn's gravitational parameter
    saturn_equatorial = 60268  # km
    saturn_polar = 54364  # km
    r_titan = 1.2*10**6  # km
    r_encel = 238000  # km
    v_titan = 5.57  # km/s

    # flight path angle array
    gamma = np.linspace(1, 360, 360)
    # inclination angle array
    incl = np.linspace(1, 360, 360)
    
    # calculate stuff if data not present
    # calc_new_v(gamma, incl)
    
    ################################ PLOTTING ################################

    # NOTES:
    # first input for each function is inclination second is flight path angle
    # both can range from 0 - 360 degrees
    # each function call will generate a new plot

    # plot a single inclination curve
    plot_single_inclination(45, 360)

    # plot the inclination curves in multiples of 10
    inclination = [1, 10, 20, 30, 40, 50, 60, 70, 80, 89]
    plot_multiple_inclinations(inclination, 360)

    # plot the entire first quadrant
    plot_inclination_range(90, 90)
