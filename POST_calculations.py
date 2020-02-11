# -*- coding: utf-8 -*-
'''
Created on Mon Feb 10 11:29:33 2020
@author: sam
'''


import numpy as np


def calculate_stuff(vsc_inf, dec, th_titan):
    '''
    Converts from incoming relative to Saturn to Titan at encounter altitude
    Args:
        vsc_inf - v infinity of the spacecraft
        dec - inbound declination in degrees
        th_titan - titan intercept location in degrees
    Returns:
        declination - relative to titan at the encounter altitude
        v_inf - relative to titan at the encounter altitude
    '''
    
    # spacecraft velocity at arrival to Titan
    vsc_1 = np.sqrt(2*((vsc_inf**2/2) + mu_saturn/r_titan))
    # velocity of spacecraft in cartesian coordinates
    vsc_x = 0
    vsc_y = vsc_1*np.cos(dec*(np.pi/180))
    vsc_z = vsc_1*np.sin(dec*(np.pi/180))
    # velocity of titan in cartesian coordinates
    v_titanx = v_titan*np.sin(th_titan*(np.pi/180))
    v_titany = -v_titan*np.cos(th_titan*(np.pi/180))
    v_titanz = 0
    # relative velocity between the spacecraft and titan
    v_relx = vsc_x - v_titanx
    v_rely = vsc_y - v_titany
    v_relz = vsc_z - v_titanz
    # convert back to spherical coordinates
    v_mag = np.sqrt(v_relx**2 + v_rely**2 + v_relz**2)
    declination = (180/np.pi)*np.arctan(v_relz/np.sqrt(v_relx**2 + v_rely**2))
    # v spacecraft - titan at intercept
    v_inf = np.sqrt(2*((v_mag**2/2) + mu_titan/r_encounter))
    
    return declination, v_inf
    
  
if __name__ == '__main__':
    
    # define some constants
    v_titan = 5.57 #km/s
    r_titan = 1221865 #km
    mu_saturn = 37931187.9
    mu_titan = 0.0225*(3.986*10**5)
    # radius from Titan encounter
    r_encounter = 49000 #km
    
    # test case
    vsc_inf = 8.32
    dec = 2.5
    th_titan = 180
    declination, v_inf = calculate_stuff(8.32, 2.5, 180)
    print('Reference Frame Transformation \
          \nRelative to Saturn: \
          \nv_inf: {} km/s \
          \ndeclination: {} deg \
          \ntheta: {} deg \
          \n \
          \nRelative to Titan: \
          \nv_inf: {} km/s\
          \ndeclination: {} deg\
          \ntheta: {} deg \
          '.format(vsc_inf, dec, th_titan, v_inf, declination, th_titan))
