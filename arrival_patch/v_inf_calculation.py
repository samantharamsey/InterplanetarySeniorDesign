# -*- coding: utf-8 -*-
'''
Created on Sat Apr 11 20:21:52 2020

@author: sam
'''

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def calc_vinf_titan(v_inf, dec, intercept):
    '''
    calculates v_infinity wrt Titan
    Assumptions:
        coordinate system - 0 deg is at W
    Args:
        v_inf - v_infinity wrt Saturn im km/s
        dec - arrival declination in deg
        intercept - intercept position in deg
    Returns:
        v_inf_mag - v_infinity wrt Titan magnitude
    '''
    mu  = 3.795*10**7
    r = 1.22*10**6
    v_titan = 5.57
    dec = dec*(np.pi/180)
    inter = intercept*(np.pi/180)
    
    # calculate probe velocity at titans orbital distance
    v1 = np.sqrt(2*((v_inf**2/2) + (mu/r)))
    
    # separate into components
    vx1 = 0
    vy1 = v1*np.cos(dec)
    vz1 = v1*np.sin(dec)
    
    # titans velocity components based on quadrant
    if 0 < intercept < 90:
        vxt = v_titan*np.sin(inter)
        vyt = v_titan*np.cos(inter)
        vzt = 0
    elif intercept == 90:
        vxt = v_titan
        vyt = 0
        vzt = 0 
    elif 90 < intercept < 180:
        phi = 180 - inter
        vxt = v_titan*np.sin(phi)
        vyt = -v_titan*np.cos(phi)
        vzt = 0
    elif intercept == 180:
        vxt = 0
        vyt = -v_titan
        vzt = 0
    elif 180 < intercept < 270:
        phi = 270 - inter
        vxt = -v_titan*np.cos(phi)
        vyt = -v_titan*np.sin(phi)
        vzt = 0
    elif intercept == 270:
        vxt = -v_titan
        vyt = 0
        vzt = 0 
    elif 270 < intercept < 360:
        phi = 360 - inter
        vxt = -v_titan*np.sin(phi)
        vyt = v_titan*np.cos(phi)
        vzt = 0
    elif intercept == 360 or intercept == 0:
        vxt = 0
        vyt = v_titan
        vzt = 0
    else:
        print('error: intercept location')
        breakpoint()
        
    # find v_infinity wrt Titian
    v_infx = vx1 + vxt
    v_infy = vy1 + vyt
    v_infz = vz1 + vzt
    
    # v_infinity magnitude
    v_inf_mag = np.sqrt(v_infx**2 + v_infy**2 + v_infz**2)
    
    return v_infx, v_infy, v_infz, v_inf_mag
        
    

if __name__ == '__main__':
    
    filepath = r'C:\Users\saman\OneDrive\Desktop\InterplanetarySeniorDesign\arrival_patch'
    filename = r'\10_declination'
    data = pd.DataFrame([])
    v_infx = []
    v_infy = []
    v_infz = []
    for i in range(360):
        resultx, resulty, resultz, resultmag = calc_vinf_titan(6, 10, i)
        v_infx.append(resultx)
        v_infy.append(resulty)
        v_infz.append(resultz)
       
        data = data.append(pd.DataFrame({'intercept location (deg)': i,
                                         'v_infinity wrt Titian magnitude (km/s)': resultmag,
                                         'v_infinity x component': resultx,
                                         'v_infinity y component': resulty,
                                         'v_infinity z component': resultz},
                                         index=[0]), ignore_index=True)
    # send results to excel
    data.to_csv(filepath + filename + r'.csv', index = False)
    # send results to HDF5 - faster loading
    data.to_hdf(filepath + filename + r'.hdf', key = 'df')
    
    # plot stuff
    fig1 = plt.figure()
    plt.plot(v_infx)
    plt.plot(v_infy)
    plt.plot(v_infz)
    plt.ylabel('velocity (km/s)')
    plt.xlabel('intercept location (deg)')
    plt.legend(['x component', 'y component', 'z component'])
    plt.title('V infinity wrt Titan: v_inf wrt Saturn = 6, dec = 10')
    