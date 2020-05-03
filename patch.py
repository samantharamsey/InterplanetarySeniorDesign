# -*- coding: utf-8 -*-
'''
Created on Sun May  3 12:22:44 2020

@author: sam
'''


import pandas as pd
import numpy as np
import math as m


def load_data(filepath, post, script):
    '''
    loads data from POST and script into dataframes
    Args:
        fileparth - path to files
        post - post data filename and extension
        script - script data filename and extension
    '''
    
    post_data = pd.read_excel(filepath + post)
    script_data = pd.read_hdf(filepath + script)
    
    return post_data, script_data
    

def referenceframe_transformation(data, intercept):
    '''
    Converts from Saturn centered reference frame to Titan centered
    Args:
        data - dataframe
        intercept - space craft intercept location in degrees
    '''
    
    # wrt Saturn
    vxs = data['Saturn v1 max']*np.cos(data['Saturn inclination (deg)']*
                 (np.pi/180))*np.sin(data['Saturn fpa (deg)']*(np.pi/180))
    vys = data['Saturn v1 max']*np.cos(data['Saturn inclination (deg)']*
                 (np.pi/180))*np.cos(data['Saturn fpa (deg)']*(np.pi/180))
    vzs = data['Saturn v1 max']*np.sin(data['Saturn inclination (deg)']*
                 (np.pi/180))
    
    # wrt Titan
    vxt = 0
    vyt = 0
    vzt = 0
    
    if intercept == 90:
        vxt = v_titan
    elif intercept == 180:
        vyt = -v_titan
    elif intercept == 270:
        vxt = -v_titan
    elif intercept == 360 or intercept == 0:
        vyt = v_titan
    else:
        vxt = v_titan*np.sin(intercept*(np.pi/180))
        vyt = v_titan*np.cos(intercept*(np.pi/180))

    vx = vxs + vxt
    vy = vys + vyt
    vz = vzs + vzt
    vmag = []
    for i in range(len(vx)):
        norm = np.linalg.norm([vx[i], vy[i], vz[i]])
        vmag.append(norm)
    
    # Create new Dataframe
    data = pd.DataFrame({'vxi': vx,
                         'vyi': vy,
                         'vzi': vz,
                         'vmag': vmag})       
               
    script_file2 = r'\3D_potato_extended'
    data.to_hdf(filepath + script_file2 + r'.hdf', key = 'df')  
    
    
def get_comp_data(post_data, script_data, intercept):
    '''
    Loads POST and script data for the purpose of comparison
    Args:
        post_data - dataframe of POST file
        script_data - dataframe of script file
        decimal - number of decimals to include
        intercept - intercept location of probe at Titan
    Returns:
        post_comp - only final xyz velocity components converted to km/s
        script_comp - only final xyz velocity components in km/s
    '''
    
    post_comp = pd.concat([post_data.iloc[-1,7::]/1000], axis = 1).T
    vmag = np.linalg.norm(post_comp[['vxi', 'vyi', 'vzi']].values,axis=1)
    post_comp = post_comp.append(pd.DataFrame({'vmag': vmag}))
    referenceframe_transformation(script_data, intercept)
    script_comp = pd.read_hdf(filepath + script_file2)
    
    return post_comp, script_comp
    
    
if __name__ == '__main__':
    
    # some constants
    mu = 37.931*10**6  #saturn mu
    r_titan = 1.2*10**6  #km
    r_encel = 238000  #km
    v_titan = 5.57  #km/s

    filepath = r'C:\Spice_Kernels'
    post_file = r'\tin7_2.xlsx'
    script_file = r'\3D_potato.hdf'
    script_file2 = r'\3D_potato_extended.hdf'
    post_data, script_data = load_data(filepath, post_file, script_file)
    
    post_comp, script_comp = get_comp_data(post_data, script_data, 60)
    result = script_comp.merge(post_comp, how='inner')
    
    
    vmag_post = post_comp['vmag'][0]
    vmag_script = script_comp['vmag']
    vx_post = post_comp['vxi'][2479]
    vy_post = post_comp['vyi'][2479]
    vz_post = post_comp['vzi'][2479]
    vx_script = script_comp['vxi']
    vy_script = script_comp['vyi']
    vz_script = script_comp['vzi']
    
    # compare the output of POST to the output of the script
    v_mag_count = 0
    v_comp_count = 0
    x_comp_count = 0
    y_comp_count = 0
    z_comp_count = 0
    good_indexes = []
    tol = .14

    good = pd.DataFrame()
    
    for i in range(len(vmag_script)):
    # check velocities after truncating decimal places for the script output
        if m.fabs(vmag_post - vmag_script[i]) < tol:
            v_mag_count += 1

            if m.fabs(vx_script[i] - vx_post) < tol and m.fabs(vy_script[i] - vy_post) < tol and m.fabs(vz_script[i] - vz_post) < tol:
                v_comp_count += 1
                good_indexes.append(i)
    
                good = good.append(pd.DataFrame({'vx_script': vx_script[i],
                             'vy_script': vy_script[i],
                             'vz_script': vz_script[i],
                             'vx_post': vx_post,
                             'vy_post': vy_post,
                             'vz_post': vz_post},
                                        index=[0]), ignore_index=True)
            if m.fabs(vx_script[i] - vx_post) < tol:
                x_comp_count += 1
            if m.fabs(vy_script[i] - vy_post) < tol:
                y_comp_count += 1
            if m.fabs(vz_script[i] - vz_post) < tol:
                z_comp_count += 1
                
    print(good)
    print("number of cases with matching velocity magnitudes:", v_mag_count)
    print("number of cases with all matching velocity components:", v_comp_count)
    print("number of cases with matching x velocity components:", x_comp_count)
    print("number of cases with matching y velocity components:", y_comp_count)
    print("number of cases with matching z velocity components:", z_comp_count)
    print(good_indexes)