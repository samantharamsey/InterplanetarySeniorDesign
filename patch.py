# -*- coding: utf-8 -*-
'''
Created on Sun May  3 12:22:44 2020

@author: sam
'''


import pandas as pd
import numpy as np


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
    

def referenceframe_transformation(data):
    '''
    Converts from Saturn centered reference frame to Titan centered
    '''
    
    # wrt Saturn
    vx_max = data['Saturn v1 max']*np.cos(data['Saturn inclination (deg)']*
                 (np.pi/180))*np.sin(data['Saturn fpa (deg)']*(np.pi/180))
    vy_max = data['Saturn v1 max']*np.cos(data['Saturn inclination (deg)']*
                 (np.pi/180))*np.cos(data['Saturn fpa (deg)']*(np.pi/180))
    vz_max = data['Saturn v1 max']*np.sin(data['Saturn inclination (deg)']*
                 (np.pi/180))
    
    # wrt Titan
    new_vymax = vy_max - v_titan
    
    # Create new Dataframe
    data = pd.DataFrame({'vxi': vx_max,
                         'vyi': new_vymax,
                         'vzi': vz_max})                      
    script_file2 = r'\3D_potato_extended'
    data.to_hdf(filepath + script_file2 + r'.hdf', key = 'df')  
    
    
def get_comp_data(post_data, script_data, decimal):
    '''
    Loads POST and script data for the purpose of comparison
    Args:
        post_data - dataframe of POST file
        script_data - dataframe of script file
        decimal - number of decimals to include
    Returns:
        post_comp - only final xyz velocity components converted to km/s
        script_comp - only final xyz velocity components in km/s
    '''
    
    post_comp = pd.concat([post_data.iloc[-1,7::]/1000], axis = 1).T
    post_comp = post_comp.round(decimals = decimal)
    referenceframe_transformation(script_data)
    script_comp = pd.read_hdf(filepath + script_file2)
    script_comp = script_comp.round(decimals = decimal)
    
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
    
    post_comp, script_comp = get_comp_data(post_data, script_data, 2)
    result = script_comp.merge(post_comp, how='inner')
