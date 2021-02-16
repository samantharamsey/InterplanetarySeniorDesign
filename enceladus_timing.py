# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 11:44:44 2021

@author: sam
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def period(a, mu):
    ''' calculates period of the orbit '''
    return 2*np.pi*np.sqrt((a**3)/mu)


if __name__ == '__main__':
    
    # define some constants
    mu = 37.931*10**6  # saturn's gravitational parameter
    saturn_equatorial = 60268  # km
    saturn_polar = 54364  # km
    r_titan = 1.2*10**6  # km
    r_encel = 238000  # km
    v_titan = 5.57  # km/s
    
    # load pre-existing file into the dataframe
    filepath = r'C:\Senior_Design\TitanAGAMission\Data\PotatoData'
    filename = r'\3D_potato_small'
    data = pd.read_hdf(filepath + filename + r'.hdf')  
    
    T = period(data['Semimajor Axis'], mu)
    data.insert(12, 'Period (Saturn)', T, True)
    
    # remove bugged data
    data.drop(data.loc[data['Saturn inclination (deg)']==180].index, inplace=True)
    data.drop(data.loc[data['Saturn fpa (deg)']== 90].index, inplace=True)
    data.drop(data.loc[data['Saturn fpa (deg)']==180].index, inplace=True)
    data.drop(data.loc[data['Saturn fpa (deg)']==270].index, inplace=True)
    data.drop(data.loc[data['Saturn fpa (deg)']==360].index, inplace=True)
    data.drop(data.loc[data['Saturn v1 max'] > 7].index, inplace=True)
    
    index = []
    for i in range(len(data)):
        index.append(i)
    
    plt.scatter(data['Saturn fpa (deg)'], data['Period (Saturn)'])
            