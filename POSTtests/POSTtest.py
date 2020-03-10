# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:45:09 2020

@author: sam, hannah
"""
import pandas as pd
import math as m


def chop(num, digits) -> float:
    step = 10.0 ** digits
    return m.trunc(step*num) / step


if __name__ == '__main__':
    
    # load in the data from POST
    filepath1 = r'C:\Spice_Kernels'
    filename1 = r'\tin7'
    test_data = pd.read_csv(filepath1 + filename1 + r'.hdf', delim_whitespace=True)
    # create a new DF with only the last lines
    post_data = pd.DataFrame([])
    post_data = post_data.append(test_data[-1:])
    
    # load in the final orbit data
    filepath2 = r'C:\Spice_Kernels'
    filename2 = r'\3D_potato'
    script_data = pd.read_hdf(filepath2 + filename2 + r'.hdf')

    # pull out the columns in the data from script that can match the data from post
    v1 = script_data['Titan v1']
    e = script_data['Eccentricity']
    la = script_data['Longitude of Ascending Node']
    ta = script_data['True Anomaly']

    # pull out values from the post data that we'll actually use
    v1_post = post_data['veli']
    v1_post = v1_post.iloc[0] / 1000
    e_post = post_data['eccen']
    e_post = e_post.iloc[0]
    la_post = post_data['longi']
    la_post = la_post.iloc[0]
    ta_post = post_data['truan']
    ta_post = ta_post.iloc[0]

    # compare the output of POST to the output of the script
    listy = []
    v_count = 0
    e_count = 0
    la_count = 0
    ta_count = 0
    for i in range(1, script_data.shape[0]):
        # set number of decimal places to truncate vales to for comparison
        j = 2

        # check velocities
        if chop(v1_post, j) == chop(v1.loc[i], j):
            v_count += 1
            print("velocity", i)

        # check eccentricities
        if chop(e_post, j) == chop(e.loc[i], j):
            e_count += 1
            print("eccentricity", i)

        # check longitudes of ascending node
        if chop(la_post, j) == chop(la.loc[i], j):
            la_count += 1
            print("longitude of ascending node", i)

        # check longitudes of ascending node
        if chop(ta_post, j) == chop(ta.loc[i], j):
            ta_count += 1
            print("true anomaly", i)

    print(v_count)
    print(e_count)
    print(la_count)
    print(ta_count)
