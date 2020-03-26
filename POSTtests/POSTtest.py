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
    post_data = post_data.append(test_data.loc[[1]])
    post_data = post_data.append(test_data[-1:])
    
    # load in the final orbit data
    filepath2 = r'C:\Spice_Kernels'
    filename2 = r'\3D_potato'
    script_data = pd.read_hdf(filepath2 + filename2 + r'.hdf')

    # pull out the columns in the data from script that can match the data from post
    v1 = script_data['Titan v1']
    fpa = script_data['Titan fpa (deg)']

    # pull out entry and exit values from the post data
    v1_post = post_data['veli']
    v1_en = v1_post.iloc[0] / 1000
    v1_ex = v1_post.iloc[1] / 1000
    e_post = post_data['eccen']
    e_en = e_post.iloc[0]
    e_ex = e_post.iloc[1]
    la_post = post_data['longi']
    la_en = la_post.iloc[0]
    la_ex = la_post.iloc[1]
    ta_post = post_data['truan']
    ta_en = m.radians(ta_post.iloc[0])
    ta_ex = m.radians(ta_post.iloc[1])
    lat_post = post_data['gclat']
    lat_en = lat_post.iloc[0]
    lat_ex = lat_post.iloc[1]

    # calculate the things
    phi1 = lat_en
    phi2 = lat_ex
    delta_lambda = m.fabs(la_ex - la_en)
    delta_sigma = m.acos(m.sin(phi1) * m.sin(phi2) + m.cos(phi1) * m.cos(phi2) * m.cos(delta_lambda))
    ta_inf_i = m.acos(-1/e_en)
    ta_inf_o = m.acos(-1/e_ex)
    turn_angle = m.fabs(ta_inf_i - ta_en) + m.fabs(ta_inf_o - ta_ex) + delta_sigma - 2*m.pi

    # compare the output of POST to the output of the script
    v_count = 0
    fpa_count = 0

    # set number of decimal places to truncate vales to for comparison
    j = 3
    v1_ex = chop(v1_ex, j)
    turn_angle = chop(turn_angle, j)
    tot = 0
    min = 1000
    max = 0

    for i in range(1, script_data.shape[0]):

        # check velocities
        if v1_ex == chop(v1.loc[i], j):
            v_count += 1
            print("velocity", i)

        # check flight path angle and turn angle
        if turn_angle == m.radians(chop(fpa.loc[i], j)):
            fpa_count += 1
            print("turn/flight path angle", i)
        else:
            diff = m.fabs(turn_angle - m.radians(chop(fpa.loc[i], j)))
            tot += diff
            if diff > max:
                max = diff
            elif diff < min:
                min = diff
            # print(turn_angle, m.radians(chop(fpa.loc[i], j)))
            # print(diff)

    # find out some stuff about the difference between turn angle and fpa to figure out what's going wrong
        # average difference
        # minimum difference
        # maximum difference
    avg = tot / (i + 1)
    print(avg)
    print(min)
    print(max)

    print(v_count)
    print(fpa_count)
