# -*- coding: utf-8 -*-
'''
Created on Tue Mar 10 17:23:35 2020

@author: sam
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    
    # load in the data from POST - cl off
    filepath = r'C:\Users\saman\OneDrive\Desktop\InterplanetarySeniorDesign\POSTdata'
    filename1 = r'\cl_off'
    test_data = pd.read_csv(filepath + filename1 + r'.hdf', delim_whitespace = True)
    
    # load in the data from POST part 2 - clcd off
    filename2 = r'\clcd_off'
    test_data2 = pd.read_csv(filepath + filename2 + r'.hdf', delim_whitespace = True)
    
    # load in the data from POST part 2 - clcd on
    filename3 = r'\clcd_on'
    test_data3 = pd.read_csv(filepath + filename3 + r'.hdf', delim_whitespace = True)
    
    theta = np.linspace(0, 2*np.pi, 360)
    r = np.full((360), 3575*10**3)
    
    plt.rcParams.update({'font.size': 10})
    plt.rcParams['font.family'] = 'times new roman'
    fig1 = plt.figure()
    ax = fig1.add_subplot(111, polar = True)
    ax.set_theta_zero_location("S")
    ax.set_xticklabels([])
    ax.set_yticklabels([10, 20, 30, 40, 50, 60])
    ax.set_ylabel('10$^3$ km')
    ax.set_rlabel_position(270)
    ax.yaxis.grid(linewidth = 0.1)
    ax.xaxis.grid(linewidth = 0)
    ax.plot((test_data['longi'] + 4)*(np.pi/180), test_data['gcrad'])
    ax.plot((test_data2['longi'] + 4)*(np.pi/180), test_data2['gcrad'], ls = 'dashed')
    ax.plot((test_data3['longi'] + 4)*(np.pi/180), test_data3['gcrad'], ls = 'dotted')
    ax.plot(theta, r)
    circle = plt.Circle((0.0, 0.0), 2575*10**3, transform = ax.transData._b, 
                        color = 'grey', alpha = 0.5)
    ax.add_artist(circle)
    plt.title('Titan Aerocapture')
    plt.legend(['CL off', 
                'CLCD off',
                'CLCD on',
                'Titan Atmospheric Radius'], loc = 3)
    plt.figtext(0.7, 0.8,  'Atmospheric Entry Velocity = 8.29 km/s')
    plt.figtext(0.7, 0.77, 'Coefficient of Lift = 0.359')
    plt.figtext(0.7, 0.74, 'Coefficient of Drag = 1.18')
    plt.figtext(0.7, 0.71, 'Mass = 600 kg')
    plt.figtext(0.7, 0.68, 'Reference Area = 1 m$^2$')
    plt.figtext(0.7, 0.65, 'Entry Angle = -36.77$^{\circ}$')
    plt.figtext(0.7, 0.62, 'Bank Angle = 180$^{\circ}$')
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, polar = True)
    ax2.set_theta_zero_location("S")
    ax2.set_xticklabels([])
    ax2.set_yticklabels([2, 4, 6, 8, 10])
    ax2.set_ylabel('10$^3$ km')
    ax2.set_rlabel_position(270)
    ax2.yaxis.grid(linewidth = 0.1)
    ax2.xaxis.grid(linewidth = 0)
    ax2.plot((test_data['longi'][2615:3437] + 4)*(np.pi/180), test_data['gcrad'][2615:3437])
    ax2.plot((test_data2['longi'][2615:3370] + 4)*(np.pi/180), test_data2['gcrad'][2615:3370], ls = 'dashed')
    ax2.plot((test_data3['longi'][2615:3786] + 4)*(np.pi/180), test_data3['gcrad'][2615:3786], ls = 'dotted')
    ax2.plot(theta, r)
    circle2 = plt.Circle((0.0, 0.0), 2575*10**3, transform = ax2.transData._b, 
                        color = 'grey', alpha = 0.5)
    ax2.add_artist(circle2)
    plt.title('Titan Aerogravity Assist')
    plt.legend(['CL off', 
                'CLCD off',
                'CLCD on',
                'Titan Atmospheric Radius'], loc = 3)
    plt.figtext(0.7, 0.8,  'Atmospheric Entry Velocity = 8.29 km/s')
    plt.figtext(0.7, 0.77, 'Coefficient of Lift = 0.359')
    plt.figtext(0.7, 0.74, 'Coefficient of Drag = 1.18')
    plt.figtext(0.7, 0.71, 'Mass = 600 kg')
    plt.figtext(0.7, 0.68, 'Reference Area = 1 m$^2$')
    plt.figtext(0.7, 0.65, 'Entry Angle = -36.77$^{\circ}$')
    plt.figtext(0.7, 0.62, 'Bank Angle = 180$^{\circ}$')
    
    plt.show()