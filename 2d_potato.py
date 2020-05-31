# -*- coding: utf-8 -*-
'''
Created on Sat Nov 16 11:29:58 2019
@author: sam
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# define some constants
mu = 37.931*10**6   # saturns gravitational parameter
saturn_equatorial = 60268   # km
saturn_polar = 54364    # km
r_titan = 1.2*10**6  # titans orbit radius
r_encel = 238000    # km
v_titan = 5.57  # km/s


def vprp(fpa, v1):
    '''
    Calculates the velocity and radius at perigee
    Args:
        fpa - flight path angle in degrees
        v1 - velocity post Titan aerocapture
    Returns:
        vp - subsequent velocity at perigee
        rp - subsequent radius at perigee
    '''
    
    gamma = fpa*(np.pi/180)
    mu = 3.794*10**7
    r1 = 1.2*10**6
    
    a = (1/2)
    b = -(mu/(v1*r1*np.cos(gamma)))
    c = -((v1**2/2)-(mu/r1))
    
    vp = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
    rp = (v1*r1*np.cos(gamma))/vp
    
    return vp, rp



gamma = np.linspace(1, 360, 360)  # flight path angle array
data = pd.DataFrame([])  # initialize an empty dataframe

for i in gamma:

    mu = 37.931*10**6  # saturn's gravitational parameter
    saturn_equatorial = 60268  # km
    saturn_polar = 54364  # km

    r1 = 1.2*10**6  # titans orbit radius / radius of spacecraft post flyby
    r_encel = 238000  # km
    rp_max = r_encel + 1

    r_min = 70000  # km
    rp_min = r_min + 1

    e_max = 2
    e_min = 2

    if 0 < i < 90 or 270 < i < 360:
        v1_max = 7.94  # km/s
        v1_min = 7.94  # km/s

        # calculate the maximums
        while rp_max > r_encel or e_max >= 1.0:
            v1_max = v1_max - 0.01
            vp_max, rp_max = vprp(i, v1_max)

            E_max = (1/2)*vp_max**2 - (mu/rp_max) # energy equation
            H_max = vp_max*rp_max # specific angular momentum
            a_max = -mu/(2*E_max)
            e_max = (a_max - rp_max)/a_max

        # calculate the minimums
        while rp_min > r_min or e_min >= 1.0:
            v1_min = v1_min - 0.01
            vp_min, rp_min = vprp(i, v1_min)

            E_min = (1/2)*vp_min**2 - (mu/rp_min) # energy equation
            H_min = vp_min*rp_min # specific angular momentum
            a_min = -mu/(2*E_min)
            e_min = (a_min - rp_min)/a_min
            
        # print some stuff to see progress
        print('%-13s %-20s %-20s %-20s'
              % ('v1', 'gamma', 'vp', 'rp'))
        print('%5.1f %20.10f %20.10f %20.10f'
              % (v1_max, i, vp_max, rp_max))

    elif 90 < i < 270:
        v1_max = -7.94 # km/s
        v1_min = -7.94 # km/s

        # calculate the  maximums
        while rp_max > r_encel or e_max >= 1.0:
            v1_max = v1_max + 0.01
            vp_max, rp_max = vprp(i, v1_max)

            E_max = (1/2)*vp_max**2 - (mu/rp_max) # energy equation
            H_max = vp_max*rp_max # specific angular momentum
            a_max = -mu/(2*E_max)
            e_max = (a_max - rp_max)/a_max

        # calculate the minimums
        while rp_min > r_min or e_min >= 1.0:
            v1_min = v1_min + 0.01
            vp_min, rp_min = vprp(i, v1_min)

            E_min = (1/2)*vp_min**2 - (mu/rp_min) # energy equation
            H_min = vp_min*rp_min # specific angular momentum
            a_min = -mu/(2*E_min)
            e_min = (a_min - rp_min)/a_min

    data = data.append(pd.DataFrame({'Saturn fpa (deg)': i,
                                    'Saturn v1 max': abs(v1_max),
                                    'Saturn v1x max': v1_max*np.sin(i*(np.pi/180)),
                                    'Saturn v1y max': v1_max*np.cos(i*(np.pi/180)),
                                    'vp max': vp_max,
                                    'rp max': rp_max,
                                    'max semimajor axis': a_max,
                                    'max eccentricity': e_max,
                                    'max energy': E_max,
                                    'max momentum': H_max,
                                    'Saturn fpa (deg)': i,
                                    'Saturn v1 min': abs(v1_min),
                                    'Saturn v1x min': v1_min*np.sin(i*(np.pi/180)),
                                    'Saturn v1y min': v1_min*np.cos(i*(np.pi/180)),
                                    'vp min': vp_min,
                                    'rp min': rp_min,
                                    'min semimajor axis': a_min,
                                    'min eccentricity': e_min,
                                    'min energy': E_min,
                                    'min momentum': H_min},
                                    index=[0]), ignore_index=True)

# wrt Saturn
vx_max = data['Saturn v1 max']*np.sin(data['Saturn fpa (deg)']*(np.pi/180))
vy_max = data['Saturn v1 max']*np.cos(data['Saturn fpa (deg)']*(np.pi/180))
vz_max = [0]*len(vx_max)

vx_min = data['Saturn v1 min']*np.sin(data['Saturn fpa (deg)']*(np.pi/180))
vy_min = data['Saturn v1 min']*np.cos(data['Saturn fpa (deg)']*(np.pi/180))
vz_min = [0]*len(vx_min)

# wrt Titan
# v_inf = v1 - v_titan in the y-direction
new_vymax = vy_max - v_titan
new_vymin = vy_min - v_titan

vmax = []
vmin = []
max_fpa = []
min_fpa = []

for i in range(360):
    vmax_titan = np.linalg.norm([vx_max[i], new_vymax[i], vz_max[i]])
    vmax.append(vmax_titan)
    vmin_titan = np.linalg.norm([vx_min[i], new_vymin[i], vz_min[i]])
    vmin.append(vmin_titan)

titan_fpa_max = np.arctan2(vx_max, new_vymax)
titan_fpa_min = np.arctan2(vx_min, new_vymin)

data.insert(4, 'Titan max fpa', titan_fpa_max*(180/np.pi), True)
data.insert(14, 'Titan min fpa', titan_fpa_min*(180/np.pi), True)

data.insert(5, 'Titan v1 max magnitude', vmax, True)
data.insert(6, 'Titan v1x max', vx_max, True)
data.insert(7, 'Titan v1y max', new_vymax, True)

data.insert(18, 'Titan v1 min magnitude', vmin, True)
data.insert(19, 'Titan v1x min', vx_min, True)
data.insert(20, 'Titan v1y min', new_vymin, True)

data.to_csv(r'C:\Spice_Kernels\potato.csv', index=False)

plt.rcParams.update({'font.size': 12})
plt.rcParams['font.family'] = 'times new roman'
f1 = plt.figure()
ax1 = f1.add_subplot(111, polar=True)
ax1.set_theta_zero_location("N")
ax1.set_theta_direction(-1)
ax1.plot(data['Titan max fpa'][:64]*(np.pi/180), data['Titan v1 max magnitude'][:64], c='green')
ax1.plot(data['Titan max fpa'][63:77]*(np.pi/180), data['Titan v1 max magnitude'][63:77], c='orange')
ax1.plot(data['Titan min fpa'][:76]*(np.pi/180), data['Titan v1 min magnitude'][:76], c='blue')
ax1.plot(data['Titan max fpa'][76:104]*(np.pi/180), data['Titan v1 max magnitude'][76:104], c='red')
ax1.plot(data['Titan min fpa'][75:104]*(np.pi/180), data['Titan v1 min magnitude'][75:104], c='red')
ax1.plot(data['Titan max fpa'][103:116]*(np.pi/180), data['Titan v1 max magnitude'][103:116], c='orange')
ax1.plot(data['Titan max fpa'][115:244]*(np.pi/180), data['Titan v1 max magnitude'][115:244], c='green')
ax1.plot(data['Titan min fpa'][103:256]*(np.pi/180), data['Titan v1 min magnitude'][103:256], c='blue')
ax1.plot(data['Titan max fpa'][243:257]*(np.pi/180), data['Titan v1 max magnitude'][243:257], c='orange')
ax1.plot(data['Titan max fpa'][256:284]*(np.pi/180), data['Titan v1 max magnitude'][256:284], c='red')
ax1.plot(data['Titan min fpa'][255:284]*(np.pi/180), data['Titan v1 min magnitude'][255:284], c='red')
ax1.plot(data['Titan max fpa'][283:296]*(np.pi/180), data['Titan v1 max magnitude'][283:296], c='orange')
ax1.plot(data['Titan max fpa'][295:]*(np.pi/180), data['Titan v1 max magnitude'][295:], c='green')
ax1.plot(data['Titan min fpa'][283:]*(np.pi/180), data['Titan v1 min magnitude'][283:], c='blue')
ax1.fill(data['Titan max fpa'][::]*(np.pi/180), data['Titan v1 max magnitude'][::], 'darkseagreen',
         data['Titan min fpa'][::]*(np.pi/180), data['Titan v1 min magnitude'][::], 'white')
plt.title('Family of v1 Velocity Vectors wrt Titan')
plt.legend(['Maximum Allowable v1 Value',
            'Maximum Allowable v1 Value - Hitting Escape Velocity Constraint',
            'Minimum Allowable v1 Value',
            'Bad Trajectories - Perigee Radius Smaller than Saturn\'s Radius'])

f2 = plt.figure()
ax2 = f2.add_subplot(111, polar=True)
ax2.set_theta_zero_location("N")
ax2.set_theta_direction(-1)
ax2.plot(data['Saturn fpa (deg)'][:64]*(np.pi/180), data['Saturn v1 max'][:64], c='green')
ax2.plot(data['Saturn fpa (deg)'][63:77]*(np.pi/180), data['Saturn v1 max'][63:77], c='orange')
ax2.plot(data['Saturn fpa (deg)'][76:104]*(np.pi/180), data['Saturn v1 max'][76:104], c='red')
ax2.plot(data['Saturn fpa (deg)'][:76]*(np.pi/180), data['Saturn v1 min'][:76], c='blue')
ax2.plot(data['Saturn fpa (deg)'][103:116]*(np.pi/180), data['Saturn v1 max'][103:116], c='orange')
ax2.plot(data['Saturn fpa (deg)'][115:244]*(np.pi/180), data['Saturn v1 max'][115:244], c='green')
ax2.plot(data['Saturn fpa (deg)'][243:257]*(np.pi/180), data['Saturn v1 max'][243:257], c='orange')
ax2.plot(data['Saturn fpa (deg)'][256:284]*(np.pi/180), data['Saturn v1 max'][256:284], c='red')
ax2.plot(data['Saturn fpa (deg)'][283:296]*(np.pi/180), data['Saturn v1 max'][283:296], c='orange')
ax2.plot(data['Saturn fpa (deg)'][295:]*(np.pi/180), data['Saturn v1 max'][295:], c='green')
ax2.plot(data['Saturn fpa (deg)'][:76]*(np.pi/180), data['Saturn v1 min'][:76], c='blue')
ax2.plot(data['Saturn fpa (deg)'][75:104]*(np.pi/180), data['Saturn v1 min'][75:104], c='red')
ax2.plot(data['Saturn fpa (deg)'][103:256]*(np.pi/180), data['Saturn v1 min'][103:256], c='blue')
ax2.plot(data['Saturn fpa (deg)'][255:284]*(np.pi/180), data['Saturn v1 min'][255:284], c='red')
ax2.plot(data['Saturn fpa (deg)'][283:]*(np.pi/180), data['Saturn v1 min'][283:], c='blue')
ax2.fill_between(data['Saturn fpa (deg)'][::]*(np.pi/180), data['Saturn v1 max'][::],
                 data['Saturn v1 min'][::], facecolor='darkseagreen')
plt.title('Family of Acceptable Post-AGA Velocity Vectors wrt Saturn')
plt.legend(['Trajectories with $r_p$ about equal to $r_{enceladus}$',
            'Escape Velocity Constraint',
            'Trajectories with Perigee Radius Smaller than Saturn\'s Radius',
            'Minimum Required Velocity'], loc=8)
ax2.set_ylabel('Velocity in km/s')
ax2.set_rlabel_position(270)
ax2.set_xlabel('Flight Path Angle in Degrees')
ax2.set_xticklabels([0, 45, 90, 135, 180, 225, ' ', 315])
plt.show()


