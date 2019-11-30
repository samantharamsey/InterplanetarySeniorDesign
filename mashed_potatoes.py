# -*- coding: utf-8 -*-
'''
Created on Sat Nov 16 11:29:58 2019
@author: sam
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


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


gamma = np.linspace(1, 359, 359) # flight path angle array
data = pd.DataFrame([]) # initialize an empty dataframe

for i in gamma:
    
    mu = 37.931*10**6 # saturns gravitational parameter
    saturn_equatorial = 60268 # km
    saturn_polar = 54364 # km
    
    r1 = 1.2*10**6 # titans orbit radius / radius of spacecraft post flyby
    r_encel = 238000 #km
    rp_max = r_encel + 1
    
    r_min = 70000 # km
    rp_min = r_min + 1
    
    e_max = 2
    e_min = 2
    
    if 0 < i < 90 or 270 < i < 360:
        v1_max = 7.94 #km/s 
        v1_min = 7.94 #km/s  
        
        # calculate the maxiumums
        while rp_max > r_encel or e_max >= 1.0:
            v1_max = v1_max - 0.01
            vp_max, rp_max = vprp(i, v1_max)
            
            # print some stuff to see progress
            print('%-13s %-20s %-20s %-20s'  
              %('v1', 'gamma', 'vp', 'rp'))
            print('%5.1f %20.10f %20.10f %20.10f' 
                  %(v1_max, i, vp_max, rp_max))
            
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
            
    elif 90 < i < 270:
        v1_max = -7.94 #km/s 
        v1_min = -7.94 #km/s 
        
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
    
    
    data = data.append(pd.DataFrame({'gamma (deg)': i, 
                                     'v1 max': abs(v1_max),
                                     'vp max': vp_max, 
                                     'rp max': rp_max, 
                                     'max semimajor axis': a_max,
                                     'max eccentricity': e_max, 
                                     'max energy': E_max,
                                     'max momentum': H_max,
                                     'v1 min': abs(v1_min),
                                     'vp min': vp_min, 
                                     'rp min': rp_min, 
                                     'min semimajor axis': a_min,
                                     'min eccentricity': e_min, 
                                     'min energy': E_min,
                                     'min momentum': H_min},
                                     index = [0]), ignore_index = True)
    
data.to_csv(r'C:\Users\saman\OneDrive\Desktop\potato.csv', index = False)

ax = plt.subplot(111, polar = True)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.plot(data['gamma (deg)'][:64]*(np.pi/180), data['v1 max'][:64], c = 'green')
ax.plot(data['gamma (deg)'][63:77]*(np.pi/180), data['v1 max'][63:77], c = 'orange')
ax.plot(data['gamma (deg)'][76:104]*(np.pi/180), data['v1 max'][76:104], c = 'red')
ax.plot(data['gamma (deg)'][:76]*(np.pi/180), data['v1 min'][:76], c = 'blue')
ax.plot(data['gamma (deg)'][103:116]*(np.pi/180), data['v1 max'][103:116], c = 'orange')
ax.plot(data['gamma (deg)'][115:244]*(np.pi/180), data['v1 max'][115:244], c = 'green')
ax.plot(data['gamma (deg)'][243:257]*(np.pi/180), data['v1 max'][243:257], c = 'orange')
ax.plot(data['gamma (deg)'][256:284]*(np.pi/180), data['v1 max'][256:284], c = 'red')
ax.plot(data['gamma (deg)'][283:296]*(np.pi/180), data['v1 max'][283:296], c = 'orange')
ax.plot(data['gamma (deg)'][295:]*(np.pi/180), data['v1 max'][295:], c = 'green')
ax.plot(data['gamma (deg)'][:76]*(np.pi/180), data['v1 min'][:76], c = 'blue')
ax.plot(data['gamma (deg)'][75:104]*(np.pi/180), data['v1 min'][75:104], c = 'red')
ax.plot(data['gamma (deg)'][103:256]*(np.pi/180), data['v1 min'][103:256], c = 'blue')
ax.plot(data['gamma (deg)'][255:284]*(np.pi/180), data['v1 min'][255:284], c = 'red')
ax.plot(data['gamma (deg)'][283:]*(np.pi/180), data['v1 min'][283:], c = 'blue')
plt.title('Family of v1 Velocity Vectors wrt Saturn')
plt.legend(['Good Trajectories - r_p about equal to r_enceladus', 
            'OK Trajectories - Hitting Escape Velocity Constraint', 
            'Bad Trajectories - Perigee Radius Smaller than Saturns Radius',
            'Minimum Required Velocity'], loc = 8)
plt.show()