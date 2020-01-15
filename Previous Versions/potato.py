# -*- coding: utf-8 -*-
'''
Created on Sat Nov 16 11:29:58 2019
@author: sam
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from Orbital_Functions import RV2COE

gamma = np.linspace(0.0, 2*np.pi, 72) # flight path angle array
data = pd.DataFrame([]) # initialize an empty dataframe

for i in gamma:
    
    mu = 37.931*10**6 # saturns gravitational parameter
    r1 = 1186780.3668940281 # titans orbit radius / radius of spacecraft post flyby
    r_encel = 238037 #km
    rp = r_encel + 1
    v1 = 5 #km/s  

    E = (1/2)*v1**2 - (mu/r1) # energy equation
    H = v1*r1*np.cos(i) # specific angular momentum
    ecc = 2
        
    while rp > r_encel or E >= 0.0:
        
        v1 = v1 - 0.01
        E = (1/2)*v1**2 - (mu/r1) # energy equation
        H = v1*r1*np.cos(i) # specific angular momentum
        
        def equations(p):
            '''
            Defining simultaneous equations with multiple 
            unknowns to solve for rp and vp
            Args:
                p - single input for unknown variables
            '''
            vp, rp = p
            eq1 = (1/2)*(vp**2) - (mu/rp) - E
            eq2 = vp*rp - H
            return eq1, eq2

        x, y =  fsolve(equations, (15, 238037))
        print('%-13s %-20s %-20s %-20s'  
          %('v1', 'gamma', 'vp', 'rp'))
        print('%5.1f %20.10f %20.10f %20.10f' 
              %(v1, i, x, y))

        rp = y
        vp = x
        
        # Conversions
        fpa = i*(180/np.pi)
        vpx, vpy = vp*np.sin(i), vp*np.cos(i)
        rpx, rpy = rp*np.sin(i), rp*np.cos(i)
        fpa = i*(180/np.pi)
        state = np.array([rpx, rpy, 0.0, vpx, vpy, 0.0])
        hvec, hmag, nvec, nmag, evec, emag, sme, a, p, inc, omega, omtrue, an = RV2COE(mu, state)
        ecc = emag
        
#        emag = emag
        
    data = data.append(pd.DataFrame({'v1': v1, 'gamma (deg)': fpa, 
                                     'vp': vp, 'rp': rp, 
                                     'angular momentum': hmag,
                                     'eccentricity': emag, 
                                     'specific mechanical energy': sme,
                                     'energy': E,
                                     'semi-major axis': a, 
                                     'semiparameter': p, 
                                     'inclination': inc, 
                                     'RAAN': omega,  
                                     'true anomaly': an,
                                     'gamma (rad)': i}, index = [0]), ignore_index = True)

ax = plt.subplot(111, polar = True)
ax.plot(data['gamma (rad)'][::], data['v1'][::])
plt.show()