# -*- coding: utf-8 -*-
"""
11/25/19
@author: hannah
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main():

    # initializing these because Pycharm didn't believe they were real outside of the if statements
    v1 = 0
    vp = 0
    a = 0
    E = 0
    H = 0

    def vprp(fpa, v1):
        """
        Calculates the velocity and radius at perigee
        Args:
            fpa - flight path angle in degrees
            v1 - velocity post Titan aerocapture
        Returns:
            vp - subsequent velocity at perigee
            rp - subsequent radius at perigee
        """

        gamma = fpa*(np.pi/180)
        mu = 3.794*10**7
        r1 = 1.2*10**6

        a = (1/2)
        b = -(mu/(v1*r1*np.cos(gamma)))
        c = -((v1**2/2)-(mu/r1))

        vp = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
        rp = (v1*r1*np.cos(gamma))/vp

        return vp, rp

    gamma = np.linspace(0, 359, 360)    # flight path angle array
    data = pd.DataFrame([])     # initialize an empty dataframe

    for i in gamma:
        mu = 37.931*10**6   # saturn's gravitational parameter
        r1 = 1.2*10**6      # titans orbit radius / radius of spacecraft post flyby
        r_encel = 238000    # km
        r_saturn = 58232    # km
        r_min = r_saturn + 1000
        rp = r_saturn + 1
        e = 2

        if 0 < i < 90 or 270 < i < 360:
            v1 = 0   # km/s
            while e >= 1.0:
                v1 += 1
                vp, rp = vprp(i, v1)

                print('%-13s %-20s %-20s %-20s'
                      % ('v1', 'gamma', 'vp', 'rp'))
                print('%5.1f %20.10f %20.10f %20.10f'
                      % (v1, i, vp, rp))

                E = (1/2)*vp**2 - (mu/rp)   # energy equation
                H = vp*rp   # specific angular momentum
                a = -mu/(2*E)
                e = (a - rp)/a

        elif 90 < i < 270:
            v1 = 0
            while e >= 1.0:
                v1 -= 1
                vp, rp = vprp(i, v1)

                print('%-13s %-20s %-20s %-20s'
                      % ('v1', 'gamma', 'vp', 'rp'))
                print('%5.1f %20.10f %20.10f %20.10f'
                      % (v1, i, vp, rp))

                E = (1/2)*vp**2 - (mu/rp)   # energy equation
                H = vp*rp   # specific angular momentum
                a = -mu/(2*E)
                e = (a - rp)/a

        data = data.append(pd.DataFrame({'v1': abs(v1),
                                         'gamma (deg)': i,
                                         'vp': vp,
                                         'rp': rp,
                                         'semimajor axis': a,
                                         'eccentricity': e,
                                         'energy': E,
                                         'momentum': H}, index=[0]), ignore_index=True)

    data.to_csv(r'C:\Spice_Kernels\potato.csv', index=False)

    ax = plt.subplot(111, polar=True)
    ax.set_theta_zero_location("N")
    ax.plot(data['gamma (deg)'][:65]*(np.pi/180), data['v1'][:65], c='green')
    ax.plot(data['gamma (deg)'][64:78]*(np.pi/180), data['v1'][64:78], c='orange')
    ax.plot(data['gamma (deg)'][77:105]*(np.pi/180), data['v1'][77:105], c='red')
    ax.plot(data['gamma (deg)'][104:117]*(np.pi/180), data['v1'][104:117], c='orange')
    ax.plot(data['gamma (deg)'][116:245]*(np.pi/180), data['v1'][116:245], c='green')
    ax.plot(data['gamma (deg)'][244:258]*(np.pi/180), data['v1'][244:258], c='orange')
    ax.plot(data['gamma (deg)'][257:285]*(np.pi/180), data['v1'][257:285], c='red')
    ax.plot(data['gamma (deg)'][284:297]*(np.pi/180), data['v1'][284:297], c='orange')
    ax.plot(data['gamma (deg)'][296:]*(np.pi/180), data['v1'][296:], c='green')
    plt.title('Family of Minimum v1 Velocity Vectors')
    plt.legend(['Good Trajectories', 'OK Trajectories - Hitting Escape Velocity Constraint',
               'Bad Trajectories - Perigee Radius Smaller than Saturn\'s Radius'], loc=8)
    plt.show()


if __name__ == '__main__':
    main()
