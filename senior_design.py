# -*- coding: utf-8 -*-
'''
Created on Mon Nov  4 21:06:57 2019

@author: sam
'''

import numpy as np
import spiceypy as spice
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve
from matplotlib import animation


############################# FURNSH THE KERNELS ##############################

spice.tkvrsn('TOOLKIT')

path = r'C:/Users/saman/OneDrive/Desktop/spice_kernels/'
spice.furnsh(path + r'de430.bsp')
spice.furnsh(path + r'naif0009.tls')
spice.furnsh(path + r'sat425.bsp')
spice.furnsh(path + r'cpck05Mar2004.tpc')
spice.furnsh(path + r'020514_SE_SAT105.bsp')
spice.furnsh(path + r'981005_PLTEPH-DE405S.bsp')
spice.furnsh(path + r'030201AP_SK_SM546_T45.bsp')
spice.furnsh(path + r'04135_04171pc_psiv2.bc')
spice.furnsh(path + r'sat425.inp')
spice.furnsh(path + r'sat425.cmt')

############################### FIND THE MOONS ################################

step = 4000
utc = ['Jun 20, 2022', 'Dec 1, 2022']

etOne = spice.str2et(utc[0])
etTwo = spice.str2et(utc[1])
print('ET One: {}, ET Two: {}'.format(etOne, etTwo))
times = [x*(etTwo-etOne)/step + etOne for x in range(step)]
print(times[0:3])


encel_pos, encel_lightTimes = spice.spkpos('Enceladus', times, 'J2000', 
                                     'NONE', 'SATURN BARYCENTER')
titan_pos, titan_lightTimes = spice.spkpos('Titan', times, 'J2000', 
                                     'NONE', 'SATURN BARYCENTER')
titan_posa = np.asarray(titan_pos).T
encel_posa = np.asarray(encel_pos).T
fig = plt.figure(figsize=(9, 9))
ax  = fig.add_subplot(111, projection='3d')
ax.plot(titan_posa[0], titan_posa[1], titan_posa[2])
ax.plot(encel_posa[0], encel_posa[1], encel_posa[2])

u = np.linspace(0, 2 * np.pi, 39)
v = np.linspace(0, np.pi, 21)
x = 60268 * np.outer(np.cos(u), np.sin(v))
y = 60268 * np.outer(np.sin(u), np.sin(v))
z = 5436.4 * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, rstride=3, cstride=3, color='blue', shade=0)

def rotate(angle):
    ax.view_init(azim = angle)

rot_animation = animation.FuncAnimation(fig, rotate, 
                                        frames=np.arange(0, 362, 2), 
                                        interval=100)

################################# DO THE MATH #################################

# Saturn gravitational parameter
mu = 37.931*10**6 # km^3/s^2

# Random value for velocity of the spacecraft when it reaches Titan
vp = 5.393 #km/s

# Enceladus' orbital distance from Saturn
vi = 238037 #km
encel_ecc = 0.0047

# Compute r when the spacecraft reaches Titan
# Not necessarily the same as Titans perigee radius
# Dependant on first value of utc defined above
# Equal to the center of Titan, need to adjust to Titan's SOI
titan_rp = np.linalg.norm(titan_pos[0]) #km
titan_ecc = 0.02

h = vp*titan_rp 
E = (1/2)*vp**2 - (mu/titan_rp)

def f(vars):
    vi, ri, gamma = vars
    eq1 = (1/2)*vi**2 - (mu/ri)
    eq2 = vi*ri*np.cos(gamma)
    return [eq1, eq2]

gamma = np.arange(0, 2*np.pi, 5*(180/np.pi))
sol = []
for i in gamma:
    x, y =  fsolve(f, (1, 1, i))
    sol.append(x, y)  
print(sol)

# Want the spacecraft to cross Enceladus' orbit when it is actually there
# Worry about redevous timing later
# rp - will not change - based on arrival date
# ri - will vary, but is still constrained to Enceladus' orbit
# vp - variable
# vi - variable
# fpa - range from 0 to 2pi in 5 deg increments

# Need to take into consideration the time of flight between rp and ri and then
# find the orbit that allows for most frequent flybys of Enceladus
# Save redevous for later




# clock angle - what is the angle between the line between saturn and titand 
# periapse around titan at the flyby
# take inbound vinf vec do everything planar, rotate that angle around 360 deg
# two body OM to propagate into periapse and back then used patched conics
# take hyperbolic excess vel vec, R portion of spherical coord

# Assume circles
# Circular restriced 3 body problem

# pick periapse r
# pick hyp excess vel mag
# find flyby
# rotate the vinf vec
# vary RA of vinf vec, and vary vinf, vary periapse r
# turning angle func of ecc in flyby orbit
# semi major axis through energy
# chaper 12.2
# equation 12.10
# tabgent to fpa

# anderson slides

# adam.r.harden@nasa.gov







