# -*- coding: utf-8 -*-
'''
Created on Thu Dec  5 12:26:07 2019
@author: sam
'''

import numpy as np


def example1():
    '''
    Example 9.1
    A spacecraft is to be launched from Earth to Venus using a Hohmann transfer 
    (assume the orbits of Earth and Venus are circular). It is to do a flyby of 
    Venus that will then take it toward Mercury. The Venus flyby will be planned 
    so that the heliocentric trajectory after the flyby is an ellipse with a 
    period equal to 5/2 of the orbital period of Mercury.
    
    a) What semimajor axis is required of the heliocentric orbit after the Venus
       flyby? 
       Solution: a = 10.67*10^7 km = 0.7129 au
    b) What is the required turn angle, delta, at Venus that will yeild a post-
       flyby heliocentric orbit having the semimajor axis determined in part (a)?
       Solution: v = 34.76 km/s; delta = 97.72 deg
    c) What flyby altitude at Venus is required to yield this delta? Is it a
       feasible flyby altitude? Solution: rp = 6053 km; yes because rp > r Venus.
    d) What will be the eccentricity of the heliocentric, post-flyby orbit?
       Solution: e = 0.0774
    e) With the post-flyby semimajor axcis and eccentricity of parts (a) and (d),
       is it possible for the spcecraft to cross the highly eccentric orbit of 
       Mercury? 
       Solution: In order for a crossing to be possible, the post-flyby
       orbit perihelion must be smaller than the apehelion of Mercury's orbit 
       (0.4667 AU). This is not true for this flyby.
    '''
    
    # some constants
    mu_sun = 1.327*(10**11)
    T_merc = 87.96*24*(60**2)
    a_venus = 1.081608*(10**8)
    mu_venus = 3.24859*(10**5)
    v_venus = 35.02

    print(' ')
    print('Example 9.1 Solutions:')
    
    # part a
    T_postflyby = (5/2)*T_merc
    a_merc = (mu_sun*((T_postflyby/(2*np.pi))**2))**(1/3)
    print('a) {} km/s'.format(a_merc))
    
    # part b
    # post flyby velocity using the vis-visa equation
    v_0 = np.sqrt(mu_sun*((2/a_venus) - (1/a_merc))) # km/s
    # v perihelion post hohman transfer from earth
    v_perihelion = 37.73
    v_inf = v_perihelion - v_venus
    gamma = np.arccos((v_venus**2 + v_inf**2 - v_0**2)/(2*v_venus*v_inf))
    delta = np.pi - gamma
    print('b) {} degrees'.format(delta*(180/np.pi)))
    
    # part c
    rp = ((np.sin(delta/2))**-1 - 1)/(v_inf**2/mu_venus)
    print('c) {} km - feasible'.format(rp))
    
    # part d
    beta = v_inf/(v_0/np.sin(gamma))
    h = a_venus*v_0*np.sin(np.pi/2 - beta)
    e = np.sqrt(-1*(h**2/(mu_sun*a_merc) - 1))
    print('c) {}'.format(e))
    
    # part e
    rp_postflyby = a_merc*(1 - e)
    print('e) rp = {} > r Mercury - not possible'.format(rp_postflyby))


def example2():
    '''
    Example 9.2
    A spacecraft is being sent tfrom Earth to the outer solar system. A flyby of
    Saturn is to be used to add energy to the orbit and change its inclination.
    The spacecraft travels from Earth orbit about the sun (no hyperbolic escape
    needed) to Saturn on a Hohmann transfer trajectory. Assume that Saturn's orbit
    lies on the ecliptic plane, and assume Saturn's orbit is circular.
    
    a) If the flyby is directly over the south pole os Saturn, at an altitude of 
       150,000 km, the spacecraft will be sent out of the ecliptic plane. That is, 
       the rotation of the v infinity vector will add a component to the final 
       (post flyby) velocity that is normal to the ecliptic place. What will be
       the post flyby inclination of the spacecrafts heliocentric orbit? 
       Solution: i = 21.34 degrees
    b) What will be the semimajor axis of the spacecraft orbut after the flyby?
       Solution: v = 13.13 km/s
    c) What will be the eccentricity of the spacecraft orbit after the flyby?
       Solution: e = 0.854
    '''
    
    # some constants
    r_saturn = 9.49*6378.12
    a_saturn_au = 9.554909595
    au2km = 1.4959965*10**8
    tu2s = 5.0226757*10**6
    autu2kms = au2km/tu2s
    kms2autu = tu2s/au2km
    mu_sun = 1
    mu_saturn = 3.794*10**7
    print(' ')
    print('Example 9.2 Solutions:')
    
    # part a
    rp = 150000 + r_saturn
    a_h = (a_saturn_au + mu_sun)/2 # au
    v_A = np.sqrt(mu_sun*((2/a_saturn_au) - (1/a_h)))*autu2kms
    v_saturn = np.sqrt(mu_sun/a_saturn_au)*autu2kms
    v_inf = v_saturn - v_A
    delta = 2*(np.arcsin(1/(1 + ((rp*v_inf**2)/(mu_saturn)))))
    e = 1/(1 + ((rp*v_inf**2)/(mu_saturn)))
    # law of cosines to find v0
    v_0 = np.sqrt(v_saturn**2 + v_inf**2 - 2*v_saturn*v_inf*np.cos(delta))
    i = np.arcsin(v_inf/(v_0/np.sin(delta)))*(180/np.pi)
    print('a) {} degrees'.format(i))
    
    # part b
    a = a_saturn_au/(2 - (a_saturn_au*((v_0*(kms2autu))**2)/mu_sun))
    print('b) {} au'.format(a))
    
    # part c
    e = -1*(a_saturn_au/a - 1)
    print('c) {}'.format(e))


def example3():
    '''
    Determining the Perilune Altitude for an Apollo Free_Return Trajectory.
    Suppose the radial component of the geocentric spacecraft velocity at lunar 
    arrival is 0.75 km/s. (This number is chosen as being representative of the 
    actual Apollo missions.) What is the perilune altitude that will yield a 
    free-return trajectory?
    '''
    
    # some constants
    a = 384000 #km
    mu_earth = 3.986004415*10**5 # km^3/s^2
    mu_moon = 0.00490*10**6 # km^3/s^2
    r_moon = 1738 #km
    v_r = 0.75 # km/s
    v_moon = 1.023 #km/s
    e = 0.983
    print(' ')
    print('Example 9.3 Solution:')
    
    v_0 = np.sqrt((mu_earth*a*(1 - e)*(1 + e)))/(384400 + 300)
    v_mag = np.sqrt(v_r**2 + v_0**2)
    fpa = np.arccos(v_0/v_mag)
    v_inf_moon = np.sqrt(v_mag**2 + v_moon**2 - 2*v_mag*v_moon*np.cos(fpa))
    new_e = 1.499
    rp_moon = ((new_e - 1)*mu_moon)/v_inf_moon**2
    perilune = rp_moon - r_moon
    print('Required perilune altitude is {} km'.format(perilune))


if __name__ == '__main__':
    example1()
    example2()
    example3()
