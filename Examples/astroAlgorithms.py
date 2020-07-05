# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 17:03:28 2019

@author: sam
"""

import numpy as np
import math
         

def Julian_Date(year, month, day, hour, mins, sec):
    ''' Converts a calendar date to Julian date '''
    JD = 367*year - ((7*(year + (int(month + 9)/12))))/4 + int(275*month)/9 + day + 1721013.5 + ((((sec/60) + mins)/60) + hour)/24
    return JD


def Modified_Julian_Date(year, month, day, hour, mins, sec):
    ''' Converts a calendar date to modified Julian date '''
    MJD = Julian_Date(year, month, day, hour, mins, sec) - 2400000.5
    return MJD


def Relative_Velocity_Vector(state, omega):
    ''' Determines the relative velocity vector from initial state conditions
        and rotational velocity '''
    vvec = state[3:]
    rvec = state[:3]
    vrel = np.subtract(vvec, np.cross(omega, rvec))
    vrelI = vrel[0]
    vrelJ = vrel[1]
    vrelK = vrel[2]
    return np.array([vrelI, vrelJ, vrelK])


def Drag_Acceleration_Vector(Cd, amratio, state, omega, radius):
    ''' Determines acceleration vector due to drag on an orbiting body
        Cd: drag coefficient (usually 2.2)
        amratio: area to mass ratio of the spacecraft
        state: initial state vector
        omega: rotational velocity of the orbited body
        radius: radius of orbited body '''
    # relative velocity vector
    vrel = Relative_Velocity_Vector(state, omega)
    vrelmag = np.linalg.norm(vrel)
    vrelunit = vrel/vrelmag
    # determine altitude
    rvec = state[:3]
    rmag = np.linalg.norm(rvec)
    alt = rmag - radius
    # determine atmosphereic density
    if alt < 500:
        rho = 3.725*(10**-12)
    elif alt < 600:
        rho = 6.967*(10**-13)
    elif alt < 700:
        rho = 1.454*(10**-13)
    elif alt < 800:
        rho = 3.614*(10**-14)
    elif alt < 900:
        rho = 1.170*(10**-14)
    elif alt < 1000:
        rho = 5.245*(10**-15)
    else:
        rho = 3.019*(10**-15)
    # calculate drag acceleration
    acc = (-1/2)*(Cd*amratio)*rho*(vrelmag**2)*vrelunit
    dragaccI = acc[0]
    dragaccJ = acc[1]
    dragaccK = acc[2]
    return np.array([dragaccI, dragaccJ, dragaccK])


def Third_Body_Acceleration(sat_rvec, body_rvec, earth_rvec, earthmu, bodymu):
    ''' Determines the acceleration vector due to a third bodies gravitational force
        assuming the mass of the orbiting object is negligible
        sat_rvec: initial position vector of the orbiting object
        body_rvec: initial position vector of the third body
        earth_rvec: initial position vector of the earth or orbited body
        earthmu: gravitational constant of the earth or orbited body
        bodymu: gravitational constant of the third body '''
    # earth to satellite
    esrvec = sat_rvec - earth_rvec
    esrmag = np.linalg.norm(esrvec)
    # satellite to third body
    s3rvec = sat_rvec - body_rvec
    s3rmag = np.linalg.norm(s3rvec)
    # earth to third body
    e3rvec = earth_rvec - body_rvec
    e3rmag = np.linalg.norm(e3rvec)
    # acceleration vector due to third body 
    acc = -((earthmu*esrvec)/(esrmag**3)) + bodymu*((s3rvec/(s3rmag**3)) - (e3rvec/(e3rmag**3)))
    return np.array([acc[0], acc[1], acc[2]])


def Third_Body_Acceleration(sat_rvec, body_rvec, bodymu):
    ''' Determines the acceleration vector due to a third bodies gravitational force
        assuming Solar-System Barycenter and the mass of the orbiting object is negligible
        sat_rvec: initial position vector of the orbiting object
        body_rvec: initial position vector of the third body
        bodymu: gravitational constant of the third body '''
    # third body to satellite
    bsrvec = sat_rvec - body_rvec
    bsrmag = np.linalg.norm(bsrvec)

    # acceleration vector due to third body 
    acc = -((bodymu*bsrvec)/(bsrmag**3))
    return np.array([acc[0], acc[1], acc[2]])


def Solar_Radiation_Pressure(cr, amratio, sat_rvec, sun_rvec, earth_rvec):
    ''' Determines acceleration vector due to solar radiation pressure assuming
        the force of solar pressure is a constant 4.57*(10**-6)*(1000**2) N/km^2
        cr: reflectivity - a value between 0.0 and 2.0 which indicates how the 
            satellite reflects incoming radiation
        amratio: the exposed area to the sun / mass the satellite
        sat_rvec: initial position vector of the satellite
        sun_rvec: initial position vector of the sun
        earth_rvec: initial position vector of the earth '''
    # if not eclipsed by earth
    psrp = 4.57*(10**-6)*(1000**2) # N/km^2
    # if eclipsed by earth
    psrp = 0
    rvec = sat_rvec - sun_rvec
    rmag = np.linalg.norm(rvec)
    runit = rvec/rmag
    acc = -psrp*cr*amratio*runit
    return np.array([acc[0], acc[1], acc[2]])


def J2(j2value, state, req, mu):
    rvec = state[:3] # r vector
    rmag = np.linalg.norm(rvec) # r magnitude
    multiplier = (-3/2)*j2value*(mu/(rmag**2))*((req/rmag)**2)
    J2I = multiplier*(1 - 5*((rvec[2]/rmag)**2))*(rvec[0]/rmag)
    J2J = multiplier*(1 - 5*((rvec[2]/rmag)**2))*(rvec[1]/rmag)
    J2K = multiplier*(1 - 5*((rvec[2]/rmag)**2))*(rvec[2]/rmag)
    return np.array([J2I, J2J, J2K])


def J3(j3value, state, req, mu):
    rvec = state[:3] # r vector
    rmag = np.linalg.norm(rvec) # r magnitude
    multiplier = (1/2)*j3value*(mu/(rmag**2))*((req/rmag)**3)
    J3I = multiplier*(5*((7*(rvec[2]/rmag)**3) - (3*(rvec[2]/rmag)))*(rvec[0]/rmag))
    J3J = multiplier*(5*((7*(rvec[2]/rmag)**3) - (3*(rvec[2]/rmag)))*(rvec[1]/rmag))
    J3K = multiplier*(3*(1 - 10*((rvec[2]/rmag)**2) + (35/3)*((rvec[2]/rmag)**4)))
    return np.array([J3I, J3J, J3K])


def J4(j4value, state, req, mu):
    rvec = state[:3]
    rmag = np.linalg.norm(rvec)
    multiplier = (5/8)*j4value*(mu/(rmag**2))*((req/rmag)**4)
    J4I = multiplier*((3 - 42*((rvec[2]/rmag)**3) - 63*((rvec[2]/rmag)**4))*(rvec[0]/rmag))
    J4J = multiplier*((3 - 42*((rvec[2]/rmag)**3) - 63*((rvec[2]/rmag)**4))*(rvec[1]/rmag))
    J4K = multiplier*((15 - 70*((rvec[2]/rmag)**3) - 63*((rvec[2]/rmag)**4))*(rvec[2]/rmag))
    return np.array([J4I, J4J, J4K])


def J5(j5value, state, req, mu):
    rvec = state[:3]
    rmag = np.linalg.norm(rvec)
    multiplier = (j5value/8)*(mu/(rmag**2))*((req/rmag)**5)
    J5I = multiplier*((3*(35*(rvec[2]/rmag) - 210*((rvec[2]/rmag)**3) + 231*((rvec[2]/rmag)**5)))*(rvec[0]/rmag))
    J5J = multiplier*((3*(35*(rvec[2]/rmag) - 210*((rvec[2]/rmag)**3) + 231*((rvec[2]/rmag)**5)))*(rvec[1]/rmag))
    J5K = multiplier*(693*((rvec[2]/rmag)**6) - 945*((rvec[2]/rmag)**4) + 315*((rvec[2]/rmag)**2) - 15)
    return np.array([J5I, J5J, J5K])


def J6(j6value, state, req, mu):
    rvec = state[:3]
    rmag = np.linalg.norm(rvec)
    multiplier = -(j6value/16)*(mu/(rmag**2))*((req/rmag)**6)
    J6I = multiplier*((35 - 945*((rvec[2]/rmag)**2) + 3465*((rvec[2]/rmag)**4) - 3003*((rvec[2]/rmag)**6))*(rvec[0]/rmag))
    J6J = multiplier*((35 - 945*((rvec[2]/rmag)**2) + 3465*((rvec[2]/rmag)**4) - 3003*((rvec[2]/rmag)**6))*(rvec[1]/rmag))
    J6K = multiplier*((245 - 2205*((rvec[2]/rmag)**2) + 4851*((rvec[2]/rmag)**4) - 3003*((rvec[2]/rmag)**6))*(rvec[2]/rmag))
    return np.array([J6I, J6J, J6K])


def RV2COE(mu, state):
    ''' Converts from an initial state of position and velocity vectors to
        the classical orbital elements '''
    # tolerances
    zero = 0.00000001
    tol = 0.0000001
    # classical orbital elements
    rvec = state[:3] # r vector
    rmag = np.linalg.norm(rvec) # r magnitude
    vvec = state[3:] # velocity vector
    vmag = np.linalg.norm(vvec) # velocity magnitude
    hvec = np.cross(rvec, vvec) # angular momentum vector
    hmag = np.linalg.norm(hvec) # angular momentum magnitude
    k = np.array([0, 0, 1]) # unit vector in k direction
    nvec = np.cross(k, hvec) # node vector
    nmag = np.linalg.norm(nvec) # node magnitude
    evec = ((vmag**2 - (mu/rmag))*rvec - np.dot(rvec, vvec)*vvec) / mu # eccentricity vector
    emag = np.linalg.norm(evec) # eccentricity magnitude
    if emag < zero:
        emag = 0
    sme = vmag**2/2 - mu/rmag # specific mechanical energy
    if emag != 1:
        a = -(mu/(2*sme)) # semi major axis
        p = a*(1 - emag**2) # semiparameter: semi-latus rectum or ellipse
    else:
        a = math.inf # a = infinity
        p = hmag**2/mu
    i = np.arccos(hvec[2]/hmag) # inclination
    omega = np.arccos(nvec[0]/nmag)
    if nvec[1] < 0:
        omega = 2*np.pi - omega
    om = np.arccos((np.dot(nvec, evec))/(nmag*emag))
    if evec[2] < 0:
        om = 2*np.pi - om
    an = np.arccos((np.dot(evec, rvec))/(emag*rmag)) # true anomaly
    if np.dot(rvec, vvec) < 0:
        an = 2*np.pi - an
    # convert to degrees
    i = i*(180/np.pi)
    omega = omega*(180/np.pi)
    om = om*(180/np.pi)
    an = an*(180/np.pi)     
    # Special Cases  
    # Elliptical Equatorial
    if (zero < emag < 1) & (i == 0 or i == 180): # tolerance based
        omtrue = np.arccos(evec[0]/emag)
        if evec[1] < 0:
            omtrue = 2*np.pi - omtrue
        return {'angular momentum vector': hvec, 'angular momentum magnitude': hmag, 'node vector': nvec, 
            'node magnitude': nmag, 'eccentricity vector': evec, 'eccentricity magnitude': emag, 
            'specific mechanical energy': sme, 'semi-major axis': a, 'semi-latus rectum': p, 'inclination': i, 
            'right ascension of the ascending node': omega, 'argument of perigee': omtrue, 'true anomaly': an}
    # Circular Inclined
    elif (emag < zero) & (0-tol > i > tol):
        umag = np.arccos(np.dot(nvec, rvec)/(nmag*rmag))
        if rvec[2] < zero:
            umag = 2*np.pi - umag
        return {'angular momentum vector': hvec, 'angular momentum magnitude': hmag, 'node vector': nvec, 
            'node magnitude': nmag, 'eccentricity vector': evec, 'eccentricity magnitude': emag, 
            'specific mechanical energy': sme, 'semi-major axis': a, 'semi-latus rectum': p, 'inclination': i, 
            'right ascension of the ascending node': omega, 'argument of perigee': om, 'true an': umag}
    # Circular Equatorial
    elif (emag == 0) & (0-tol > i < tol): # +/- tolerance
        lamtrue = np.arccos(rvec[0]/rmag)
        if rvec[1] < zero:
            lamtrue = 2*np.pi - lamtrue
        return {'angular momentum vector': hvec, 'angular momentum magnitude': hmag, 'node vector': nvec, 
            'node magnitude': nmag, 'eccentricity vector': evec, 'eccentricity magnitude': emag, 
            'specific mechanical energy': sme, 'semi-major axis': a, 'semi-latus rectum': p, 'inclination': i, 
            'right ascension of the ascending node': omega, 'argument of perigee': om, 'lamda true': lamtrue}
    # if not special case
    else:
        return {'angular momentum vector': hvec, 'angular momentum magnitude': hmag, 'node vector': nvec, 
            'node magnitude': nmag, 'eccentricity vector': evec, 'eccentricity magnitude': emag, 
            'specific mechanical energy': sme, 'semi-major axis': a, 'semi-latus rectum': p, 'inclination': i, 
            'right ascension of the ascending node': omega, 'argument of perigee': om, 'true anomaly': an}


def COE2RV(mu, semip, truean, ecc, inc, OM, om):
    ''' Converts from classical orbital elements to an initial state vector 
        containing position and velocity vectors '''
    # tolerances
    zero = 0.00000001
    tol = 0.0000001
    # special cases
    # Circular Equatorial
    if (ecc == 0) & (0-tol > inc < tol): # +/- tolerance
        om =0
        OM = 0
        truean = lamtrue
    # Circular Inclined
    elif (ecc < zero) & (0-tol > inc > tol): # i > tolerance
        om = 0
        truean = arglat
    # Elliptical Equatorial
    elif (zero < ecc < 1) & (inc == 0 or inc == 180):
        OM = 0
        om = omtrue
    # Position vector in perifocal coordinate system
    rI = (semip*np.cos(truean))/(1 * ecc*np.cos(truean))
    rJ = (semip*np.cos(truean))/(1 * ecc*np.cos(truean))
    rK = 0
    rpqw = np.array([rI, rJ, rK])
    # Velocity vector in perifocal coordinate system
    vI = (-np.sqrt(mu/semip)*np.sin(truean))
    vJ = (np.sqrt(mu/semip)(ecc + np.cos(truean)))
    vK = 0
    vpqw = np.array([vI, vJ, vK])
    # r v in IJK coordinates
    
    
        


def sysEOM(t, state):
    ''' Generic form of a system of equations of motion for a satellite orbiting Earth '''
    mu = 398600.4415
    x = state[0]
    y = state[1]
    z = state[2]
    dx = state[3]
    dy = state[4]
    dz = state[5]
    rad = np. linalg.norm([x, y, z])
    ddx = (-mu/(rad**2))*(x/rad)
    ddy = (-mu/(rad**2))*(y/rad)
    ddz = (-mu/(rad**2))*(z/rad)
    
    return np.array([dx, dy, dz, ddx, ddy, ddz])


if __name__ == '__main__':
    mu = 398600.4415
    case1 = np.array([-3.85888212e+02,  9.46360783e+03,  0.00000000e+00, -6.51512124e+00, 1.36447314e+00,  0.00000000e+00])
    case2 = np.array([-6.70556108e+03, -7.41858818e+03,  0.00000000e+00,  3.31188428e+00, -2.99356721e+00,  4.46430533e+00])
    case3 = np.array([1.00000000e+04, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 6.31348114e+00, 0.00000000e+00])
    result1 = RV2COE(mu, case1)
    result2 = RV2COE(mu, case2)
    result3 = RV2COE(mu, case3)
    statevec = np.array([1282.09095, 6936.412744, -964.5288, -3.831107, -0.125767, -6.353187])
    result = RV2COE(mu, statevec)
    