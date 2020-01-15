import pandas as pd
import numpy as np
from numpy import linalg as la
import matplotlib.pyplot as plt
from mayavi import mlab


def main():

    # constants
    r_encel = 238000
    mu = 3.794*10**7

    # initialize
    r = 238001
    v = 7944

    for gamma in np.linspace(0, 2*np.pi):
        for i in np.linspace(0, 2*np.pi):
            while r > r_encel:
                v -= 0.1
                """
                a = mu / (2 * (v**2/2 - mu/r))
                if a < 0:
                    break
                print("a", a)

                h = r*v*np.cos(gamma)
                p = h**2 / mu
                print("p", p)
                e = np.sqrt(1 - p/a)
                theta = np.arccos((p/r - 1) / e)

                # get v & r vectors in PQW
                v_P = np.sqrt(mu/p) * np.sin(theta)
                v_Q = np.sqrt(mu/p) * (e + np.cos(theta))
                r_P = r * np.cos(theta)
                r_Q = r * np.sin(theta)

                # calculate w

                # calculate nu
                # nu = np.pi - w

                # calculate r
                # r = (a * (1 - e**2)) / (1 + e * (np.cos(nu)))
                r = np.sqrt(r_P**2 + r_Q**2)
                print(r)
                """

                vv = v * np.array([np.sin(gamma), np.cos(gamma), np.sin(i)])
                rr = r * np.array([1, 0, 0])
                [p, e, new_i, big_omega, little_omega, nu] = rv_to_coe(rr, vv, mu)
                if e > 1:
                    break

                print("p =", p, "e =", e, "i =", np.degrees(new_i), "W = ", np.degrees(big_omega), "w =",
                      np.degrees(little_omega), "v =", np.degrees(nu))
                r = p / (1 + e * np.cos(nu))

    # Examples to verify functions
    """
    rv_to_coe(np.array([2, 0, 0]), np.array([0, 1, 0]), 1)
    rv_to_coe(np.array([3*np.sqrt(3)/4, 3/4, 0]), np.array([-1/(2*np.sqrt(2)), np.sqrt(3)/(2*np.sqrt(2)), 1/np.sqrt(2)]), 1)
    coe_to_rv(2.25, 0.5, np.radians(45), np.radians(30), 0, 0, 1)
    """


# BMW 2.5, 2.6
def coe_to_rv(p, e, i, big_omega, little_omega, nu, mu):
    r = p / (1 + e*np.cos(nu))

    rr_pqw = np.array([r * np.cos(nu), r * np.sin(nu), 0])
    vv_pqw = np.sqrt(mu/p) * np.array([-np.sin(nu), (e + np.cos(nu)), 0])
    print("r PQW =", rr_pqw, "v PQW =", vv_pqw)

    cW = np.cos(big_omega)
    cw = np.cos(little_omega)
    ci = np.cos(i)
    sW = np.sin(big_omega)
    sw = np.sin(little_omega)
    si = np.sin(i)

    r11 = cW*cw - sW*sw*ci
    r12 = -cW*sw - sW*cw*ci
    r13 = sW*si
    r21 = sW*cw + cW*sw*ci
    r22 = -sW*sw + cW*cw*ci
    r23 = -cW*si
    r31 = sw*si
    r32 = cw*si
    r33 = ci

    tran_matrix = np.array([[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]])

    rr = np.matmul(tran_matrix, rr_pqw)
    vv = np.matmul(tran_matrix, vv_pqw)

    print("r =", rr, "\nv =", vv)
    return rr, vv


# BMW 2.4
def rv_to_coe(rr, vv, mu):
    r = la.norm(rr)
    v = la.norm(vv)

    ii = np.array([1, 0, 0])
    kk = np.array([0, 0, 1])

    hh = np.cross(rr, vv)
    h = la.norm(hh)
    nn = np.cross(kk, hh)
    n = la.norm(nn)
    ee = 1/mu * ((v**2 - mu/r)*rr - np.dot(rr, vv)*vv)

    # a = -mu / (2*(v**2/2 - mu/r))
    p = h**2/mu
    e = la.norm(ee)
    i = np.arccos(np.dot(hh, kk) / h)

    big_omega = np.arccos(np.dot(ii, nn) / n)
    if np.sin(big_omega) < 0:
        big_omega = 2*np.pi - big_omega

    little_omega = np.arccos(np.dot(nn, ee) / (n*e))
    if ee[2] < 0:
        little_omega = 2*np.pi - little_omega

    nu = np.pi - little_omega
    """
    nu = np.arccos(np.dot(ee, rr) / (e*r))
    if np.dot(rr, vv) < 0:
        nu = 2*np.pi - nu

    u = np.arccos(np.dot(nn, rr) / (n*r))
    if rr[2] < 0:
        u = 2*np.pi - u
    """

    return p, e, i, big_omega, little_omega, nu


if __name__ == '__main__':
    main()
