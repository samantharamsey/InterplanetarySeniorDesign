import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import fsolve
from matplotlib import animation

mu = 37.931 * 10 ** 6  # km^3/s^2
r_encel = 238037
r_titan = 1186780.3668940281
r_saturn = 58232


def main():

    # calculate_v1(r_encel, 360)
    calculate_v1(r_saturn + 1000, 360)


def calculate_v1(rp, n):

    data = pd.DataFrame([])
    ax = plt.subplot(111, polar=True)
    v1 = [0]*n
    bad = []
    good = []

    r1 = r_titan
    a = (r1 + rp)/2
    e = 1 + rp/a
    E = -mu / (2 * a)
    vp = np.sqrt(2*(E + mu/rp))
    h = vp * rp
    v_esc_p = np.sqrt(2 * (mu/rp))
    v_esc_1 = np.sqrt(2 * (mu/r1))
    # print("escape v at r1", v_esc_1)

    for i in range(n):
        gamma = np.radians(i)
        v1[i] = h / (r1 * np.cos(gamma))
        if abs(v1[i]) > v_esc_1:
            print("bad v", v1[i], "at", np.degrees(gamma))
            bad.append(np.degrees(gamma))
        else:
            good.append(np.degrees(gamma))

        if (i != 90) & (i != 270):
            data = data.append(pd.DataFrame({'v1': v1[i], 'gamma (rad)': gamma}, index=[0]), ignore_index=True)
        # else:
            # print("bad", i)

    ax.plot(data['gamma (rad)'][::], data['v1'][::])

    print(v1)
    print(len(good))
    print(max(v1))
    #print(bad)
    plt.show()


if __name__ == '__main__':
    main()
