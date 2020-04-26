import pandas as pd
import math as m
import numpy as np
import matplotlib.pyplot as plt


mu = 3.795*10**7
r = 1.22*10**6
v_titan = 5.57


def chop(num, digits) -> float:
    step = 10.0 ** digits
    return m.trunc(step*num) / step


def calc_titan_vel(vx1, vy1, vz1, intercept):
    vxt = 0
    vyt = 0
    vzt = 0

    if 0 < intercept < 90 or 90 < intercept < 180 or 180 < intercept < 270 or 270 < intercept < 360:
        vxt = v_titan * m.sin(m.radians(intercept))
        vyt = v_titan * m.cos(m.radians(intercept))
    elif intercept == 90:
        vxt = v_titan
    elif intercept == 180:
        vyt = -v_titan
    elif intercept == 270:
        vxt = -v_titan
    elif intercept == 360 or intercept == 0:
        vyt = v_titan
    else:
        print('error: intercept location')
        breakpoint()

    v_inf_x = vx1 + vxt
    v_inf_y = vy1 + vyt
    v_inf_z = vz1 + vzt
    v_inf_mag = m.sqrt(v_inf_x**2 + v_inf_y**2 + v_inf_z**2)

    return v_inf_x, v_inf_y, v_inf_z, v_inf_mag


def calc_vinf_for_post(v_inf, dec, intercept):
    # calculate probe velocity at titans orbital distance
    v1 = np.sqrt(2*((v_inf**2/2) + (mu/r)))

    # separate into components
    vx1 = 0
    vy1 = v1*np.cos(m.radians(dec))
    vz1 = v1*np.sin(m.radians(dec))

    return calc_titan_vel(vx1, vy1, vz1, intercept)


def post_to_patch(POST_filepath, POST_filename, patch_filepath, patch_filename):
    """
    POST -> Patch
    Takes output from POST and compares it to
    """

    # load in the data from POST
    post_data = pd.read_excel(POST_filepath + POST_filename + r'.xlsx')

    vx1 = post_data['vxi'].iloc[-1]/1000
    vy1 = post_data['vyi'].iloc[-1]/1000
    vz1 = post_data['vzi'].iloc[-1]/1000
    intercept = post_data['trunmx'].iloc[-1]

    # load in the final orbit data from the script calculations
    script_data = pd.read_hdf(patch_filepath + patch_filename + r'.hdf')
    # pull out the titan-centered velocity
    v1 = script_data['Titan v1']
    # print(script_data.columns)

    [v_inf_x, v_inf_y, v_inf_z, v_inf_mag] = calc_titan_vel(vx1, vy1, vz1, intercept)

    # compare the output of POST to the output of the script
    v_count = 0
    tot = 0
    min = 1000
    max = 0
    j = 3
    v_inf_mag = chop(v_inf_mag, j)

    for i in range(1, script_data.shape[0]):
        # print(v_inf_mag, chop(v1.loc[i], j))
        # check velocities after truncating decimal places for the script output
        if v_inf_mag == chop(v1.loc[i], j):
            v_count += 1
        # this isn't relevant, just a check
        else:
            diff = m.fabs(v_inf_mag - chop(v1.loc[i], j))
            tot += diff
        if diff > max:
            max = diff
        elif diff < min:
            min = diff

    avg = tot / (i + 1)
    print("average difference in velocities", avg)
    print("minimum difference in velocities", min)
    print("maximum difference in velocities", max)

    # print out the stats on the POST-script comparison
    print("number of cases with matching velocities:", v_count)


def mane_to_post(filepath, filename):
    """
    MAnE -> POST
    Turn MAnE output to corresponding input for POST
    Takes filepath and filename for the output CSV
    Converts declination and Saturn-centered velocity to Titan-centered velocity
    Uses calc_vinf_for_post
    """

    data = pd.DataFrame([])
    v_infx = []
    v_infy = []
    v_infz = []
    for i in range(360):
        resultx, resulty, resultz, resultmag = calc_vinf_for_post(5.7622, 14.2513, i)
        v_infx.append(resultx)
        v_infy.append(resulty)
        v_infz.append(resultz)

        data = data.append(pd.DataFrame({'intercept location (deg)': i,
                                         'v_infinity wrt Titian magnitude (km/s)': resultmag,
                                         'v_infinity x component': resultx,
                                         'v_infinity y component': resulty,
                                         'v_infinity z component': resultz},
                                        index=[0]), ignore_index=True)
    # send results to excel
    data.to_csv(filepath + filename + r'.csv', index=False)
    # send results to HDF5 - faster loading
    data.to_hdf(filepath + filename + r'.hdf', key='df')

    # plot stuff
    fig1 = plt.figure()
    plt.plot(v_infx)
    plt.plot(v_infy)
    plt.plot(v_infz)
    plt.ylabel('velocity (km/s)')
    plt.xlabel('intercept location (deg)')
    plt.legend(['x component', 'y component', 'z component'])
    plt.title('V infinity wrt Titan: v_inf wrt Saturn = 6, dec = 10')
    plt.show()


"""
Main Function
Uncomment what you want to run
"""
if __name__ == '__main__':
    # mane_to_post( r'C:\Spice_Kernels', r'\10_declination')
    post_to_patch(r'C:\Spice_Kernels', r'\tin7_2', r'C:\Spice_Kernels',  r'\3D_potato')
