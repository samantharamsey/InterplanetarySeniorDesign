import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# load pre-existing file into the dataframe
# UPDATE FILEPATH BEFORE RUNNING
filepath = r'C:\Users\saman\OneDrive\Desktop\InterplanetarySeniorDesign'
filename = r'\3D_potato'
data = pd.read_hdf(filepath + filename + r'.hdf')

# define some constants
mu = 37.931 * 10 ** 6  # saturn's gravitational parameter
saturn_equatorial = 60268  # km
saturn_polar = 54364  # km
r_titan = 1.2 * 10 ** 6  # km
r_encel = 238000  # km
v_titan = 5.57  # km/s


def convert_reference_frame(v1_max, inc, fpa, v_body2, body2_code):
    """
    Converts from Saturn centered reference frame to Titan centered
    *** Appends to DataFrame - will not overwrite columns ***
    """

    # wrt Saturn
    vx_max = v1_max * np.cos(inc * (np.pi / 180)) * np.sin(fpa * (np.pi / 180))
    vy_max = v1_max * np.cos(inc * (np.pi / 180)) * np.cos(fpa * (np.pi / 180))
    vz_max = v1_max * np.sin(inc * (np.pi / 180))

    # wrt Titan
    # v_inf = v1 - v_titan in the y-direction
    vy_max_converted = vy_max - v_body2

    for i in range(v1_max.shape[0]):
        v_max_converted = np.linalg.norm([vx_max[i], vy_max_converted[i], vz_max[i]])
    converted_fpa_max = np.arctan2(vx_max, vy_max_converted)

    if body2_code == 0:
        data.insert(12, 'Titan fpa (deg)', converted_fpa_max * (180 / np.pi), True)
        data.insert(13, 'Titan v1', v_max_converted, True)

    elif body2_code == 1:
        data.insert(12, 'Saturn fpa (deg)', converted_fpa_max * (180 / np.pi), True)
        data.insert(13, 'Saturn v1', v_max_converted, True)

    else:
        print("Invalid body inputted")

    data.to_hdf(filepath + filename + r'.hdf', key='df')
