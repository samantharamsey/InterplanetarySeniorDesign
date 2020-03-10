import pandas as pd
import numpy as np


# data is the dataframe, filepath is the path to the hdf, and filename is the name of the hdf
# body2_code is a signifier for which body is being converted to
#   0 means Saturn -> Titan
#   1 means Titan -> Saturn
def convert_reference_frame_dataframe(data, filepath, filename, body2_code):
    """
    Converts from Saturn centered reference frame to Titan centered
    *** Appends to DataFrame - will not overwrite columns ***
    """

    # pull necessary columns out of dataframe
    if body2_code == 0:
        v1_max_o = data['Saturn v1 max']
        inc_o = data['Saturn fpa (deg']
        fpa_o = data['Saturn inclination (deg)']
        v_body2 = 5.57

    elif body2_code == 1:
        v1_max_o = data['Titan v1 max']
        inc_o = data['Titan fpa (deg']
        fpa_o = data['Titan inclination (deg)']
        v_body2 = 9.68

    else:
        print("Unknown input for body code")
        return

    # wrt Saturn
    vx_max = v1_max_o * np.cos(inc_o * (np.pi / 180)) * np.sin(fpa_o * (np.pi / 180))
    vy_max = v1_max_o * np.cos(inc_o * (np.pi / 180)) * np.cos(fpa_o * (np.pi / 180))
    vz_max = v1_max_o * np.sin(inc_o * (np.pi / 180))

    # wrt Titan
    # v_inf = v1 - v_titan in the y-direction
    vy_max_converted = vy_max - v_body2

    # convert v_max and fpa_max
    for i in range(v1_max_o.shape[0]):
        v_max_converted = np.linalg.norm([vx_max[i], vy_max_converted[i], vz_max[i]])
    converted_fpa_max = np.arctan2(vx_max, vy_max_converted)

    # add to dataframe
    if body2_code == 0:
        data.insert(12, 'Titan fpa (deg)', converted_fpa_max * (180 / np.pi), True)
        data.insert(13, 'Titan v1', v_max_converted, True)

    elif body2_code == 1:
        data.insert(12, 'Saturn fpa (deg)', converted_fpa_max * (180 / np.pi), True)
        data.insert(13, 'Saturn v1', v_max_converted, True)

    # save to file
    data.to_hdf(filepath + filename + r'.hdf', key='df')
    return data


# v1_max_o is the original v1_max, inc_o is the original inclination, and fpa_o is the original fpa
# body2_code is a signifier for which body is being converted to
#   0 means Saturn -> Titan
#   1 means Titan -> Saturn
def convert_reference_frame_single_input(v1_max_o, inc_o, fpa_o, body2_code):
    """
    Converts from Saturn centered reference frame to Titan centered
    *** Appends to DataFrame - will not overwrite columns ***
    """

    # pull necessary columns out of dataframe
    if body2_code == 0:
        v_body2 = 5.57

    elif body2_code == 1:
        v_body2 = 9.68

    else:
        print("Unknown input for body code")
        return

    # wrt Saturn
    vx_max = v1_max_o * np.cos(inc_o * (np.pi / 180)) * np.sin(fpa_o * (np.pi / 180))
    vy_max = v1_max_o * np.cos(inc_o * (np.pi / 180)) * np.cos(fpa_o * (np.pi / 180))
    vz_max = v1_max_o * np.sin(inc_o * (np.pi / 180))

    # wrt Titan
    # v_inf = v1 - v_titan in the y-direction
    vy_max_converted = vy_max - v_body2

    # convert v_max and fpa_max
    v_max_converted = np.linalg.norm([vx_max, vy_max_converted, vz_max])
    converted_fpa_max = np.arctan2(vx_max, vy_max_converted)

    # add to dataframe
    if body2_code == 0:
        print('Titan fpa (deg)', converted_fpa_max * (180 / np.pi))
        print('Titan v1', v_max_converted)

    elif body2_code == 1:
        print('Saturn fpa (deg)', converted_fpa_max * (180 / np.pi))
        print('Saturn v1', v_max_converted)
