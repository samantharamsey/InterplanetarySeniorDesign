"""
author: hannah

converts the output of MAnE to the input of POST
"""
import math as m


def main():
    # constants
    mu_s = 3.793E7                  # gravitational parameter of saturn, km3/s2
    r_t = 1.222E6                   # orbital distance of titan, km

    # inputs from mane - modify to pull in multiple
    v_inf_mane = 5.7622             # excess speed form mane, km/s
    dec_mane = m.radians(14.2513)   # declination from mane, convert deg to rad

    # calculations
    v_inf = m.sqrt(2*(v_inf_mane**2/2 + mu_s/r_t))
    vv_inf = [0, v_inf * m.cos(dec_mane), v_inf * m.sin(dec_mane)]

    # find titan's velocity somehow
    v_titan = [1, 1, 0]     # placeholder - put real values in
    total_v_t = vv_inf + v_titan


if __name__ == '__main__':
    main()
