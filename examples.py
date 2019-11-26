import numpy as n


def main():

    mu_sun = 1.327 * 10**11
    period_mercury = 87.96 * 24 * 60**2
    a_venus = 1.081608 * 10**8
    mu_venus = 3.24859 * 10**5
    v_venus = 35
    a_earth = 1.496 * 10**8

    """example 1: Earth -> Hohmann transfer -> Venus flyby -> orbit around Mercury"""
    # find the semi major axis for the heliocentric orbit after Venus flyby
    period_post_flyby = 5/2 * period_mercury
    a_post_flyby = n.cbrt(mu_sun * (period_post_flyby / (2 * n.pi))**2)
    print("post flyby a =", a_post_flyby)

    # find the turn angle at Venus necessary for the heliocentric orbit after Venus flyby
    v_post_flyby = n.sqrt(mu_sun * (2/a_venus - 1/a_post_flyby))
    print("post flyby v =", v_post_flyby)
    v_perihelion = n.sqrt(mu_sun * (2/a_venus - 1/(a_venus + a_earth)))
    print("perihelion v =", v_perihelion)
    v_inf = v_perihelion - v_venus
    print("v infinity =", v_inf)
    delta = n.pi - n.arccos((v_post_flyby**2 - v_venus**2 - v_inf**2)/(-2 * v_venus * v_inf))
    print("turn angle =", n.degrees(delta))

    # find the flyby altitude to get this turn angle
    rp = mu_venus/v_inf**2 * (1/(n.sin(delta/2)) - 1)
    print("rp =", rp)


if __name__ == '__main__':
    main()
