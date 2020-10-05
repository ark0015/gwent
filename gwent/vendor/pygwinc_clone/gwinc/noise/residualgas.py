from __future__ import division
from numpy import sqrt, log, pi

from .. import const


def gas(f, ifo):
    """Residual gas noise in arm cavities

    The following function models the noise spectrum caused by the
    passage of residual gas molecules through the laser beams in the
    arms of the interferometer.

    The method used here is presented by Rainer Weiss, Micheal
    E. Zucker, and Stanley E. Whitcomb in their paper Optical
    Pathlength Noise in Sensitive Interferometers Due to Residual Gas.

    Added to Bench by Zhigang Pan, Summer 2006
    Cleaned up by PF, Apr 07
    Eliminated numerical integration and substituted first order
    expansion of exp, to speed it up.

    """
    L      = ifo.Infrastructure.Length
    Lambda = ifo.Laser.Wavelength
    k      = const.kB
    T      = ifo.Infrastructure.Temp
    P      = ifo.Infrastructure.ResidualGas.pressure            # Pressure inside the vacuum
    M      = ifo.Infrastructure.ResidualGas.mass
    R1     = ifo.Optics.Curvature.ITM                           # Radius of curvature of ITM
    R2     = ifo.Optics.Curvature.ETM                           # Radius of curvature of ETM
    alpha  = ifo.Infrastructure.ResidualGas.polarizability      #

    rho = P /(k*T)                                           # number density of Gas
    v0  = sqrt(2*k*T / M)                                    # Mean speed of Gas

    g1 = 1 - L/R1                                 #  first resonator g-parameter of the ITM
    g2 = 1 - L/R2                                 # second resonator g-parameter of the ETM
    waist = L * Lambda / pi
    waist = waist * sqrt(((g1*g2)*(1-g1*g2))/((g1+g2-2*g1*g2)**2))
    waist = sqrt(waist)                             # Gaussian beam waist size
    zr = pi*waist**2 / Lambda                       # Rayleigh range
    z1 = -((g2*(1-g1))/(g1+g2-2*g1*g2))*L           # location of ITM relative to the waist
    z2 =  ((g1*(1-g2))/(g1+g2-2*g1*g2))*L           # location of ETM relative to the waist

    # The exponential of Eq. 1 of P940008 is expanded to first order; this
    # can be integrated analytically
    zint = log(z2 + sqrt(z2**2 + zr**2)) - log(z1 + sqrt(z1**2 + zr**2))
    zint = zint * zr/waist
    zint = zint - 2*pi*L*f/v0
    # optical path length for one arm
    zint = zint*((4*rho*(2*pi*alpha)**2)/v0)
    # eliminate any negative values due to first order approx.
    zint[zint < 0] = 0
    # account for both arms & turn into strain noise power
    n = 2 * zint / L**2
    return n
