from __future__ import division

import numpy as np
import scipy.integrate as scint

from numpy import pi, sqrt, exp
from .seismic import seisNLNM
from .. import const


def gravg(f, ifo):
    """Return estimate of newtonian noise contribribution to |h(f)|^2

    N = GRAVG(F, IFO) returns the gravity gradient (Newtonian) noise
    contribution in strain^2/Hz for four mirrors.

    References:

     Saulson 1984,           http://dx.doi.org/10.1103/PhysRevD.30.732
     Hughes and Thorne 1998, http://dx.doi.org/10.1103/PhysRevD.58.122002

     Driggers and Harms 2011, ``Results of Phase 1 Newtonian Noise
     Measurements at the LIGO Sites,'' February-March 2011.  T1100237.
     https://dcc.ligo.org/cgi-bin/private/DocDB/ShowDocument?docid=60064

    Written by Enrico Camagna (?)

    added to Bench by Gregg Harry 8/27/03
    seismic spectrum modified by Jan Harms 05/11/2010
    Calculates gravity gradient noise for four mirrors

    """

    fk = ifo.Seismic.KneeFrequency
    a = ifo.Seismic.LowFrequencyLevel
    L = ifo.Infrastructure.Length
    gamma = ifo.Seismic.Gamma
    ggcst = const.G
    rho = ifo.Seismic.Rho
    # factor to account for correlation between masses
    # and the height of the mirror above the ground
    beta = ifo.Seismic.Beta
    h = ifo.Seismic.TestMassHeight
    c_rayleigh = ifo.Seismic.RayleighWaveSpeed

    if 'Omicron' in ifo.Seismic:
        omicron = ifo.Seismic.Omicron
    else:
        omicron = 1

    # a sort of theta function (Fermi distr.)
    coeff = 3**(-gamma*f)/(3**(-gamma*f) + 3**(-gamma*fk))

    # modelization of seismic noise (vertical)
    ground = a*coeff + a*(1-coeff)*(fk/f)**2
    if 'Site' in ifo.Seismic and ifo.Seismic.Site == 'LLO':
        ground = a*coeff*(fk/f) + a*(1-coeff)*(fk/f)**2

    # effective GG spring frequency, with G gravitational
    fgg = sqrt(ggcst * rho) / (2*pi)

    # fixed numerical factors, 5/9/06, PF
    n = (beta*4*pi/L*(fgg**2/f**2)*ground)**2

    # The following two terms are corrections due to Jan Harms
    # https://git.ligo.org/rana-adhikari/CryogenicLIGO/issues/45
    # (1) projection of NN force onto the direction of the arm
    n = n * 1/2
    # (2) exponential cutoff at frequency (seismic speed)/(test mass height)
    n = n * exp(-4*pi*f*h/c_rayleigh)

    # Feedforward cancellation
    n /= (omicron**2)

    return n * ifo.gwinc.sinc_sqr


def gravg_rayleigh(f, ifo):
    """Gravity gradient noise for seismic Rayleigh waves
    Following Harms LRR: https://doi.org/10.1007/lrr-2015-3

    """
    fk = ifo.Seismic.KneeFrequency
    a = ifo.Seismic.LowFrequencyLevel
    L = ifo.Infrastructure.Length
    gamma = ifo.Seismic.Gamma
    ggcst = const.G
    rho = ifo.Seismic.Rho
    h = ifo.Seismic.TestMassHeight
    c_rayleigh = ifo.Seismic.RayleighWaveSpeed

    if 'Omicron' in ifo.Seismic:
        omicron = ifo.Seismic.Omicron
    else:
        omicron = 1

    # a sort of theta function (Fermi distr.)
    coeff = 3**(-gamma*f)/(3**(-gamma*f) + 3**(-gamma*fk))

    # modelization of seismic noise (vertical)
    ground = a*coeff + a*(1-coeff)*(fk/f)**2
    if 'Site' in ifo.Seismic and ifo.Seismic.Site == 'LLO':
        ground = a*coeff*(fk/f) + a*(1-coeff)*(fk/f)**2

    # Harms LRR eqs. 35, 96, and 98
    w = 2 * pi * f
    k = w / c_rayleigh
    kP = w / ifo.Seismic.pWaveSpeed
    kS = w / ifo.Seismic.sWaveSpeed
    qzP = sqrt(k**2 - kP**2)
    qzS = sqrt(k**2 - kS**2)
    zeta = sqrt(qzP / qzS)

    gnu = k * (1 - zeta) / (qzP - k * zeta)

    n = 2 * (2 * pi * ggcst * rho * exp(-h * k) * gnu)**2 * ground**2 / w**4

    n /= omicron**2

    return n * ifo.gwinc.dhdl_sqr


def gravg_pwave(f, ifo):
    """Gravity gradient noise for seismic p-waves
    Following Harms LRR: https://doi.org/10.1007/lrr-2015-3

    """
    ggcst = const.G
    cP = ifo.Seismic.pWaveSpeed
    levelP = ifo.Seismic.pWaveLevel
    kP = (2 * pi * f) / cP

    rho_ground = ifo.Seismic.Rho
    psd_ground_pwave = (levelP * seisNLNM(f))**2

    tmheight = ifo.Seismic.TestMassHeight 
    xP = np.abs(kP * tmheight)

    if tmheight >= 0:
        # Surface facility
        # The P-S conversion at the surface is not implemented
        height_supp_power = (3 / 2) * np.array([scint.quad(lambda th, x: np.sin(th)**3
                * np.exp(-2 * x * np.sin(th)), 0, pi / 2, args=(x,))[0]
                for x in xP])
    else:
        # Underground facility
        # The cavity effect is not included
        height_supp_power = (3 / 4) * np.array([scint.quad(lambda th, x: np.sin(th)**3
                * (2 - np.exp(-x * np.sin(th)))**2, 0, pi, args=(x,))[0]
                for x in xP])
    psd_gravg_pwave = ((2 * pi * ggcst * rho_ground)**2
            * psd_ground_pwave * height_supp_power)
    psd_gravg_pwave *= 4 / (2 * pi * f)**4
    return psd_gravg_pwave * ifo.gwinc.dhdl_sqr


def gravg_swave(f, ifo):
    """Gravity gradient noise for seismic s-waves
    Following Harms LRR: https://doi.org/10.1007/lrr-2015-3

    """
    ggcst = const.G
    cS = ifo.Seismic.sWaveSpeed
    levelS = ifo.Seismic.sWaveLevel
    kS = (2 * pi * f) / cS

    rho_ground = ifo.Seismic.Rho
    psd_ground_swave = (levelS * seisNLNM(f))**2

    tmheight = ifo.Seismic.TestMassHeight 
    xS = np.abs(kS * tmheight)

    # For both surface and underground facilities
    height_supp_power = (3 / 2) * np.array([scint.quad(lambda th, x: np.sin(th)**3
            * np.exp(-2 * x * np.sin(th)), 0, pi / 2, args=(x,))[0]
            for x in xS])
    psd_gravg_swave = ((2 * pi * ggcst * rho_ground)**2
            * psd_ground_swave * height_supp_power)
    psd_gravg_swave *= 4 / (2 * pi * f)**4

    return psd_gravg_swave * ifo.gwinc.dhdl_sqr


def atmois(f, ifo):
    """
    Atmospheric infrasound calculation following Harms LRR.

    """
    p_air = ifo.Atmospheric.AirPressure
    rho_air = ifo.Atmospheric.AirDensity
    ai_air = ifo.Atmospheric.AdiabaticIndex
    c_sound = ifo.Atmospheric.SoundSpeed

    L = ifo.Infrastructure.Length
    ggcst = const.G
    h = ifo.Seismic.TestMassHeight

    w = 2 * pi * f
    k = w / c_sound

    # Pressure spectrum
    try:
        a_if = ifo.Atmospheric.InfrasoundLevel1Hz
        e_if = ifo.Atmospheric.InfrasoundExponent
        psd_if = (a_if * f**e_if)**2
    except AttributeError:
        psd_if = atmoBowman(f)**2

    # Harms LRR (2015), eq. 172
    # https://doi.org/10.1007/lrr-2015-3
    # With an extra factor 2 for two arms
    # And with the Bessel terms ignored... for 4 km this amounts to a 10%
    # correction at 10 Hz and a 30% correction at 1 Hz
    coupling_if = 4./3 * (4 * pi / (k * w**2) * ggcst * rho_air / (ai_air * p_air))**2

    n_if = coupling_if * psd_if

    return n_if * ifo.gwinc.dhdl_sqr


def atmoBowman(f):
    """The Bowman infrasound model
   
    """ 
    freq = np.array([
                0.01, 0.0155, 0.0239, 0.0367, 0.0567,
                0.0874, 0.1345, 0.2075, 0.32, 0.5,
                0.76, 1.17, 1.8, 2.79, 4.3,
                6.64, 10, 100])
    pressure_asd = np.sqrt([22.8, 4, 0.7, 0.14, 0.027, 0.004,
                            0.0029, 0.0039, 7e-4, 1.44e-4, 0.37e-4,
                            0.12e-4, 0.56e-5, 0.35e-5, 0.26e-5, 0.24e-5,
                            2e-6, 2e-6])
    return 10**(np.interp(np.log10(f), np.log10(freq), np.log10(pressure_asd)))
