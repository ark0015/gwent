from __future__ import division
import numpy as np
from numpy import log10
from scipy.interpolate import PchipInterpolator as interp1d


def seismic(f, ifo):
    """Seismic noise.
    """
    return seismicAll(f, ifo)[0]


def seismicAll(f, ifo):
    """Seismic noise.

    Return (noise, noise_vertical, noise_horizontal)

    """
    hTable = ifo.Suspension.hTable
    vTable = ifo.Suspension.vTable

    theta = ifo.Suspension.VHCoupling.theta

    # noise input, horizontal and vertical
    if 'PlatformMotion' in ifo.Seismic:
        if ifo.Seismic.PlatformMotion == 'BSC':
            nt, nr = seisBSC(f)
        elif ifo.Seismic.PlatformMotion == '6D':
            nt, nr = seis6D(f)
        else:
            nt, nr = seisBSC(f)
    else:
        nt, nr = seisBSC(f)

    # horizontal noise total
    nh = (abs(hTable)**2) * nt**2

    # vertical noise total
    nv = (abs(theta * vTable)**2) * nt**2

    # new total noise
    n = nv + nh

    # Convert into Strain PSD (4 TMs)
    nh *= 4 * ifo.gwinc.dhdl_sqr
    nv *= 4 * ifo.gwinc.dhdl_sqr
    n *= 4 * ifo.gwinc.dhdl_sqr

    return n, nh, nv


def seisBSC(f):
    """Rough ISI noise source spectra.

    Returns ISI (translational, rotational) DOFs

    """
    SEI_F = np.array([0.01, 0.03, 0.1, 0.2, 0.5, 1, 10, 30, 300])

    # translational DOFs
    SEI_T = np.array([3e-6, 1e-6, 2e-7, 2e-7, 8e-10, 1e-11, 3e-13, 3e-14, 3e-14])
    nt = 10**(interp1d(SEI_F, log10(SEI_T))(f))

    # rotational DOFs
    SEI_R = np.array([1e-8, 3e-8, 2e-8, 1e-8, 4e-10, 1e-11, 3e-13, 3e-14, 3e-14])
    nr = 10**(interp1d(SEI_F, log10(SEI_R))(f))

    return nt, nr


def seis6D(f):
    """ISI noise source spectra with a 6D seismometer.

    This largely follows Mow-Lowry and Martynov, arXiv:1801.01468.

    """

    SEI_F = np.array([0.01, 0.03, 0.1, 0.2, 0.5, 1, 10, 100, 300])

    SEI_T_self = np.array([1e-7, 1e-9, 3e-11, 6e-12, 3e-13, 1e-13, 3e-14, 1e-14, 1e-14])
    nt_self = 10**(interp1d(SEI_F, log10(SEI_T_self))(f))
    nt_gnd = 10*seisNLNM(f)
    blend_t = np.abs(100/(1+1j*f/0.01)**4)
    nt = np.sqrt(nt_self**2 + (blend_t * nt_gnd)**2)

    SEI_R_self = np.array([2e-11, 5e-12, 1e-12, 6e-13, 3e-13, 2e-13, 6e-14, 2e-14, 2e-14])
    nr_self = 10**(interp1d(SEI_F, log10(SEI_R_self))(f))
    nr_gnd = np.abs(1e-7/(1+1j*f/0.001))
    blend_r = np.abs(100/(1+1j*f/0.01)**4)
    nr = np.sqrt(nr_self**2 + (blend_r * nr_gnd)**2)

    return nt, nr


def seisNLNM(f):
    """The Peterson New Low-Noise Model.

    Returns a displacement ASD.

    """
    Pl = np.array([
       1.00e-02, 1.00e-01, 1.70e-01, 4.00e-01, 8.00e-01, 1.24e+00,
       2.40e+00, 4.30e+00, 5.00e+00, 6.00e+00, 1.00e+01, 1.20e+01,
       1.56e+01, 2.19e+01, 3.16e+01, 4.50e+01, 7.00e+01, 1.01e+02,
       1.54e+02, 3.28e+02, 6.00e+02, 1.00e+04])
    Al = np.array([
       -156.72, -162.36, -166.7 , -170.  , -166.4 , -168.6 , -159.98,
       -141.1 ,  -71.36,  -97.26, -132.18, -205.27,  -37.65, -114.37,
       -160.58, -187.5 , -216.47, -185.  , -168.34, -217.43, -258.28,
       -346.88])
    Bl = np.array([
          5.64,    5.64,    0.  ,   -8.3 ,   28.9 ,   52.48,   29.81,
          0.  ,  -99.77,  -66.49,  -31.57,   36.16, -104.33,  -47.1 ,
        -16.28,    0.  ,   15.7 ,    0.  ,   -7.61,   11.9 ,   26.6 ,
         48.75])
    nlnm = 10**(np.interp(1/f, Pl, Al+Bl*np.log10(Pl))/20) / (2 * np.pi * f)**2
    return nlnm


def seisNHNM(f):
    """The Peterson New High-Noise Model.

    Returns a displacement ASD.

    """

    Pl = np.array([
        1.00e-01, 2.20e-01, 3.20e-01, 8.00e-01, 3.80e+00,
        4.60e+00, 6.30e+00, 7.90e+00, 1.54e+01, 2.00e+01,
        3.54e+02,
        ])  
    Al = np.array([
        -108.73, -150.34, -122.31, -116.85, -108.48,
         -74.66,    0.66,  -93.37,   73.54, -151.52,
        -206.66,
        ])
    Bl = np.array([
         -17.23,  -80.50,  -23.87,   32.51,   18.08,
         -32.95, -127.18,  -22.42, -162.98,   10.01,
          31.63,
        ])
    nhnm = 10**(np.interp(1/f, Pl, Al+Bl*np.log10(Pl))/20) / (2 * np.pi * f)**2
    return nhnm
