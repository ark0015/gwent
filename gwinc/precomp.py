from __future__ import division
from numpy import pi, sqrt, sin, exp
from scipy.io import loadmat
import logging
import copy

from . import const
from . import suspension
from .struct import Struct
from .noise.coatingthermal import getCoatDopt


def precompIFO(f, ifoin, PRfixed=True):
    """Add precomputed data to the IFO model.

    To prevent recomputation of these precomputed data, if the
    ifo argument contains ifo.gwinc.PRfixed, and this matches
    the argument PRfixed, no changes are made.

    This function DOES NOT modify IFO. Its return value is a copy of
    IFO with precomputations filled out.

    """
    ifo = copy.deepcopy(ifoin)

    if 'gwinc' not in ifo:
        ifo.gwinc = Struct()

    ifo.gwinc.PRfixed = PRfixed

    ##############################
    # derived temp

    if 'Temp' not in ifo.Materials.Substrate:
        ifo.Materials.Substrate.Temp = ifo.Infrastructure.Temp

    ##############################
    # suspension vertical-horizontal coupling

    if 'VHCoupling' not in ifo.Suspension:
        ifo.Suspension.VHCoupling = Struct()
        ifo.Suspension.VHCoupling.theta = ifo.Infrastructure.Length / const.R_earth

    ##############################
    # optics values

    # calculate optics' parameters
    ifo.Materials.MirrorVolume = pi*ifo.Materials.MassRadius**2 * \
                                 ifo.Materials.MassThickness
    ifo.Materials.MirrorMass = ifo.Materials.MirrorVolume * \
                               ifo.Materials.Substrate.MassDensity
    ifo.Optics.ITM.Thickness = ifo.Materials.MassThickness

    # coating layer optical thicknesses - mevans 2 May 2008
    if 'CoatLayerOpticalThickness' not in ifo.Optics.ITM:
        T = ifo.Optics.ITM.Transmittance
        dL = ifo.Optics.ITM.CoatingThicknessLown
        dCap = ifo.Optics.ITM.CoatingThicknessCap
        ifo.Optics.ITM.CoatLayerOpticalThickness = getCoatDopt(ifo, T, dL, dCap=dCap)
        T = ifo.Optics.ETM.Transmittance
        dL = ifo.Optics.ETM.CoatingThicknessLown
        dCap = ifo.Optics.ETM.CoatingThicknessCap
        ifo.Optics.ETM.CoatLayerOpticalThickness = getCoatDopt(ifo, T, dL, dCap=dCap)

    ##############################
    # beam parameters

    armlen = ifo.Infrastructure.Length

    # g-factors
    g1 = 1 - armlen / ifo.Optics.Curvature.ITM
    g2 = 1 - armlen / ifo.Optics.Curvature.ETM
    gcav = sqrt(g1 * g2 * (1 - g1 * g2))
    gden = g1 - 2 * g1 * g2 + g2

    if (g1 * g2 * (1 - g1 * g2)) <= 0:
        raise Exception('Unstable arm cavity g-factors.  Change ifo.Optics.Curvature')
    elif gcav < 1e-3:
        logging.warning('Nearly unstable arm cavity g-factors.  Reconsider ifo.Optics.Curvature')

    ws = sqrt(armlen * ifo.Laser.Wavelength / pi)
    w1 = ws * sqrt(abs(g2) / gcav)
    w2 = ws * sqrt(abs(g1) / gcav)

    # waist size
    w0 = ws * sqrt(gcav / abs(gden))
    zr = pi * w0**2 / ifo.Laser.Wavelength
    z1 = armlen * g2 * (1 - g1) / gden
    z2 = armlen * g1 * (1 - g2) / gden

    ifo.Optics.ITM.BeamRadius = w1
    ifo.Optics.ETM.BeamRadius = w2

    ##############################
    # calc power and IFO parameters

    pbs, parm, finesse, prfactor, Tpr = precompPower(ifo, PRfixed)

    ifo.gwinc.pbs = pbs
    ifo.gwinc.parm = parm
    ifo.gwinc.finesse = finesse
    ifo.gwinc.prfactor = prfactor
    ifo.gwinc.gITM = g1
    ifo.gwinc.gETM = g2
    ifo.gwinc.BeamWaist = w0
    ifo.gwinc.BeamRayleighRange = zr
    ifo.gwinc.BeamWaistToITM = z1
    ifo.gwinc.BeamWaistToETM = z2
    ifo.Optics.PRM.Transmittance = Tpr

    ##############################
    # saved seismic spectrum

    if 'darmSeiSusFile' in ifo.Seismic and ifo.Seismic.darmSeiSusFile:
        darmsei = loadmat(ifo.Seismic.darmSeiSusFile)
        ifo.Seismic.darmseis_f = darmsei['darmseis_f'][0]
        ifo.Seismic.darmseis_x = darmsei['darmseis_x'][0]

    # --------------------------------------------------------
    # Suspension
    # if the suspension code supports different temps for the stages
    fname = eval('suspension.susp{}'.format(ifo.Suspension.Type))
    hForce, vForce, hTable, vTable = fname(f, ifo)

    try:
        # full TF (conventional)
        ifo.Suspension.hForce = hForce.fullylossy
        ifo.Suspension.vForce = vForce.fullylossy
        # TFs with each stage lossy
        ifo.Suspension.hForce_singlylossy = hForce.singlylossy
        ifo.Suspension.vForce_singlylossy = vForce.singlylossy
    except:
        ifo.Suspension.hForce = hForce
        ifo.Suspension.vForce = vForce

    ifo.Suspension.hTable = hTable
    ifo.Suspension.vTable = vTable

    ##############################
    # strain to length conversion
    #
    # used for converting displacement to strain in all noise calculations

    ifo.gwinc.dhdl_sqr, ifo.gwinc.sinc_sqr = dhdl(f, ifo.Infrastructure.Length)

    return ifo


def precompPower(ifo, PRfixed=True):
    """Compute power on beamsplitter, finesse, and power recycling factor.

    """
    c       = const.c
    pin     = ifo.Laser.Power
    t1      = sqrt(ifo.Optics.ITM.Transmittance)
    r1      = sqrt(1 - ifo.Optics.ITM.Transmittance)
    r2      = sqrt(1 - ifo.Optics.ETM.Transmittance)
    t5      = sqrt(ifo.Optics.PRM.Transmittance)
    r5      = sqrt(1 - ifo.Optics.PRM.Transmittance)
    loss    = ifo.Optics.Loss                          # single TM loss
    bsloss  = ifo.Optics.BSLoss
    acoat   = ifo.Optics.ITM.CoatingAbsorption
    pcrit   = ifo.Optics.pcrit

    # Finesse, effective number of bounces in cavity, power recycling factor
    finesse = 2*pi / (t1**2 + 2*loss)        # arm cavity finesse
    neff    = 2 * finesse / pi

    # Arm cavity reflectivity with finite loss
    garm = t1 / (1 - r1*r2*sqrt(1-2*loss))  # amplitude gain wrt input field
    rarm = r1 - t1 * r2 * sqrt(1-2*loss) * garm

    if PRfixed:
        Tpr = ifo.Optics.PRM.Transmittance  # use given value
    else:
        Tpr = 1-(rarm*sqrt(1-bsloss))**2 # optimal recycling mirror transmission
        t5 = sqrt(Tpr)
        r5 = sqrt(1 - Tpr)
    prfactor = t5**2 / (1 + r5 * rarm * sqrt(1-bsloss))**2

    pbs  = pin * prfactor          # BS power from input power
    parm = pbs * garm**2 / 2       # arm power from BS power

    asub = 1.3*2*ifo.Optics.ITM.Thickness*ifo.Optics.SubstrateAbsorption
    pbsl = 2*pcrit/(asub+acoat*neff) # bs power limited by thermal lensing

    if pbs > pbsl:
        logging.warning('P_BS exceeds BS Thermal limit!')

    return pbs, parm, finesse, prfactor, Tpr


def dhdl(f, armlen):
    """Strain to length conversion for noise power spetra

    This takes into account the GW wavelength and is only important
    when this is comparable to the detector arm length.

    From R. Schilling, CQG 14 (1997) 1513-1519, equation 5,
    with n = 1, nu = 0.05, ignoring overall phase and cos(nu)^2.
    A small value of nu is used instead of zero to avoid infinities.

    Returns the square of the dh/dL function, and the same divided by
    the arm length squared.

    """
    c = const.c
    nu_small = 15*pi/180
    omega_arm = pi * f * armlen / c
    omega_arm_f = (1 - sin(nu_small)) * pi * f * armlen / c
    omega_arm_b = (1 + sin(nu_small)) * pi * f * armlen / c
    sinc_sqr = 4 / abs(sin(omega_arm_f) * exp(-1j * omega_arm) / omega_arm_f
                       + sin(omega_arm_b) * exp(1j * omega_arm) / omega_arm_b)**2
    # keep DC value equal to 1
    sinc_sqr /= sinc_sqr[0]
    dhdl_sqr = sinc_sqr / armlen**2
    return dhdl_sqr, sinc_sqr
