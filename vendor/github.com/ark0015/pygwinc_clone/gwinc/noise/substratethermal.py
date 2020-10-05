from __future__ import division, print_function
from numpy import exp, inf, pi, sqrt
import numpy as np
import scipy.special
import scipy.integrate

from .. import const
from ..const import BESSEL_ZEROS as zeta
from ..const import J0M as j0m


def carrierdensity_adiabatic(f, ifo):
    """strain noise psd arising from charge carrier density
    fluctuations in ITM substrate (for semiconductor substrates)."""

    Omega = 2*pi*f
    H = ifo.Materials.MassThickness
    gammaElec = ifo.Materials.Substrate.ElectronIndexGamma
    gammaHole = ifo.Materials.Substrate.HoleIndexGamma
    diffElec = ifo.Materials.Substrate.ElectronDiffusion
    diffHole = ifo.Materials.Substrate.HoleDiffusion
    cdDens = ifo.Materials.Substrate.CarrierDensity
    r0 = ifo.Optics.ITM.BeamRadius/np.sqrt(2)
    L = ifo.Infrastructure.Length
    gPhase = ifo.gwinc.finesse*2/pi

    psdElec = 4*H*gammaElec**2*cdDens*diffElec/(pi*r0**4*Omega**2) # units are meters
    psdHole = 4*H*gammaHole**2*cdDens*diffHole/(pi*r0**4*Omega**2) # units are meters
    psdMeters = 2 * (psdElec + psdHole) # electrons and holes for two ITMs
    n = psdMeters / (gPhase*L)**2
    return n


def carrierdensity_exact(f, ifo):
    """Strain noise arising from charge carrier density fluctuations in ITM substrate

    For semiconductor substrates

    """
    w = ifo.Optics.ITM.BeamRadius
    L = ifo.Infrastructure.Length
    H = ifo.Materials.MassThickness
    kBT = const.kB * ifo.Materials.Substrate.Temp
    hbar = const.hbar
    c = const.c
    
    diffElec = ifo.Materials.Substrate.ElectronDiffusion
    diffHole = ifo.Materials.Substrate.HoleDiffusion
    mElec = ifo.Materials.Substrate.ElectronEffMass
    mHole = ifo.Materials.Substrate.HoleEffMass
    cdDens = ifo.Materials.Substrate.CarrierDensity
    gammaElec = ifo.Materials.Substrate.ElectronIndexGamma
    gammaHole = ifo.Materials.Substrate.HoleIndexGamma
    gPhase = ifo.gwinc.finesse*2/pi

    omega = 2*pi*f

    def integrand(k, om, D):
        return D * k**3 * exp(-k**2 * w**2/4) / (D**2 * k**4 + om**2)
    
    integralElec = np.array([scipy.integrate.quad(lambda k: integrand(k, om, diffElec), 0, inf)[0] for om in omega])
    integralHole = np.array([scipy.integrate.quad(lambda k: integrand(k, om, diffHole), 0, inf)[0] for om in omega])

    # From P1400084 Heinert et al. Eq. 15 
    #psdCD = @(gamma,m,int) 2*(3/pi^7)^(1/3)*kBT*H*gamma^2*m/hbar^2*cdDens^(1/3)*int; %units are meters
    def psdCD(gamma,m,int_):
        return 2/pi * H * gamma**2 * cdDens * int_ #units are meters
    
    psdElec = psdCD(gammaElec, mElec, integralElec)
    psdHole = psdCD(gammaHole, mHole, integralHole)
    
    psdMeters = 2 * (psdElec + psdHole)
    
    n = psdMeters / (gPhase*L)**2

    return n

carrierdensity = carrierdensity_adiabatic


def thermorefractiveITM_adiabatic(f, ifo):
    """strain noise psd arising from thermorefractive
    fluctuations in ITM substrate (for semiconductor substrates)."""

    Omega = 2*pi*f
    H = ifo.Materials.MassThickness
    beta = ifo.Materials.Substrate.dndT
    kappa = ifo.Materials.Substrate.MassKappa
    rho = ifo.Materials.Substrate.MassDensity
    C = ifo.Materials.Substrate.MassCM
    Temp = ifo.Materials.Substrate.Temp
    kBT = const.kB * Temp
    r0 = ifo.Optics.ITM.BeamRadius/np.sqrt(2)
    L = ifo.Infrastructure.Length
    gPhase = ifo.gwinc.finesse*2/pi

    psd = 4*H*beta**2*kappa*kBT*Temp/(pi*r0**4*Omega**2*(rho*C)**2) # units are meters
    psdMeters = 2*psd # two ITMs
    n = psdMeters / (gPhase*L)**2
    return n


def thermorefractiveITM_exact(f, ifo):
    """Strain noise from thermorefractive fluctuations in ITM substrate

    For semiconductor substrates.

    """
    
    w = ifo.Optics.ITM.BeamRadius
    L = ifo.Infrastructure.Length
    H = ifo.Materials.MassThickness
    kBT = const.kB * ifo.Materials.Substrate.Temp
    Temp = ifo.Materials.Substrate.Temp
    c = const.c
    
    rho = ifo.Materials.Substrate.MassDensity
    beta = ifo.Materials.Substrate.dndT
    C = ifo.Materials.Substrate.MassCM
    kappa = ifo.Materials.Substrate.MassKappa

    gPhase = ifo.gwinc.finesse*2/pi

    omega = 2*pi*f
    
    def integrand(k,om,D):
        return D * k**3 * exp(-k**2 * w**2/4) / (D**2 * k**4 + om**2)
    
    inte = np.array([scipy.integrate.quad(lambda k: integrand(k, om, kappa/(rho*C)), 0, inf)[0] for om in omega])
    
    # From P1400084 Heinert et al. Eq. 15 
    #psdCD = @(gamma,m,int) 2*(3/pi^7)^(1/3)*kBT*H*gamma^2*m/hbar^2*cdDens^(1/3)*int; %units are meters
    psdTR = lambda int_: 2/pi * H * beta**2 * kBT * Temp / (rho*C) * int_; #units are meters
    
    
    psd = psdTR(inte)
    
    psdMeters = 2*psd # two itms
    
    n = psdMeters / (gPhase*L)**2

    return n

thermorefractiveITM = thermorefractiveITM_adiabatic


def subbrownian(f, ifo):
    """Strain noise from the Brownian thermal noise due to substrate mechanical loss

    """
    wITM = ifo.Optics.ITM.BeamRadius
    wETM = ifo.Optics.ETM.BeamRadius
    Y = ifo.Materials.Substrate.MirrorY
    sigma = ifo.Materials.Substrate.MirrorSigma

    c2 = ifo.Materials.Substrate.c2
    n = ifo.Materials.Substrate.MechanicalLossExponent
    alphas = ifo.Materials.Substrate.Alphas
    L = ifo.Infrastructure.Length
    kBT = const.kB * ifo.Materials.Substrate.Temp

    # Bulk substrate contribution
    phibulk = c2 * f**n

    cITM, aITM = subbrownianFiniteCorr(ifo, 'ITM')
    cETM, aETM = subbrownianFiniteCorr(ifo, 'ETM')
    cbulk = 8 * kBT * (aITM + aETM) * phibulk / (2 * pi * f)

    # Surface loss contribution
    # csurfETM = alphas/(Y*pi*wETM^2);
    # csurfITM = alphas/(Y*pi*wITM^2);

    csurfETM = alphas*(1-2*sigma)/((1-sigma)*Y*pi*wETM**2)
    csurfITM = alphas*(1-2*sigma)/((1-sigma)*Y*pi*wITM**2)
    csurf = 8 * kBT * (csurfITM + csurfETM) / (2 * pi * f)

    # account for 2 ITM and 2 ETM, and convert to strain whith 1/L^2
    n = 2 * (csurf + cbulk) * ifo.gwinc.dhdl_sqr

    return n


def subbrownianFiniteCorr(ifo, opticName):
    """Amplitude coefficient of mirror thermal noise

    Contribution for finite-size test masses.
    
    [cftm, aftm] = subbrownianFiniteCorr(ifo, opticName)
    cftm = finite mirror correction factor
    aftm = amplitude coefficient for thermal noise:
           thermal noise contribution to displacement noise is
           S_x(f) = (8 * kB * T / (2*pi*f)) * Phi(f) * aftm
    
    Equation references to Bondu, et al. Physics Letters A 246 (1998)
    227-236 (hereafter BHV) and Liu and Thorne gr-qc/0002055 (hereafter LT)

    """
    # get some numbers
    a = ifo.Materials.MassRadius
    h = ifo.Materials.MassThickness
    w = ifo.Optics[opticName].BeamRadius
    Y = ifo.Materials.Substrate.MirrorY
    sigma = ifo.Materials.Substrate.MirrorSigma

    # do the work
    r0 = w / sqrt(2)                     # LT uses e-folding of power
    km = zeta/a

    Qm = exp(-2*km*h)                    # LT eq. 35a

    Um = (1-Qm)*(1+Qm)+4*h*km*Qm
    Um = Um/((1-Qm)**2-4*(km*h)**2*Qm)   # LT 53 (BHV eq. btwn 29 & 30)

    x = exp(-(zeta*r0/a)**2/4)
    s = sum(x/(zeta**2*j0m))             # LT 57

    x2 = x*x
    U0 = sum(Um*x2/(zeta*j0m**2))
    U0 = U0*(1-sigma)*(1+sigma)/(pi*a*Y) # LT 56 (BHV eq. 3)

    p0 = 1/(pi*a**2)                     # LT 28
    DeltaU = (pi*h**2*p0)**2
    DeltaU = DeltaU + 12*pi*h**2*p0*sigma*s
    DeltaU = DeltaU + 72*(1-sigma)*s**2
    DeltaU = DeltaU*a**2/(6*pi*h**3*Y)   # LT 54

    aftm = DeltaU + U0                   # LT 58 (eq. following BHV 31)

    # amplitude coef for infinite TM
    #   factored out: (8 * kB * T * Phi) / (2 * pi * f)
    aitm = (1 - sigma**2) / (2 * sqrt(2 * pi) * Y * r0) # LT 59

    # finite mirror correction
    cftm = aftm / aitm

    return cftm, aftm


def subtherm(f, ifo):
    """Noise from thermoelastic fluctuations in mirror

    """
    wITM = ifo.Optics.ITM.BeamRadius
    wETM = ifo.Optics.ETM.BeamRadius
    sigma = ifo.Materials.Substrate.MirrorSigma

    L = ifo.Infrastructure.Length
    kBT = const.kB * ifo.Materials.Substrate.Temp

    rho = ifo.Materials.Substrate.MassDensity
    kappa = ifo.Materials.Substrate.MassKappa # thermal conductivity
    alpha = ifo.Materials.Substrate.MassAlpha # thermal expansion
    CM = ifo.Materials.Substrate.MassCM # heat capacity @ constant mass
    Temp = ifo.Materials.Substrate.Temp # temperature

    S0 = 8*(1+sigma)**2*kappa*alpha**2*Temp*kBT # note kBT has factor Temp
    S0 = S0/(sqrt(2*pi)*(CM*rho)**2)
    SITM = S0/(wITM/sqrt(2))**3 # LT 18 less factor 1/omega^2
    SETM = S0/(wETM/sqrt(2))**3 # LT 18 less factor 1/omega^2

    # Corrections for finite test masses:
    SITM = SITM * subthermFiniteCorr(ifo, 'ITM')
    SETM = SETM * subthermFiniteCorr(ifo, 'ETM')

    # 2 mirrors each type, factor omega^2, dimensionless strain
    n = 2 * (SITM + SETM)/(2*pi*f)**2 * ifo.gwinc.dhdl_sqr

    return n


def subthermFiniteCorr(ifo, opticName):
    """Finite size test mass correction to noise amplitude coefficient

    (Liu & Thorne gr-qc/0002055 equation 46)
    
    Equation references to Bondu, et al. Physics Letters A 246 (1998)
    227-236 (hereafter BHV) or Liu and Thorne gr-qc/0002055 (hereafter LT)

    """
    # extract some numbers
    a = ifo.Materials.MassRadius
    h = ifo.Materials.MassThickness
    w = ifo.Optics[opticName].BeamRadius
    sigma = ifo.Materials.Substrate.MirrorSigma

    # do the work
    r0 = w/sqrt(2) # LT uses power e-folding
    km = zeta/a

    Qm = exp(-2*km*h) # LT eq. 35a

    pm = exp(-(km*r0)**2/4)/(pi*(a*j0m)**2) # LT 37

    c0 = 6*(a/h)**2*sum(j0m*pm/zeta**2) # LT 32
    c1 = -2*c0/h # LT 32
    p0 = 1/(pi*a**2) # LT 28
    c1 = c1 + p0/(2*h) # LT 40

    coeff = (1-Qm)*((1-Qm)*(1+Qm)+8*h*km*Qm)
    coeff = coeff + 4*(h*km)**2*Qm*(1+Qm)
    coeff = coeff*km*(pm*j0m)**2*(1-Qm)
    coeff = coeff/((1-Qm)**2-4*(h*km)**2*Qm)**2
    coeff = sum(coeff) + h*c1**2/(1+sigma)**2
    coeff = coeff*(sqrt(2*pi)*r0)**3*a**2 # LT 46
    return coeff
