from numpy import pi, sqrt

from . import const


def sql(ifo):
    """Computer standard quantum limit (SQL) for IFO"""
    c = const.c
    Parm = ifo.gwinc.parm
    w0 = 2 * pi * c / ifo.Laser.Wavelength
    m = ifo.Materials.MirrorMass
    Titm = ifo.Optics.ITM.Transmittance
    Tsrm = ifo.Optics.SRM.Transmittance
    tSR = sqrt(Tsrm)
    rSR = sqrt(1 - Tsrm)
    fSQL = (1/(2*pi))*(8/c)*sqrt((Parm*w0)/(m*Titm))*(tSR/(1+rSR))
    return fSQL


def computeFCParams(ifo, fcParams):
    """Compute ideal filter cavity Tin, detuning [Hz] and bandwidth [Hz]

    """
    # FC parameters
    c = const.c
    fsrFC = c / (2 * fcParams.L)
    lossFC = fcParams.Lrt + fcParams.Te

    fSQL = sql(ifo)

    # detuning and cavity bandwidth (D&D paper P1400018 and/or PRD)
    eps = 4 / (2 + sqrt(2 + 2 * sqrt(1 + (4 * pi * fSQL / (fsrFC * lossFC))**4)))
    s1eps = sqrt(1 - eps)

    # cavity bandwidth [Hz]
    gammaFC = fSQL / sqrt(s1eps + s1eps**3)
    # cavity detuning [Hz]
    detuneFC = s1eps * gammaFC

    # input mirror transmission
    TinFC = 4 * pi * gammaFC / fsrFC - lossFC
    if TinFC < lossFC:
        raise RuntimeError('IFC: Losses are too high! %.1f ppm max.' % 1e6 * gammaFC / fsrFC)

    # Add to fcParams structure
    fcParams.Ti = TinFC
    fcParams.fdetune = -detuneFC
    fcParams.gammaFC = gammaFC
    fcParams.fsrFC = fsrFC

    return fcParams
