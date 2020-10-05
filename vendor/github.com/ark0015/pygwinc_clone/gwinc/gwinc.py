from __future__ import division
import numpy as np
from numpy import pi, sqrt
from collections import OrderedDict
import logging

from . import load_ifo
from .precomp import precompIFO
from . import noise
from .plot import plot_noise


def gwinc(freq, ifoin, source=None, plot=False, PRfixed=True):
    """Calculate strain noise budget for a specified interferometer model.

    Argument `freq` is the frequency array for which the noises will
    be calculated, and `ifoin` is the IFO model (see the `load_ifo()`
    function).

    If `source` structure provided, so evaluates the sensitivity of
    the detector to several potential gravitational wave
    sources.

    If `plot` is True a plot of the budget will be created.

    Returns tuple of (score, noises, ifo)

    """
    # assume generic aLIGO configuration
    # FIXME: how do we allow adding arbitrary addtional noise sources
    # from just ifo description, without having to specify full budget
    Budget, ifo_, freq_, plot_style = load_ifo('aLIGO')

    # add some precomputed info to the ifo struct
    #this implicitly deepcopies and the return value is the copy
    ifo = precompIFO(freq, ifoin, PRfixed)

    traces = Budget(freq, ifo=ifo).calc_trace()

    noises = {n:d[0] for n, d in traces.items()}
    noises['Freq'] = freq

    pbs      = ifo.gwinc.pbs
    parm     = ifo.gwinc.parm
    finesse  = ifo.gwinc.finesse
    prfactor = ifo.gwinc.prfactor
    if ifo.Laser.Power * prfactor != pbs:
        pass
        #warning(['Thermal lensing limits input power to ' num2str(pbs/prfactor, 3) ' W']);

    #TODO decide if all this below this should remain, since it is already inside of __main__

    # report astrophysical scores if so desired
    score = None
    if source:
        score = int73(freq, noises['Total'], ifo, source)
        score.Omega = intStoch(freq, noises['Total'], 0, ifo, source)

    # --------------------------------------------------------
    # output graphics

    if plot:
        # Report input parameters
        if ifo.Optics.Type == 'DualCarrier_new':     #include the case for Dual carrier
            finesseA = 2*pi/ifo.Optics.ITM.TransmittanceD1
            finesseB = 2*pi/ifo.Optics.ITM.TransmittanceD2
            pbsA = ifo.Laser.PBSD1
            pbsB = ifo.Laser.PBSD2
            logging.info('Finesse for carrier A:  %7.2f' % finesseA)
            logging.info('Finesse for carrier B:  %7.2f' % finesseB)
            logging.info('Power Recycling Factor: %7.2f' % ifo.PRCgain)
            logging.info('Arm power for carrier A:%7.2f kW' % (finesseA*2/pi*pbsA/2/1000))
            logging.info('Arm power for carrier B:%7.2f kW' % (finesseB*2/pi*pbsB/2/1000))
            logging.info('Power on beam splitter for carrier A: %7.2f W' % pbsA)
            logging.info('Power on beam splitter for carrier B: %7.2f W' % pbsB)
            logging.info('Laser Power for Carrier A:     %7.2f Watt' % ifo.LP1)
            logging.info('Laser Power for Carrier B:     %7.2f Watt' % ifo.LP2)
            logging.info('SRM Detuning for Carrier A:    %7.2f degree' % (ifo.Optics.SRM.TunephaseD1*180/pi))
            logging.info('SRM Detuning for Carrier B:    %7.2f degree' % (ifo.Optics.SRM.TunephaseD2*180/pi))
            logging.info('SRM transmission for Carrier A:%9.4f' % ifo.Optics.SRM.TransmittanceD1)
            logging.info('SRM transmission for Carrier B:%9.4f' % ifo.Optics.SRM.TransmittanceD2)
            logging.info('ITM transmission for Carrier A:%9.4f' % ifo.Optics.ITM.TransmittanceD1)
            logging.info('ITM transmission for Carrier B:%9.4f' % ifo.Optics.ITM.TransmittanceD2)
            logging.info('PRM transmission for both:     %9.4f' % ifo.Optics.PRM.Transmittance)
        else:
            logging.info('Laser Power:            %7.2f Watt' % ifo.Laser.Power)
            logging.info('SRM Detuning:           %7.2f degree' % (ifo.Optics.SRM.Tunephase*180/pi))
            logging.info('SRM transmission:       %9.4f' % ifo.Optics.SRM.Transmittance)
            logging.info('ITM transmission:       %9.4f' % ifo.Optics.ITM.Transmittance)
            logging.info('PRM transmission:       %9.4f' % ifo.Optics.PRM.Transmittance)
            logging.info('Finesse:                %7.2f' % finesse)
            logging.info('Power Recycling Gain:   %7.2f' % prfactor)
            logging.info('Arm Power:              %7.2f kW' % (parm/1000))
            logging.info('Power on BS:            %7.2f W' % pbs)

        # coating and substrate thermal load on the ITM
        PowAbsITM = (pbs/2) * \
                    np.hstack([(finesse*2/pi) * ifo.Optics.ITM.CoatingAbsorption,
                               (2 * ifo.Materials.MassThickness) * ifo.Optics.ITM.SubstrateAbsorption])

        # Stefan's Mysterious TCS Section ~~~~ `~~ ` ` 324@#%@#$ !
        M = np.array([[ifo.TCS.s_cc, ifo.TCS.s_cs], [ifo.TCS.s_cs, ifo.TCS.s_ss]])
        S_uncorr = PowAbsITM.T*M*PowAbsITM
        TCSeff = 1-sqrt(ifo.TCS.SRCloss/S_uncorr)

        logging.info('Thermal load on ITM:    %8.3f W' % sum(PowAbsITM))
        logging.info('Thermal load on BS:     %8.3f W' %
                     (ifo.Materials.MassThickness*ifo.Optics.SubstrateAbsorption*pbs))
        if (ifo.Laser.Power*prfactor != pbs):
            logging.info('Lensing limited input power: %7.2f W' % (pbs/prfactor))

        if source:
            logging.info('BNS Inspiral Range:     ' + str(score.effr0ns) + ' Mpc/ z = ' + str(score.zHorizonNS))
            logging.info('BBH Inspiral Range:     ' + str(score.effr0bh) + ' Mpc/ z = ' + str(score.zHorizonBH))
            logging.info('Stochastic Omega: %4.1g Universes' % score.Omega)

        plot_noise(ifo, traces, **plot_style)
    return score, noises, ifo
