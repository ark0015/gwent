from __future__ import division
from numpy import pi, imag, zeros
import numpy as np

from .. import const


def susptherm(f, ifo):
    """Suspention thermal noise.

    'Temp' must either be set for each stage individually, or globally
    in ifo.Suspension.Temp.  The latter will be preferred if
    specified, so if you wish to use per-stage tempurature you must
    remove Suspension.Temp.

    Assumes suspension transfer functions and V-H coupling have been
    pre-calculated and populated into the relevant `ifo` struct
    fields.

    """
    # Assign Physical Constants
    kB = const.kB

    # and vertical to beamline coupling angle
    theta = ifo.Suspension.VHCoupling.theta

    noise = zeros((1, f.size))

    # if the temperature is uniform along the suspension
    if 'Temp' in ifo.Suspension:
        ##########################################################
        # Suspension TFs
        ##########################################################

        hForce = ifo.Suspension.hForce
        vForce = ifo.Suspension.vForce

        ##########################################################
        # Thermal Noise Calculation
        ##########################################################

        # convert to beam line motion
        #  theta is squared because we rotate by theta into the suspension
        #  basis, and by theta to rotate back to the beam line basis
        dxdF = hForce + theta**2 * vForce

        # thermal noise (m^2/Hz) for one suspension
        w = 2*pi*f
        noise = 4 * kB * ifo.Suspension.Temp * abs(imag(dxdF)) / w

    # if the temperature is set for each suspension stage
    else:
        ##########################################################
        # Suspension TFs
        ##########################################################

        hForce = ifo.Suspension.hForce_singlylossy[:, :]
        vForce = ifo.Suspension.vForce_singlylossy[:, :]

        ##########################################################
        # Thermal Noise Calculation
        ##########################################################

        dxdF = zeros(hForce.shape, dtype=complex)
        for n, stage in enumerate(reversed(ifo.Suspension.Stage)):
            # add up the contribution from each stage

            # convert to beam line motion.  theta is squared because
            # we rotate by theta into the suspension basis, and by
            # theta to rotate back to the beam line basis
            dxdF[n, :] = hForce[n, :] + theta**2 * vForce[n, :]

            # thermal noise (m^2/Hz) for one suspension
            w = 2*pi*f
            noise += 4 * kB * stage.Temp * abs(imag(dxdF[n, :])) / w

    # 4 masses, turn into gravitational wave strain
    noise *= 4 * ifo.gwinc.dhdl_sqr

    return np.squeeze(noise)
