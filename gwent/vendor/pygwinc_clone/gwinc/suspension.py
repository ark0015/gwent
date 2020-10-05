from __future__ import division
from numpy import pi, sqrt, sin, cos, tan, real, imag, zeros
import numpy as np
import logging

from . import const
from .struct import Struct


# supported fiber geometries
FIBER_TYPES = [
    'Round',
    'Ribbon',
    'Tapered',
]


def generate_symbolic_tfs(stages=4):
    import sympy as sp

    # construct quad pendulum equation of motion matrix
    ksyms = sp.numbered_symbols('k')
    msyms = sp.numbered_symbols('m')
    w = sp.symbols('w')
    k = [next(ksyms) for n in range(stages)]
    m = [next(msyms) for n in range(stages)]
    A = sp.zeros(stages)
    for n in range(stages-1):
        # mass and restoring forces (diagonal elements)
        A[n, n] = k[n] + k[n+1] - m[n] * w**2
        # couplings to stages above and below
        A[n, n+1] = -k[n+1]
        A[n+1, n] = -k[n+1]
    # mass and restoring force of bottom stage
    A[-1, -1] = k[-1] - m[-1] * w**2

    # want TM equations of motion, so index 4
    b = sp.zeros(stages, 1)
    b[-1] = 1

    # solve linear system
    xsyms = sp.numbered_symbols('x')
    x = [next(xsyms) for n in range(stages)]
    ans = sp.linsolve((A, b), x)
    return ans


def tst_force_to_tst_displ(k, m, f):
    """transfer function for quad pendulum
    
    """
    k0, k1, k2, k3 = k
    m0, m1, m2, m3 = m
    w = 2*pi*f
    X3 = (k2**2*(k0 + k1 - m0*w**2) + (k1**2 - (k0 + k1 - m0*w**2)*(k1 + k2 - m1*w**2))*(k2 + k3 - m2*w**2))/(-k3**2*(k1**2 - (k0 + k1 - m0*w**2)*(k1 + k2 - m1*w**2)) + (k3 - m3*w**2)*(k2**2*(k0 + k1 - m0*w**2) - (-k1**2 + (k0 + k1 - m0*w**2)*(k1 + k2 - m1*w**2))*(k2 + k3 - m2*w**2)))
    return X3


def top_displ_to_tst_displ(k, m, f):
    """transfer function for quad pendulum
    
    """
    k0, k1, k2, k3 = k
    m0, m1, m2, m3 = m
    w = 2*pi*f
    X0 = k1*k2*k3/(k3**2*(k1**2 - (k0 + k1 - m0*w**2)*(k1 + k2 - m1*w**2)) - (k3 - m3*w**2)*(k2**2*(k0 + k1 - m0*w**2) + (k1**2 - (k0 + k1 - m0*w**2)*(k1 + k2 - m1*w**2))*(k2 + k3 - m2*w**2)))
    return X0 * k0


def suspQuad(f, ifo, material='Silica'):
    """Suspension for quadruple pendulum

    `f` is frequency vector, `ifo` is IFO model.  `material` specifies
    material used for test mass suspension stage.  steel used for all
    other stages.  Violin modes are included.

    ifo.Suspension.FiberType should be: 0=round, 1=ribbons.

    hForce, vForce = transfer functions from the force on the TM to TM
    motion these should have the correct losses for the mechanical
    system such that the thermal noise is:

    dxdF = force on TM along beam line to position of TM along beam line
         = hForce + theta^2 * vForce
         = admittance / (i * w)

    where theta = ifo.Suspension.VHCoupling.theta.

    Since this is just suspension thermal noise, the TM internal modes
    and coating properties should not be included.
    
    hTable, vTable = TFs from support motion to TM motion
    
    Ah = horizontal equations of motion
    Av = vertical equations of motion
    
    Adapted from code by Morag Casey (Matlab) and Geppo Cagnoli
    (Maple).  Modification for the different temperatures between the
    stages by K.Arai.

    """
    g = const.g

    sus = ifo.Suspension

    # bottom stage fiber Type
    FiberType = sus.FiberType
    assert FiberType in FIBER_TYPES

    ##############################
    # Parameter Assignment

    def scarray():
        return np.zeros(len(sus.Stage))

    ds_w  = scarray()
    dil0 = scarray()
    mass = scarray()
    length = scarray()
    kv0 = scarray()
    r_w = scarray()
    t_b = scarray()
    N_w = scarray()
    Temp = scarray()
    # wire material properties
    alpha_w = scarray()
    beta_w = scarray()
    rho_w = scarray()
    C_w = scarray()
    K_w = scarray()
    Y_w = scarray()
    phi_w = scarray()
    # blade material properties
    alpha_b = scarray()
    beta_b = scarray()
    rho_b = scarray()
    C_b = scarray()
    K_b = scarray()
    Y_b = scarray()
    phi_b = scarray()

    # NOTE: reverse suspension list so that the last stage in the list
    # is the test mass
    last_stage = len(sus.Stage) - 1
    for n, stage in enumerate(reversed(sus.Stage)):
        ##############################
        # main suspension parameters

        mass[n] = stage.Mass
        length[n] = stage.Length
        dil0[n] = stage.Dilution
        kv0[n] = stage.K

        if np.isnan(stage.WireRadius):
            if 'FiberRadius' in stage:
                r_w[n] = stage.FiberRadius
            elif FiberType == 'Ribbon':
                r_w[n] = sus.Ribbon.Width
            else:
                r_w[n] = sus.Fiber.Radius
        else:
            r_w[n] = stage.WireRadius

        # blade thickness
        t_b[n] = stage.Blade
        # number of support wires
        N_w[n] = stage.NWires

        if 'Temp' in stage:
            Temp[n] = stage.Temp
        else:
            Temp[n] = sus.Temp

        ##############################
        # support wire material parameters

        if 'WireMaterial' in stage:
            WireMaterial = stage.WireMaterial
        elif n == last_stage:
            WireMaterial = material
        else:
            WireMaterial = 'C70Steel'

        logging.debug('stage {} wires: {}'.format(n, WireMaterial))
        wireMat = sus[WireMaterial]

        alpha_w[n] = wireMat.Alpha  # coeff. thermal expansion
        beta_w[n] = wireMat.dlnEdT  # temp. dependence Youngs modulus
        rho_w[n] = wireMat.Rho      # mass density
        C_w[n] = wireMat.C          # heat capacity
        K_w[n] = wireMat.K          # W/(m kg)
        Y_w[n] = wireMat.Y          # Young's modulus
        phi_w[n] = wireMat.Phi      # loss angle

        # surface loss dissipation depth
        if 'Dissdepth' in wireMat:
            ds_w[n] = wireMat.Dissdepth
        else:
            # otherwise ignore surface effects
            ds_w[n] = 0

        ##############################
        # support blade material parameters

        if 'BladeMaterial' in stage:
            BladeMaterial = stage.BladeMaterial
        else:
            BladeMaterial = 'MaragingSteel'

        logging.debug('stage {} blades: {}'.format(n, BladeMaterial))
        bladeMat = sus[BladeMaterial]

        alpha_b[n] = bladeMat.Alpha   # coeff. thermal expansion
        beta_b[n] = bladeMat.dlnEdT   # temp. dependence Youngs modulus
        rho_b[n] = bladeMat.Rho       # mass density
        C_b[n] = bladeMat.C           # heat capacity
        K_b[n] = bladeMat.K           # W/(m kg)
        Y_b[n] = bladeMat.Y           # Young's modulus
        phi_b[n] = bladeMat.Phi       # loss angle

    # weight support by lower stages
    Mg = g * np.flipud(np.cumsum(np.flipud(mass)))

    # Correction for the pendulum restoring force
    kh0 = Mg / length              # N/m, horiz. spring constant, stage n

    ##############################
    # Thermoelastic Calculations for wires and blades

    # wire geometry
    tension = Mg / N_w           # Tension
    xsect = pi * r_w**2          # cross-sectional area
    xII = r_w**4 * pi / 4        # x-sectional moment of inertia
    mu_h = 4 / r_w               # surface to volume ratio, horizontal
    mu_v = 2 / r_w               # surface to volume ratio, vertical (wire)

    # horizontal TE time constant, wires
    # The constant 7.37e-2 is 1/(4*q0^2) from eq 12, C. Zener 10.1103/PhysRev.53.90 (1938)
    tau_h = 7.37e-2 * 4 * (rho_w * C_w * xsect) / (pi * K_w)

    # vertical TE time constant, blades
    tau_v = (rho_b * C_b * t_b**2) / (K_b * pi**2)

    # vertical delta, blades
    delta_v = Y_b * alpha_b**2 * Temp / (rho_b * C_b)

    # deal with ribbon geometry for last stage
    if FiberType == 'Ribbon':
        W = sus.Ribbon.Width
        t = sus.Ribbon.Thickness
        xsect[-1] = W * t                   # cross-sectional area
        xII[-1] = (W * t**3)/12             # x-sectional moment of inertia
        mu_v[-1] = 2 * (W + t) / (W * t)
        mu_h[-1] = mu_v[-1] * (3 * N_w[-1] * W + t) / (N_w[-1] * W + t)
        tau_h[-1] = (rho_w[-1] * C_w[-1] * t**2) / (K_w[-1] * pi**2)

    # horizontal delta, wires
    delta_h = (alpha_w - tension * beta_w / (xsect * Y_w))**2 * Y_w * Temp / (rho_w * C_w)

    # deal with tapered geometry for last stage
    if FiberType == 'Tapered':
        r_end = sus.Fiber.EndRadius

        # recompute these for
        xsectEnd = pi * r_end**2      # cross-sectional area (for delta_h)
        xII[-1] = pi * r_end**4 / 4    # x-sectional moment of inertia
        mu_h[-1] = 4 / r_end           # surface to volume ratio, horizontal

        # use this xsect for thermo-elastic noise
        delta_h[-1] = (alpha_w[-1] - tension[-1] * beta_w[-1] / (xsectEnd * Y_w[-1]))**2 * Y_w[-1] * Temp[-1] / (rho_w[-1] * C_w[-1])

    # bending length, and dilution factors
    d_bend = sqrt(Y_w * xII / tension)
    dil = length / d_bend
    # dil(~isnan(dil0)) = dil0(~isnan(dil0))
    dil = np.where(~np.isnan(dil0), dil0, dil)

    ##############################
    # Loss Calculations for wires and blades

    # these calculations use the frequency vector
    w = 2 * pi * f

    phih = np.zeros([len(sus.Stage), len(w)], dtype=complex)
    kh = np.zeros([len(sus.Stage), len(w)], dtype=complex)
    phiv = np.zeros([len(sus.Stage), len(w)], dtype=complex)
    kv = np.zeros([len(sus.Stage), len(w)], dtype=complex)

    for n, stage in enumerate(sus.Stage):
        # horizontal loss factor, wires
        phih[n, :] = phi_w[n] * (1 + mu_h[n] * ds_w[n]) + delta_h[n] * tau_h[n] * w / (1 + w**2 * tau_h[n]**2)

        # complex spring constant, horizontal
        kh[n, :] = kh0[n] * (1 + 1j * phih[n, :] / dil[n])

        # vertical loss factor, blades
        phiv[n, :] = phi_b[n] + delta_v[n] * tau_v[n] * w / (1 + w**2 * tau_v[n]**2)

        # complex spring constant, vertical
        kv[n, :] = kv0[n] * (1 + 1j * phiv[n, :])

    ##############################
    # last suspension stage
    # Equations from "GG" (maybe?)
    #   Suspensions thermal noise in the LIGO gravitational wave detector
    #   Gabriela Gonzalez, Class. Quantum Grav. 17 (2000) 4409?4435
    #
    # Note:
    #  I in GG = xII
    #  rho in GG = rho_w * xsect
    #  delta in GG = d_bend
    ##############################

    ### Vertical (bounce) ###
    # loss factor, last stage suspension, vertical (no blades)
    phiv[-1, :] = phi_w[-1] * (1 + mu_v[-1] * ds_w[-1])

    # vertical Young's modulus
    Y_v = Y_w[-1] * (1 + 1j * phiv[-1, :])

    # vertical spring constant, last stage
    k_z = sqrt(rho_w[-1] / Y_v) * w
    kv4 = N_w[-1] * xsect[-1] * Y_v * k_z / (tan(k_z * length[-1]))

    # deal with tapered geometry for last stage
    if FiberType == 'Tapered' and 'EndLength' in sus.Fiber:
        l_end = 2 * sus.Fiber.EndLength
        l_mid = length[-1] - l_end

        kv_mid = N_w[-1] * xsect[-1] * Y_v * k_z / (tan(k_z * l_mid))
        kv_end = N_w[-1] * xsectEnd * Y_v * k_z / (tan(k_z * l_end))
        kv4 = kv_mid * kv_end / (kv_mid + kv_end)

    if np.isnan(kv0[-1]):
        kv[-1, :] = kv4 # no blades
    else:
        kv[-1, :] = kv[-1, :] * kv4 / (kv[-1, :] + kv4) # with blades

    ### Horizontal (pendulum and violins) ###
    # horizontal Young's modulus
    Y_h  = Y_w[-1] * (1 + 1j * phih[-1, :])

    # simplification factors for later calculations
    ten4 = tension[-1]                          # T in GG
    k4 = sqrt(rho_w[-1] * xsect[-1] / ten4) * w	# k in GG
    d_bend4 = sqrt(Y_h * xII[-1] / ten4)        # complex d_bend(4)
    dk4 = k4 * d_bend4

    # simp3a is inherited from the previous version of this suspension
    # thermal noise calculation (part of simp3).
    #
    # I'm not sure where this comes from, but it differs from
    # 1 by less than 1e-6 for frequencies below 100Hz (and less than
    # 1e-3 up to 10kHz).  What's more, the units appear to be wrong.
    # (Missing factor of rho * length?)
    # [mevans June 2015]
    #
    # simp3a = sqrt(1 + d_bend4 .* xsect(4) .* w.^2 / ten4)
    simp3a = 1

    coskl = simp3a * cos(k4 * length[-1])
    sinkl = sin(k4 * length[-1])

    # numerator, horiz spring constant, last stage
    #   numerator of K_xx in eq 9 of GG
    #     = T k (cos(k L) + k delta sin(k L))
    #   for w -> 0, this reduces to N_w * T * k
    kh4num  = N_w[-1] * ten4 * k4 * simp3a * (simp3a**2 + dk4**2) * (coskl + dk4 * sinkl)

    # denominator, horiz spring constant, last stage
    #   D after equation 8 in GG
    #   D = sin(k L) - 2 k delta cos(k L)
    #   for w -> 0, this reduces to k (L - 2 delta)
    kh4den = ((simp3a**2 - dk4**2) * sinkl - 2 * dk4 * coskl)

    # horizontal spring constant, last stage
    #   K_xx in eq 9 of GG
    kh[-1, :] = kh4num / kh4den

    ###############################################################
    # Equations of motion for the system
    ###############################################################

    # Calculate TFs turning on the loss of each stage one by one
    hForce = Struct()
    vForce = Struct()
    hForce.singlylossy = np.zeros([len(sus.Stage), len(w)], dtype=complex)
    vForce.singlylossy = np.zeros([len(sus.Stage), len(w)], dtype=complex)

    for n in range(len(mass)):
        # horizontal
        # only the imaginary part of the specified stage is used.
        k = real(kh) + 1j*imag([kh[0,:]*(n==0),
                                kh[1,:]*(n==1),
                                kh[2,:]*(n==2),
                                kh[3,:]*(n==3)])
        # calculate TFs
        hForce.singlylossy[n,:] = tst_force_to_tst_displ(k, mass, f)

        # vertical
        # only the imaginary part of the specified stage is used
        k = real(kv) + 1j*imag([kv[0,:]*(n==0),
                                kv[1,:]*(n==1),
                                kv[2,:]*(n==2),
                                kv[3,:]*(n==3)])
        # calculate TFs
        vForce.singlylossy[n,:] = tst_force_to_tst_displ(k, mass, f)

    # calculate horizontal TFs with all losses on
    hForce.fullylossy = tst_force_to_tst_displ(kh, mass, f)
    hTable = top_displ_to_tst_displ(kh, mass, f)

    # calculate vertical TFs with all losses on
    vForce.fullylossy = tst_force_to_tst_displ(kv, mass, f)
    vTable = top_displ_to_tst_displ(kv, mass, f)

    return hForce, vForce, hTable, vTable


def suspBQuad(f, ifo):
    """Wrapper of suspQuad to use Silicon for final stage

    FIXME: material should be specified in ifo.Suspension.Stage

    """
    return suspQuad(f, ifo, material='Silicon')
