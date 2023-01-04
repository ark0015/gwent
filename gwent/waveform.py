import os

import astropy.constants as const
import astropy.units as u
import lalsimulation
import numpy as np
from astropy.cosmology import WMAP9 as cosmo

import gwent

from . import binary

current_path = os.path.abspath(gwent.__path__[0])
load_directory = os.path.join(current_path, "LoadFiles/")


def Get_Waveform(
    source,
    approximant="pyPhenomD",
    pct_of_peak=0.01,
    lalsuite_kwargs={},
    out_frame="observer",
):
    """Gets the frequency domain waveform of a particular source using approximant type in approximant.
    this can either be the python implementation of IMRPhenomD given below, or a waveform modelled in
    LIGO's lalsuite package. If a lalsuite approximant is used, source will be given plus and cross
    properties with their phase and amplitude for user convenience.

    Parameters
    ----------
    source: object
            source object that contains all necessary source parameters
    approximant: str, optional
            the approximant used to calculate the frequency domain waveform of the source.
            Can either be the python implementation of IMRPhenomD (``'pyPhenomD'``, the default) given below,
            or a waveform modelled in LIGO's ``lalsuite``'s ``lalsimulation`` package.
    pct_of_peak: float, optional
            ``pyPhenomD`` kwarg, the percentange of the strain at merger that dictates the maximum frequency the waveform is calculated at in geometrized units (G=c=1)
    lalsuite_kwargs: dict, optional
            More specific user-defined kwargs for the different lalsuite waveforms
    out_frame: str, {'observer','source'}
        Determines whether the returned frequency is in the source or observer frame.

    Returns
    -------
    freqs: numpy array of floats
        the waveform frequencies in Hertz
    strain: numpy array of floats
        the waveform amplitude strain in dimensionless units
    """
    if approximant == "pyPhenomD":
        if not hasattr(source, "_fitcoeffs"):
            Get_Fitcoeffs(source)
        if not all([hasattr(source, "_phenomD_f"), hasattr(source, "_phenomD_h")]):
            [source._phenomD_f, source._phenomD_h] = Get_PyPhenomD(
                source, pct_of_peak=pct_of_peak
            )

        return Strain_Conv(
            source, source._phenomD_f, source._phenomD_h, out_frame=out_frame
        )
    elif approximant in dir(lalsimulation) or isinstance(approximant, int):
        if isinstance(approximant, int):
            approx = approximant
        else:
            approx = getattr(lalsimulation, approximant)

        waveform_dict = {}

        M_time = (
            source.M.to("kg") * const.G / const.c**3
        )  # Converts M = [M] to M = [sec]

        # Used for finding the starting frequency based on our SNR calculation assumption of
        # wanting to integrate from observed frequency (f(T_obs_before_merger)) till merger
        if not hasattr(source, "f_T_obs"):
            if not hasattr(source, "instrument"):
                # No instrument? We set some time observed
                binary.Check_Freq_Evol(source, T_evol=4.0 * u.yr)
            else:
                binary.Check_Freq_Evol(source)
        # Uses source f_T_obs (observed frequency T_obs before merger) in detector frame so we convert to source frame
        default_f_min = source.f_T_obs * (1 + source.z) / 2.0
        # Uses Mf space as the end frequency (might be a little arbitrary, but lalsuite cuts the waveform off after f_cut anyway), then convert to Hz for lalsuite
        default_f_max = 1.0 / M_time
        # Find even spacing given nfreqs of the source (default is 1e3)
        # default_deltaF = (default_f_max - default_f_min) / source.nfreqs

        # Spacing required to include defaut_f_min in frequency
        # WARNING: THIS CAN RESULT IN AN EXTREMELY LONG ARRAY
        default_deltaF = default_f_min

        waveform_dict.update(
            {
                "m1": source.M.to("kg").value / (1 + source.q),
                "m2": source.M.to("kg").value * source.q / (1 + source.q),
                "S1x": 0,
                "S1y": 0,
                "S1z": 0,
                "S2x": 0,
                "S2y": 0,
                "S2z": 0,
                "distance": cosmo.luminosity_distance(source.z).to("m").value,
                "inclination": 0,
                "phiRef": 0,
                "longAscNodes": 0,
                "eccentricity": 0,
                "meanPerAno": 0,
                "deltaF": default_deltaF.value,
                "f_min": default_f_min.value,
                "f_max": default_f_max.value,
                "f_ref": 0,
                "LALpars": {},
                "approximant": approx,
            }
        )

        if "m1" in lalsuite_kwargs.keys():
            waveform_dict["m1"] = lalsuite_kwargs["m1"]
        if "m2" in lalsuite_kwargs.keys():
            waveform_dict["m2"] = lalsuite_kwargs["m2"]
        if "S1x" in lalsuite_kwargs.keys():
            waveform_dict["S1x"] = lalsuite_kwargs["S1x"]
        if "S1y" in lalsuite_kwargs.keys():
            waveform_dict["S1y"] = lalsuite_kwargs["S1y"]

        if "S1z" in lalsuite_kwargs.keys():
            waveform_dict["S1z"] = lalsuite_kwargs["S1z"]
        elif "chi1" in lalsuite_kwargs.keys():
            waveform_dict["S1z"] = lalsuite_kwargs["chi1"]
        elif hasattr(source, "chi1"):
            waveform_dict["S1z"] = source.chi1

        if "S2x" in lalsuite_kwargs.keys():
            waveform_dict["S2x"] = lalsuite_kwargs["S2x"]
        if "S2y" in lalsuite_kwargs.keys():
            waveform_dict["S2y"] = lalsuite_kwargs["S2y"]

        if "S2z" in lalsuite_kwargs.keys():
            waveform_dict["S2z"] = lalsuite_kwargs["S2z"]
        elif "chi2" in lalsuite_kwargs.keys():
            waveform_dict["S2z"] = lalsuite_kwargs["chi2"]
        elif hasattr(source, "chi2"):
            waveform_dict["S2z"] = source.chi2

        if "distance" in lalsuite_kwargs.keys():
            waveform_dict["distance"] = lalsuite_kwargs["distance"]
        if "inclination" in lalsuite_kwargs.keys():
            waveform_dict["inclination"] = lalsuite_kwargs["inclination"]
        if "phiRef" in lalsuite_kwargs.keys():
            waveform_dict["phiRef"] = lalsuite_kwargs["phiRef"]
        if "longAscNodes" in lalsuite_kwargs.keys():
            waveform_dict["longAscNodes"] = lalsuite_kwargs["longAscNodes"]
        if "eccentricity" in lalsuite_kwargs.keys():
            waveform_dict["eccentricity"] = lalsuite_kwargs["eccentricity"]
        if "meanPerAno" in lalsuite_kwargs.keys():
            waveform_dict["meanPerAno"] = lalsuite_kwargs["meanPerAno"]

        if "deltaF" in lalsuite_kwargs.keys():
            waveform_dict["deltaF"] = lalsuite_kwargs["deltaF"]
        if "f_min" in lalsuite_kwargs.keys():
            if isinstance(lalsuite_kwargs["f_min"], u.Quantity):
                waveform_dict["f_min"] = lalsuite_kwargs["f_min"].value
            else:
                waveform_dict["f_min"] = lalsuite_kwargs["f_min"]
        elif hasattr(source, "f_min"):
            new_f_min = source.f_min / M_time
            waveform_dict["f_min"] = new_f_min.to("Hz").value

        if "f_max" in lalsuite_kwargs.keys():
            if isinstance(lalsuite_kwargs["f_max"], u.Quantity):
                waveform_dict["f_max"] = lalsuite_kwargs["f_max"].value
            else:
                waveform_dict["f_max"] = lalsuite_kwargs["f_max"]
        elif hasattr(source, "f_max"):
            new_f_max = source.f_max / M_time
            waveform_dict["f_max"] = new_f_max.to("Hz").value

        if "f_ref" in lalsuite_kwargs.keys():
            if isinstance(lalsuite_kwargs["f_ref"], u.Quantity):
                waveform_dict["f_ref"] = lalsuite_kwargs["f_ref"].value
            else:
                waveform_dict["f_ref"] = lalsuite_kwargs["f_ref"]
        if "LALpars" in lalsuite_kwargs.keys():
            waveform_dict["LALpars"] = lalsuite_kwargs["LALpars"]

        Get_Proper_Freq_Params(waveform_dict)

        return Get_LALSuite_Waveform(source, waveform_dict, out_frame=out_frame)
    else:
        raise ValueError(
            f"{approximant}, not an available waveform. Must select either pyPhenomD or one from lalsuite's FD waveform approximants"
        )


def Get_Proper_Freq_Params(waveform_dict):
    """Used to check if the calculated array is too big, and to assure f_min is included by adjusting ``deltaF``"""
    in_deltaF = waveform_dict["deltaF"]
    in_f_min = waveform_dict["f_min"]
    in_f_max = waveform_dict["f_max"]

    # Smallest deltaF before lalsuite can't shift t_coalescence to 0 (?)
    deltaF_min = 5e-10
    # Maximum array size before we found trouble with memory
    array_max = 2e5

    if in_deltaF != in_f_min:
        out_deltaF = in_f_min
        out_f_min = in_f_min
        if out_deltaF < deltaF_min:
            print("deltaF is too small, setting to min value we checked: ", deltaF_min)
            out_deltaF = deltaF_min
            out_f_min = deltaF_min
    elif in_deltaF < deltaF_min:
        print("deltaF is too small, setting to min value we checked: ", deltaF_min)
        out_deltaF = deltaF_min
        out_f_min = deltaF_min
    else:
        out_deltaF = in_deltaF
        out_f_min = in_f_min

    # Size based on np arange size
    arr_size = np.ceil((in_f_max - (out_f_min)) / (out_deltaF))
    if arr_size > array_max:
        errstr_1 = f"predicted frequency array size is very large: {arr_size}. "
        errstr_1 += "This will probably cause memory errors..."
        errstr_1 += "We will attempt to shrink it by reducing f_max."
        # raise ValueError(errstr_1)
        # print(errstr_1)
        new_arr_size = arr_size
        out_f_max = in_f_max
        scale = 2.0
        while new_arr_size > array_max:
            out_f_max = out_f_max / scale
            new_arr_size = np.ceil((out_f_max - (out_f_min)) / (out_deltaF))
            scale += 1
    else:
        out_f_max = in_f_max

    waveform_dict.update({"deltaF": out_deltaF, "f_min": out_f_min, "f_max": out_f_max})
    return waveform_dict


def Get_Fitcoeffs(source):
    """Loads Quasi-Normal Mode fitting files for speed later."""
    fit_coeffs_filedirectory = os.path.join(
        load_directory, "PhenomDFiles/fitcoeffsWEB.dat"
    )
    source._fitcoeffs = np.loadtxt(fit_coeffs_filedirectory)


def Get_Amp_Phase(h):
    """Separates the amplitude and phase from complex strain polarizations (``hcross``,``hplus``)."""
    amp = np.abs(h)
    phase = np.unwrap(np.angle(h))
    return amp, phase


def Get_Full_Amp(h_plus_f, h_cross_f):
    """Gets the raw amplitude from the plus and cross GW polarizations."""
    return np.sqrt((np.abs(h_cross_f)) ** 2 + (np.abs(h_plus_f)) ** 2)


def Strain_Conv(
    source, freqs, strain, inverse=False, in_frame="source", out_frame="observer"
):
    """Converts frequency and strain in natural units (G=c=1) to Hertz and raw Fourier strain amplitude (1/Hertz), respectively.
    If inverse is true, it does the reverse and assumes the strain and frequency are given in the detector frame.

    Parameters
    ----------
    source
        Instance of gravitational wave source class
    freqs: array
        the frequency of the source in either natural units (G=c=1) or Hertz
    strain: array
        the strain of the source in natural units (G=c=1) or raw Fourier strain amplitude (1/Hertz)
    inverse: bool, optional
        Converts non-naturalized (Hertz and dimensionless) frequency and strain to ``Mf`` and strain in G=c=1 units
    in_frame: str, {'source','observer'}
        If inverse is true, determines whether the source frequency ``f_gw`` is in the source or observer frame.
    out_frame: str, {'observer','source'}
        Determines whether the returned frequency is in the source or observer frame.
    """
    DL = cosmo.luminosity_distance(source.z)
    DL = DL.to("m")

    m_conv = const.G / const.c**3  # Converts M = [M] to M = [sec]
    M_time = source.M.to("kg") * m_conv
    # M_redshifted_time = source.M.to("kg") * (1 + source.z) * m_conv

    # frequency and strain of source in source frame
    freq_conv = 1 / M_time
    # Normalized factor to match Stationary phase approx at low frequencies
    strain_conv = np.sqrt(5 / 24 / np.pi) * (const.c / DL) * M_time**2
    if inverse:
        if in_frame == "source":
            # converted to source frame natural units
            conv_freqs = freqs / freq_conv
            conv_strain = strain / strain_conv
        elif in_frame == "observer":
            # converted to source frame natural units
            conv_freqs = (freqs / freq_conv) * (1 + source.z)
            conv_strain = strain / strain_conv / (1 + source.z) ** 2
        else:
            raise ValueError("The reference frame can only be observer or source.")
    else:
        if out_frame == "source":
            # converted to source frame physical units
            conv_freqs = freqs * freq_conv
            conv_strain = strain * strain_conv
        elif out_frame == "observer":
            # converted to detector frame physical units
            conv_freqs = (freqs * freq_conv) / (1 + source.z)
            conv_strain = strain * strain_conv * (1 + source.z) ** 2
        else:
            raise ValueError("The reference frame can only be observer or source.")

    return [conv_freqs, conv_strain]


def Get_LALSuite_Waveform(source, waveform_dict, out_frame="observer"):
    """Gets the frequency domain waveform of a particular source using a waveform modelled in
    LIGO's ``lalsuite`` package. The source is given plus and cross
    properties with their phase and amplitude for user convenience.

    Parameters
    ----------
    source: object
            source object from ``binary``, contains all source parameters
    waveform_dict: dictionary
            The dictionary is comprised of necessities for the ``SimInspiralChooseFDWaveform`` call in ``lalsimulation`` comprised of:
            ``m1``    mass of companion 1 (kg)
            ``m2``    mass of companion 2 (kg)
            ``S1x``    x-component of the dimensionless spin of object 1
            ``S1y``    y-component of the dimensionless spin of object 1
            ``S1z``    z-component of the dimensionless spin of object 1
            ``S2x``    x-component of the dimensionless spin of object 2
            ``S2y``    y-component of the dimensionless spin of object 2
            ``S2z``    z-component of the dimensionless spin of object 2
            ``distance``    distance of source (m)
            ``inclination``    inclination of source (rad)
            ``phiRef``    reference orbital phase (rad)
            ``longAscNodes``    longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation
            ``eccentricity``    eccentricity at reference epoch
            ``meanPerAno``    mean anomaly of periastron
            ``deltaF``    sampling interval (Hz)
            ``f_min``    starting GW frequency (Hz)
            ``f_max``    ending GW frequency (Hz)
            ``f_ref``    Reference frequency (Hz)
            ``LALparams``    LAL dictionary containing accessory parameters
            ``approximant``    post-Newtonian approximant to use for waveform production
    out_frame: str, {'observer','source'}
        Determines whether the returned frequency is in the source or observer frame.
    """
    h_f_plus, h_f_cross = lalsimulation.SimInspiralChooseFDWaveform(**waveform_dict)
    h_f_plus_amp, h_f_plus_phase = Get_Amp_Phase(h_f_plus.data.data)
    h_f_cross_amp, h_f_cross_phase = Get_Amp_Phase(h_f_cross.data.data)

    source.h_f_plus_amp = h_f_plus_amp
    source.h_f_plus_phase = h_f_plus_phase
    source.h_f_cross_amp = h_f_cross_amp
    source.h_f_cross_phase = h_f_cross_phase

    full_amp = Get_Full_Amp(h_f_plus_amp, h_f_cross_amp)
    # Need to trim because SimInspiralChooseFDWaveform returns nans and zeros outside of ranges (f_min,f_cutoff)
    trimmed_full_amp = []
    for amp in full_amp:
        if not np.isnan(amp) and amp != 0.0:
            trimmed_full_amp.append(amp)
    trimmed_full_amp = np.asarray(trimmed_full_amp)
    lin_freqs = np.arange(0, waveform_dict["f_max"], waveform_dict["deltaF"])
    trimmed_freqs = lin_freqs[1 : trimmed_full_amp.shape[0] + 1]

    # The difference between the two is `gwent` has a factor of :math:`\sqrt{\frac{1}{24}}(1+z)^{2}`
    # and `LALSuite` has a factor of :math:`2\sqrt{\frac{1}{64}}`.
    # Thus, :math:`h_{\mathrm{gwent}} = \sqrt{\frac{2}{3}}(1+z)^{2}h_{\mathrm{LAL}}`,
    # this factor is reduced by :math:`\sqrt{\frac{1}{2}}` if using the cross and plus polarizations to get the total Fourier strain amplitude.

    if out_frame == "observer":
        # frequency and strain of source in detector frame and physical units: Hertz and raw Fourier strain amplitude (1/Hertz)
        freqs = trimmed_freqs / (1 + source.z) * u.Hz
        strain = np.sqrt(1 / 3) * (1 + source.z) ** 2 * trimmed_full_amp / u.Hz
    elif out_frame == "source":
        # frequency and strain of source in source frame and physical units: Hertz and raw Fourier strain amplitude (1/Hertz)
        freqs = trimmed_freqs * u.Hz
        strain = np.sqrt(1 / 3) * trimmed_full_amp / u.Hz
    else:
        raise ValueError("The reference frame can only be observer or source.")

    return [freqs, strain]


def Get_PyPhenomD(source, pct_of_peak=0.01):
    """Uses Mass Ratio (``q`` <= 18), aligned spins (abs(a/m)~0.85 or when q=1 abs(a/m)<0.98),
    fitting coefficients for QNM type, and sampling rate
    Returns the frequency, the Phenom amplitude of the inspiral-merger-ringdown
    Uses methods found in <https://arxiv.org/abs/1508.07253> and <https://arxiv.org/abs/1508.07250>

    Parameters
    ----------
    source: object
            source object from ``binary``, contains all source parameters
    pct_of_peak: float, optional
            the percentange of the strain at merger that dictates the maximum frequency the waveform is calculated at in geometrized units (G=c=1)

    Returns
    -------
    Mf: numpy array of floats
        the waveform frequencies in geometrized units (G=c=1)
    fullwaveform: numpy array of floats
        the waveform strain in geometrized units (G=c=1)

    """
    f_min = source.f_min
    N = source.nfreqs
    q = source.q
    x1 = source.chi1
    x2 = source.chi2
    fitcoeffs = source._fitcoeffs

    # M = m1+m2 #Total Mass
    # q = m2/m1 #Mass Ratio: Paper tested up to 18
    # eta = m1*m2/M**2 reduced mass: Paper tested up to 0.05 (q=18)
    eta = q / (q + 1) ** 2
    x_PN = chi_PN(eta, x1, x2)  # PN reduced spin parameter
    a_f = a_final(x1, x2, q, eta)  # dimensionless spin

    ##################
    # Finds f_ringdown and f_damp from fit taken from <https://arxiv.org/abs/gr-qc/0512160>
    n = 0  # QNM indices
    ell = 2
    m = 2
    numn = 3  # number of n's included in the table

    index = (ell - 2) * (2 * ell + 1) * numn + (ell - m) * numn + n
    f_fit = fitcoeffs[index][3:6]
    q_fit = fitcoeffs[index][6:9]

    omega_RD = f_fit[0] + f_fit[1] * (1 - a_f) ** f_fit[2]  # M omega_{lmn}
    tau = (
        2 * (q_fit[0] + q_fit[1] * (1 - a_f) ** q_fit[2]) / omega_RD
    )  # tau_{lmn}/M = 2 Q_{lmn}/(M omega_{lmn})
    ########################
    f_RD = omega_RD / 2 / np.pi
    f_damp = 1 / tau / 2 / np.pi

    Gamma1 = Lambda(eta, x_PN, 4)
    Gamma2 = Lambda(eta, x_PN, 5)
    Gamma3 = Lambda(eta, x_PN, 6)

    f_peak = Calc_f_peak(f_RD, f_damp, [Gamma1, Gamma2, Gamma3])

    f1 = 0.014
    f3 = f_peak
    f2 = (f1 + f3) / 2

    cutoffFreq = Find_Cutoff_Freq(
        f_RD, f_damp, [Gamma1, Gamma2, Gamma3], pct_of_peak=pct_of_peak
    )

    # If lowest frequency is greater than cutoffFreq, then raise error.
    if f_min >= cutoffFreq:
        raise ValueError(
            "Lower frequency bound (ie. f_min) must be lower than that of the merger ringdown."
        )

    Mf = np.logspace(np.log10(f_min), np.log10(cutoffFreq), N)

    v1 = A_insp(f1, eta, x1, x2, x_PN)
    v2 = Lambda(eta, x_PN, 3)
    v3 = A_MR(f3, f_RD, f_damp, [Gamma1, Gamma2, Gamma3])
    fund1 = DA_insp(f1, eta, x1, x2, x_PN)
    fund3 = DA_MR(f3, f_RD, f_damp, [Gamma1, Gamma2, Gamma3])

    #############################
    # Calculate Solutions to eqn 21 in intermediate region
    Del_solns = A_intermediate(
        f1, f2, f3, v1, v2, v3, fund1, fund3
    )  # Solutions to eqn 21

    ##############################
    # Calculate all sections of waveform and Paste together
    indxf1 = np.argmin(np.abs(Mf - f1))
    indxfpeak = np.argmin(np.abs(Mf - f_peak))

    tmpinspiral = A_norm(Mf[0 : indxf1 + 1], eta) * A_insp(
        Mf[0 : indxf1 + 1], eta, x1, x2, x_PN
    )
    tmpintermediate = A_norm(Mf[indxf1 + 1 : indxfpeak], eta) * A_int(
        Mf[indxf1 + 1 : indxfpeak], Del_solns
    )
    tmpmergerringdown = A_norm(Mf[indxfpeak:], eta) * A_MR(
        Mf[indxfpeak:], f_RD, f_damp, [Gamma1, Gamma2, Gamma3]
    )
    fullwaveform = np.hstack((tmpinspiral, tmpintermediate, tmpmergerringdown))

    return [Mf, fullwaveform]


def A_norm(freqs, eta):
    """Calculates the constant scaling factor A_0

    Parameters
    ----------
    freqs: array
        The frequencies in Natural units (Mf, G=c=1) of the waveform
    eta: float
        The reduced mass ratio

    """
    const = np.sqrt(2 * eta / 3 / np.pi ** (1 / 3))
    return const * freqs ** -(7 / 6)


def A_insp(freqs, eta, x1, x2, X_PN):
    """Calculates the Inspiral Amplitude

    Parameters
    ----------
    freqs: array
        The frequencies in Natural units (``Mf``, G=c=1) of the waveform
    eta: float
        The reduced mass ratio
    x1: float
        The dimensionless spin parameter abs(a/m) for black hole ``m1``.
    x2: float
        The dimensionless spin parameter abs(a/m) for black hole ``m2``.
    x_PN: float
        The PN reduced spin parameter

    """
    A_PN = 0.0
    A_higher = 0.0
    for i in range(7):
        A_PN = A_PN + PN_coeffs(eta, x1, x2, i) * (np.pi * freqs) ** (i / 3)
        if i >= 1 and i <= 3:
            A_higher = A_higher + Lambda(eta, X_PN, i - 1) * freqs ** ((6 + i) / 3)
    return A_PN + A_higher


def DA_insp(freqs, eta, x1, x2, X_PN):
    """Calculates Derivative of the inspiral amplitude.

    Parameters
    ----------
    freqs: array
        The frequencies in Natural units (``Mf``, G=c=1) of the waveform
    eta: float
        The reduced mass ratio
    x1: float
        The dimensionless spin parameter abs(a/m) for black hole ``m1``.
    x2: float
        The dimensionless spin parameter abs(a/m) for black hole ``m2``.
    x_PN: float
        The PN reduced spin parameter

    """
    DA_PN = 0.0
    DA_higher = 0.0
    for i in range(7):
        PN_const = np.pi ** (i / 3) * (i / 3) * PN_coeffs(eta, x1, x2, i)
        DA_PN = DA_PN + PN_const * (freqs) ** ((i - 3) / 3)
        if i >= 1 and i <= 3:
            higher_const = ((6 + i) / 3) * Lambda(eta, X_PN, i - 1)
            DA_higher = DA_higher + higher_const * freqs ** ((i + 3) / 3)

    return DA_PN + DA_higher


def A_MR(freqs, f_RD, f_damp, Gammas):
    """Calculates the Normalized Merger-Ringdown Amplitude

    Parameters
    ----------
    freqs: array
        The frequencies in Natural units (``Mf``, G=c=1) of the waveform
    f_RD: float
        Frequency of the Ringdown transition
    f_damp: float
        Damping frequency
    Gammas: array-like
        Normalizes lorentzian to correct shape

    """
    varf = freqs - f_RD
    fg_d = Gammas[2] * f_damp
    return (
        (Gammas[0] * fg_d)
        / (varf**2 + fg_d**2)
        * np.exp(-(Gammas[1] / fg_d) * varf)
    )


def DA_MR(freqs, f_RD, f_damp, Gammas):
    """Calculates Derivative of the Merger-Ringdown Amplitude

    Parameters
    ----------
    freqs: array
        The frequencies in Natural units (``Mf``, G=c=1) of the waveform
    f_RD: float
        Frequency of the Ringdown transition
    f_damp: float
        Damping frequency
    Gammas: array-like
        Normalizes lorentzian to correct shape

    """
    varf = freqs - f_RD
    fg_d = Gammas[2] * f_damp
    A_MR_0 = A_MR(freqs, f_RD, f_damp, Gammas)
    return -A_MR_0 * (2 * varf / (varf**2 + fg_d**2) + Gammas[1] / fg_d)


def A_intermediate(f1, f2, f3, v1, v2, v3, d1, d3):
    """Solves system of equations for intermediate amplitude matching"""
    Mat = np.array(
        [
            [1.0, f1, f1**2, f1**3, f1**4],
            [1.0, f2, f2**2, f2**3, f2**4],
            [1.0, f3, f3**2, f3**3, f3**4],
            [0.0, 1.0, 2 * f1, 3 * f1**2, 4 * f1**3],
            [0.0, 1.0, 2 * f3, 3 * f3**2, 4 * f3**3],
        ],
        dtype="float",
    )
    a = np.array([v1, v2, v3, d1, d3], dtype="float")
    return np.linalg.solve(Mat, a)


def A_int(freqs, delt):
    """Calculates the Intermediate Amplitude

    Parameters
    ----------
    freqs: array
        The frequencies in Natural units (``Mf``, G=c=1) of the waveform
    delt: array
        Coefficient solutions to match the inspiral to the merger-ringdown portion of the waveform

    """
    return (
        delt[0]
        + delt[1] * freqs
        + delt[2] * freqs**2
        + delt[3] * freqs**3
        + delt[4] * freqs**4
    )


def Lambda(eta, x_PN, lmbda):
    """Gets the Lambdas from Eqn 31 in <https://arxiv.org/abs/1508.07253>

    Parameters
    ----------
    eta: float
        The reduced mass ratio
    x_PN: float
        The PN reduced spin parameter
    lmbda: int
        Iterator for different Lambda variables using the zeta function

    """
    xi = x_PN - 1
    xi2 = xi * xi
    xi3 = xi2 * xi
    eta2 = eta * eta
    if lmbda == 0:  # rho1
        coeffs = zeta(0)
    elif lmbda == 1:  # rho2
        coeffs = zeta(1)
    elif lmbda == 2:  # rho3
        coeffs = zeta(2)
    elif lmbda == 3:  # v2
        coeffs = zeta(3)
    elif lmbda == 4:  # gamma1
        coeffs = zeta(4)
    elif lmbda == 5:  # gamma2
        coeffs = zeta(5)
    elif lmbda == 6:  # gamma3
        coeffs = zeta(6)

    return (
        coeffs[0]
        + coeffs[1] * eta
        + (coeffs[2] + coeffs[3] * eta + coeffs[4] * eta2) * xi
        + (coeffs[5] + coeffs[6] * eta + coeffs[7] * eta2) * xi2
        + (coeffs[8] + coeffs[9] * eta + coeffs[10] * eta2) * xi3
    )


def zeta(k):
    """Coefficients in table 5 of <https://arxiv.org/abs/1508.07253>"""
    if k == 0:  # rho 1
        coeffs = [
            3931.9,
            -17395.8,
            3132.38,
            343966.0,
            -1.21626e6,
            -70698.0,
            1.38391e6,
            -3.96628e6,
            -60017.5,
            803515.0,
            -2.09171e6,
        ]
    elif k == 1:  # rho 2
        coeffs = [
            -40105.5,
            112253.0,
            23561.7,
            -3.47618e6,
            1.13759e7,
            754313.0,
            -1.30848e7,
            3.64446e7,
            596227.0,
            -7.42779e6,
            1.8929e7,
        ]
    elif k == 2:  # rho 3
        coeffs = [
            83208.4,
            -191238.0,
            -210916.0,
            8.71798e6,
            -2.69149e7,
            -1.98898e6,
            3.0888e7,
            -8.39087e7,
            -1.4535e6,
            1.70635e7,
            -4.27487e7,
        ]
    elif k == 3:  # v 2
        coeffs = [
            0.814984,
            2.57476,
            1.16102,
            -2.36278,
            6.77104,
            0.757078,
            -2.72569,
            7.11404,
            0.176693,
            -0.797869,
            2.11624,
        ]
    elif k == 4:  # gamma 1
        coeffs = [
            0.0069274,
            0.0302047,
            0.00630802,
            -0.120741,
            0.262716,
            0.00341518,
            -0.107793,
            0.27099,
            0.000737419,
            -0.0274962,
            0.0733151,
        ]
    elif k == 5:  # gamma 2
        coeffs = [
            1.01034,
            0.000899312,
            0.283949,
            -4.04975,
            13.2078,
            0.103963,
            -7.02506,
            24.7849,
            0.030932,
            -2.6924,
            9.60937,
        ]
    elif k == 6:  # gamma 3
        coeffs = [
            1.30816,
            -0.00553773,
            -0.0678292,
            -0.668983,
            3.40315,
            -0.0529658,
            -0.992379,
            4.82068,
            -0.00613414,
            -0.384293,
            1.75618,
        ]
    return coeffs


def PN_coeffs(eta, x1, x2, i):
    """Gets the PN Amplitude coefficients

    Parameters
    ----------
    eta: float
        The reduced mass ratio
    x1: float
        The dimensionless spin parameter abs(a/m) for black hole ``m1``.
    x2: float
        The dimensionless spin parameter abs(a/m) for black hole ``m2``.
    i: int
        iterator to dictate which PN Amplitude to use

    Notes
    -----
    Coefficients in appendix B (eqns B14-B20) of <https://arxiv.org/abs/1508.07253>

    """
    delta = np.sqrt(1.0 - 4.0 * eta)
    chi_s = (x1 + x2) / 2.0
    chi_a = (x1 - x2) / 2.0
    if i == 0:
        A_i = 1
    elif i == 1:
        A_i = 0
    elif i == 2:
        A_i = -323 / 224 + (451 / 168) * eta
    elif i == 3:
        A_i = (27 / 8) * delta * chi_a + (27 / 8 - (11 / 6) * eta) * chi_s
    elif i == 4:
        A_i = (
            -27312085 / 8128512
            - (1975055 / 338688) * eta
            + (105271 / 24192) * eta**2
            + (-81 / 32 + 8 * eta) * chi_a**2
            - 81 / 16 * delta * chi_a * chi_s
            + (-81 / 32 + 17 / 8 * eta) * chi_s**2
        )
    elif i == 5:
        A_i = (
            -85 * np.pi / 64
            + 85 * np.pi / 16 * eta
            + (285197 / 16128 - 1579 / 4032 * eta) * delta * chi_a
            + (285197 / 16128 - 15317 / 672 * eta - 2227 / 1008 * eta**2) * chi_s
        )
    elif i == 6:
        A_i = (
            -177520268561 / 8583708672
            + (545384828789 / 5007163392 - 205 * np.pi**2 / 48) * eta
            - 3248849057 / 178827264 * eta**2
            + 34473079 / 6386688 * eta**3
            + (1614569 / 64512 - 1873643 / 16128 * eta + 2167 / 42 * eta**2)
            * chi_a**2
            + (31 * np.pi / 12 - 7 * np.pi / 3 * eta) * chi_s
            + (1614569 / 64512 - 61391 / 1344 * eta + 57451 / 4032 * eta**2)
            * chi_s**2
            + delta
            * chi_a
            * (31 * np.pi / 12 + (1614569 / 32256 - 165961 / 2688 * eta) * chi_s)
        )
    return A_i


def Calc_f_peak(f_RD, f_damp, Gammas):
    """Calculates the frequency at the peak of the merger

    Parameters
    ----------
    f_RD: float
        Frequency of the Ringdown transition
    f_damp: float
        Damping frequency
    Gammas: array-like
        Normalizes lorentzian to correct shape

    Notes
    -----
    There is a problem with this expression from the paper becoming imaginary if ``gamma2`` >= 1
    so if ``gamma2`` >= 1 then set the square root term to zero.

    """
    if Gammas[1] <= 1:
        f_max = np.abs(
            f_RD + f_damp * Gammas[2] * (np.sqrt(1 - Gammas[1] ** 2) - 1) / Gammas[1]
        )
    else:
        f_max = np.abs(f_RD + (f_damp * Gammas[2] * -1) / Gammas[1])
    return f_max


def Find_Cutoff_Freq(f_RD, f_damp, Gammas, pct_of_peak=0.0001):
    """Cutoff signal when the amplitude is a factor of 10 below the value at f_RD

    Parameters
    ----------
    f_RD: float
        Frequency of the Ringdown transition
    f_damp: float
        Damping frequency
    Gammas: array-like
        Normalizes lorentzian to correct shape
    pct_of_peak: float, optional
        the percentange of the strain at merger that dictates the maximum
        frequency the waveform is calculated at in geometrized units (G=c=1)

    """
    tempfreqs = np.logspace(np.log10(f_RD), np.log10(10 * f_RD), 100)
    cutoffAmp = pct_of_peak * A_MR(
        f_RD, f_RD, f_damp, [Gammas[0], Gammas[1], Gammas[2]]
    )
    merger_ringdown_Amp = A_MR(
        tempfreqs, f_RD, f_damp, [Gammas[0], Gammas[1], Gammas[2]]
    )
    cutoffindex = np.argmin(np.abs(cutoffAmp - merger_ringdown_Amp))
    return tempfreqs[cutoffindex]


def a_final(x1, x2, q, eta):
    """The Final spin of the binary remnant black hole

    Parameters
    ----------
    x1: float
        The dimensionless spin parameter abs(a/m) for black hole ``m1``.
    x2: float
        The dimensionless spin parameter abs(a/m) for black hole ``m2``.
    q: float
        The mass ratio ``m1/m2``, ``m1<=m2``
    eta: float
        The reduced mass ratio

    Notes
    -----
    Uses eq. 3 in <https://arxiv.org/abs/0904.2577>, changed to match our q convention
    :math:`a=J/M^{2}` where ``J = x1*m1**2 + x2*m2**2``

    """
    a = (q**2 * x1 + x2) / (q**2 + 1)
    s4 = -0.1229
    s5 = 0.4537
    t0 = -2.8904
    t2 = -3.5171
    t3 = 2.5763
    return (
        a
        + s4 * a**2 * eta
        + s5 * a * eta**2
        + t0 * a * eta
        + 2 * np.sqrt(3) * eta
        + t2 * eta**2
        + t3 * eta**3
    )


def chi_PN(eta, x1, x2):
    """Calculates the PN reduced spin parameter

    Parameters
    ----------
    eta: float
        The reduced mass ratio
    x1: float
        The dimensionless spin parameter abs(a/m) for black hole ``m1``.
    x2: float
        The dimensionless spin parameter abs(a/m) for black hole ``m2``.

    Notes
    -----
    See Eq 5.9 in <https://arxiv.org/abs/1107.1267v2>

    """
    delta = np.sqrt(1.0 - 4.0 * eta)
    chi_s = (x1 + x2) / 2.0
    chi_a = (x1 - x2) / 2.0
    return chi_s * (1.0 - eta * 76.0 / 113.0) + delta * chi_a
