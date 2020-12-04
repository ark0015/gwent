import numpy as np
import os
import astropy.constants as const
import astropy.units as u
import scipy.interpolate as interp
from astropy.cosmology import z_at_value
from astropy.cosmology import WMAP9 as cosmo

import gwent
from . import waveform
from . import utils

current_path = os.path.abspath(gwent.__path__[0])
load_directory = os.path.join(current_path, "LoadFiles/")


class BinaryBlackHole:
    """Base Class for frequency domain strains from Binary Black Holes.

    Parameters
    ----------
    M : float
        Total mass of the black hole binary (m1+m2)
    q : float
        Mass ratio of the black hole binary (m1/m2, m1<m2)
    z : float
        Redshift of the black hole binary

    load_location : string, optional
        the directory of the loaded file, (ie. '/path/to/file')

    Notes
    -----
    IMRPhenomD waveforms calibrated for q = m1/m2 < 18

    """

    def __init__(self, *args, **kwargs):
        if len(args) == 3:
            [M, q, z] = args
        elif len(args) == 5:
            [M, q, z, _, _] = args
        else:
            raise ValueError(
                "args must be a list of 3 ([M,q,z]) or 5 ([M,q,z,chi1,chi2])"
            )
        self.M = M
        self.q = q
        self.z = z

        if "load_location" in kwargs.keys():
            self.load_location = kwargs["load_location"]

        if hasattr(self, "load_location"):
            self.Load_Data()

    @property
    def M(self):
        self._M = utils.make_quant(self._M, "M_sun")
        return self._M

    @M.setter
    def M(self, value):
        self.var_dict = ["M", value]
        self._M = self._return_value

    @property
    def q(self):
        return self._q

    @q.setter
    def q(self, value):
        self.var_dict = ["q", value]
        self._q = self._return_value

    @property
    def z(self):
        return self._z

    @z.setter
    def z(self, value):
        self.var_dict = ["z", value]
        self._z = self._return_value

    @property
    def h_f(self):
        if not hasattr(self, "_h_f"):
            raise NotImplementedError(
                "The strain must be defined inside BBHFrequencyDomain or BBHTimeDomain classes."
            )
        return self._h_f

    @h_f.setter
    def h_f(self, value):
        self._h_f = value

    @property
    def f(self):
        if not hasattr(self, "_f"):
            raise NotImplementedError(
                "The frequency must be defined inside BBHFrequencyDomain or BBHTimeDomain classes."
            )
        return self._f

    @f.setter
    def f(self, value):
        self._f = value

    @property
    def var_dict(self):
        return self._var_dict

    @var_dict.setter
    def var_dict(self, value):
        utils.Get_Var_Dict(self, value)

    def Load_Data(self):
        if hasattr(self, "load_location"):
            if os.path.exists(self.load_location):
                self._load_data = np.loadtxt(self.load_location)
            else:
                raise IOError(
                    "File %s does not exist, please assign load_location a correct filepath."
                    % self.load_location
                )
        else:
            raise ValueError(
                'load_location is not assigned, please set with name_of_BBH.load_location="path/to/file".'
            )


class BBHFrequencyDomain(BinaryBlackHole):
    """Subclass of BinaryBlackHole for BBH GWs generated in the frequency domain.

    Parameters
    ----------
    chi1 : float
        The dimensionless spin parameter abs(a/m) for black hole m1.
    chi2 : float
        The dimensionless spin parameter abs(a/m) for black hole m2

    f_low : float, optional
        The lowest frequency in natural units (Mf, G=c=1) at which the BBH waveform is calculated
    nfreqs : int, optional
        The number of frequencies at which the BBH waveform is calculated
    instrument: object, optional
        If assigned, the optimal frequency (ie. most sensitive frequency) of the detector is used as
        the binary's GW frequency

    Notes
    -----
    IMRPhenomD waveforms calibrated for aligned spins chi_1, chi_2 = abs(a/m) <= 0.85 or if q=1 abs(a/m)<0.98

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if len(args) == 5:
            [_, _, _, chi1, chi2] = args
            self.chi1 = chi1
            self.chi2 = chi2

        for keys, value in kwargs.items():
            if keys == "f_low":
                self.f_low = value
            elif keys == "f_high":
                self.f_high = value
            elif keys == "nfreqs":
                self.nfreqs = value
            elif keys == "instrument":
                self.instrument = value
                self.Check_Freq_Evol()
            elif keys == "f_gw":
                self.f_gw = value
            elif keys == "chi1":
                self.chi1 = value
            elif keys == "chi2":
                self.chi2 = value
            elif keys == "approximant":
                self.approximant = value
            elif keys == "lalsuite_kwargs":
                self.lalsuite_kwargs = value

        # If self has either chi1 or chi2, assign the other zero
        if any([hasattr(self, "chi1"), hasattr(self, "chi2")]) and not all(
            [hasattr(self, "chi1"), hasattr(self, "chi2")]
        ):
            if hasattr(self, "chi1"):
                self.chi2 = 0.0
            elif hasattr(self, "chi2"):
                self.chi1 = 0.0

        if not hasattr(self, "nfreqs"):
            self.nfreqs = int(1e3)
        if not hasattr(self, "f_low"):
            self.f_low = 1e-6
        if not hasattr(self, "approximant"):
            self.approximant = "pyPhenomD"
        if not hasattr(self, "lalsuite_kwargs"):
            self.lalsuite_kwargs = {}

        self.Get_Fitcoeffs()

    @property
    def chi1(self):
        return self._chi1

    @chi1.setter
    def chi1(self, value):
        self.var_dict = ["chi1", value]
        self._chi1 = self._return_value

    @property
    def chi2(self):
        return self._chi2

    @chi2.setter
    def chi2(self, value):
        self.var_dict = ["chi2", value]
        self._chi2 = self._return_value

    @property
    def instrument(self):
        return self._instrument

    @instrument.setter
    def instrument(self, value):
        self._instrument = value

    @property
    def h_gw(self):
        if not hasattr(self, "_h_gw"):
            if hasattr(self, "_instrument"):
                self._h_gw = Get_Mono_Strain(self, in_frame="source").to("")
            else:
                raise ValueError(
                    "No instrument assigned, please fix it. "
                    'Try: "source.instrument = instrument".'
                )
        return self._h_gw

    @h_gw.setter
    def h_gw(self, value):
        self._h_gw = value

    @h_gw.deleter
    def h_gw(self):
        del self._h_gw

    @property
    def f_gw(self):
        if not hasattr(self, "_f_gw"):
            if hasattr(self, "_instrument"):
                self._f_gw = self.instrument.f_opt
            else:
                raise ValueError(
                    "No GW frequency or instrument assigned, please fix it. "
                    "Try: assigning source.f_gw or source.instrument."
                )
        return self._f_gw

    @f_gw.setter
    def f_gw(self, value):
        self._f_gw = value

    @f_gw.deleter
    def f_gw(self):
        del self._f_gw

    @property
    def h_f(self):
        if not hasattr(self, "_h_f"):
            if self.approximant == "pyPhenomD":
                if not all([hasattr(self, "_phenomD_f"), hasattr(self, "_phenomD_h")]):
                    self.Get_PhenomD_Strain()
                [_, self._h_f] = waveform.Strain_Conv(
                    self, self._phenomD_f, self._phenomD_h
                )
            else:
                if not hasattr(self, "_h_f"):
                    [self._f, self._h_f] = waveform.Get_Waveform(
                        self,
                        approximant=self.approximant,
                        lalsuite_kwargs=self.lalsuite_kwargs,
                    )
                else:
                    [_, self._h_f] = waveform.Get_Waveform(
                        self,
                        approximant=self.approximant,
                        lalsuite_kwargs=self.lalsuite_kwargs,
                    )
        return self._h_f

    @h_f.deleter
    def h_f(self):
        del self._h_f

    @property
    def f(self):
        if not hasattr(self, "_f"):
            if self.approximant == "pyPhenomD":
                if not all([hasattr(self, "_phenomD_f"), hasattr(self, "_phenomD_h")]):
                    self.Get_PhenomD_Strain()
                [self._f, _] = waveform.Strain_Conv(
                    self, self._phenomD_f, self._phenomD_h
                )
            else:
                if not hasattr(self, "_h_f"):
                    [self._f, self._h_f] = waveform.Get_Waveform(
                        self,
                        approximant=self.approximant,
                        lalsuite_kwargs=self.lalsuite_kwargs,
                    )
                else:
                    [self._f, _] = waveform.Get_Waveform(
                        self,
                        approximant=self.approximant,
                        lalsuite_kwargs=self.lalsuite_kwargs,
                    )
        return self._f

    @f.deleter
    def f(self):
        del self._f

    def Get_Fitcoeffs(self):
        """Loads Quasi-Normal Mode fitting files for speed later."""
        fit_coeffs_filedirectory = os.path.join(
            load_directory, "PhenomDFiles/fitcoeffsWEB.dat"
        )
        self._fitcoeffs = np.loadtxt(fit_coeffs_filedirectory)

    def Get_PhenomD_Strain(self):
        """Gets the BBH's frequency and waveform from IMRPhenomD."""
        if not hasattr(self, "_fitcoeffs"):
            self.Get_Fitcoeffs()
        [self._phenomD_f, self._phenomD_h] = waveform.Get_PyPhenomD(self)


class BBHTimeDomain(BinaryBlackHole):
    """Subclass of BinaryBlackHole for input in the time domain"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.Get_hf_from_hcross_hplus()

    @property
    def t(self):
        if not hasattr(self, "_t"):
            self._t = self._load_data[:, 0]
        self._t = utils.make_quant(self._t, "s")
        return self._t

    @property
    def h_plus_t(self):
        if not hasattr(self, "_h_plus_t"):
            self._h_plus_t = self._load_data[:, 1]
        return self._h_plus_t

    @property
    def h_cross_t(self):
        if not hasattr(self, "_h_cross_t"):
            self._h_cross_t = self._load_data[:, 1]
        return self._h_cross_t

    @property
    def h_f(self):
        if not hasattr(self, "_h_f"):
            [natural_f, natural_h] = self.Get_hf_from_hcross_hplus()
            [_, self._h_f] = waveform.Strain_Conv(self, natural_f, natural_h)
        return self._h_f

    @h_f.deleter
    def h_f(self):
        del self._h_f

    @property
    def f(self):
        if not hasattr(self, "_f"):
            [natural_f, natural_h] = self.Get_hf_from_hcross_hplus()
            [self._f, _] = waveform.Strain_Conv(self, natural_f, natural_h)
        return self._f

    @f.deleter
    def f(self):
        del self._f

    def Get_hf_from_hcross_hplus(self, interp_res="coarse", windowing="left"):
        """Converts dimensionless, time domain strain to frequency space using a windowed fft

        Parameters
        ----------
        interp_res : {'coarse','fine'}, optional
            'coarse' uses maximum difference between subsequent time steps for interpolation
            'fine' uses minimum difference between subsequent time steps for interpolation
        windowing : {'left','right','all'}, optional
            'left' windows the left side of the time data
            'right' windows the right side of the time data
            'all' windows the both the left and right side of the time data

        Returns
        -------
        natural_f : array
            The frequency of the input source in natural units (G=c=1)
        natural_h : array
            The strain of the input source in natural units (G=c=1)

        """

        # Interpolate time to evenly sampled data, can be fine or coarse
        diff_t = np.diff(self.t.value)
        if interp_res == "fine":
            dt = min(diff_t)
        elif interp_res == "coarse":
            dt = max(diff_t)

        interp_t = np.arange(self.t[0].value, self.t[-1].value, dt)
        # interpolate strain to evenly sampled data for FFT
        h_cross_t = interp.interp1d(self.t, self.h_cross_t, kind="cubic")
        h_plus_t = interp.interp1d(self.t, self.h_plus_t, kind="cubic")
        interp_h_cross_t = h_cross_t(interp_t)
        interp_h_plus_t = h_plus_t(interp_t)

        # Filter/Window
        hann_window = np.hanning(len(interp_t))  # Two sided
        if windowing == "left":
            #########################
            """Applies window to first (left) half"""
            first_half = hann_window[
                : int(len(interp_t) / 2)
            ]  # Only need tapering on first half of waveform
            second_half = np.ones(
                len(interp_t) - len(first_half)
            )  # no windowing on second half of waveform
            #########################
            window = np.append(
                first_half, second_half
            )  # Only apply window to first half of waveform
        elif windowing == "right":
            #########################
            """Applies window to second (right) half"""
            second_half = hann_window[
                int(len(interp_t) / 2) :
            ]  # Only need tapering on second half of waveform
            first_half = np.ones(
                len(interp_t) - len(second_half)
            )  # no windowing on first half of waveform
            #########################
            window = np.append(first_half, second_half)
        elif windowing == "all":
            window = hann_window
        # Window!
        win_h_cross_t = np.multiply(interp_h_cross_t, window)
        win_h_plus_t = np.multiply(interp_h_plus_t, window)

        # FFT the two polarizations
        h_cross_f = np.fft.fft(win_h_cross_t)
        h_plus_f = np.fft.fft(win_h_plus_t)
        freqs = np.fft.fftfreq(len(interp_t), d=dt)

        # cut = np.abs(freqs).argmax() #Cut off the negative frequencies
        f_cut_low = 3e-3  # Low Cutoff frequency
        f_cut_high = 1.5e-1  # High Cutoff frequency
        cut_low = np.abs(
            freqs - f_cut_low
        ).argmin()  # Cut off frequencies lower than a frequency
        cut_high = np.abs(
            freqs - f_cut_high
        ).argmin()  # Cut off frequencies higher than a frequency
        # cut=int(len(freqs)*0.9) #Cut off percentage of frequencies
        h_cross_f = h_cross_f[cut_low:cut_high]
        h_plus_f = h_plus_f[cut_low:cut_high]
        natural_f = freqs[cut_low:cut_high]

        # Combine them for raw spectral power
        natural_h_f = np.sqrt((np.abs(h_cross_f)) ** 2 + (np.abs(h_plus_f)) ** 2)
        return [natural_f, natural_h_f]


def Get_Char_Strain(source):
    """Converts source strain to characteristic strain

    Parameters
    ----------
    source : object
        Instance of gravitational wave source class

    """
    h_char = np.sqrt(4 * source.f ** 2 * source.h_f ** 2)
    return h_char


def Get_Mono_Char_Strain(
    source, freq, in_frame="observer", out_frame="observer", inc=None
):
    """Calculates the characteristic strain from a slowly evolving (monochromatic) binary black hole

    Parameters
    ----------
    source : object
        Instance of gravitational wave source class
    in_frame : str, {'observer','source'}
        Determines whether the source frequency f_gw is in the source or observer frame. 
    out_frame : str, {'observer','source'}
        Determines whether the returned frequency is in the source or observer frame.
    inc : float, optional
        The inclination of the source in radians. If inc is None, the strain is \
        sky and inclination averaged strain from Robson et al. 2019 (eqn 27) <https://arxiv.org/pdf/1803.01944.pdf> \

    Returns
    -------
    float
        The characteristic strain of a monochromatic source in the source frame.
    """

    strain_amp = Get_Mono_Strain(
        source, inc=inc, in_frame=in_frame, out_frame=out_frame
    )
    f_dot = source.Get_F_Dot(freq, in_frame=in_frame, out_frame=out_frame)

    if out_frame == "observer":
        if in_frame == "source":
            return source.f_gw / (1 + source.z) * np.sqrt(2 / f_dot) * strain_amp
        else:
            return source.f_gw * np.sqrt(2 / f_dot) * strain_amp
    else:
        if in_frame == "source":
            return source.f_gw * np.sqrt(2 / f_dot) * strain_amp
        else:
            return source.f_gw * (1 + source.z) * np.sqrt(2 / f_dot) * strain_amp


def Get_Mono_Strain(source, freq, inc=None, in_frame="source", out_frame="source"):
    """Calculates the fourier domain strain from a monochromatic binary black hole.

    Parameters
    ----------
    source : object
        Instance of gravitational wave source class
    freq : float
        the binary GW frequency in the corresponding in frame.
    inc : float, optional
        The inclination of the source in radians. If inc is None, the strain is \
        sky and inclination averaged strain from Robson et al. 2019 (eqn 27) <https://arxiv.org/pdf/1803.01944.pdf> \
    in_frame : str, {'source','observer'}
        Determines whether the frequency is in the source or observer frame.
    out_frame : str, {'observer','source'}
        Determines whether the returned frequency is in the source or observer frame.

    Returns
    -------
    float
        The strain of a monochromatic source in the source frame.

    """
    f_gw = utils.make_quant(freq, "Hz")

    if in_frame == "observer":
        f_obs_source = f_gw * (1 + source.z)
    elif in_frame == "source":
        f_obs_source = f_gw
    else:
        raise ValueError("The reference frame can only be observer or source.")

    DL = cosmo.luminosity_distance(source.z)
    DL = DL.to("m")

    # Converts M = [M] to M = [sec]
    m_conv = const.G / const.c ** 3

    eta = source.q / (1 + source.q) ** 2
    M_time = source.M.to("kg") * m_conv
    M_chirp = eta ** (3 / 5) * M_time

    if inc is not None:
        if inc > np.pi or inc < -np.pi:
            raise ValueError("Inclination must be between -pi and pi.")
        a = (1 + np.cos(inc) ** 2) / 2
        b = np.cos(inc)
        const_val = 4.0 * np.sqrt(a ** 2 + b ** 2)
    else:
        const_val = 8 / np.sqrt(5)

    if out_frame == "observer":
        return (
            const_val
            * (const.c / DL)
            * (np.pi * f_obs_source / (1 + source.z)) ** (2.0 / 3.0)
            * (M_chirp * (1 + source.z)) ** (5.0 / 3.0)
        )
    elif out_frame == "source":
        return (
            const_val
            * (const.c / DL)
            * (np.pi * f_obs_source) ** (2.0 / 3.0)
            * M_chirp ** (5.0 / 3.0)
        )
    else:
        raise ValueError("The reference frame can only be observer or source.")


def Get_F_Dot(source, freq, in_frame="observer", out_frame="observer"):
    """Calculates the change in frequency of a binary black hole at a given frequency.

    Parameters
    ----------
    source : object
        Instance of gravitational wave source class
    freq : float
        the binary GW frequency in the corresponding in frame.
    in_frame : str, {'observer','source'}
        Determines whether the given frequency is in the source or observer frame.
    out_frame : str, {'observer','source'}
        Determines whether the returned frequency is in the source or observer frame.

    """
    freq = utils.make_quant(freq, "Hz")
    m_conv = const.G / const.c ** 3  # Converts M = [M] to M = [sec]
    eta = source.q / (1 + source.q) ** 2

    # Always assume the mass is in the source frame
    M_time = source.M.to("kg") * m_conv
    M_chirp = eta ** (3 / 5) * M_time

    if in_frame == "observer":
        f_obs_source = freq * (1 + source.z)
    elif in_frame == "source":
        f_obs_source = freq
    else:
        raise ValueError("The reference frame can only be observer or source.")

    if out_frame == "observer":
        return (
            96
            / 5
            * np.pi ** (8 / 3)
            * ((M_chirp * (1 + source.z)) ** (5 / 3))
            * ((f_obs_source / (1 + source.z)) ** (11 / 3))
        )
    elif out_frame == "source":
        return (
            96
            / 5
            * np.pi ** (8 / 3)
            * (M_chirp ** (5 / 3))
            * (f_obs_source ** (11 / 3))
        )
    else:
        raise ValueError("The reference frame can only be observer or source.")


def Get_Time_From_Merger(source, freq, in_frame="observer", out_frame="source"):
    """Calculates the time from merger of a binary black hole given a frequency.

    Parameters
    ----------
    source : object
        Instance of gravitational wave source class
    freq : float
        the binary GW frequency in the corresponding in frame.
    in_frame : str, {'observer','source'}
        Determines whether the given frequency is in the source or observer frame.
    out_frame : str, {'source','observer'}
        Determines whether the returned frequency is in the source or observer frame.

    """
    freq = utils.make_quant(freq, "1/s")
    m_conv = const.G / const.c ** 3  # Converts M = [M] to M = [sec]
    eta = source.q / (1 + source.q) ** 2

    # Always assume the mass is in the source frame
    M_time = source.M.to("kg") * m_conv
    M_chirp = eta ** (3 / 5) * M_time

    if in_frame == "observer":
        f_obs_source = freq * (1 + source.z)
    elif in_frame == "source":
        f_obs_source = freq
    else:
        raise ValueError("The reference frame can only be observer or source.")

    if out_frame == "observer":
        return (
            5
            / 256
            * M_chirp
            * (1 + source.z)
            * (np.pi * M_chirp * f_obs_source) ** (-8 / 3)
        )
    elif out_frame == "source":
        return 5 / 256 * M_chirp * (np.pi * M_chirp * f_obs_source) ** (-8 / 3)
    else:
        raise ValueError("The reference frame can only be observer or source.")


def Get_Source_Freq(source, tau, in_frame="observer"):
    """Calculates the binary black hole's gravitational wave frequency given a time from merger in the source frame

    Parameters
    ----------
    source : object
        Instance of gravitational wave source class
    tau : int, float, Quantity
        the time from merger in the respective frame.
        If not an astropy quantity, assumed to be a time in seconds
    in_frame : str, {'observer','source'}
        Determines whether the given frequency is in the source or observer frame.

    """
    tau = utils.make_quant(tau, "s")
    if in_frame == "observer":
        tau_source = tau / (1 + source.z)
    elif in_frame == "source":
        tau_source = tau
    else:
        raise ValueError("The reference frame can only be observer or source.")

    m_conv = const.G / const.c ** 3  # Converts M = [M] to M = [sec]
    eta = source.q / (1 + source.q) ** 2

    M_time = source.M.to("kg") * m_conv
    M_chirp = eta ** (3 / 5) * M_time

    return 1.0 / 8.0 / np.pi / M_chirp * (5 * M_chirp / tau_source) ** (3.0 / 8.0)


def Check_Freq_Evol(source, T_evol=None, T_evol_frame="observer", f_gw_frame="source"):
    """Checks the frequency evolution of the black hole binary.

    Parameters
    ----------
    T_evol : int,float, Quantity, optional
        The length of time the binary may evolve, if not provided it is assumed source
        has an instrument assigned and the instruments observation time will be used.
    T_evol_frame : str, {'observer','source'}
        Determines whether the given T_evol is in the source or observer frame.
    f_gw_frame : str, {'observer','source'}
        Determines whether the frequency is in the source or observer frame.
        May not be used if source has an instrument assigned.

    Notes
    -----
    If the frequency of the binary does evolve over more than one bin,
    (ie f(T_evol)-f(t_init) = delf_obs < 1/T_evol), it is monochromatic, so we set the frequency
    to the optimal frequency of the detector

    Otherwise it is chirping and evolves over the observation and we
    set the starting frequency we observe it at to f(Tobs), which is the
    frequency at an observation time before merger

    To get the change in frequency, we use eqn 41 from Hazboun,Romano, and Smith (2019) <https://arxiv.org/abs/1907.04341>
    which uses binomial expansion of f_T_evol_inst - f_init_inst and thus will never be imaginary

    """

    m_conv = const.G / const.c ** 3  # Converts M = [M] to M = [sec]
    eta = source.q / (1 + source.q) ** 2

    M_time = source.M.to("kg") * m_conv
    M_chirp_source = eta ** (3 / 5) * M_time

    if T_evol is not None:
        T_evol = utils.make_quant(T_evol, "s")
        if hasattr(source, "_instrument"):
            if hasattr(source, "_f_gw"):
                t_init_source = Get_Time_From_Merger(
                    source, source.f_gw, in_frame=f_gw_frame
                )
            else:
                # Assumes f_init is the optimal frequency in the instrument frame to get t_init_source
                t_init_source = Get_Time_From_Merger(
                    source, source.instrument.f_opt, frame="observer"
                )
            # f(T_evol), the frequency of the source at T_evol before merger
            f_T_evol_source = Get_Source_Freq(source, T_evol, in_frame=T_evol_frame)
        elif hasattr(source, "_f_gw") and not hasattr(source, "_instrument"):
            t_init_source = Get_Time_From_Merger(
                source, source.f_gw, in_frame=f_gw_frame
            )
            # f(T_evol), the frequency of the source at T_evol before merger
            f_T_evol_source = Get_Source_Freq(source, T_evol, in_frame=T_evol_frame)
        else:
            raise ValueError(
                "Must either assign T_evol a value, or assign the source an instrument."
            )
    else:
        if hasattr(source, "instrument"):
            T_evol_frame == "observer"
            T_evol = utils.make_quant(np.max(np.unique(source.instrument.T_obs)), "s")
            # Assumes f_init is the optimal frequency in the instrument frame to get t_init_source
            t_init_source = Get_Time_From_Merger(
                source, source.instrument.f_opt, in_frame="observer"
            )
            # f(T_evol), the frequency of the source at T_evol before merger
            f_T_evol_source = Get_Source_Freq(source, T_evol, in_frame="observer")
        else:
            raise ValueError(
                "Must either assign T_evol a value, or assign the source an instrument."
            )

    # Assumes t_init is in source frame, can either be randomly drawn
    # t_init_source = np.random.uniform(0,100)*u.yr

    if T_evol_frame == "observer":
        T_obs = T_evol
        T_evol_source = T_evol / (1 + source.z)
    elif T_evol_frame == "source":
        T_obs = T_evol * (1 + source.z)
        T_evol_source = T_evol
    else:
        raise ValueError("The reference frame can only be observer or source.")

    # f(T_evol) in the instrument frame (ie. the observed frequency)
    source.f_T_obs = f_T_evol_source / (1 + source.z)

    delf_source = (
        1.0
        / 8.0
        / np.pi
        / M_chirp_source
        * (5 * M_chirp_source / t_init_source) ** (3.0 / 8.0)
        * (3 * T_evol_source / 8 / t_init_source)
    )
    delf_obs = delf_source / (1 + source.z)

    if delf_obs < (1 / T_obs):
        source.ismono = True
    else:
        source.ismono = False
