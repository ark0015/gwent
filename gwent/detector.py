import copy
import os
import sys

import astropy.constants as const
import astropy.units as u
import hasasia.sensitivity as hassens
import hasasia.sim as hassim
import numpy as np
import scipy.stats as stats

import gwent

from . import utils

current_path = os.path.abspath(gwent.__path__[0])
# sys.path.insert(0, current_path + "/vendor/pygwinc_clone")
import gwinc

load_directory = os.path.join(current_path, "LoadFiles/")


class PTA:
    r"""
    Class to make a PTA instrument using the methods of Hazboun, Romano, Smith (2019) <https://arxiv.org/abs/1907.04341>

    Parameters
    ----------

    name: string
        name of the instrument

    n_p: int
        the number of pulsars in the PTA

    T_obs: float,array, list, optional
        the observation time of the PTA in [years]
        If ``T_obs`` is a float, the observation time is used as the observation time for all pulsars
        If ``T_obs`` is a list of length of ``n_p``, the observation times are used as the corresponding pulsar observation times.
        If ``T_obs`` is a list of length 2, it is assumed the values are the minimum and maximum observation time values
        (ie. [min,max]) from which individual pulsar observation times are uniformly sampled.
    sigma: float, array, list, optional
        the rms error on the pulsar TOAs in [sec]
        If ``sigma`` is a float, the given rms error is used for all pulsars.
        If ``sigma`` is a list of length of n_p, the amplitudes are used as the corresponding pulsar rms error.
        If ``sigma`` is a list of length 2, it is assumed the values are the minimum and maximum rms errors values
        (ie. [min,max]) from which individual pulsar rms errors are uniformly sampled.
    cadence: float, array, list, optional
        How often the pulsars are observed in [num/year]
        If ``cadence`` is a float, the given cadence is used for all pulsars.
        If ``cadence`` is a list of length of n_p, the amplitudes are used as the corresponding pulsar cadence.
        If ``cadence`` is a list of length 2, it is assumed the values are the minimum and maximum cadence values
        (ie. [min,max]) from which individual pulsar cadences are uniformly sampled.
    sb_amp: float, optional
        Amplitude of the Stochastic background added as red noise
    sb_alpha: float, optional
        the Stochastic background power law, if empty and sb_amp is set, it is assumed to be -2/3
    rn_amp: array, list, optional
        Individual pulsar red noise amplitude.
        If ``rn_amp`` is a list of length of ``n_p``, the amplitudes are used as the corresponding pulsar RN injection.
        If ``rn_amp`` is a list of length 2, it is assumed the values are the minimum and maximum RN amplitude values
        (ie. [min,max]) from which individual pulsar RN amplitudes are uniformly sampled.
    rn_alpha: array, list, optional
        Individual pulsar red noise alpha (power law spectral index).
        If ``rn_alpha`` is a list of length of ``n_p``, the alpha indices are used as the corresponding pulsar RN injection.
        If ``rn_alpha`` is a list of length 2, it is assumed the values are the minimum and maximum RN alpha values
        (ie. [min,max]) from which individual pulsar RN alphas are uniformly sampled.
    phi: list, array, optional
        Individual pulsar longitude in ecliptic coordinates.
        If not defined, NANOGrav 11yr pulsar locations are used.
        If ``n_p`` > 34 (the number of pulsars in the 11yr dataset),
        it draws more pulsars from distributions based on the NANOGrav 11yr pulsars.
    theta: array, list, optional
        Individual pulsar colatitude in ecliptic coordinates.
        If not defined, NANOGrav 11yr pulsar locations are used.
        If ``n_p`` > 34 (the number of pulsars in the 11yr dataset),
        it draws more pulsars from distributions based on the NANOGrav 11yr pulsars.
    use_11yr: bool, optional
        Uses the NANOGrav 11yr noise as the individual pulsar noises,
        if ``n_p`` > 34 (the number of pulsars in the 11yr dataset),
        it draws more pulsars from distributions based on the NANOGrav 11yr pulsar noise
    use_rn: bool, optional
        If no rn_amp assigned, uses the NANOGrav 11yr noise as the individual pulsar RN noises,
        if ``n_p`` > 34 (the number of pulsars in the 11yr dataset),
        it draws more pulsars from distributions based on the NANOGrav 11yr pulsar noise
    load_location: string, optional
        If you want to load a PTA curve from a file, it's the file path
    I_type: string, {'E','A','h'}
        Sets the type of input data.
        'E' is the effective strain spectral density :math:`S_{n}(f)` ('ENSD'),
        'A' is the amplitude spectral density, :math:`\sqrt{S_{n}(f)}` ('ASD'),
        'h' is the characteristic strain :math:`h_{n}(f)` ('h')
    f_min: float, optional
        Assigned lowest frequency of PTA (default assigns 1/(5*T_obs))
    f_max: float, optional
        Assigned highest frequency of PTA (default is Nyquist freq cadence/2)
    nfreqs: int, optional
        Number of frequencies in logspace the sensitivity is calculated
    nbins: int, optional
        Used to add values to every bin for sampled parameters. Default is 8 for smooth, non-zero distributions.
        Changing this could change distribution, so be wary, not sure how much it affects anything.
    """

    def __init__(self, name, *args, **kwargs):
        self.name = name
        if len(args) != 0:
            if len(args) == 1:
                self.n_p = args[0]
            else:
                raise ValueError("Too many args, not enough kwargs!")

        for keys, value in kwargs.items():
            if keys == "T_obs":
                self.T_obs = value
            elif keys == "sigma":
                self.sigma = value
            elif keys == "cadence":
                self.cadence = value
            elif keys == "rn_amp":
                self.rn_amp = value
            elif keys == "rn_alpha":
                self.rn_alpha = value
            elif keys == "phi":
                self.phi = value
            elif keys == "theta":
                self.theta = value
            elif keys == "sb_amp":
                self.sb_amp = value
            elif keys == "sb_alpha":
                self.sb_alpha = value
            elif keys == "use_11yr":
                self.use_11yr = value
            elif keys == "use_rn":
                self.use_rn = value
            elif keys == "load_location":
                self.load_location = value
            elif keys == "I_type":
                self.I_type = value
            elif keys == "f_min":
                self.f_min = utils.make_quant(value, "Hz")
            elif keys == "f_max":
                self.f_max = utils.make_quant(value, "Hz")
            elif keys == "nfreqs":
                self.nfreqs = value
            elif keys == "nbins":
                self.nbins = value
            else:
                raise ValueError("%s is not an accepted input option." % keys)

        if not hasattr(self, "nfreqs"):
            self.nfreqs = int(1e3)
        if not hasattr(self, "nbins"):
            self.nbins = 8
        if hasattr(self, "load_location"):
            Load_Data(self)

        if not hasattr(self, "use_11yr"):
            self.use_11yr = False
        if not hasattr(self, "use_rn"):
            self.use_rn = False

        if hasattr(self, "f_min") and hasattr(self, "f_max"):
            self.fT = (
                np.logspace(
                    np.log10(self.f_min.value), np.log10(self.f_max.value), self.nfreqs
                )
                * u.Hz
            )

    @property
    def n_p(self):
        return self._n_p

    @n_p.setter
    def n_p(self, value):
        self.var_dict = ["n_p", value]
        self._n_p = self._return_value

    @property
    def T_obs(self):
        return self._T_obs

    @T_obs.setter
    def T_obs(self, value):
        self.var_dict = ["T_obs", value]
        if not isinstance(self._return_value, u.Quantity):
            self._return_value = utils.make_quant(self._return_value, "yr")
        self._T_obs = self._return_value

    @property
    def cadence(self):
        return self._cadence

    @cadence.setter
    def cadence(self, value):
        self.var_dict = ["cadence", value]
        if not isinstance(self._return_value, u.Quantity):
            self._return_value = utils.make_quant(self._return_value, "1/yr")
        self._cadence = self._return_value

    @property
    def sigma(self):
        return self._sigma

    @sigma.setter
    def sigma(self, value):
        self.var_dict = ["sigma", value]
        if not isinstance(self._return_value, u.Quantity):
            self._return_value = utils.make_quant(self._return_value, "s")
        self._sigma = self._return_value

    @property
    def phi(self):
        return self._phi

    @phi.setter
    def phi(self, value):
        self.var_dict = ["phi", value]
        self._phi = self._return_value

    @property
    def theta(self):
        return self._theta

    @theta.setter
    def theta(self, value):
        self.var_dict = ["theta", value]
        self._theta = self._return_value

    @property
    def rn_amp(self):
        return self._rn_amp

    @rn_amp.setter
    def rn_amp(self, value):
        self.var_dict = ["rn_amp", value]
        self._rn_amp = self._return_value

    @property
    def rn_alpha(self):
        return self._rn_alpha

    @rn_alpha.setter
    def rn_alpha(self, value):
        self.var_dict = ["rn_alpha", value]
        self._rn_alpha = self._return_value

    @property
    def var_dict(self):
        return self._var_dict

    @var_dict.setter
    def var_dict(self, value):
        utils.Get_Var_Dict(self, value)

    @property
    def fT(self):
        if not hasattr(self, "_fT"):
            # frequency sampled from 1/observation time to nyquist frequency (c/2)
            # 5 is the default value for now (from Hazboun et al. 2019)
            if not hasattr(self, "_T_obs") or not hasattr(self, "_cadence"):
                self.Init_PTA()
            T_obs_sec = np.max(self.T_obs).to("s").value
            cadence_sec = np.max(self.cadence).to("1/s").value
            self._fT = (
                np.logspace(
                    np.log10(1 / (5 * T_obs_sec)),
                    np.log10(cadence_sec / 2),
                    self.nfreqs,
                )
                * u.Hz
            )
        return self._fT

    @fT.setter
    def fT(self, value):
        self._fT = value

    @fT.deleter
    def fT(self):
        del self._fT

    @property
    def h_n_f(self):
        """Effective Strain Noise Amplitude"""
        if not hasattr(self, "_h_n_f"):
            if hasattr(self, "_I_data"):
                if self._I_Type == "h":
                    self._h_n_f = self._I_data[:, 1]
                elif self._I_Type == "ENSD":
                    self._h_n_f = np.sqrt(self.S_n_f * self.fT)
                elif self._I_Type == "ASD":
                    S_n_f_sqrt = self._I_data[:, 1]
                    self._h_n_f = S_n_f_sqrt * np.sqrt(self.fT.value)
            else:
                if not hasattr(self, "_sensitivitycurve"):
                    self.Init_PTA()
                self._h_n_f = self._sensitivitycurve.h_c
        return self._h_n_f

    @h_n_f.setter
    def h_n_f(self, value):
        self._h_n_f = value

    @h_n_f.deleter
    def h_n_f(self):
        del self._h_n_f

    @property
    def S_n_f(self):
        """Effective noise power amplitude"""
        if not hasattr(self, "_S_n_f"):
            if hasattr(self, "_I_data"):
                if self._I_Type == "ASD":
                    S_n_f_sqrt = self._I_data[:, 1]
                    self._S_n_f = S_n_f_sqrt**2 / u.Hz
                elif self._I_Type == "ENSD":
                    self._S_n_f = self._I_data[:, 1] / u.Hz
                elif self._I_Type == "h":
                    self._S_n_f = self.h_n_f**2 / self.fT
            else:
                if not hasattr(self, "_sensitivitycurve"):
                    self.Init_PTA()
                self._S_n_f = self._sensitivitycurve.S_eff
                self._S_n_f = utils.make_quant(self._S_n_f, "1/Hz")
        return self._S_n_f

    @S_n_f.setter
    def S_n_f(self, value):
        self._S_n_f = value

    @S_n_f.deleter
    def S_n_f(self):
        del self._S_n_f

    @property
    def f_opt(self):
        """The optimal frequency of the instrument ie. the frequecy at the lowest strain noise PSD"""
        self._f_opt = self.fT[np.argmin(self.S_n_f)]
        return self._f_opt

    def Load_NANOGrav_11yr_Params(self):
        """Loads in NANOGrav 11yr data

        Notes
        -----
        The file is in the form of observation times (``T_obs``) in the first column,
        sky locations (``phi``,``theta``) in the second and third columns,
        Individual Pulsar cadences and WN RMS (sigmas) in the fourth and fifth,
        RN Amplitudes, and RN Alphas in the last two columns.
        """
        NANOGrav_11yr_params_filedirectory = os.path.join(
            load_directory, "InstrumentFiles/NANOGrav/NANOGrav_11yr_params.txt"
        )
        self._NANOGrav_11yr_params = np.loadtxt(NANOGrav_11yr_params_filedirectory)

    def Get_Param_Distributions(self, var_name, NG_11yr_idx):
        """Gets the noise parameter values (``sigma``, ``rn_amp``, ``rn_alpha``) and sky locations (``phis``, ``thetas``)
        and generates populated arrays from which distributions can be made. If no user values for a param are given,
        it uses the NANOGrav 11yr parameters.

        var_name: string
            The name of the noise parameter to assign sampled parameters
        NG_11yr_idx: int
            Index of corresponding value in NANOGrav 11yr params
        """

        if not hasattr(self, var_name):
            samp_var = self._NANOGrav_11yr_params[NG_11yr_idx]
            if var_name == "phi":
                # Add non-zero probability of picking 0 and 2pi
                return np.append(samp_var, np.linspace(0.0, 2 * np.pi, self.nbins))
            elif var_name == "theta":
                # Add non-zero probability of picking 0 and pi
                return np.append(samp_var, np.linspace(0.0, np.pi, self.nbins))
            elif var_name == "rn_amp":
                return np.append(
                    samp_var,
                    np.logspace(
                        min(np.log10(samp_var)),
                        max(np.log10(samp_var)),
                        self.nbins,
                    ),
                )
            else:
                return np.append(
                    samp_var, np.linspace(min(samp_var), max(samp_var), self.nbins)
                )
        else:
            var = getattr(self, var_name)
            if isinstance(var, u.Quantity):
                var = var.value
            if isinstance(var, (list, np.ndarray)):
                if not self.var_dict[var_name]["sampled"]:
                    if len(var) == self.n_p:
                        return var
                    elif len(var) == 1:
                        return np.ones(self.n_p) * var
                    else:
                        if self.var_dict["n_p"]["sampled"]:
                            unique_vals = np.unique(var)
                            if len(unique_vals) == 1:
                                return unique_vals[0]
                            else:
                                if var_name == "rn_amp":
                                    return np.append(
                                        var,
                                        np.logspace(
                                            min(np.log10(unique_vals)),
                                            max(np.log10(unique_vals)),
                                            self.nbins,
                                        ),
                                    )
                                else:
                                    return np.append(
                                        var,
                                        np.linspace(
                                            min(unique_vals),
                                            max(unique_vals),
                                            self.nbins,
                                        ),
                                    )
                        else:
                            raise ValueError(
                                f"{var_name} must be a single value, or the same length as n_p: {self.n_p}"
                            )
                else:
                    if len(var) == 2:
                        # Uniformly sample in logspace
                        if var_name == "rn_amp":
                            samp_var = np.random.uniform(
                                np.log10(var[0]), np.log10(var[1]), size=self.n_p
                            )
                        else:
                            # Uniformly Sample in linspace
                            samp_var = np.random.uniform(var[0], var[1], size=self.n_p)
                    elif len(var) == self.n_p:
                        samp_var = var
                    else:
                        raise ValueError(
                            f"To sample {var_name}, it must be either [min,max], or an array of individual pulsar {var_name} of length n_p: {self.n_p}"
                        )

                    if var_name == "rn_amp":
                        add_var = np.logspace(min(samp_var), max(samp_var), self.nbins)
                        samp_var = 10**samp_var
                        return np.append(samp_var, add_var)
                    else:
                        return np.append(
                            samp_var,
                            np.linspace(min(samp_var), max(samp_var), self.nbins),
                        )
            else:
                if var_name in self.var_dict.keys():
                    if not self.var_dict[var_name]["sampled"]:
                        return np.ones(self.n_p) * var
                else:
                    self.var_dict[var_name]["sampled"] = True
                    samp_var = self._NANOGrav_11yr_params[NG_11yr_idx]
                    if var_name == "phi":
                        # Add non-zero probability of picking 0 and 2pi
                        return np.append(
                            samp_var, np.linspace(0.0, 2 * np.pi, self.nbins)
                        )
                    elif var_name == "theta":
                        # Add non-zero probability of picking 0 and pi
                        return np.append(samp_var, np.linspace(0.0, np.pi, self.nbins))
                    elif var_name == "rn_amp":
                        return np.append(
                            samp_var,
                            np.logspace(
                                min(np.log10(samp_var)),
                                max(np.log10(samp_var)),
                                self.nbins,
                            ),
                        )
                    else:
                        return np.append(
                            samp_var,
                            np.linspace(min(samp_var), max(samp_var), self.nbins),
                        )

    def Get_Sample_Draws(self, var_name, num_draws):
        """For observation times, all noise parameters (``sigma``, ``rn_amp``, ``rn_alpha``), ``cadence``, and sky locations (``phis``, ``thetas``),
        uses the individual parameter value ranges and generates distributions from which to draw new parameters.

        Parameters
        ----------
        var_name: string
            The name of the noise parameter to assign sampled parameters
        num_draws: int,float
            The number of draws to return

        Notes
        -----
        To draw from the generated distributions, one does ``draws = self._distribution.rvs(size=sample_size)``
        """
        var_list = ["T_obs", "phi", "theta", "cadence", "sigma", "rn_amp", "rn_alpha"]

        for i in range(len(var_list)):
            if var_name == var_list[i]:
                NG_11yr_idx = i

        samp_var = self.Get_Param_Distributions(var_name, NG_11yr_idx)
        if isinstance(samp_var, (list, np.ndarray)):
            if len(samp_var) > 1:
                if var_name in ["rn_amp"]:
                    var_hist = np.histogram(
                        samp_var,
                        bins=np.logspace(
                            min(np.log10(samp_var)), max(np.log10(samp_var)), self.nbins
                        ),
                        density=True,
                    )
                else:
                    var_hist = np.histogram(samp_var, bins=self.nbins, density=True)
                var_dist = stats.rv_histogram(var_hist)
                return var_dist.rvs(size=num_draws)
            else:
                return np.ones(num_draws) * samp_var[0]
        else:
            return np.ones(num_draws) * samp_var

    def Init_PTA(self):
        """Initializes a PTA in hasasia

        Notes
        -----
        Assigns pulsar parameters based on what the initial values were given per parameter.
        If necessary parameters are left unassigned, it uses 11yr values for ``n_p`` <= 34, and samples from the 11yr parameters if ``n_p`` > 34
        If a range of values were given for a parameter, the per pulsar parameters are drawn from a uniform distribution
        assigns the new pulsar parameters to the corresponding PTA class parameter.
        See Hazboun, Romano, Smith (2019) <https://arxiv.org/abs/1907.04341> for details

        """
        if not hasattr(self, "_NANOGrav_11yr_params"):
            self.Load_NANOGrav_11yr_Params()

        [
            NG_T_obs,
            NG_phis,
            NG_thetas,
            NG_cadences,
            NG_sigmas,
            NG_rn_amps,
            NG_rn_alphas,
        ] = self._NANOGrav_11yr_params
        var_list = ["T_obs", "phi", "theta", "cadence", "sigma", "rn_amp", "rn_alpha"]

        for i, var in enumerate(var_list):
            if var in self.var_dict.keys():
                if self.var_dict[var]["sampled"]:
                    setattr(self, var, self.Get_Sample_Draws(var, self.n_p))
                else:
                    if self.var_dict["n_p"]["sampled"]:
                        prev_var = getattr(self, var)
                        if isinstance(prev_var, u.Quantity):
                            prev_var = prev_var.value
                        if isinstance(prev_var, (list, np.ndarray)):
                            if len(prev_var) > self.n_p:
                                setattr(self, var, prev_var[: self.n_p])
                            elif len(prev_var) < self.n_p:
                                n_added_p = self.n_p - len(prev_var)
                                var_draw = self.Get_Sample_Draws(var, n_added_p)
                                setattr(
                                    self,
                                    var,
                                    np.append(prev_var, var_draw),
                                )
                            else:
                                pass
                        else:
                            # Constant values for all pulsars
                            setattr(self, var, self.Get_Param_Distributions(var, i))
                    else:
                        # Constant values for all pulsars
                        setattr(self, var, self.Get_Param_Distributions(var, i))
            else:
                # Assign/sample values for values needed to make a sensitivity curve
                if var in ["T_obs", "phi", "theta", "cadence", "sigma"]:
                    # 34 pulsars in the 11yr dataset (ie. len(phis))
                    if self.use_11yr:
                        if self.n_p <= len(self._NANOGrav_11yr_params[i]):
                            setattr(
                                self, var, self._NANOGrav_11yr_params[i][: self.n_p]
                            )
                        else:
                            n_added_p = self.n_p - len(self._NANOGrav_11yr_params[i])
                            var_draw = self.Get_Sample_Draws(var, n_added_p)
                            setattr(
                                self,
                                var,
                                np.append(self._NANOGrav_11yr_params[i], var_draw),
                            )
                    else:
                        setattr(self, var, self.Get_Sample_Draws(var, self.n_p))
                else:
                    if self.use_rn:
                        if self.use_11yr:
                            if self.n_p <= len(self._NANOGrav_11yr_params[i]):
                                setattr(
                                    self, var, self._NANOGrav_11yr_params[i][: self.n_p]
                                )
                            else:
                                n_added_p = self.n_p - len(
                                    self._NANOGrav_11yr_params[i]
                                )
                                var_draw = self.Get_Sample_Draws(var, n_added_p)
                                setattr(
                                    self,
                                    var,
                                    np.append(self._NANOGrav_11yr_params[i], var_draw),
                                )
                        else:
                            setattr(self, var, self.Get_Sample_Draws(var, self.n_p))

        if hasattr(self, "rn_amp"):
            if hasattr(self, "sb_amp"):
                psrs = hassim.sim_pta(
                    timespan=self.T_obs.value,
                    cad=self.cadence.value,
                    sigma=self.sigma.value,
                    phi=self.phi,
                    theta=self.theta,
                    Npsrs=self.n_p,
                    A_rn=self.rn_amp,
                    alpha=self.rn_alpha,
                    A_gwb=self.sb_amp,
                    freqs=self.fT.value,
                )
            else:
                psrs = hassim.sim_pta(
                    timespan=self.T_obs.value,
                    cad=self.cadence.value,
                    sigma=self.sigma.value,
                    phi=self.phi,
                    theta=self.theta,
                    Npsrs=self.n_p,
                    A_rn=self.rn_amp,
                    alpha=self.rn_alpha,
                    freqs=self.fT.value,
                )
        elif hasattr(self, "sb_amp"):
            if not hasattr(self, "sb_alpha"):
                self.sb_alpha = -2 / 3.0
            # Make a set of psrs with the same parameters with a sb as red noise
            psrs = hassim.sim_pta(
                timespan=self.T_obs.value,
                cad=self.cadence.value,
                sigma=self.sigma.value,
                phi=self.phi,
                theta=self.theta,
                Npsrs=self.n_p,
                A_rn=self.sb_amp,
                alpha=self.sb_alpha,
                freqs=self.fT.value,
            )
        else:
            psrs = hassim.sim_pta(
                timespan=self.T_obs.value,
                cad=self.cadence.value,
                sigma=self.sigma.value,
                phi=self.phi,
                theta=self.theta,
                Npsrs=self.n_p,
                freqs=self.fT.value,
            )
        # Turn of sampling for an initialized PTA (except if sampling n_p)
        for key, value in self.var_dict.items():
            if key == "n_p":
                pass
            else:
                value["sampled"] = False
        # Get Spectra of pulsars
        spectra = []
        for p in psrs:
            sp = hassens.Spectrum(p, freqs=self.fT.value)
            spectra.append(sp)

        self._sensitivitycurve = hassens.DeterSensitivityCurve(spectra)


class Interferometer:
    r"""
    Class to make an interferometer

    Parameters
    ----------

    name: string
        name of the instrument
    T_obs: float
        the observation time of the PTA in [years]
    load_location: string, optional
        If you want to load an instrument curve from a file, it's the file path
    I_type: string, {'E','A','h'}
        Sets the type of input data.
        ``'E'`` is the effective strain spectral density :math:`S_{n}(f)` (``'ENSD'``),
        ``'A'`` is the amplitude spectral density, :math:`\sqrt{S_{n}(f)}` (``'ASD'``),
        ``'h'`` is the characteristic strain :math:`h_{n}(f)` (``'h'``)
    f_min: float, optional
        Assigned lowest frequency of instrument (default is assigned in particular child classes)
    f_max: float, optional
        Assigned highest frequency of instrument (default is assigned in particular child classes)
    nfreqs: int, optional
        Number of frequencies in logspace the sensitivity is calculated (default is 1e3)

    """

    def __init__(self, name, T_obs, **kwargs):
        self.name = name
        self.T_obs = T_obs
        for keys, value in kwargs.items():
            if keys == "load_location":
                self.load_location = value
            elif keys == "I_type":
                self.I_type = value
            elif keys == "f_min":
                self.f_min = utils.make_quant(value, "Hz")
            elif keys == "f_max":
                self.f_max = utils.make_quant(value, "Hz")
            elif keys == "nfreqs":
                self.nfreqs = value

        if hasattr(self, "load_location"):
            Load_Data(self)

    @property
    def T_obs(self):
        return self._T_obs

    @T_obs.setter
    def T_obs(self, value):
        self.var_dict = ["T_obs", value]
        if not isinstance(self._return_value, u.Quantity):
            self._return_value = utils.make_quant(self._return_value, "yr")
        self._T_obs = self._return_value

    @property
    def var_dict(self):
        return self._var_dict

    @var_dict.setter
    def var_dict(self, value):
        utils.Get_Var_Dict(self, value)

    @property
    def fT(self):
        if not hasattr(self, "_fT"):
            if hasattr(self, "_I_data"):
                self._fT = self._I_data[:, 0] * u.Hz
            if isinstance(self, SpaceBased):
                self.Set_T_Function_Type()
            if isinstance(self, GroundBased):
                self._fT = (
                    np.logspace(
                        np.log10(self.f_min.value),
                        np.log10(self.f_max.value),
                        self.nfreqs,
                    )
                    * u.Hz
                )
        return self._fT

    @fT.setter
    def fT(self, value):
        self._fT = value

    @fT.deleter
    def fT(self):
        del self._fT

    @property
    def f_opt(self):
        """The optimal frequency of the instrument ie. the frequecy at the lowest strain noise PSD"""
        self._f_opt = self.fT[np.argmin(self.S_n_f)]
        return self._f_opt

    @property
    def P_n_f(self):
        """Strain power sensitivity."""
        raise NotImplementedError(
            "Power Spectral Density method must be defined inside SpaceBased or GroundBased classes."
        )

    @property
    def S_n_f(self):
        """Effective Noise Power Specral Density"""
        if not hasattr(self, "_S_n_f"):
            if hasattr(self, "_I_data"):
                if self._I_Type == "ASD":
                    S_n_f_sqrt = self._I_data[:, 1]
                    self._S_n_f = S_n_f_sqrt**2 / u.Hz
                elif self._I_Type == "ENSD":
                    self._S_n_f = self._I_data[:, 1] / u.Hz
                elif self._I_Type == "h":
                    self._S_n_f = self.h_n_f**2 / self.fT
            else:
                raise NotImplementedError(
                    "Effective Noise Power Spectral Density method must be defined inside SpaceBased or GroundBased classes."
                )
        return self._S_n_f

    @S_n_f.deleter
    def S_n_f(self):
        del self._S_n_f

    @property
    def h_n_f(self):
        """Characteristic Strain/effective strain noise amplitude"""
        if not hasattr(self, "_h_n_f"):
            if hasattr(self, "_I_data") and self._I_Type == "h":
                self._h_n_f = self._I_data[:, 1]
            else:
                self._h_n_f = np.sqrt(self.fT * self.S_n_f)
        return self._h_n_f

    @h_n_f.deleter
    def h_n_f(self):
        del self._h_n_f


class GroundBased(Interferometer):
    """
    Class to make a Ground Based Instrument using the Interferometer base class

    Parameters
    ----------
    noise_dict: dictionary, optional
        A nested noise dictionary that has the main variable parameter name(s) in the top level,
        the next level of the dictionary contains the subparameter variable name(s) and the desired value
        to which the subparameter will be changed. The subparameter value can also be an array/list of the
        [value,min,max] if one wishes to vary the instrument over then min/max range.

    """

    def __init__(self, name, T_obs, **kwargs):
        super().__init__(name, T_obs, **kwargs)

        for keys, value in kwargs.items():
            if keys == "noise_dict":
                if isinstance(value, dict):
                    self.noise_dict = value
                else:
                    raise ValueError(keys + " must be a dictionary of noise sources.")

        if not hasattr(self, "nfreqs"):
            self.nfreqs = int(1e3)
        if not hasattr(self, "f_min"):
            self.f_min = 1.0 * u.Hz
        if not hasattr(self, "f_max"):
            self.f_max = 1e4 * u.Hz

        if not hasattr(self, "load_location"):
            if not hasattr(self, "noise_dict"):
                self.Init_GroundBased()
            else:
                self.Set_Noise_Dict(self.noise_dict)

    @property
    def P_n_f(self):
        """Power Spectral Density."""
        return (self._trace.psd / u.Hz) ** 2

    @property
    def S_n_f(self):
        """Effective Noise Power Spectral Density"""
        if not hasattr(self, "_S_n_f"):
            if hasattr(self, "_I_data"):
                if self._I_Type == "ASD":
                    S_n_f_sqrt = self._I_data[:, 1]
                    self._S_n_f = S_n_f_sqrt**2 / u.Hz
                elif self._I_Type == "ENSD":
                    self._S_n_f = self._I_data[:, 1] / u.Hz
                elif self._I_Type == "h":
                    self._S_n_f = self.h_n_f**2 / self.fT
            else:
                if not any(
                    hasattr(self, attr) for attr in ["_noise_budget", "_base_inst"]
                ):
                    self.Init_GroundBased()
                self._S_n_f = self._trace.psd / u.Hz
                # self._S_n_f = self._noise_budget(self.fT.value,ifo=self._ifo).calc()/u.Hz
        return self._S_n_f

    @S_n_f.deleter
    def S_n_f(self):
        del self._S_n_f

    def Init_GroundBased(self):
        """Initialized the Ground Based detector in gwinc"""

        base_inst = [name for name in self.name.split() if name in gwinc.ifo.IFOS]
        if len(base_inst) == 1:
            self._base_inst = base_inst[0]
        else:
            print(f"You must select a base instrument model from {gwinc.ifo.IFOS}")
            print(
                "Setting base instrument to aLIGO. To change base instrument, include different model in class name and reinitialize."
            )
            self._base_inst = "aLIGO"

        self._noise_budget = gwinc.load_budget(self._base_inst)
        self._ifo = copy.deepcopy(self._noise_budget.ifo)
        self._trace = self._noise_budget.run(freq=self.fT.value, ifo=self._ifo)

    def Set_Noise_Dict(self, noise_dict):
        """Sets new values in the nested dictionary of variable noise values

        Parameters
        ----------
        noise_dict: dictionary
            A nested noise dictionary that has the main variable parameter name(s) in the top level,
            the next level of the dictionary contains the subparameter variable name(s) and the desired value
            to which the subparameter will be changed. The subparameter value can also be an array/list of the
            [value,min,max] if one wishes to vary the instrument over then min/max range.

        Examples
        --------
        ``obj.Set_Noise_Dict({'Infrastructure':{'Length':[3000,1000,5000],'Temp':500},'Laser':{'Wavelength':1e-5,'Power':130}})``

        """
        if not hasattr(self, "_noise_budget"):
            self.Init_GroundBased()
        if isinstance(noise_dict, dict):
            for base_noise, inner_noise_dict in noise_dict.items():
                if base_noise in self._ifo.keys():
                    for sub_noise, sub_noise_val in inner_noise_dict.items():
                        if sub_noise in self._ifo[base_noise].keys():
                            if isinstance(sub_noise_val, dict):
                                for (
                                    sub_sub_noise,
                                    sub_sub_noise_val,
                                ) in sub_noise_val.items():
                                    self.var_dict = [
                                        base_noise
                                        + " "
                                        + sub_noise
                                        + " "
                                        + sub_sub_noise,
                                        sub_sub_noise_val,
                                    ]
                                    setattr(
                                        getattr(self._ifo, base_noise)[sub_noise],
                                        sub_sub_noise,
                                        self._return_value,
                                    )
                                    self._trace = self._noise_budget.run(
                                        freq=self.fT.value, ifo=self._ifo
                                    )
                            else:
                                self.var_dict = [
                                    base_noise + " " + sub_noise,
                                    sub_noise_val,
                                ]
                                setattr(
                                    getattr(self._ifo, base_noise),
                                    sub_noise,
                                    self._return_value,
                                )
                                self._trace = self._noise_budget.run(
                                    freq=self.fT.value, ifo=self._ifo
                                )
                        else:
                            raise ValueError(
                                sub_noise
                                + " is not a subparameter variable noise source.\
                                Try calling Get_Noise_Dict on your GroundBased object to find acceptable variables."
                            )
                else:
                    err_mssg = (
                        base_noise
                        + " is not a valid parameter variable noise source.\n"
                    )
                    err_mssg += "Try calling Get_Noise_Dict on your GroundBased object to find acceptable variables."
                    raise ValueError(err_mssg)
        else:
            raise ValueError("Input must be a dictionary of noise sources.")
        # Overwrite old S_n_f with newly set parameters
        if hasattr(self, "S_n_f"):
            del self.S_n_f

    def Get_Noise_Dict(self):
        """Gets and prints the available variable noises in the detector design"""
        i = 0
        for key_1, item_1 in self._ifo.items():
            print(key_1, "Parameters:")
            if hasattr(item_1, "items"):
                for key_2, item_2 in item_1.items():
                    if isinstance(item_2, np.ndarray):
                        i += 1
                        print("    ", key_2, ": array of shape", item_2.shape)
                    elif isinstance(item_2, list):
                        i += 1
                        print("    ", key_2, ": array of shape", len(item_2))
                    elif isinstance(item_2, (int, float)):
                        i += 1
                        print("    ", key_2, ":", item_2)
                    elif isinstance(item_2, gwinc.struct.Struct):
                        print("    ", key_2, "Subparameters:")
                        for key_3, item_3 in item_2.items():
                            if isinstance(item_3, np.ndarray):
                                i += 1
                                print(
                                    "    ",
                                    "    ",
                                    key_3,
                                    ": array of shape",
                                    item_3.shape,
                                )
                            elif isinstance(item_3, list):
                                i += 1
                                print(
                                    "    ",
                                    "    ",
                                    key_3,
                                    ": array of shape",
                                    len(item_3),
                                )
                            elif isinstance(item_3, (int, float)):
                                i += 1
                                print("    ", "    ", key_3, ":", item_3)
                    else:
                        i += 1
                        print("    ", key_2, ":", item_2)

            print(" ")
            print("Number of Variables: ", i)


class SpaceBased(Interferometer):
    """
    Class to make a Space Based Instrument using the Interferometer base class

    Parameters
    ----------
    L: float
        the armlength the of detector in [meters]
    A_acc: float
        the Amplitude of the Acceleration Noise in [meters/second^2]
    f_acc_break_low: float
        the lower break frequency of the acceleration noise in [Hz]
    f_acc_break_high: float
        the higher break frequency of the acceleration noise in [Hz]
    A_IFO: float
        the amplitude of the interferometer
    T_type: string, {`'N'`,`'A'`}
        Picks the transfer function generation method
        ``'N'`` uses the numerically approximated method in Robson, Cornish, and Liu, 2019
        ``'A'`` uses the analytic fit in Larson, Hiscock, and Hellings, 2000
    Background: Boolean, optional
        Add in a Galactic Binary Confusion Noise
    Background_model: int, {0,1}
        Used to selectGalactic Binary Confusion Noise model: 0 is the Cornish and Robson 2017 model
        with corrections from Schmitz, Kai 2020 (https://www.mdpi.com/2073-8994/12/9/1477),
        and 1 is a model of the Galactic confusion noise parameterized as a function of T_obs
    openingangle: int,float,Quantity, optional
        Determines the opening angle of the instrument, if unassigned use the value in in Robson, Cornish, and Liu, 2019.

    """

    def __init__(self, name, T_obs, *args, **kwargs):
        super().__init__(name, T_obs, **kwargs)

        for keys, value in kwargs.items():
            if keys == "T_type":
                self.T_type = value
            elif keys == "Background":
                self.Background = value
            elif keys == "Background_model":
                self.Background = True
                self.Background_model = value
            elif keys == "openingangle":
                self.openingangle = value

        if not hasattr(self, "nfreqs"):
            self.nfreqs = int(1e3)
        if not hasattr(self, "f_min"):
            self.f_min = 1e-5 * u.Hz
        if not hasattr(self, "f_max"):
            self.f_max = 1.0 * u.Hz
        if hasattr(self, "Background"):
            if not hasattr(self, "Background_model"):
                self.Background_model = 0
        else:
            self.Background = False

        if len(args) != 0:
            [L, A_acc, f_acc_break_low, f_acc_break_high, A_IFO, f_IFO_break] = args
            self.L = L
            self.A_acc = A_acc
            self.f_acc_break_low = f_acc_break_low
            self.f_acc_break_high = f_acc_break_high
            self.A_IFO = A_IFO
            self.f_IFO_break = f_IFO_break

        if not hasattr(self, "load_location"):
            if not hasattr(self, "T_type"):
                self.T_type = "N"
            self.Set_T_Function_Type()

    @property
    def L(self):
        return self._L

    @L.setter
    def L(self, value):
        self.var_dict = ["L", value]
        if not isinstance(self._return_value, u.Quantity):
            self._return_value = utils.make_quant(self._return_value, "m")
        self._L = self._return_value

    @property
    def A_acc(self):
        return self._A_acc

    @A_acc.setter
    def A_acc(self, value):
        self.var_dict = ["A_acc", value]
        if not isinstance(self._return_value, u.Quantity):
            self._return_value = utils.make_quant(self._return_value, "m/s2")
        self._A_acc = self._return_value

    @property
    def f_acc_break_low(self):
        return self._f_acc_break_low

    @f_acc_break_low.setter
    def f_acc_break_low(self, value):
        self.var_dict = ["f_acc_break_low", value]
        if not isinstance(self._return_value, u.Quantity):
            self._return_value = utils.make_quant(self._return_value, "Hz")
        self._f_acc_break_low = self._return_value

    @property
    def f_acc_break_high(self):
        return self._f_acc_break_high

    @f_acc_break_high.setter
    def f_acc_break_high(self, value):
        self.var_dict = ["f_acc_break_high", value]
        if not isinstance(self._return_value, u.Quantity):
            self._return_value = utils.make_quant(self._return_value, "Hz")
        self._f_acc_break_high = self._return_value

    @property
    def A_IFO(self):
        return self._A_IFO

    @A_IFO.setter
    def A_IFO(self, value):
        self.var_dict = ["A_IFO", value]
        if not isinstance(self._return_value, u.Quantity):
            self._return_value = utils.make_quant(self._return_value, "m")
        self._A_IFO = self._return_value

    @property
    def f_IFO_break(self):
        return self._f_IFO_break

    @f_IFO_break.setter
    def f_IFO_break(self, value):
        self.var_dict = ["f_IFO_break", value]
        if not isinstance(self._return_value, u.Quantity):
            self._return_value = utils.make_quant(self._return_value, "Hz")
        self._f_IFO_break = self._return_value

    @property
    def P_n_f(self):
        """Power Spectral Density"""
        if not hasattr(self, "_P_n_f"):
            if not hasattr(self, "_T_Function_Type"):
                self.Set_T_Function_Type()

            P_acc = (
                self.A_acc**2
                * (1 + (self.f_acc_break_low / self.fT) ** 2)
                * (1 + (self.fT / (self.f_acc_break_high)) ** 4)
                / (2 * np.pi * self.fT) ** 4
            )  # Acceleration Noise
            P_IMS = self.A_IFO**2 * (
                1 + (self.f_IFO_break / self.fT) ** 4
            )  # Displacement noise of the interferometric TM--to-TM

            f_trans = const.c / 2 / np.pi / self.L  # Transfer frequency
            self._P_n_f = (
                (P_IMS + 2 * (1 + np.cos(self.fT.value / f_trans.value) ** 2) * P_acc)
                / self.L**2
                / u.Hz
            )
        return self._P_n_f

    @P_n_f.deleter
    def P_n_f(self):
        del self._P_n_f

    @property
    def S_n_f(self):
        """Effective Noise Power Specral Density"""
        if not hasattr(self, "_S_n_f"):
            if hasattr(self, "_I_data"):
                if self._I_Type == "ASD":
                    S_n_f_sqrt = self._I_data[:, 1]
                    self._S_n_f = S_n_f_sqrt**2 / u.Hz
                elif self._I_Type == "ENSD":
                    self._S_n_f = self._I_data[:, 1] / u.Hz
                elif self._I_Type == "h":
                    self._S_n_f = self.h_n_f**2 / self.fT
            else:
                if self.Background:
                    self._S_n_f = (
                        self.P_n_f + self.Add_Background()
                    ) / self.transferfunction**2
                else:
                    self._S_n_f = self.P_n_f / self.transferfunction**2

        return self._S_n_f

    @S_n_f.deleter
    def S_n_f(self):
        del self._S_n_f

    def Load_Transfer_Function(self):
        # Numerical transfer function
        Numerical_Transfer_Function_filedirectory = os.path.join(
            load_directory, "NumericalTransferFunction/transfer.dat"
        )
        Numerical_Transfer_Function_data = np.loadtxt(
            Numerical_Transfer_Function_filedirectory
        )
        self._transferfunctiondata = Numerical_Transfer_Function_data

    def Get_Numeric_Transfer_Function(self):
        if not hasattr(self, "_transferfunctiondata"):
            self.Load_Transfer_Function()

        fc = const.c / (2 * self.L)  # light round trip freq
        LISA_Transfer_Function_f = fc * self._transferfunctiondata[:, 0]

        idx_f_5 = np.abs(LISA_Transfer_Function_f - self.f_min).argmin()
        idx_f_1 = np.abs(LISA_Transfer_Function_f - self.f_max).argmin()

        # 3/10 is normalization 2/5sin(openingangle)
        # Some papers use 3/20, not summing over 2 independent low-freq data channels
        self.transferfunction = (
            np.sqrt(3 / 10) * self._transferfunctiondata[idx_f_5:idx_f_1, 1]
        )
        self.fT = LISA_Transfer_Function_f[idx_f_5:idx_f_1]

    def Get_Analytic_Transfer_Function(self):
        """Response function approximation from Calculation described by Cornish, Robson, Liu 2019
        Uses ``openingangle`` property to determine the opening angle of the instrument, if unassigned use the value in the above paper.
        """
        if hasattr(self, "openingangle"):
            openingangle = self.openingangle
        else:
            openingangle = None
        self.fT = (
            np.logspace(
                np.log10(self.f_min.value), np.log10(self.f_max.value), self.nfreqs
            )
            * u.Hz
        )
        f_L = const.c / 2 / np.pi / self.L  # Transfer frequency
        # 3/10 is normalization 2/5sin(openingangle)
        if isinstance(openingangle, (int, float, u.Quantity)):
            R_f = (2 * np.sin(openingangle) ** 2) / 5 / (1 + 0.6 * (self.fT / f_L) ** 2)
        else:
            R_f = 3 / 10 / (1 + 0.6 * (self.fT / f_L) ** 2)
        self.transferfunction = np.sqrt(R_f)

    def Set_T_Function_Type(self):
        if self.T_type == "n" or self.T_type == "N":
            self._T_type = "numeric"
        elif self.T_type == "a" or self.T_type == "A":
            self._T_type = "analytic"
        else:
            print("\nYou can get the transfer function via 2 methods:")
            print(
                ' *To use the numerically approximated method in Robson, Cornish, and Liu, 2019, input "N".'
            )
            print(
                ' *To use the analytic fit in Larson, Hiscock, and Hellings, 2000, input "A".'
            )
            calc_type = input("Please select the calculation type: ")
            self.Set_T_Function_Type(calc_type)

        if self._T_type == "numeric":
            self.Get_Numeric_Transfer_Function()
        if self._T_type == "analytic":
            self.Get_Analytic_Transfer_Function()

    def Add_Background(self):
        """
        Galactic confusions noise parameters for 6months, 1yr, 2yr, and 4yr
        corresponding to array index 0,1,2,3 respectively
        Uses Background_model to select between Galactic Binary Confusion Noise models: 0 is the Cornish and Robson 2017 model
        with corrections from Schmitz, Kai 2020 (https://www.mdpi.com/2073-8994/12/9/1477),
        and 1 is a model of the Galactic confusion noise parameterized as a function of ``T_obs``
        """
        if self.Background_model == 0:
            """
            Galactic confusions noise parameters for 6months, 1yr, 2yr, and 4yr
            corresponding to array index 0,1,2,3 respectively
            From Cornish and Robson 2017 with corrections from Schmitz, Kai 2020 (https://www.mdpi.com/2073-8994/12/9/1477)
            """
            f = self.fT.to("Hz")
            A = 9e-38 / u.Hz
            a = [0.133, 0.171, 0.165, 0.138]
            b = [x / u.mHz for x in [0.243, 0.292, 0.299, -0.221]]
            k = [x / u.mHz for x in [0.482, 1.020, 0.611, 0.521]]
            g = [x / u.mHz for x in [0.917, 1.680, 1.340, 1.680]]
            f_k = [x * u.mHz for x in [2.58, 2.15, 1.73, 1.13]]
            f_ref = 1000 * u.mHz
            if self.T_obs < 1.0 * u.yr:
                index = 0
            elif self.T_obs >= 1.0 * u.yr and self.T_obs < 2.0 * u.yr:
                index = 1
            elif self.T_obs >= 2.0 * u.yr and self.T_obs < 4.0 * u.yr:
                index = 2
            else:
                index = 3

            gfkf = g[index].to("1/Hz") * (f_k[index].to("Hz") - f)
            kf = k[index].to("1/Hz") * f
            S_c_f = (
                A
                * ((1.0 * u.mHz).to("Hz") / f) ** (7 / 3)
                * np.exp(
                    -((f / f_ref.to("Hz")) ** a[index])
                    - b[index].to("1/Hz") * f * np.sin(kf.to("").value)
                )
                * (1 + np.tanh(gfkf.to("").value))
            )  # White Dwarf Background Noise
        elif self.Background_model == 1:
            """
            Galactic confusions noise parameters as a function of T_obs
            """
            A = 1.4e-44
            agam = 1100
            bgam = 3 / 10
            afk = 0.0016
            bfk = -2 / 9
            f_k = afk * self.T_obs**bfk
            gamma = agam * self.T_obs**bgam

            f = self.fT.value
            S_c_f = (
                A
                * (f ** (-7 / 3))
                * (1 + np.tanh(gamma.value * (f_k.value - f)))
                * (1 / u.Hz)
            )  # White Dwarf Background Noise
        return S_c_f


def Load_Data(detector):
    """
    Function to load in a file to initialize any detector.

    Parameters
    ----------
    detector: object
        Instance of a detector class

    """
    if not hasattr(detector, "I_type"):
        print("Is the data:")
        print(' *Effective Noise Spectral Density - "E"')
        print(' *Amplitude Spectral Density- "A"')
        print(' *Effective Strain - "h"')
        detector.I_type = input("Please enter one of the answers in quotations: ")
        Load_Data(detector)

    if detector.I_type == "E" or detector.I_type == "e":
        detector._I_Type = "ENSD"
    elif detector.I_type == "A" or detector.I_type == "a":
        detector._I_Type = "ASD"
    elif detector.I_type == "h" or detector.I_type == "H":
        detector._I_Type = "h"
    else:
        print("Is the data:")
        print(' *Effective Noise Spectral Density - "E"')
        print(' *Amplitude Spectral Density- "A"')
        print(' *Effective Strain - "h"')
        detector.I_type = input("Please enter one of the answers in quotations: ")
        Load_Data(detector)

    detector._I_data = np.loadtxt(detector.load_location)
    detector.fT = detector._I_data[:, 0] * u.Hz
