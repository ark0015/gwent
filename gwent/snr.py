import astropy.units as u
import numpy as np
import scipy.interpolate as interp

from . import binary, detector


def Get_SNR_Matrix(
    source, instrument, var_x, sample_rate_x, var_y, sample_rate_y, **kwargs
):
    """Calculates SNR Matrix

    Parameters
    ----------
    source: object
        Instance of a gravitational wave source class
    instrument: object
        Instance of a gravitational wave detector class
    var_x: str
        x-axis variable
    sample_rate_x: int
        Number of samples at which ``SNRMatrix`` is calculated corresponding to the x-axis variable
    var_y: str
        y-axis variable
    sample_rate_y: array
        samples at which ``SNRMatrix`` was calculated corresponding to the y-axis variable
    inc: int, float, Quantity, Optional
        The inclination of the source in degrees
    integral_consts: int, float, Optional
        Used to adjust the SNR scaling in ``Calc_Chirp_SNR``
    method: str, {'SPA','PN'}
        Switches between methods of calculating the monochromatic strain based on the stationary phase approximation,
        or a rescaling of the source waveform in the low frequency regime (Post-Newtonian approximation)

    Returns
    -------
    sample_x: array
        samples at which SNRMatrix was calculated corresponding to the x-axis variable
    sample_y: array
        samples at which SNRMatrix was calculated corresponding to the y-axis variable
    SNRMatrix: array-like
        the ``sample_rate_y`` X ``sample_rate_x`` matrix at which the SNR was calculated corresponding to the particular x and y-axis variable choices

    Notes
    -----
    Uses the variable given and the data range to sample the space either logrithmically or linearly based on the
    selection of variables. Then it computes the SNR for each value.
    Returns the variable ranges used to calculate the SNR for each matrix, then returns the SNRs with size of the ``sample_y``X``sample_x``
    """
    if "inc" in kwargs.keys():
        inc = kwargs["inc"]
    else:
        inc = None

    if "integral_consts" in kwargs.keys():
        integral_consts = kwargs["integral_consts"]
    else:
        integral_consts = None

    if "method" in kwargs.keys():
        method = kwargs["method"]
    else:
        method = "SPA"

    source.instrument = instrument
    # Get Samples for variables
    [sample_x, sample_y, recalculate_strain, recalculate_noise] = Get_Samples(
        source, instrument, var_x, sample_rate_x, var_y, sample_rate_y
    )

    switch = False
    if recalculate_noise == "y":
        """Calculation runs much faster when it doesn't recalculate
        the noise every time."""
        switch = True
        recalculate_noise = "x"
        original_sample_x = sample_x
        original_sample_y = sample_y
        original_var_x = var_x
        original_var_y = var_y
        var_x = original_var_y
        var_y = original_var_x
        sample_x = original_sample_y
        sample_y = original_sample_x

    sampleSize_x = len(sample_x)
    sampleSize_y = len(sample_y)
    SNRMatrix = np.zeros((sampleSize_y, sampleSize_x))

    for i in range(sampleSize_x):
        if recalculate_noise in ["x", "both"]:
            # Update Attribute (also updates dictionary)
            if isinstance(instrument, detector.GroundBased):
                var_x_names = var_x.split()
                if len(var_x_names) == 2:
                    updated_dict_x = {var_x_names[0]: {var_x_names[1]: sample_x[i]}}
                elif len(var_x_names) == 3:
                    updated_dict_x = {
                        var_x_names[0]: {var_x_names[1]: {var_x_names[2]: sample_x[i]}}
                    }
                instrument.Set_Noise_Dict(updated_dict_x)
            else:
                setattr(instrument, var_x, sample_x[i])
            Recalculate_Noise(source, instrument)
        elif recalculate_noise in ["neither"]:
            if var_x == "chii":
                # Used to change both spins simultaneously
                # Update Attribute (also updates dictionary)
                setattr(source, "chi1", sample_x[i])
                setattr(source, "chi2", sample_x[i])
            else:
                # Update Attribute (also updates dictionary)
                setattr(source, var_x, sample_x[i])

        for j in range(sampleSize_y):
            if recalculate_noise in ["x", "neither"]:
                if var_y == "chii":
                    # Used to change both spins simultaneously
                    # Update Attribute (also updates dictionary)
                    setattr(source, "chi1", sample_y[j])
                    setattr(source, "chi2", sample_y[j])
                else:
                    # Update Attribute (also updates dictionary)
                    setattr(source, var_y, sample_y[j])
            elif recalculate_noise in ["both"]:
                # Update Attribute (also updates dictionary)
                if isinstance(instrument, detector.GroundBased):
                    var_y_names = var_y.split()
                    if len(var_y_names) == 2:
                        updated_dict_y = {var_y_names[0]: {var_y_names[1]: sample_y[i]}}
                    elif len(var_y_names) == 3:
                        updated_dict_y = {
                            var_y_names[0]: {
                                var_y_names[1]: {var_y_names[2]: sample_y[i]}
                            }
                        }
                    instrument.Set_Noise_Dict(updated_dict_y)
                else:
                    setattr(instrument, var_y, sample_y[j])
                Recalculate_Noise(source, instrument)

            binary.Check_Freq_Evol(
                source, T_evol=None, T_evol_frame="observer", f_gw_frame="observer"
            )
            if source.ismono:
                # Monochromatic Source and not diff EOB SNR
                if hasattr(source, "h_gw"):
                    del source.h_gw
                if method == "PN":
                    if recalculate_strain:
                        # If we need to calculate the waveform everytime
                        # Delete old PhenomD waveform
                        if hasattr(source, "_phenomD_f"):
                            del source._phenomD_f
                        if hasattr(source, "_phenomD_h"):
                            del source._phenomD_h
                    if hasattr(source, "f"):
                        del source.f
                    if hasattr(source, "h_f"):
                        del source.h_f
                SNRMatrix[j, i] = Calc_Mono_SNR(
                    source, instrument, inc=inc, method=method
                )
            else:  # Chirping Source
                if recalculate_strain:
                    # If we need to calculate the waveform everytime
                    # Delete old PhenomD waveform
                    if hasattr(source, "_phenomD_f"):
                        del source._phenomD_f
                    if hasattr(source, "_phenomD_h"):
                        del source._phenomD_h
                if hasattr(source, "f"):
                    del source.f
                if hasattr(source, "h_f"):
                    del source.h_f
                SNRMatrix[j, i] = Calc_Chirp_SNR(
                    source, instrument, integral_consts=integral_consts
                )
    if switch:
        return [original_sample_x, original_sample_y, SNRMatrix.T]
    else:
        return [sample_x, sample_y, SNRMatrix]


def Get_Samples(source, instrument, var_x, sample_rate_x, var_y, sample_rate_y):
    """Gets the x and y-axis samples

    Parameters
    ----------
    source: object
        Instance of a gravitational wave source class
    instrument: object
        Instance of a gravitational wave detector class
    var_x: str
        x-axis variable
    sample_rate_x: int
        Number of samples at which ``SNRMatrix`` is calculated corresponding to the x-axis variable
    var_y: str
        y-axis variable
    sample_rate_y: array
        samples at which ``SNRMatrix`` was calculated corresponding to the y-axis variable

    Returns
    -------
    sample_x: array
        samples at which ``SNRMatrix`` was calculated corresponding to the x-axis variable
    sample_y: array
        samples at which ``SNRMatrix`` was calculated corresponding to the y-axis variable

    Notes
    -----
        The function uses that to create a
        sample space for the variable either in linear space or logspace for ``M, z, L, A_acc``
        for everything else.
    """
    sample_x = []
    sample_y = []
    recalculate_strain = False
    recalculate_noise = "neither"

    # Used to change both spins simultaneously, should be arbitary if one uses chi1 or chi2 since they are set to the same value in `Get_SNR_Matrix`
    if var_x == "chii":
        var_x = "chi1"
    elif var_y == "chii":
        var_y = "chi1"

    if var_x in source.var_dict.keys():
        if isinstance(source.var_dict[var_x]["min"], u.Quantity):
            var_x_min = source.var_dict[var_x]["min"].value
            var_x_max = source.var_dict[var_x]["max"].value
        else:
            var_x_min = source.var_dict[var_x]["min"]
            var_x_max = source.var_dict[var_x]["max"]
    elif var_x in instrument.var_dict.keys():
        recalculate_noise = "x"
        if isinstance(instrument.var_dict[var_x]["min"], u.Quantity):
            var_x_min = instrument.var_dict[var_x]["min"].value
            var_x_max = instrument.var_dict[var_x]["max"].value
        else:
            var_x_min = instrument.var_dict[var_x]["min"]
            var_x_max = instrument.var_dict[var_x]["max"]
    else:
        raise ValueError(var_x + " is not a variable in the source or the instrument.")

    if var_y in source.var_dict.keys():
        if isinstance(source.var_dict[var_y]["min"], u.Quantity):
            var_y_min = source.var_dict[var_y]["min"].value
            var_y_max = source.var_dict[var_y]["max"].value
        else:
            var_y_min = source.var_dict[var_y]["min"]
            var_y_max = source.var_dict[var_y]["max"]
    elif var_y in instrument.var_dict.keys():
        if recalculate_noise == "x":
            recalculate_noise = "both"
        else:
            recalculate_noise = "y"
        if isinstance(instrument.var_dict[var_y]["min"], u.Quantity):
            var_y_min = instrument.var_dict[var_y]["min"].value
            var_y_max = instrument.var_dict[var_y]["max"].value
        else:
            var_y_min = instrument.var_dict[var_y]["min"]
            var_y_max = instrument.var_dict[var_y]["max"]
    else:
        raise ValueError(var_y + " is not a variable in the source or the instrument.")

    if var_x in ["q", "chi1", "chi2"] or var_y in ["q", "chi1", "chi2"]:
        recalculate_strain = True  # Must recalculate the waveform at each point

    # order of magnitude cut
    oom_cut = 2.0
    if (
        var_x_min is not None and var_x_max is not None
    ):  # If the variable has non-None 'min',and 'max' dictionary attributes
        if var_x == "n_p":
            instrument.var_dict[var_x]["sampled"] = True
            # sample in integer steps
            sample_range = var_x_max - var_x_min
            if sample_range > 10:
                sample_rate = max(2, int(sample_range / 10))
                sample_x = np.arange(var_x_min, var_x_max, sample_rate)
                if var_x_max not in sample_x:
                    sample_x = np.append(sample_x, var_x_max)
            else:
                sample_x = np.arange(var_x_min, var_x_max + 1)
        else:
            # Any other variables get sorted to linear if max-min < order of magnitude cut (oom_cut)
            # Otherwse the samples are in logspace
            if var_x_max <= 0.0 or var_x_min <= 0.0:
                sample_x = np.linspace(var_x_min, var_x_max, sample_rate_x)
            else:
                scale = np.log10(var_x_max) - np.log10(var_x_min)
                if scale >= oom_cut:
                    sample_x = np.logspace(
                        np.log10(var_x_min), np.log10(var_x_max), sample_rate_x
                    )
                else:
                    sample_x = np.linspace(var_x_min, var_x_max, sample_rate_x)
    else:
        raise ValueError(var_x + " does not have an assigned min and/or max.")

    if (
        var_y_min is not None and var_y_max is not None
    ):  # If the variable has non-None 'min',and 'max' dictionary attributes
        if var_y == "n_p":
            instrument.var_dict[var_y]["sampled"] = True
            # sample in integer steps
            sample_range = var_y_max - var_y_min
            if sample_range > 10:
                sample_rate = max(2, int(sample_range / 10))
                sample_y = np.arange(var_y_min, var_y_max, sample_rate)
                if var_y_max not in sample_y:
                    sample_y = np.append(sample_y, var_y_max)
            else:
                sample_y = np.arange(var_y_min, var_y_max + 1)
        else:
            # Any other variables get sorted to linear if max-min < order of magnitude cut (oom_cut)
            # Otherwse the samples are in logspace
            if var_y_max <= 0.0 or var_y_min <= 0.0:
                sample_y = np.linspace(var_y_min, var_y_max, sample_rate_y)
            else:
                scale = np.log10(var_y_max) - np.log10(var_y_min)
                if scale >= oom_cut:
                    sample_y = np.logspace(
                        np.log10(var_y_min), np.log10(var_y_max), sample_rate_y
                    )
                else:
                    sample_y = np.linspace(var_y_min, var_y_max, sample_rate_y)
    else:
        raise ValueError(var_y + " does not have an assigned min and/or max value.")

    return sample_x, sample_y, recalculate_strain, recalculate_noise


def Recalculate_Noise(source, instrument):
    """Recalculate noise curves if something is varied

    Parameters
    ----------
    source: object
        Instance of a gravitational wave source class
    instrument: object
        Instance of a gravitational wave detector class
    """
    if hasattr(instrument, "I_type") or hasattr(instrument, "load_location"):
        raise ValueError("Cannot vary a loaded instrument's parameters")

    if not isinstance(instrument, detector.GroundBased):
        if hasattr(instrument, "P_n_f"):
            del instrument.P_n_f

    if hasattr(instrument, "fT"):
        del instrument.fT
    if hasattr(instrument, "S_n_f"):
        del instrument.S_n_f
    if hasattr(instrument, "h_n_f"):
        del instrument.h_n_f

    if isinstance(instrument, detector.PTA) and hasattr(
        instrument, "_sensitivitycurve"
    ):
        del instrument._sensitivitycurve
    if hasattr(source, "instrument"):
        source.instrument = instrument


def Calc_Mono_SNR(source, instrument, inc=None, method="SPA"):
    r"""Calculates the SNR for a monochromatic source

    Parameters
    ----------
    source: object
        Instance of a gravitational wave source class
    instrument: object
        Instance of a gravitational wave detector class
    inc: None,float,int, optional
        The inclination of the monochromatic source in radians.
    method: str, {'SPA','PN'}
        Switches between methods of calculating the monochromatic strain based on the stationary phase approximation,
        or a rescaling of the source waveform in the low frequency regime (Post-Newtonian approximation)

    Notes
    -----
    When comparing :math:`h_{0}` from ``Get_Mono_Strain`` (i.e., the typical monochromatic stationary phase approximation with the IMRPhenomD :math:`h_{0}` from :math:`|\tilde{h}(f)|\sqrt{2\dot{f}}`,
    the former is a factor of :math:`\frac{\pi}{2}` larger. We scale this out to make the transition between monochromatic and chirping if ``'SPA'`` is used.
    """
    if not hasattr(source, "instrument"):
        source.instrument = instrument

    # Assumes mass and frequency in source class are in the source frame and observer frame, respectively
    source.h_gw = binary.Get_Mono_Strain(
        source,
        inc=inc,
        freq=source.f_gw,
        f_gw_frame="observer",
        pn_frame="observer",
        out_frame="observer",
        method=method,
    )
    if method == "SPA":
        scale = 2 / np.pi
    else:
        scale = 1.0

    indxfgw = np.abs(instrument.fT - source.f_gw).argmin()
    if indxfgw == 0 or indxfgw >= len(instrument.fT) - 1:
        # The source frequency is assumed to be outside the instrument's frequency, thus the SNR is ~0.
        # print(f"Your assigned source GW frequency is {source.f_gw} and the instrument frequency range is [{np.unique(np.min(instrument.fT))[0]:.1e},{np.unique(np.max(instrument.fT))[0]:.1e}]")
        return 1e-30
    else:
        return (
            scale
            * source.h_gw
            * np.sqrt(
                np.max(np.unique(instrument.T_obs.to("s"))) / instrument.S_n_f[indxfgw]
            )
        )


def Calc_Chirp_SNR(source, instrument, integral_consts=None):
    """Calculates the SNR for an evolving source

    Parameters
    ----------
    source: object
        Instance of a gravitational wave source class
    instrument: object
        Instance of a gravitational wave detector class
    integral_consts: int, float, Optional
        Used to adjust the SNR scaling

    Notes
    -----
    Uses an interpolated method to align waveform and instrument noise, then integrates
    over the overlapping region. See eqn 18 from Robson,Cornish,and Liu 2018 <https://arxiv.org/abs/1803.01944>
    Values outside of the sensitivity curve are arbitrarily set to 1e30 so the SNR is effectively 0

    """

    # Previously, it was designed to integrate from initial observed frequency f(t_init) to f(t_init-T_obs)
    # Does not work unless t_init is randomly sampled, which we don't do
    # indxfgw_start = np.abs(source.f-source.f_init).argmin()
    # indxfgw_end = np.abs(source.f-source.f_T_obs).argmin()

    if not hasattr(source, "instrument"):
        source.instrument = instrument
    if not hasattr(source, "f_T_obs"):
        binary.Check_Freq_Evol(source)

    # Only want to integrate from observed frequency (f(T_obs_before_merger)) till merger
    indxfgw_start = np.abs(source.f - source.f_T_obs).argmin()
    if indxfgw_start == 0:
        statement_1 = "Uh, you probably should set your source f_min to lower. "
        statement_1 += f"Your minimum calculated frequency is {source.f[0]} and f(T_obs) is {source.f_T_obs}"
        print(statement_1)
    indxfgw_end = len(source.f)
    if indxfgw_start >= indxfgw_end - 1:
        # If the SMBH has already merged set the SNR to ~0
        return 1e-30
    else:
        f_cut = source.f[indxfgw_start:indxfgw_end]
        h_cut = source.h_f[indxfgw_start:indxfgw_end]

    #################################
    # Interpolate the Strain Noise Spectral Density to only the frequencies the
    # strain runs over
    # Set Noise to 1e30 outside of signal frequencies
    S_n_f_interp_old = interp.interp1d(
        np.log10(instrument.fT.value),
        np.log10(instrument.S_n_f.value),
        kind="cubic",
        fill_value=30.0,
        bounds_error=False,
    )
    S_n_f_interp = 10 ** S_n_f_interp_old(np.log10(f_cut.value))

    if not isinstance(integral_consts, (int, float)):
        integral_consts = 16.0 / 5.0

    # CALCULATE SNR FOR BOTH NOISE CURVES
    denom = S_n_f_interp  # Sky Averaged Noise Spectral Density
    numer = h_cut**2
    integrand = numer / denom
    if isinstance(integrand, u.Quantity) and isinstance(f_cut, u.Quantity):
        SNRsqrd = integral_consts * np.trapz(
            integrand.value, f_cut.value, axis=0
        )  # SNR**2
    else:
        SNRsqrd = integral_consts * np.trapz(integrand, f_cut, axis=0)  # SNR**2

    return np.sqrt(SNRsqrd)
