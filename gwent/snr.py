import numpy as np
import scipy.interpolate as interp
import scipy.integrate as integrate

import astropy.units as u

from . import detector
from . import binary
from .utils import make_quant


def Get_SNR_Matrix(source, instrument, var_x, sample_rate_x, var_y, sample_rate_y,**kwargs):
    """Calculates SNR Matrix

    Parameters
    ----------
    source : object
        Instance of a gravitational wave source class
    instrument : object
        Instance of a gravitational wave detector class
    var_x : str
        x-axis variable
    sample_rate_x : int
        Number of samples at which SNRMatrix is calculated corresponding to the x-axis variable
    var_y : str
        y-axis variable
    sample_rate_y : array
        samples at which SNRMatrix was calculated corresponding to the y-axis variable

    Returns
    -------
    sample_x : array
        samples at which SNRMatrix was calculated corresponding to the x-axis variable
    sample_y : array
        samples at which SNRMatrix was calculated corresponding to the y-axis variable
    SNRMatrix : array-like
        the sample_rate_y X sample_rate_x matrix at which the SNR was calculated corresponding to the particular x and y-axis variable choices

    Notes
    -----
    Uses the variable given and the data range to sample the space either logrithmically or linearly based on the
    selection of variables. Then it computes the SNR for each value.
    Returns the variable ranges used to calculate the SNR for each matrix, then returns the SNRs with size of the sample_yXsample_x

    """
    for keys,value in kwargs.items():
        if keys == 'inc':
            inc = value
        elif keys == 'integral_consts':
            integral_consts = value
        else:
            raise ValueError("%s is not an accepted input option." % keys)

    if 'inc' not in locals():
        inc = None
    if 'integral_consts' not in locals():
        integral_consts = None


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
                updated_dict_x = {var_x_names[0]: {var_x_names[1]: sample_x[i]}}
                instrument.Set_Noise_Dict(updated_dict_x)
            else:
                setattr(instrument, var_x, sample_x[i])
            Recalculate_Noise(source, instrument)
        elif recalculate_noise in ["neither"]:
            # Update Attribute (also updates dictionary)
            setattr(source, var_x, sample_x[i])

        for j in range(sampleSize_y):
            if recalculate_noise in ["x", "neither"]:
                # Update Attribute (also updates dictionary)
                setattr(source, var_y, sample_y[j])
            elif recalculate_noise in ["both"]:
                # Update Attribute (also updates dictionary)
                if isinstance(instrument, detector.GroundBased):
                    var_y_names = var_y.split()
                    updated_dict_y = {var_y_names[0]: {var_y_names[1]: sample_y[i]}}
                    instrument.Set_Noise_Dict(updated_dict_y)
                else:
                    setattr(instrument, var_y, sample_y[j])
                Recalculate_Noise(source, instrument)

            source.Check_Freq_Evol()
            if source.ismono:  # Monochromatic Source and not diff EOB SNR
                if hasattr(source, "h_gw"):
                    del source.h_gw
                SNRMatrix[j, i] = Calc_Mono_SNR(source, instrument,inc=inc)
            else:  # Chirping Source
                if (
                    recalculate_strain == True
                ):  # If we need to calculate the waveform everytime
                    # Delete old PhenomD waveform
                    if hasattr(source, "_phenomD_f"):
                        del source._phenomD_f
                    if hasattr(source, "_phenomD_h"):
                        del source._phenomD_h
                if hasattr(source, "f"):
                    del source.f
                if hasattr(source, "h_f"):
                    del source.h_f
                SNRMatrix[j, i] = Calc_Chirp_SNR(source, instrument,integral_consts=integral_consts)

    if switch:
        return [original_sample_x, original_sample_y, SNRMatrix.T]
    else:
        return [sample_x, sample_y, SNRMatrix]


def Get_Samples(source, instrument, var_x, sample_rate_x, var_y, sample_rate_y):
    """Gets the x and y-axis samples

    Parameters
    ----------
    source : object
        Instance of a gravitational wave source class
    instrument : object
        Instance of a gravitational wave detector class
    var_x : str
        x-axis variable
    sample_rate_x : int
        Number of samples at which SNRMatrix is calculated corresponding to the x-axis variable
    var_y : str
        y-axis variable
    sample_rate_y : array
        samples at which SNRMatrix was calculated corresponding to the y-axis variable

    Returns
    -------
    sample_x : array
        samples at which SNRMatrix was calculated corresponding to the x-axis variable
    sample_y : array
        samples at which SNRMatrix was calculated corresponding to the y-axis variable

    Notes
    -----
        The function uses that to create a
        sample space for the variable either in linear space or logspace for M,z,L,A_acc
        for everything else.

    """
    sample_x = []
    sample_y = []
    recalculate_strain = False
    recalculate_noise = "neither"

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
        var_x_min != None and var_x_max != None
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
        var_y_min != None and var_y_max != None
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
    source : object
        Instance of a gravitational wave source class
    instrument : object
        Instance of a gravitational wave detector class
    """
    if hasattr(instrument, "I_type") or hasattr(instrument, "load_location"):
        raise ValueError("Cannot vary a loaded instrument's parameters")

    if not isinstance(instrument, detector.GroundBased):
        if hasattr(instrument, "P_n_f"):
            del instrument.P_n_f
    if hasattr(instrument, "S_n_f"):
        del instrument.S_n_f
    if hasattr(instrument, "h_n_f"):
        del instrument.h_n_f
    if hasattr(instrument, "fT"):
        del instrument.fT

    if isinstance(instrument, detector.PTA) and hasattr(
        instrument, "_sensitivitycurve"
    ):
        del instrument._sensitivitycurve
    if hasattr(source, "instrument"):
        source.instrument = instrument


def Calc_Mono_SNR(source, instrument,inc=None):
    """Calculates the SNR for a monochromatic source

    Parameters
    ----------
    source : object
        Instance of a gravitational wave source class
    instrument : object
        Instance of a gravitational wave detector class
    inc : None,float,int, optional
        The inclination of the monochromatic source in radians.

    """
    if not hasattr(source,'instrument'):
        source.instrument = instrument

    if not isinstance(instrument, detector.PTA):
        source.h_gw = binary.Get_Mono_Strain(source, inc=inc)
    else:
        if inc is not None:
            source.h_gw = binary.Get_Mono_Strain(source, inc=inc)
        else:
            source.h_gw = binary.Get_Mono_Strain(source, inc=0.0)
    indxfgw = np.abs(instrument.fT - source.f_gw).argmin()

    return source.h_gw * np.sqrt(
        np.max(instrument.T_obs.to("s")) / instrument.S_n_f[indxfgw]
    )


def Calc_Chirp_SNR(source, instrument,integral_consts=None):
    """Calculates the SNR for an evolving source

    Parameters
    ----------
    source : object
        Instance of a gravitational wave source class
    instrument : object
        Instance of a gravitational wave detector class

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

    if not hasattr(source,'instrument'):
        source.instrument = instrument
    if not hasattr(source,'f_T_obs'):
        source.Check_Freq_Evol()
        
    # Only want to integrate from observed frequency (f(T_obs_before_merger)) till merger
    indxfgw_start = np.abs(source.f - source.f_T_obs).argmin()
    indxfgw_end = len(source.f)
    if indxfgw_start >= len(source.f) - 1:
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
    S_n_f_interp_new = S_n_f_interp_old(np.log10(f_cut.value))
    S_n_f_interp = 10 ** S_n_f_interp_new

    if isinstance(instrument, detector.PTA):
        # Rescaled by 1.5 to make SNR plots match...
        integral_consts = 4.0 * 1.5

    elif isinstance(instrument, detector.SpaceBased) or isinstance(
        instrument, detector.GroundBased
    ):
        integral_consts = 16.0 / 5.0

    # CALCULATE SNR FOR BOTH NOISE CURVES
    denom = S_n_f_interp  # Sky Averaged Noise Spectral Density
    numer = h_cut ** 2
    integrand = numer / denom
    if isinstance(integrand, u.Quantity) and isinstance(f_cut, u.Quantity):
        SNRsqrd = integral_consts * np.trapz(
            integrand.value, f_cut.value, axis=0
        )  # SNR**2
    else:
        SNRsqrd = integral_consts * np.trapz(integrand, f_cut, axis=0)  # SNR**2

    return np.sqrt(SNRsqrd)


def Save_SNR(
    sample_x, sample_y, SNRMatrix, save_location, SNR_filename, sample_filename
):
    """Saves SNR Matrix

    Parameters
    ----------
    sample_x : array
        samples at which SNRMatrix was calculated corresponding to the x-axis variable
    sample_y : array
        samples at which SNRMatrix was calculated corresponding to the y-axis variable
    SNRMatrix : array-like
        the matrix at which the SNR was calculated corresponding to the particular x and y-axis variable choices
    save_location : str
        the directory to which the Samples and SNR are saved
    SNR_filename : str
        the name of the SNR file
    sample_filename : str
        the name of the sample file

    """
    np.savetxt(save_location + SNR_filename, SNRMatrix)
    np.savetxt(save_location + sample_filename, np.transpose([sample_x, sample_y]))
