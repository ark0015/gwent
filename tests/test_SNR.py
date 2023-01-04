#!/usr/bin/env python
# coding: utf-8

# # Using `gwent` to Calculate Signal-to-Noise Ratios

import astropy.units as u
import matplotlib.pyplot as plt
import pytest

import gwent
from gwent import binary, detector, snr, snrplot

# We need to get the file directories to load in the instrument files.
load_directory = gwent.__path__[0] + "/LoadFiles/InstrumentFiles/"

# Number of SNRMatrix rows
sampleRate_y = 5
# Number of SNRMatrix columns
sampleRate_x = 5


# q = m2/m1 reduced mass
q = 1.0
q_min = 1.0
q_max = 18.0

# Chi = S_i*L/m_i**2, spins of each mass i
chi1 = 0.0  # spin of m1
chi2 = 0.0  # spin of m2
chi_min = -0.85  # Limits of PhenomD for unaligned spins
chi_max = 0.85

z = 0.1  # Redshift
z_min = 1e-2
z_max = 1e2


@pytest.fixture
def source_pta():
    # M = m1+m2 Total Mass
    M = 1e8
    M_min = 1e8
    M_max = 1e11

    source_pta = binary.BBHFrequencyDomain(M, q, z, chi1, chi2)
    source_pta.M = [M, M_min, M_max]
    source_pta.q = [q, q_min, q_max]
    source_pta.chi1 = [chi1, chi_min, chi_max]
    source_pta.chi2 = [chi2, chi_min, chi_max]
    source_pta.z = [z, z_min, z_max]
    source_pta.f_min = 1e-12

    return source_pta


@pytest.fixture
def source_space_based():
    # M = m1+m2 Total Mass
    M = 1e6
    M_min = 1e1
    M_max = 1e10

    source_space_based = binary.BBHFrequencyDomain(M, q, z, chi1=chi1, chi2=chi2)
    source_space_based.M = [M, M_min, M_max]
    source_space_based.q = [q, q_min, q_max]
    source_space_based.chi1 = [chi1, chi_min, chi_max]
    source_space_based.chi2 = [chi2, chi_min, chi_max]
    source_space_based.z = [z, z_min, z_max]

    return source_space_based


@pytest.fixture
def source_ground_based():
    # M = m1+m2 Total Mass
    M = 10.0
    M_min = 1e0
    M_max = 1e5

    source_ground_based = binary.BBHFrequencyDomain(M, q, z, chi1, chi2)
    source_ground_based.M = [M, M_min, M_max]
    source_ground_based.q = [q, q_min, q_max]
    source_ground_based.chi1 = [chi1, chi_min, chi_max]
    source_ground_based.chi2 = [chi2, chi_min, chi_max]
    source_ground_based.z = [z, z_min, z_max]

    return source_ground_based


# aLIGO calculation using pygwinc
T_obs = 4.0 * u.yr  # Observing time in years
T_obs_min = 1.0 * u.yr
T_obs_max = 10.0 * u.yr
noise_dict = {
    "Infrastructure": {"Length": [3995, 2250, 4160], "Temp": [295, 100, 1e3]},
    "Laser": {"Power": [125, 10, 1e4]},
    "Materials": {"Substrate": {"Temp": [295, 10, 300]}},
    "Seismic": {"Gamma": [0.8, 0.1, 1.0], "Rho": [1.8e3, 1500, 2500]},
}


@pytest.fixture
def aLIGO_gwinc():
    aLIGO_gwinc = detector.GroundBased(
        "aLIGO gwinc",
        T_obs,
        noise_dict=noise_dict,
        f_min=1.0,
        f_max=1e4,
        nfreqs=int(1e3),
    )
    aLIGO_gwinc.T_obs = [T_obs, T_obs_min, T_obs_max]
    return aLIGO_gwinc


# NANOGrav calculation using 11.5yr parameters https://arxiv.org/abs/1801.01837
T_obs = 15.0 * u.yr  # Observing time in years
T_obs_min = 5.0 * u.yr
T_obs_max = 30.0 * u.yr

sigma = 100.0 * u.ns.to("s") * u.s  # rms timing residuals in seconds
sigma_min = 100.0 * u.ns.to("s") * u.s
sigma_max = 500.0 * u.ns.to("s") * u.s

N_p = 18  # Number of pulsars
N_p_min = 18
N_p_max = 22

cadence = 1.0 / (
    2 * u.wk.to("yr") * u.yr
)  # Avg observation cadence of 1 every 2 weeks in num/year
cadence_min = 2.0 / u.yr
cadence_max = 1.0 / (u.wk.to("yr") * u.yr)


@pytest.fixture
def NANOGrav_WN():
    NANOGrav_WN = detector.PTA(
        "NANOGrav", N_p, T_obs=T_obs, sigma=sigma, cadence=cadence
    )
    NANOGrav_WN.T_obs = [T_obs, T_obs_min, T_obs_max]
    NANOGrav_WN.sigma = [sigma, sigma_min, sigma_max]
    NANOGrav_WN.n_p = [N_p, N_p_min, N_p_max]
    NANOGrav_WN.cadence = [cadence, cadence_min, cadence_max]
    return NANOGrav_WN


# L3 proposal
# Default Params from https://arxiv.org/abs/1702.00786
T_obs = 4.0 * u.yr  # Observing time in years
T_obs_min = 1.0 * u.yr
T_obs_max = 10.0 * u.yr

L = 2.5e9 * u.m  # armlength in meters
L_min = 1.0e7 * u.m
L_max = 1.0e11 * u.m

A_acc = 3e-15 * u.m / u.s / u.s
A_acc_min = 1e-16 * u.m / u.s / u.s
A_acc_max = 1e-14 * u.m / u.s / u.s

f_acc_break_low = 0.4 * u.mHz.to("Hz") * u.Hz
f_acc_break_low_min = 0.1 * u.mHz.to("Hz") * u.Hz
f_acc_break_low_max = 1.0 * u.mHz.to("Hz") * u.Hz

f_acc_break_high = 8.0 * u.mHz.to("Hz") * u.Hz
f_acc_break_high_min = 1.0 * u.mHz.to("Hz") * u.Hz
f_acc_break_high_max = 10.0 * u.mHz.to("Hz") * u.Hz

f_IFO_break = 2.0 * u.mHz.to("Hz") * u.Hz
f_IFO_break_min = 1.0 * u.mHz.to("Hz") * u.Hz
f_IFO_break_max = 5.0 * u.mHz.to("Hz") * u.Hz

A_IFO = 10e-12 * u.m
A_IFO_min = 1.0e-12 * u.m
A_IFO_max = 2.0e-11 * u.m

Background = False
T_type = "N"


@pytest.fixture
def LISA_ESA():
    LISA_ESA = detector.SpaceBased(
        "LISA_ESA",
        T_obs,
        L,
        A_acc,
        f_acc_break_low,
        f_acc_break_high,
        A_IFO,
        f_IFO_break,
        Background=Background,
        T_type=T_type,
    )
    LISA_ESA.T_obs = [T_obs, T_obs_min, T_obs_max]
    LISA_ESA.L = [L, L_min, L_max]
    LISA_ESA.A_acc = [A_acc, A_acc_min, A_acc_max]
    LISA_ESA.f_acc_break_low = [
        f_acc_break_low,
        f_acc_break_low_min,
        f_acc_break_low_max,
    ]
    LISA_ESA.f_acc_break_high = [
        f_acc_break_high,
        f_acc_break_high_min,
        f_acc_break_high_max,
    ]
    LISA_ESA.A_IFO = [A_IFO, A_IFO_min, A_IFO_max]
    LISA_ESA.f_IFO_break = [f_IFO_break, f_IFO_break_min, f_IFO_break_max]
    return LISA_ESA


# ## SNR Calculation
# ### Global Source Params for all Fiducial Detectors
#
# * 'M' - Mass (Solar Units)
# * 'q' - Mass Ratio
# * 'chi1' - Dimensionless Spin of Black Hole 1
# * 'chi2' - Dimensionless Spin of Black Hole 2
# * 'z' - Redshift

# ### Global Detector Params
# * 'T_obs' - Detector Observation Time

# ## Create of SNR Matrices and Samples for all models

# ###GroundBased ONLY:
#
# * Any single valued variable in list of params given by:
#     instrument_GroundBased.Get_Noise_Dict()
# * To make variable in SNR, declare the main variable, then the subparameter variable
#     as a string e.g. var_x = 'Infrastructure Length', the case matters.


def test_aLIGO_params_MvIL(source_ground_based, aLIGO_gwinc):
    # Variable on x-axis
    var_x = "M"
    # Variable on y-axis
    var_y = "Infrastructure Temp"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_ground_based, aLIGO_gwinc, var_x, sampleRate_x, var_y, sampleRate_y
    )
    plt.show()
    fig, ax = fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        display_cbar=True,
        y_axis_label=False,
        smooth_contours=False,
        logLevels_min=-1.0,
        logLevels_max=5.0,
        y_axis_line=295,
        yticklabels_kwargs={"rotation": 70, "y": 0.02},
        xlabels_kwargs={"labelpad": 0.45},
        display=False,
        return_plt=True,
    )
    plt.close(fig)


def test_aLIGO_params_ILvIT(source_ground_based, aLIGO_gwinc):
    # Variable on x-axis
    var_x = "Infrastructure Length"
    # Variable on y-axis
    var_y = "Infrastructure Temp"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_ground_based, aLIGO_gwinc, var_x, sampleRate_x, var_y, sampleRate_y
    )
    plt.show()
    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        cfill=False,
        display=False,
        return_plt=True,
        x_axis_line=3995,
    )
    plt.close(fig)

    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        smooth_contours=False,
        cfill=True,
        display=False,
        return_plt=True,
        x_axis_label=False,
        y_axis_label=False,
    )
    plt.close(fig)


def test_aLIGO_params_LPvz(source_ground_based, aLIGO_gwinc):
    # Variable on x-axis
    var_x = "Laser Power"
    # Variable on y-axis
    var_y = "z"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_ground_based, aLIGO_gwinc, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        cfill=False,
        display=False,
        return_plt=True,
        x_axis_line=125,
    )
    plt.close(fig)

    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        smooth_contours=False,
        cfill=True,
        display=False,
        return_plt=True,
        x_axis_label=False,
        y_axis_label=False,
    )
    plt.close(fig)


def test_aLIGO_params_SGvMST(source_ground_based, aLIGO_gwinc):
    # Variable on x-axis
    var_x = "Seismic Gamma"
    # Variable on y-axis
    var_y = "Materials Substrate Temp"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_ground_based, aLIGO_gwinc, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        cfill=False,
        display=False,
        return_plt=True,
    )
    plt.close(fig)

    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        smooth_contours=False,
        cfill=True,
        display=False,
        return_plt=True,
        x_axis_label=False,
        y_axis_label=False,
    )
    plt.close(fig)


# ### PTA Only Params
#
# * 'N_p' - Number of Pulsars
# * 'sigma' - Root-Mean-Squared Timing Error
# * 'cadence' - Observation Cadence


def test_NANOGrav_WN_params_Mvq(source_pta, NANOGrav_WN):
    # Variable on x-axis
    var_x = "M"
    # Variable on y-axis
    var_y = "q"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_pta, NANOGrav_WN, var_x, sampleRate_x, var_y, sampleRate_y
    )
    [_, _, _] = snr.Get_SNR_Matrix(
        source_pta, NANOGrav_WN, var_x, sampleRate_x, var_y, sampleRate_y, method="PN"
    )
    fig, ax = snrplot.Plot_SNR(
        var_x, sample_x, var_y, sample_y, SNRMatrix, display=False, return_plt=True
    )
    plt.close(fig)


def test_NANOGrav_WN_params_Mvchi1(source_pta, NANOGrav_WN):
    # Variable on x-axis
    var_x = "M"
    # Variable on y-axis
    var_y = "chi1"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_pta, NANOGrav_WN, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        cfill=False,
        display=False,
        return_plt=True,
    )
    plt.close(fig)

    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        smooth_contours=False,
        cfill=True,
        display=False,
        return_plt=True,
    )
    plt.close(fig)


def test_NANOGrav_WN_params_Mvchii(source_pta, NANOGrav_WN):
    # Variable on x-axis
    var_x = "M"
    # Variable on y-axis
    var_y = "chii"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_pta, NANOGrav_WN, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        cfill=False,
        display=False,
        return_plt=True,
    )
    plt.close(fig)


def test_NANOGrav_WN_params_chiivM(source_pta, NANOGrav_WN):
    # Variable on x-axis
    var_x = "chii"
    # Variable on y-axis
    var_y = "M"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_pta, NANOGrav_WN, var_x, sampleRate_x, var_y, sampleRate_y, method="PN"
    )
    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        cfill=False,
        display=False,
        return_plt=True,
    )
    plt.close(fig)


def test_NANOGrav_WN_params_MvTobs(source_pta, NANOGrav_WN):
    # Variable on x-axis
    var_x = "M"
    # Variable on y-axis
    var_y = "T_obs"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_pta, NANOGrav_WN, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        display=False,
        return_plt=True,
        xticklabels_kwargs={"rotation": 70, "y": 0.02},
        ylabels_kwargs={"labelpad": -5},
    )
    plt.close(fig)


def test_NANOGrav_WN_params_MvNp(source_pta, NANOGrav_WN):
    # Variable on x-axis
    var_x = "M"
    # Variable on y-axis
    var_y = "n_p"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_pta, NANOGrav_WN, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x, sample_x, var_y, sample_y, SNRMatrix, display=False, return_plt=True
    )
    plt.close(fig)


def test_NANOGrav_WN_params_Mvsigma(source_pta, NANOGrav_WN):
    # Variable on x-axis
    var_x = "M"
    # Variable on y-axis
    var_y = "sigma"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_pta, NANOGrav_WN, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x, sample_x, var_y, sample_y, SNRMatrix, display=False, return_plt=True
    )
    plt.close(fig)


def test_NANOGrav_WN_params_Mvcadence(source_pta, NANOGrav_WN):
    source_pta.q = 1.0
    source_pta.chi1 = 0.0
    source_pta.chi2 = 0.0
    source_pta.z = 0.1
    source_pta.f_min = 1e-9
    T_obs = 15.0 * u.yr  # Observing time in years
    T_obs_min = 5.0 * u.yr
    T_obs_max = 30.0 * u.yr
    NANOGrav_WN.T_obs = [T_obs, T_obs_min, T_obs_max]
    NANOGrav_WN.sigma = [sigma, sigma_min, sigma_max]
    NANOGrav_WN.n_p = [N_p, N_p_min, N_p_max]
    NANOGrav_WN.cadence = [cadence, cadence_min, cadence_max]

    # Variable on x-axis
    var_x = "M"
    # Variable on y-axis
    var_y = "cadence"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_pta, NANOGrav_WN, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x, sample_x, var_y, sample_y, SNRMatrix, display=False, return_plt=True
    )
    plt.close(fig)


# ### LISA Only Params
#
# * 'L' - Detector Armlength
# * 'A_acc' - Detector Acceleration Noise
# * 'A_IFO' - Detector Optical Metrology Noise
# * 'f_acc_break_low' - The Low Acceleration Noise Break Frequency
# * 'f_acc_break_high' - The High Acceleration Noise Break Frequency
# * 'f_IFO_break' - The Optical Metrology Noise Break Frequency


def test_LISA_params_LvM(source_space_based, LISA_ESA):
    # Variable on x-axis
    var_x = "L"
    # Variable on y-axis
    var_y = "M"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_space_based, LISA_ESA, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        display=False,
        return_plt=True,
        smooth_contours=False,
    )
    plt.close(fig)


def test_LISA_params_Aaccvz(source_space_based, LISA_ESA):
    # Variable on x-axis
    var_x = "A_acc"
    # Variable on y-axis
    var_y = "z"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_space_based, LISA_ESA, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        display=False,
        return_plt=True,
        dl_axis=True,
    )
    plt.close(fig)

    fig, ax = snrplot.Plot_SNR(
        var_x,
        sample_x,
        var_y,
        sample_y,
        SNRMatrix,
        display=False,
        return_plt=True,
        lb_axis=True,
        smooth_contours=False,
    )
    plt.close(fig)


def test_LISA_params_AIFOvq(source_space_based, LISA_ESA):
    # Variable on x-axis
    var_x = "A_IFO"
    # Variable on y-axis
    var_y = "q"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_space_based, LISA_ESA, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x, sample_x, var_y, sample_y, SNRMatrix, display=False, return_plt=True
    )
    plt.close(fig)


def test_LISA_params_faccbreaklowvchi1(source_space_based, LISA_ESA):
    # Variable on x-axis
    var_x = "f_acc_break_low"
    # Variable on y-axis
    var_y = "chi1"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_space_based, LISA_ESA, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x, sample_x, var_y, sample_y, SNRMatrix, display=False, return_plt=True
    )
    plt.close(fig)


def test_LISA_params_faccbreakhighvfaccbreaklow(source_space_based, LISA_ESA):
    # Variable on x-axis
    var_x = "f_acc_break_high"
    # Variable on y-axis
    var_y = "f_acc_break_low"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_space_based, LISA_ESA, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x, sample_x, var_y, sample_y, SNRMatrix, display=False, return_plt=True
    )
    plt.close(fig)


def test_LISA_params_Tobsvfaccbreakhigh(source_space_based, LISA_ESA):
    # Variable on x-axis
    var_x = "T_obs"
    # Variable on y-axis
    var_y = "f_acc_break_high"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_space_based, LISA_ESA, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x, sample_x, var_y, sample_y, SNRMatrix, display=False, return_plt=True
    )
    plt.close(fig)


def test_LISA_params_fIFObreakvL(source_space_based, LISA_ESA):
    # Variable on x-axis
    var_x = "f_IFO_break"
    # Variable on y-axis
    var_y = "L"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_space_based, LISA_ESA, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x, sample_x, var_y, sample_y, SNRMatrix, display=False, return_plt=True
    )
    plt.close(fig)


def test_LISA_paramsMvfIFObreak(source_space_based, LISA_ESA):
    # Variable on x-axis
    var_x = "q"
    # Variable on y-axis
    var_y = "f_IFO_break"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_space_based, LISA_ESA, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x, sample_x, var_y, sample_y, SNRMatrix, display=False, return_plt=True
    )
    plt.close(fig)


def test_LISA_params_MvAacc(source_space_based, LISA_ESA):
    # Variable on x-axis
    var_x = "z"
    # Variable on y-axis
    var_y = "A_acc"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_space_based,
        LISA_ESA,
        var_x,
        sampleRate_x,
        var_y,
        sampleRate_y,
        inc=0.0,
        integral_consts=4.0,
    )
    fig, ax = snrplot.Plot_SNR(
        var_x, sample_x, var_y, sample_y, SNRMatrix, display=False, return_plt=True
    )
    plt.close(fig)


def test_LISA_params_MvAIFO(source_space_based, LISA_ESA):
    # Variable on x-axis
    var_x = "chi1"
    # Variable on y-axis
    var_y = "A_IFO"
    [sample_x, sample_y, SNRMatrix] = snr.Get_SNR_Matrix(
        source_space_based, LISA_ESA, var_x, sampleRate_x, var_y, sampleRate_y
    )
    fig, ax = snrplot.Plot_SNR(
        var_x, sample_x, var_y, sample_y, SNRMatrix, display=False, return_plt=True
    )
    plt.close(fig)
