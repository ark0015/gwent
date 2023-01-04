#!/usr/bin/env python
# coding: utf-8

import astropy.units as u
import pytest

import gwent
from gwent import binary, detector, snr

# We need to get the file directories to load in the instrument files.
load_directory = gwent.__path__[0] + "/LoadFiles/InstrumentFiles/"

# Set axes variables and sizes

# Variable on y-axis
var_y = "z"
# Number of SNRMatrix rows
sampleRate_y = 75
# Variable on x-axis
var_x = "M"
# Number of SNRMatrix columns
sampleRate_x = 75


# q = m2/m1 reduced mass
q = 1.0
q_min = 1.0
q_max = 18.0

# Chi = S_i*L/m_i**2, spins of each mass i
chi1 = 0.0  # spin of m1
chi2 = 0.0  # spin of m2
chi_min = -0.85  # Limits of PhenomD for unaligned spins
chi_max = 0.85

z = 1.0  # Redshift
z_min = 1e-2
z_max = 1e3


@pytest.fixture
def source_ground_based():
    # M = m1+m2 Total Mass
    M = 1e8
    M_min = 1e7
    M_max = 1e11

    source_ground_based = binary.BBHFrequencyDomain(M, q, z, chi1, chi2)
    source_ground_based.M = [M, M_min, M_max]
    source_ground_based.q = [q, q_min, q_max]
    source_ground_based.chi1 = [chi1, chi_min, chi_max]
    source_ground_based.chi2 = [chi2, chi_min, chi_max]
    source_ground_based.z = [z, z_min, z_max]

    return source_ground_based


@pytest.fixture
def source_pta():
    # M = m1+m2 Total Mass
    M = 1e8
    M_min = 1e7
    M_max = 1e11

    source_pta = binary.BBHFrequencyDomain(M, q, z, chi1, chi2)
    source_pta.M = [M, M_min, M_max]
    source_pta.q = [q, q_min, q_max]
    source_pta.chi1 = [chi1, chi_min, chi_max]
    source_pta.chi2 = [chi2, chi_min, chi_max]
    source_pta.z = [z, z_min, z_max]

    return source_pta


@pytest.fixture
def source_space_based():
    # M = m1+m2 Total Mass
    M = 1e6
    M_min = 1e1
    M_max = 1e10

    source_space_based = binary.BBHFrequencyDomain(M, q, z, chi1, chi2)
    source_space_based.M = [M, M_min, M_max]
    source_space_based.q = [q, q_min, q_max]
    source_space_based.chi1 = [chi1, chi_min, chi_max]
    source_space_based.chi2 = [chi2, chi_min, chi_max]
    source_space_based.z = [z, z_min, z_max]

    return source_space_based


# #### LISA Proposal 1
#
# SNR values from the ESA L3 proposal run.

# L3 proposal
# Default Params from https://arxiv.org/abs/1702.00786


def test_LISA_prop1(source_space_based):
    T_obs = 4 * u.yr  # Observing time in years
    L = 2.5e9 * u.m  # armlength in meters
    A_acc = 3e-15 * u.m / u.s / u.s
    f_acc_break_low = 0.4 * u.mHz.to("Hz") * u.Hz
    f_acc_break_high = 8.0 * u.mHz.to("Hz") * u.Hz
    f_IMS_break = 2.0 * u.mHz.to("Hz") * u.Hz
    A_IMS = 10e-12 * u.m
    Background = False
    T_type = "N"
    LISA_prop1 = detector.SpaceBased(
        "LISA_ESA",
        T_obs,
        L,
        A_acc,
        f_acc_break_low,
        f_acc_break_high,
        A_IMS,
        f_IMS_break,
        Background=Background,
        T_type=T_type,
    )
    [lisa_sample_x, lisa_sample_y, lisa_SNR] = snr.Get_SNR_Matrix(
        source_space_based, LISA_prop1, var_x, sampleRate_x, var_y, sampleRate_y
    )


# #### Einstein Telescope
# SNR values from the Einstein Telescope proposal run.
# Einstein Telescope
def test_ET(source_ground_based):
    ET_filedirectory = load_directory + "/EinsteinTelescope/"
    ET_filename = "ET_D_data.txt"
    ET_filelocation = ET_filedirectory + ET_filename
    T_obs = 4 * u.yr  # Observing time in years
    ET = detector.GroundBased("ET", T_obs, load_location=ET_filelocation, I_type="A")
    [et_sample_x, et_sample_y, et_SNR] = snr.Get_SNR_Matrix(
        source_ground_based, ET, var_x, sampleRate_x, var_y, sampleRate_y
    )


# #### aLIGO
# SNR values from the Advanced LIGO run.
# aLIGO
def test_aLIGO(source_ground_based):
    aLIGO_filedirectory = load_directory + "/aLIGO/"
    aLIGO_filename = "aLIGODesign.txt"
    aLIGO_filelocation = aLIGO_filedirectory + aLIGO_filename
    T_obs = 4 * u.yr  # Observing time in years
    aLIGO = detector.GroundBased(
        "aLIGO", T_obs, load_location=aLIGO_filelocation, I_type="A"
    )
    [aLIGO_sample_x, aLIGO_sample_y, aLIGO_SNR] = snr.Get_SNR_Matrix(
        source_ground_based, aLIGO, var_x, sampleRate_x, var_y, sampleRate_y
    )


# #### NANOGrav
# SNR values from the NANOGrav-esque run.
# NANOGrav calculation using 11.5yr parameters https://arxiv.org/abs/1801.01837
# rms timing residuals in nanoseconds to seconds
sigma_nano = 100 * u.ns.to("s") * u.s
T_nano = 15 * u.yr  # Observing time in years
N_p_nano = 18  # Number of pulsars
# Avg observation cadence of 1 every 2 weeks in number/year
cadence_nano = 1 / (2 * u.wk.to("yr") * u.yr)


def test_NANOGrav_11yr(source_pta):
    load_name = "NANOGrav_11yr_S_eff.txt"
    load_location = load_directory + "/NANOGrav/StrainFiles/" + load_name
    T_obs = 11.42 * u.yr  # Observing time in years
    nanograv = detector.PTA(
        "NANOGrav 11yr", T_obs=T_obs, load_location=load_location, I_type="E"
    )
    [nanograv_sample_x, nanograv_sample_y, nanograv_SNR] = snr.Get_SNR_Matrix(
        source_pta, nanograv, var_x, sampleRate_x, var_y, sampleRate_y
    )


# NANOGrav with White Noise only
def test_pta_NANOGrav_WN(source_pta):
    NANOGrav_WN = detector.PTA(
        "NANOGrav, WN Only",
        N_p_nano,
        T_obs=T_nano,
        sigma=sigma_nano,
        cadence=cadence_nano,
    )
    [NANOGrav_WN_sample_x, NANOGrav_WN_sample_y, NANOGrav_WN_SNR] = snr.Get_SNR_Matrix(
        source_pta, NANOGrav_WN, var_x, sampleRate_x, var_y, sampleRate_y
    )


# NANOGrav with White and Varied Red Noise
def test_pta_NANOGrav_WN_GWB(source_pta):
    NANOGrav_WN_GWB = detector.PTA(
        "NANOGrav, WN and GWB",
        N_p_nano,
        T_obs=T_nano,
        sigma=sigma_nano,
        cadence=cadence_nano,
        sb_amp=4e-16,
    )
    [
        NANOGrav_WN_GWB_sample_x,
        NANOGrav_WN_GWB_sample_y,
        NANOGrav_WN_GWB_SNR,
    ] = snr.Get_SNR_Matrix(
        source_pta, NANOGrav_WN_GWB, var_x, sampleRate_x, var_y, sampleRate_y
    )


# #### SKA
# SNR values from the SKA-esque run.
# SKA calculation using parameters and methods from https://arxiv.org/abs/0804.4476 section 7.1
def test_SKA(source_pta):
    T_obs = 15 * u.yr  # Observing time (years)
    sigma = 10 * u.ns.to("s") * u.s  # rms timing residuals in nanoseconds
    N_p = 20  # Number of pulsars
    cadence = 1 / (
        u.wk.to("yr") * u.yr
    )  # Avg observation cadence of 1 every week in num/year

    SKA = detector.PTA("SKA", N_p, T_obs=T_obs, sigma=sigma, cadence=cadence)
    [SKA_sample_x, SKA_sample_y, SKA_SNR] = snr.Get_SNR_Matrix(
        source_pta, SKA, var_x, sampleRate_x, var_y, sampleRate_y
    )
