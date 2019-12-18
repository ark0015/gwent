#!/usr/bin/env python
# coding: utf-8

# # Using `gwent` to Calculate Signal-to-Noise Ratios

# Here we present a tutorial on how to use `gwent` to calculate SNRs for the instrument models currently implemented (LISA, PTAs, aLIGO, and Einstein Telescope) with the signal being an array of coalescing Binary Black Holes.

import numpy as np
import astropy.units as u

import pytest

import os,sys
current_path = os.getcwd()
splt_path = current_path.split("/")
top_path_idx = splt_path.index('Research')
top_directory = "/".join(splt_path[0:top_path_idx+1])
gwent_path = top_directory + '/gwent/'
sys.path.insert(0,gwent_path)

import gwent
from gwent import binary
from gwent import detector
from gwent import snr
from gwent import snrplot

# We need to get the file directories to load in the instrument files.
load_directory = gwent.__path__[0] + '/LoadFiles/InstrumentFiles/'

#Number of SNRMatrix rows
sampleRate_y = 10
#Number of SNRMatrix columns
sampleRate_x = 10


# ## Source Selection Function
# Takes in a an instrument model that dictates reasonable mass ranges for the particular detector mass regime and instantiates a source with the variable ranges limited by the waveform calibration region.
# The source parameters must be set (ie. M,q,z,chi1,chi2), but one only needs to set the minima and maxima of the selected SNR axes variables.

#q = m2/m1 reduced mass
q = 1.0
q_min = 1.0
q_max = 18.0

#Chi = S_i*L/m_i**2, spins of each mass i
chi1 = 0.0 #spin of m1
chi2 = 0.0 #spin of m2
chi_min = -0.85 #Limits of PhenomD for unaligned spins
chi_max = 0.85

z = 1.0 #Redshift
z_min = 1e-2
z_max = 1e3

@pytest.fixture
def source_pta():
    #M = m1+m2 Total Mass
    M = 1e8
    M_min = 1e7
    M_max = 1e11
    
    source_pta = binary.BBHFrequencyDomain(M,q,z,chi1,chi2)
    source_pta.M = [M,M_min,M_max]
    source_pta.q = [q,q_min,q_max]
    source_pta.chi1 = [chi1,chi_min,chi_max]
    source_pta.chi2 = [chi2,chi_min,chi_max]
    source_pta.z = [z,z_min,z_max]

    return source_pta

@pytest.fixture
def source_space_based():
    #M = m1+m2 Total Mass
    M = 1e6
    M_min = 1e1
    M_max = 1e10
    
    source_space_based = binary.BBHFrequencyDomain(M,q,z,chi1,chi2)
    source_space_based.M = [M,M_min,M_max]
    source_space_based.q = [q,q_min,q_max]
    source_space_based.chi1 = [chi1,chi_min,chi_max]
    source_space_based.chi2 = [chi2,chi_min,chi_max]
    source_space_based.z = [z,z_min,z_max]

    return source_space_based

@pytest.fixture
def source_ground_based():
    #M = m1+m2 Total Mass
    M = 10
    M_min = 1e0
    M_max = 1e5
    
    source_ground_based = binary.BBHFrequencyDomain(M,q,z,chi1,chi2)
    source_ground_based.M = [M,M_min,M_max]
    source_ground_based.q = [q,q_min,q_max]
    source_ground_based.chi1 = [chi1,chi_min,chi_max]
    source_ground_based.chi2 = [chi2,chi_min,chi_max]
    source_ground_based.z = [z,z_min,z_max]

    return source_ground_based

#aLIGO calculation using pygwinc
T_obs = 4*u.yr #Observing time in years
T_obs_min = 1*u.yr
T_obs_max = 10*u.yr
noise_dict = {'Infrastructure':
                {'Length':[10000,20000,50000],'Temp':[290,200,400]}}
@pytest.fixture
def aLIGO_gwinc():
    aLIGO_gwinc = detector.GroundBased('aLIGO gwinc',T_obs,noise_dict=noise_dict,
        f_low=1.,f_high=1e4,nfreqs=1e3)        
    aLIGO_gwinc.T_obs = [T_obs,T_obs_min,T_obs_max]
    return aLIGO_gwinc

#NANOGrav calculation using 11.5yr parameters https://arxiv.org/abs/1801.01837
T_obs = 15*u.yr #Observing time in years
T_obs_min = 5*u.yr
T_obs_max = 30*u.yr

sigma = 100*u.ns.to('s')*u.s #rms timing residuals in seconds
sigma_min = 100*u.ns.to('s')*u.s
sigma_max = 500*u.ns.to('s')*u.s

N_p = 18 #Number of pulsars
N_p_min = 18
N_p_max = 22

cadence = 1/(2*u.wk.to('yr')*u.yr) #Avg observation cadence of 1 every 2 weeks in num/year
cadence_min = 2/u.yr
cadence_max = 1/(u.wk.to('yr')*u.yr)

@pytest.fixture
def NANOGrav_WN():
    NANOGrav_WN = detector.PTA('NANOGrav',T_obs,N_p,sigma,cadence)
    NANOGrav_WN.T_obs = [T_obs,T_obs_min,T_obs_max]
    NANOGrav_WN.sigma = [sigma,sigma_min,sigma_max]
    NANOGrav_WN.N_p = [N_p,N_p_min,N_p_max]
    NANOGrav_WN.cadence = [cadence,cadence_min,cadence_max]
    return NANOGrav_WN

@pytest.fixture
def NANOGrav_WN_RN():
    NANOGrav_WN_RN = detector.PTA('NANOGrav, WN and RN',T_obs,N_p,sigma,cadence,
        rn_amp=[1e-16,1e-12],rn_alpha=[-1/2,1.25])
    NANOGrav_WN_RN.T_obs = [T_obs,T_obs_min,T_obs_max]
    NANOGrav_WN_RN.sigma = [sigma,sigma_min,sigma_max]
    NANOGrav_WN_RN.N_p = [N_p,N_p_min,N_p_max]
    NANOGrav_WN_RN.cadence = [cadence,cadence_min,cadence_max]
    return NANOGrav_WN_RN

@pytest.fixture
def NANOGrav_WN_GWB():
    NANOGrav_WN_RN = detector.PTA('NANOGrav, WN and RN',T_obs,N_p,sigma,cadence,
        GWB_amp=4e-16,GWB_alpha=-2/3)
    NANOGrav_WN_RN.T_obs = [T_obs,T_obs_min,T_obs_max]
    NANOGrav_WN_RN.sigma = [sigma,sigma_min,sigma_max]
    NANOGrav_WN_RN.N_p = [N_p,N_p_min,N_p_max]
    NANOGrav_WN_RN.cadence = [cadence,cadence_min,cadence_max]
    return NANOGrav_WN_RN

        
sigma = 10*u.ns.to('s')*u.s #rms timing residuals in nanoseconds
sigma_min = 10*u.ns.to('s')*u.s
sigma_max = 100*u.ns.to('s')*u.s
N_p = 20 #Number of pulsars
cadence = 1/(u.wk.to('yr')*u.yr) #Avg observation cadence of 1 every week in num/year
        
@pytest.fixture
def SKA():
    SKA = detector.PTA('SKA',T_obs,N_p,sigma,cadence)
    SKA.T_obs = [T_obs,T_obs_min,T_obs_max]
    SKA.sigma = [sigma,sigma_min,sigma_max]
    SKA.N_p = [N_p,N_p_min,N_p_max]
    SKA.cadence = [cadence,cadence_min,cadence_max]
    return SKA
        

#L3 proposal
#Default Params from https://arxiv.org/abs/1702.00786
T_obs = 4*u.yr #Observing time in years
T_obs_min = 1*u.yr
T_obs_max = 10*u.yr

L = 2.5e9*u.m #armlength in meters
L_min = 1.0e7*u.m
L_max = 1.0e11*u.m

A_acc = 3e-15*u.m/u.s/u.s
A_acc_min = 1e-16*u.m/u.s/u.s
A_acc_max = 1e-14*u.m/u.s/u.s

f_acc_break_low = .4*u.mHz.to('Hz')*u.Hz
f_acc_break_low_min = .1*u.mHz.to('Hz')*u.Hz
f_acc_break_low_max = 1.0*u.mHz.to('Hz')*u.Hz

f_acc_break_high = 8.*u.mHz.to('Hz')*u.Hz
f_acc_break_high_min = 1.*u.mHz.to('Hz')*u.Hz
f_acc_break_high_max = 10.*u.mHz.to('Hz')*u.Hz

f_IFO_break = 2.*u.mHz.to('Hz')*u.Hz
f_IFO_break_min = 1.*u.mHz.to('Hz')*u.Hz
f_IFO_break_max = 5.*u.mHz.to('Hz')*u.Hz

A_IFO = 10e-12*u.m
A_IFO_min = 1.0e-12*u.m
A_IFO_max = 2.0e-11*u.m

Background = False
T_type = 'N'
    
@pytest.fixture
def LISA_ESA():
    LISA_ESA = detector.SpaceBased('LISA_ESA',T_obs,L,A_acc,f_acc_break_low,f_acc_break_high,A_IFO,f_IFO_break,Background=Background,T_type=T_type)
    LISA_ESA.T_obs = [T_obs,T_obs_min,T_obs_max]
    LISA_ESA.L = [L,L_min,L_max]
    LISA_ESA.A_acc = [A_acc,A_acc_min,A_acc_max]
    LISA_ESA.f_acc_break_low = [f_acc_break_low,f_acc_break_low_min,f_acc_break_low_max]
    LISA_ESA.f_acc_break_high = [f_acc_break_high,f_acc_break_high_min,f_acc_break_high_max]
    LISA_ESA.A_IFO = [A_IFO,A_IFO_min,A_IFO_max]
    LISA_ESA.f_IFO_break = [f_IFO_break,f_IFO_break_min,f_IFO_break_max]
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

#Variable on x-axis
var_x = 'M'
# ## Create of SNR Matrices and Samples for all models

# ###GroundBased ONLY:
#
# * Any single valued variable in list of params given by:
#     instrument_GroundBased.Get_Noise_Dict()
# * To make variable in SNR, declare the main variable, then the subparameter variable
#     as a string e.g. var_x = 'Infrastructure Length', the case matters.
#Variable on y-axis
var_xs_ground = ['M','Infrastructure Length']
var_ys_ground = ['Infrastructure Length','Infrastructure Temp']
def test_aLIGO_params(source_ground_based,aLIGO_gwinc):
    reset_source_ground_based = source_ground_based
    reset_aLIGO_gwinc = aLIGO_gwinc
    for var_x in var_xs_ground:
        if var_x in source_ground_based.var_dict.keys():
            reset_source_ground_based = source_ground_based
        elif var_x in aLIGO_gwinc.var_dict.keys():
            reset_aLIGO_gwinc = aLIGO_gwinc

        for var_y in var_ys_ground:
            if var_y in source_ground_based.var_dict.keys():
                reset_source_ground_based = source_ground_based
            elif var_y in aLIGO_gwinc.var_dict.keys():
                reset_aLIGO_gwinc = aLIGO_gwinc

            if var_x != var_y:
                [sample_x,sample_y,SNRMatrix] = snr.Get_SNR_Matrix(reset_source_ground_based,reset_aLIGO_gwinc,
                                                                   var_x,sampleRate_x,
                                                                   var_y,sampleRate_y)
                snrplot.Plot_SNR(source_ground_based,aLIGO_gwinc,var_x,sample_x,var_y,sample_y,
                    SNRMatrix,display=False)
# ### PTA Only Params
# 
# * 'N_p' - Number of Pulsars
# * 'sigma' - Root-Mean-Squared Timing Error
# * 'cadence' - Observation Cadence
#Variable on y-axis
var_ys_ptas = ['q','chi1','chi2','T_obs','N_p','sigma','cadence']
def test_NANOGrav_WN_params(source_pta,NANOGrav_WN):
    reset_source_pta = source_pta
    reset_NANOGrav_WN = NANOGrav_WN
    for var_y in var_ys_ptas:
        if var_y in source_pta.var_dict.keys():
            reset_source_pta = source_pta
        elif var_y in NANOGrav_WN.var_dict.keys():
            reset_NANOGrav_WN = NANOGrav_WN

        [sample_x,sample_y,SNRMatrix] = snr.Get_SNR_Matrix(reset_source_pta,reset_NANOGrav_WN,
                                                           var_x,sampleRate_x,
                                                           var_y,sampleRate_y)
        snrplot.Plot_SNR(source_pta,NANOGrav_WN,var_x,sample_x,var_y,sample_y,
            SNRMatrix,display=False)
        
def test_NANOGrav_WN_RN_params(source_pta,NANOGrav_WN_RN):
    reset_source_pta = source_pta
    reset_NANOGrav_WN_RN = NANOGrav_WN_RN
    for var_y in var_ys_ptas:
        if var_y in source_pta.var_dict.keys():
            reset_source_pta = source_pta
        elif var_y in NANOGrav_WN_RN.var_dict.keys():
            reset_NANOGrav_WN_RN = NANOGrav_WN_RN

        [sample_x,sample_y,SNRMatrix] = snr.Get_SNR_Matrix(reset_source_pta,reset_NANOGrav_WN_RN,
                                                           var_x,sampleRate_x,
                                                           var_y,sampleRate_y)
def test_NANOGrav_WN_GWB_params(source_pta,NANOGrav_WN_GWB):
    reset_source_pta = source_pta
    reset_NANOGrav_WN_GWB = NANOGrav_WN_GWB
    for var_y in var_ys_ptas:
        if var_y in source_pta.var_dict.keys():
            reset_source_pta = source_pta
        elif var_y in NANOGrav_WN_GWB.var_dict.keys():
            reset_NANOGrav_WN_GWB = NANOGrav_WN_GWB

        [sample_x,sample_y,SNRMatrix] = snr.Get_SNR_Matrix(reset_source_pta,reset_NANOGrav_WN_GWB,
                                                           var_x,sampleRate_x,
                                                           var_y,sampleRate_y)
def test_SKA_params(source_pta,SKA):
    reset_source_pta = source_pta
    reset_SKA = SKA
    for var_y in var_ys_ptas:
        if var_y in source_pta.var_dict.keys():
            reset_source_pta = source_pta
        elif var_y in SKA.var_dict.keys():
            reset_SKA = SKA

        [sample_x,sample_y,SNRMatrix] = snr.Get_SNR_Matrix(reset_source_pta,reset_SKA,
                                                           var_x,sampleRate_x,
                                                           var_y,sampleRate_y)


# ### LISA Only Params
# 
# * 'L' - Detector Armlength
# * 'A_acc' - Detector Acceleration Noise
# * 'A_IFO' - Detector Optical Metrology Noise
# * 'f_acc_break_low' - The Low Acceleration Noise Break Frequency
# * 'f_acc_break_high' - The High Acceleration Noise Break Frequency
# * 'f_IFO_break' - The Optical Metrology Noise Break Frequency
var_y = 'z'
var_xs_LISA = ['M','z','q','chi1','T_obs','L','A_acc','A_IFO','f_acc_break_low','f_acc_break_high','f_IFO_break']
var_ys_LISA = ['M','z','q','chi1','T_obs','L','A_acc','A_IFO','f_acc_break_low','f_acc_break_high','f_IFO_break']

def test_LISA_params(source_space_based,LISA_ESA):
    reset_LISA_ESA = LISA_ESA
    reset_source_space_based = source_space_based
    for var_x in var_xs_LISA:
        if var_x in source_space_based.var_dict.keys():
            reset_source_space_based = source_space_based
        elif var_x in LISA_ESA.var_dict.keys():
            reset_LISA_ESA = LISA_ESA

        for var_y in var_ys_LISA:
            if var_y in source_space_based.var_dict.keys():
                reset_source_space_based = source_space_based
            elif var_y in LISA_ESA.var_dict.keys():
                reset_LISA_ESA = LISA_ESA

            if var_x != var_y:
                [sample_x,sample_y,SNRMatrix] = snr.Get_SNR_Matrix(reset_source_space_based,reset_LISA_ESA,
                                                                   var_x,sampleRate_x,
                                                                   var_y,sampleRate_y)
                snrplot.Plot_SNR(source_space_based,LISA_ESA,var_x,sample_x,var_y,sample_y,
                    SNRMatrix,display=False)
                if var_y == 'z':
                    [sample_x,sample_y,SNRMatrix] = snr.Get_SNR_Matrix(reset_source_space_based,reset_LISA_ESA,
                                                                       var_x,sampleRate_x,
                                                                       var_y,sampleRate_y)
                    snrplot.Plot_SNR(source_space_based,LISA_ESA,var_x,sample_x,var_y,sample_y,
                        SNRMatrix,display=False,dl_axis=True)



