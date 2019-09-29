#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib import rc
import os, sys
import pytest

import astropy.constants as const
import astropy.units as u
from astropy.cosmology import z_at_value
from astropy.cosmology import WMAP9 as cosmo
from fractions import Fraction

import hasasia.sensitivity as hassens
import hasasia.sim as hassim
import hasasia.skymap as hassky

import gwent
from gwent import detector
from gwent import binary

load_directory = gwent.__path__[0] + '/LoadFiles/InstrumentFiles/'

# axissize = 14
# labelsize = 16
# legendsize = 12
# figsize = (10,8)
# colornorm = colors.Normalize(vmin=0.0, vmax=5.0)
# linesize = 3


# ####################################################################
# # Initialize different instruments

# ### aLIGO

Ground_T_obs = 4*u.yr


#aLIGO
@pytest.fixture
def aLIGO():
    aLIGO_filedirectory = load_directory + 'aLIGO/StrainFiles/'
    aLIGO_filename = 'aLIGODesign.txt'
    aLIGO_filelocation = aLIGO_filedirectory + aLIGO_filename

    aLIGO = detector.GroundBased('aLIGO',Ground_T_obs,
                                 load_location=aLIGO_filelocation,I_type='A')
    return aLIGO

#Einstein Telescope
@pytest.fixture
def ET():
    ET_filedirectory = load_directory + 'EinsteinTelescope/StrainFiles/'
    ET_filename = 'ET_B_data.txt'
    ET_filelocation = ET_filedirectory + ET_filename

    ET = detector.GroundBased('ET',Ground_T_obs,
                              load_location=ET_filelocation,I_type='A')
    return ET

SpaceBased_T_obs = 4*u.yr
LISA_Other_filedirectory = load_directory + 'LISA_Other/StrainFiles/'

#Martin data
def test_LISA_Martin():
    LISA_Martin_filename = 'LISA_Allocation_S_h_tot.txt'
    LISA_Martin_filelocation = LISA_Other_filedirectory + LISA_Martin_filename

    #Should be ENSD
    LISA_Martin = detector.SpaceBased('LISA_Martin',SpaceBased_T_obs,
                                      load_location=LISA_Martin_filelocation,
                                      I_type='E')
# ### LISA Neil Cornish data
def test_LISA_Cornish():
    #Neil Cornish data
    LISA_Neil_filedirectory = load_directory + 'LISA_Neil/StrainFiles/'
    LISA_Neil_filename = 'LISA_sensitivity.txt'
    LISA_Neil_filelocation = LISA_Neil_filedirectory + LISA_Neil_filename

    #Should be ENSD
    LISA_Neil = detector.SpaceBased('LISA_Neil',SpaceBased_T_obs,
                                    load_location=LISA_Neil_filelocation,
                                    I_type='E')

# ### LISA Larson Sensitivity Curve
def test_LISA_Larson():
    #Larson Sensitivity Curve
    LISA_Larson_filename = 'scg_6981.dat'
    LISA_Larson_filelocation = LISA_Other_filedirectory + LISA_Larson_filename

    #Should be ASD
    LISA_Larson = detector.SpaceBased('LISA_Larson',SpaceBased_T_obs,
                                      load_location=LISA_Larson_filelocation,
                                      I_type='A')

# ### NANOGrav continuous wave sensitivity

NANOGrav_filedirectory = load_directory + 'NANOGrav/StrainFiles/'

#NANOGrav continuous wave sensitivity
NANOGrav_background = 4e-16 # Unsubtracted GWB amplitude: 0,4e-16
NANOGrav_dp = 0.95 #Detection Probablility: 0.95,0.5
NANOGrav_fap = 0.0001 #False Alarm Probability: 0.05,0.003,0.001,0.0001
NANOGrav_Tobs = 15 #Observation years: 15,20,25

NANOGrav_filename = 'cw_simulation_Ared_' + str(NANOGrav_background)
NANOGrav_filename += '_dp_' + str(NANOGrav_dp) + '_fap_' + str(NANOGrav_fap)
NANOGrav_filename += '_T_' + str(NANOGrav_Tobs) + '.txt'
NANOGrav_filelocation = NANOGrav_filedirectory + NANOGrav_filename

def test_pta_Migarelli_no_GWB():
    NANOGrav_Mingarelli_no_GWB = detector.PTA('NANOGrav_Mingarelli_no_GWB',
                                              load_location=NANOGrav_filelocation)

#NANOGrav continuous wave sensitivity
NANOGrav_background_2 = 0 # Unsubtracted GWB amplitude: 0,4e-16
NANOGrav_dp_2 = 0.95 #Detection Probablility: 0.95,0.5
NANOGrav_fap_2 = 0.0001 #False Alarm Probability: 0.05,0.003,0.001,0.0001
NANOGrav_Tobs_2 = 15 #Observation years: 15,20,25

NANOGrav_filename_2 = 'cw_simulation_Ared_' + str(NANOGrav_background_2)
NANOGrav_filename_2 += '_dp_' + str(NANOGrav_dp_2)
NANOGrav_filename_2 += '_fap_' + str(NANOGrav_fap_2) + '_T_'
NANOGrav_filename_2 += str(NANOGrav_Tobs_2) + '.txt'
NANOGrav_filelocation_2 = NANOGrav_filedirectory + NANOGrav_filename_2

def test_pta_Migarelli_GWB():
    NANOGrav_Mingarelli_GWB = detector.PTA('NANOGrav_Mingarelli_GWB',
                                           load_location=NANOGrav_filelocation_2)

# ### SKA  parameters and methods from arXiv:0804.4476 section 7.1
###############################################
#SKA calculation using parameters and methods from arXiv:0804.4476 section 7.1
#sigma_rms timing residuals in nanoseconds to seconds
sigma_SKA = 10*u.ns.to('s')*u.s
T_SKA = 15*u.yr #Observing time in years
N_p_SKA = 20 #Number of pulsars
#Avg observation cadence of 1 every week in [number/yr]
cadence_SKA = 1/(u.wk.to('yr')*u.yr)

@pytest.fixture
def SKA_Hazboun():
    SKA_Hazboun = detector.PTA('SKA_Hazboun',T_SKA,N_p_SKA,
                               sigma_SKA,cadence_SKA)
    return SKA_Hazboun

def test_pta_SKA_wRN():
    SKA_Hazboun_wRN = detector.PTA('SKA_Hazboun_wRN',T_SKA,N_p_SKA,sigma_SKA,
                                   cadence_SKA,A_rn=[1e-16,1e-12],
                                   alpha_rn=[-3/4,1])
def test_pta_SKA_wGWB():
    SKA_Hazboun_wGWB = detector.PTA('SKA_Hazboun_wGWB',T_SKA,
                                    N_p_SKA,sigma_SKA,cadence_SKA,
                                    A_GWB=4e-16)

# #### Using Jeff's Methods/code https://arxiv.org/abs/1907.04341

# ### NANOGrav 11.5yr parameters https://arxiv.org/abs/1801.01837

###############################################
#NANOGrav calculation using 11.5yr parameters https://arxiv.org/abs/1801.01837
#rms timing residuals in nanoseconds to seconds
sigma_nano = 100*u.ns.to('s')*u.s
T_nano = 15*u.yr #Observing time in years
N_p_nano = 18 #Number of pulsars
#Avg observation cadence of 1 every 2 weeks in number/year
cadence_nano = 1/(2*u.wk.to('yr')*u.yr)

def test_pta_NANOGrav():
    NANOGrav_Hazboun = detector.PTA('NANOGrav_Hazboun',T_nano,
                                    N_p_nano,sigma_nano,cadence_nano)
def test_pta_NANOGrav_wRN():
    NANOGrav_Hazboun_wRN = detector.PTA('NANOGrav_Hazboun_wRN',T_nano,N_p_nano,
                                        sigma_nano,cadence_nano,A_rn=[1e-16,1e-12],
                                        alpha_rn=[-3/4,1])
def test_pta_NANOGrav_wGWB():
    NANOGrav_Hazboun_wGWB = detector.PTA('NANOGrav_Hazboun_wGWB',T_nano,N_p_nano,
                                         sigma_nano,cadence_nano,A_GWB=4e-16)

# ####################################################################
# # Calculate LISA amplitude spectral densities for various models

# In[42]:


L = 2.5*u.Gm  #armlength in Gm
L = L.to('m')
LISA_T_obs = 4*u.yr


# ### LISA Calculation from https://arxiv.org/pdf/1702.00786.pdf (Amaro-Seaone 2017)

# In[46]:


f_acc_break_low = .4*u.mHz.to('Hz')*u.Hz
f_acc_break_high = 8.*u.mHz.to('Hz')*u.Hz
f_IMS_break = 2.*u.mHz.to('Hz')*u.Hz
A_acc = 3e-15*u.m/u.s/u.s
A_IMS = 10e-12*u.m

Background = False

@pytest.fixture
def ESA_LISA():
    ESA_LISA = detector.SpaceBased('ESA_LISA', LISA_T_obs, L, A_acc,
                                    f_acc_break_low, f_acc_break_high,
                                    A_IMS, f_IMS_break, Background=Background)
    return ESA_LISA

#Neil Calculation from https://arxiv.org/pdf/1803.01944.pdf
f_acc_break_low = .4*u.mHz.to('Hz')*u.Hz
f_acc_break_high = 8.*u.mHz.to('Hz')*u.Hz
f_IMS_break = 2.*u.mHz.to('Hz')*u.Hz
A_acc = 3e-15*u.m/u.s/u.s
A_IMS = 1.5e-11*u.m
Background = False

def test_LISA_asd1():
    Neil_LISA = detector.SpaceBased('Neil_LISA', LISA_T_obs, L, A_acc,
                                    f_acc_break_low, f_acc_break_high,
                                    A_IMS, f_IMS_break, Background=Background)

# #######################################################################
# # BBH strain calculation

#Vars = [M,q,chi1,chi2,z]
M = [1e6,65.0,1e10]
q = [1.0,18.0,1.0]
x1 = [0.95,0.0,-0.95]
x2 = [0.95,0.0,-0.95]
z = [3.0,0.093,20.0]
inc = 0.0 #Doesn't really work...

Vars1 = [M[0],q[0],x1[0],x2[0],z[0]]
Vars2 = [M[1],q[1],x1[1],x2[1],z[1]]
Vars3 = [M[2],q[2],x1[2],x2[2],z[2]]
Vars4 = [M[1],q[0],x1[1],x2[1],z[1]]

def test_BBHStrain(ESA_LISA,aLIGO,SKA_Hazboun,ET):
    source_1 = binary.BBHFrequencyDomain(M[0],q[0],x1[0],x2[0],z[0],
                                         inc,instrument=ESA_LISA)
    source_2 = binary.BBHFrequencyDomain(M[1],q[1],x1[1],x2[1],z[1],
                                         inc,instrument=aLIGO)
    source_3 = binary.BBHFrequencyDomain(M[2],q[2],x1[2],x2[2],z[2],
                                         inc,instrument=SKA_Hazboun)
    source_4 = binary.BBHFrequencyDomain(M[1],q[0],x1[1],x2[1],z[1],
                                         inc,instrument=ET)


# ### Numerical Relativity from EOB subtraction

EOBdiff_filedir = gwent.__path__[0] + '/LoadFiles//DiffStrain/EOBdiff/'

def test_NR_EOB():
    diff0002 = binary.BBHTimeDomain(M[1],q[0],z[1],
                                    load_location=EOBdiff_filedir+'diff0002.dat')
    diff0114 = binary.BBHTimeDomain(M[1],q[0],z[1],
                                    load_location=EOBdiff_filedir+'diff0114.dat')
    diff0178 = binary.BBHTimeDomain(M[1],q[0],z[1],
                                    load_location=EOBdiff_filedir+'diff0178.dat')
    diff0261 = binary.BBHTimeDomain(M[1],q[0],z[1],
                                    load_location=EOBdiff_filedir+'diff0261.dat')
    diff0303 = binary.BBHTimeDomain(M[1],q[0],z[1],
                                    load_location=EOBdiff_filedir+'diff0303.dat')
