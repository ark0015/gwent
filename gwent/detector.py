import numpy as np
import os
import astropy.constants as const
import astropy.units as u
import scipy.interpolate as interp
import scipy.stats as stats

from astropy.cosmology import z_at_value
from astropy.cosmology import WMAP9 as cosmo

import gwent
from . import utils

import hasasia.sensitivity as hassens
import hasasia.sim as hassim

current_path = os.path.abspath(gwent.__path__[0])
load_directory = os.path.join(current_path,'LoadFiles/')

class PTA:
    """
    Class to make a PTA instrument using the methods of Hazboun, Romano, Smith (2019) <https://arxiv.org/abs/1907.04341>

    Parameters
    ----------

    name : string
        name of the instrument

    T_obs : float
        the observation time of the PTA in [years]
    N_p : int
        the number of pulsars in the PTA
    sigma : float
        the rms error on the pulsar TOAs in [sec]
    cadence : float
        How often the pulsars are observed in [num/year]

    load_location : string, optional
        If you want to load a PTA curve from a file, it's the file path
    I_type : string, {'E','A','h'}
        Sets the type of input data.
        'E' is the effective strain spectral density $S_{n}(f)$ ('ENSD'),
        'A' is the amplitude spectral density, $\sqrt{S_{n}(f)}$ ('ASD'),
        'h' is the characteristic strain $h_{n}(f)$ ('h')
    GWB_amp : float, optional
        Amplitude of the gravitational wave background added as red noise
    GWB_alpha : float, optional
        the GWB power law, if empty and A_GWB is set, it is assumed to be -2/3
    rn_amp : float, optional
        Individual pulsar red noise amplitude, is a list of [min,max] values from which to uniformly sample
    rn_alpha : float, optional
        Individual pulsar red noise alpha (power law), is a list of [min,max] values from which to uniformly sample
    use_11yr : bool, optional
        Uses the NANOGrav 11yr noise as the individual pulsar noises, if N_p > 34 (the number of pulsars in the 11yr), 
        it draws more pulsars from distributions based on the NANOGrav 11yr pulsar noise
    realistic_noise : float, optional
        Uses realistic noise drawn from distributions based on the NANOGrav 11yr pulsar noise
    f_low : float, optional
        Assigned lowest frequency of PTA (default assigns 1/(5*T_obs))
    f_high : float, optional
        Assigned highest frequency of PTA (default is Nyquist freq cadence/2)
    nfreqs : int, optional
        Number of frequencies in logspace the sensitivity is calculated

    """
    def __init__(self,name,*args,**kwargs):
        self.name = name
        for keys,value in kwargs.items():
            if keys == 'load_location':
                self.load_location = value
            elif keys == 'I_type':
                self.I_type = value
            elif keys == 'GWB_amp':
                self.GWB_amp = value
            elif keys == 'GWB_alpha':
                self.GWB_alpha = value
            elif keys == 'rn_amp':
                self.rn_amp_min = value[0]
                self.rn_amp_max = value[1]
            elif keys == 'rn_alpha':
                self.rn_alpha_min = value[0]
                self.rn_alpha_max = value[1]
            elif keys == 'use_11yr':
                self.use_11yr = value
            elif keys == 'realistic_noise':
                self.realistic_noise = value
            elif keys == 'f_low':
                self.f_low = utils.make_quant(value,'Hz')
            elif keys == 'f_high':
                self.f_high = utils.make_quant(value,'Hz')
            elif keys == 'nfreqs':
                self.nfreqs = value
            else:
                raise ValueError('%s is not an accepted input option.' %keys)

        if not hasattr(self,'nfreqs'):
            self.nfreqs = int(1e3)
        if hasattr(self,'load_location'):
            Load_Data(self)

        if hasattr(self,'use_11yr'):
            self.realistic_noise = True
        else:
            self.use_11yr = False
            if not hasattr(self,'realistic_noise'):
                self.realistic_noise = False

        if hasattr(self,'f_low') and hasattr(self,'f_high'):
            self.fT = np.logspace(np.log10(self.f_low.value),np.log10(self.f_high.value),self.nfreqs)

        if len(args) != 0:
            if len(args) == 1:
                T_obs = args[0]
                self.T_obs = utils.make_quant(T_obs,'yr')
            elif len(args) == 3:
                [T_obs,N_p,cadence] = args
                self.T_obs = utils.make_quant(T_obs,'yr')
                self.N_p = N_p
                self.cadence = utils.make_quant(cadence,'1/yr')
            else:
                [T_obs,N_p,sigma,cadence] = args
                self.T_obs = utils.make_quant(T_obs,'yr')
                self.N_p = N_p
                self.sigma = utils.make_quant(sigma,'s')
                self.cadence = utils.make_quant(cadence,'1/yr')

    @property
    def T_obs(self):
        return self._T_obs
    @T_obs.setter
    def T_obs(self,value):
        self.var_dict = ['T_obs',value]
        if not isinstance(self._return_value,u.Quantity):
            self._return_value = utils.make_quant(self._return_value,'yr')
        self._T_obs = self._return_value

    @property
    def N_p(self):
        return self._N_p
    @N_p.setter
    def N_p(self,value):
        self.var_dict = ['N_p',value]
        self._N_p = self._return_value

    @property
    def cadence(self):
        return self._cadence
    @cadence.setter
    def cadence(self,value):
        self.var_dict = ['cadence',value]
        if not isinstance(self._return_value,u.Quantity):
            self._return_value = utils.make_quant(self._return_value,'1/yr')
        self._cadence = self._return_value

    @property
    def sigma(self):
        return self._sigma
    @sigma.setter
    def sigma(self,value):
        self.var_dict = ['sigma',value]
        if not isinstance(self._return_value,u.Quantity):
            self._return_value = utils.make_quant(self._return_value,'s')
        self._sigma = self._return_value

    @property
    def var_dict(self):
        return self._var_dict
    @var_dict.setter
    def var_dict(self,value):
        utils.Get_Var_Dict(self,value)

    @property
    def fT(self):
        if not hasattr(self,'_fT'):
            #frequency sampled from 1/observation time to nyquist frequency (c/2)
            #5 is the default value for now (from Hazboun et al. 2019)
            T_obs_sec = self.T_obs.to('s').value
            cadence_sec = self.cadence.to('1/s').value
            self._fT = np.logspace(np.log10(1/(5*T_obs_sec)),np.log10(cadence_sec/2),self.nfreqs)*u.Hz
        return self._fT
    @fT.setter
    def fT(self,value):
        self._fT = value
    @fT.deleter
    def fT(self):
        del self._fT

    @property
    def h_n_f(self):
        """Effective Strain Noise Amplitude"""
        if not hasattr(self,'_h_n_f'):
            if hasattr(self,'_I_data'):
                if self._I_Type == 'h':
                    self._h_n_f = self._I_data[:,1]
                elif self._I_Type == 'ENSD':
                    self._h_n_f = np.sqrt(self.S_n_f*self.fT)
                elif self._I_Type == 'ASD':
                    S_n_f_sqrt = self._I_data[:,1]
                    self._h_n_f = S_n_f_sqrt*np.sqrt(self.fT.value)
            else:
                if not hasattr(self,'_sensitivitycurve'):
                    self.Init_PTA()
                self._h_n_f = self._sensitivitycurve.h_c
        return self._h_n_f
    @h_n_f.setter
    def h_n_f(self,value):
        self._h_n_f = value
    @h_n_f.deleter
    def h_n_f(self):
        del self._h_n_f

    @property
    def S_n_f(self):
        #Effective noise power amplitude
        if not hasattr(self,'_S_n_f'):
            if hasattr(self,'_I_data'):
                if self._I_Type == 'ASD':
                    S_n_f_sqrt = self._I_data[:,1]
                    self._S_n_f = S_n_f_sqrt**2/u.Hz
                elif self._I_Type == 'ENSD':
                    self._S_n_f = self._I_data[:,1]/u.Hz
                elif self._I_Type == 'h':
                    self._S_n_f = self.h_n_f**2/self.fT
            else:
                if not hasattr(self,'_sensitivitycurve'):
                    self.Init_PTA()
                self._S_n_f = self._sensitivitycurve.S_eff
                self._S_n_f = utils.make_quant(self._S_n_f,'1/Hz')
        return self._S_n_f
    @S_n_f.setter
    def S_n_f(self,value):
        self._S_n_f = value
    @S_n_f.deleter
    def S_n_f(self):
        del self._S_n_f

    @property
    def f_opt(self):
        #The optimal frequency of the instrument ie. the frequecy at the lowest strain
        if not hasattr(self,'_f_opt'):
            self._f_opt = self.fT[np.argmin(self.h_n_f)]
        return self._f_opt

    def Get_NANOGrav_Param_Distributions(self):
        """Loads in NANOGrav 11yr data

        Notes
        -----
        The file is in the form of sky locations (phi,theta) in the first two columns, 
        Individual Pulsar WN RMS (sigmas), RN Amplitudes, and RN Alphas in the last three columns, respectively.

        """
        NANOGrav_11yr_params_filedirectory = os.path.join(load_directory,'InstrumentFiles/NANOGrav/NANOGrav_11yr_params.txt')
        self._NANOGrav_11yr_params = np.loadtxt(NANOGrav_11yr_params_filedirectory)
        [phis,thetas,sigmas,rn_amps,rn_alphas] = self._NANOGrav_11yr_params

        nbins = 8
        #Add non-zero probability of picking 0 and pi
        new_thetas = np.append(thetas,np.linspace(0.,np.pi,nbins))
        new_phis = np.append(phis,np.linspace(0.,2*np.pi,nbins))
        #add non-zero probability in middle of alphas and amps
        new_rn_alphas = np.append(rn_alphas,np.linspace(min(rn_alphas),max(rn_alphas),nbins))
        new_rn_amps = np.append(rn_amps,np.logspace(min(np.log10(rn_amps)),max(np.log10(rn_amps)),nbins))
        new_sigmas = np.append(sigmas,np.linspace(min(sigmas),max(sigmas),nbins))


        phi_hist = np.histogram(new_phis, bins=nbins,density=True)
        theta_hist = np.histogram(new_thetas, bins=nbins,density=True)
        rn_alpha_hist = np.histogram(new_rn_alphas, bins=nbins,density=True)
        rn_amp_hist = np.histogram(new_rn_amps,
                                   bins=np.logspace(min(np.log10(new_rn_amps)),max(np.log10(new_rn_amps)),nbins),
                                   density=True)
        sigma_hist = np.histogram(new_sigmas, bins=nbins,density=True)

        self._phi_dist = stats.rv_histogram(phi_hist)
        self._theta_dist = stats.rv_histogram(theta_hist)
        self._rn_alpha_dist = stats.rv_histogram(rn_alpha_hist)
        self._rn_amp_dist = stats.rv_histogram(rn_amp_hist)
        self._sigma_dist = stats.rv_histogram(sigma_hist)


    def Draw_New_Pulsars(self):
        if not hasattr(self,'_NANOGrav_11yr_params'):
            self.Get_NANOGrav_Param_Distributions()

        [phis,thetas,sigmas,rn_amps,rn_alphas] = self._NANOGrav_11yr_params

        #34 pulsars in the 11yr dataset (ie. len(phis))
        if self.use_11yr:
            if self.N_p > len(phis):
                N_added_p = self.N_p - len(phis)
                theta_draws = self._theta_dist.rvs(size=N_added_p)
                phi_draws = self._phi_dist.rvs(size=N_added_p)
                rn_amp_draws = self._rn_amp_dist.rvs(size=N_added_p)
                rn_alpha_draws = self._rn_alpha_dist.rvs(size=N_added_p)
                sigma_draws = self._sigma_dist.rvs(size=N_added_p)

                new_thetas = np.append(thetas,theta_draws)
                new_phis = np.append(phis,phi_draws)
                new_rn_amps = np.append(rn_amps,rn_amp_draws)
                new_rn_alphas = np.append(rn_alphas,rn_alpha_draws)
                new_sigmas = np.append(sigmas,sigma_draws)
            else:
                new_thetas = thetas[:self.N_p]
                new_phis = phis[:self.N_p]
                new_rn_amps = rn_amps[:self.N_p]
                new_rn_alphas = rn_alphas[:self.N_p]
                new_sigmas = sigmas[:self.N_p]
        else:
            new_thetas = self._theta_dist.rvs(size=self.N_p)
            new_phis = self._phi_dist.rvs(size=self.N_p)
            new_rn_amps = self._rn_amp_dist.rvs(size=self.N_p)
            new_rn_alphas = self._rn_alpha_dist.rvs(size=self.N_p)
            new_sigmas = self._sigma_dist.rvs(size=self.N_p)

        return [new_phis,new_thetas,new_sigmas,new_rn_amps,new_rn_alphas]

    def Init_PTA(self):
        """Initializes a PTA in hasasia

        Notes
        -----
        See Hazboun, Romano, Smith (2019) <https://arxiv.org/abs/1907.04341> for details

        """
        if self.realistic_noise:
            [phis,thetas,sigmas,rn_amps,rn_alphas] = self.Draw_New_Pulsars()
            #Make a set of psrs with parameters drawn from 11yr distributaions (or real 11yr parameters if use_11yr=True)
            psrs = hassim.sim_pta(timespan=self.T_obs.value,cad=self.cadence.value,sigma=sigmas,\
                phi=phis, theta=thetas, Npsrs=self.N_p,A_rn=rn_amps,alpha=rn_alphas,freqs=self.fT.value)
        else:
            if self.N_p < 33:
                NANOGrav_11yr_params_filedirectory = os.path.join(load_directory,'InstrumentFiles/NANOGrav/NANOGrav_11yr_params.txt')
                self._NANOGrav_11yr_params = np.loadtxt(NANOGrav_11yr_params_filedirectory)
                [phis,thetas,_,_,_] = self._NANOGrav_11yr_params
                thetas = thetas[:self.N_p]
                phis = phis[:self.N_p]
            else:
                #Random Sky Locations of Pulsars
                phis = np.random.uniform(0, 2*np.pi,size=self.N_p)
                cos_theta = np.random.uniform(-1,1,size=self.N_p)
                thetas = np.arccos(cos_theta)

            if hasattr(self,'GWB_amp'):
                if not hasattr(self,'GWB_alpha'):
                    self.GWB_alpha = -2/3.
                #Make a set of psrs with the same parameters with a GWB as red noise
                psrs = hassim.sim_pta(timespan=self.T_obs.value,cad=self.cadence.value,sigma=self.sigma.value,\
                    phi=phis, theta=thetas, Npsrs=self.N_p,A_rn=self.GWB_amp,alpha=self.GWB_alpha,freqs=self.fT.value)
            elif hasattr(self,'rn_amp_min') or hasattr(self,'rn_alpha_min'):
                if not hasattr(self,'rn_amp_min'):
                    rn_amps = np.random.uniform(1e-16,1e-12,size=self.N_p)
                else:
                    rn_amps = np.random.uniform(self.rn_amp_min,self.rn_amp_max,size=self.N_p)
                if not hasattr(self,'rn_alpha_min'):
                    rn_alphas = np.random.uniform(-3/4,1,size=self.N_p)
                else:
                    rn_alphas = np.random.uniform(self.rn_alpha_min,self.rn_alpha_max,size=self.N_p)
                #Make a set of psrs with uniformly sampled red noise
                psrs = hassim.sim_pta(timespan=self.T_obs.value,cad=self.cadence.value,sigma=self.sigma.value,\
                    phi=phis, theta=thetas, Npsrs=self.N_p,A_rn=rn_amps,alpha=rn_alphas,freqs=self.fT.value)
            else:
                #Make a set of psrs with the same parameters
                psrs = hassim.sim_pta(timespan=self.T_obs.value,cad=self.cadence.value,sigma=self.sigma.value,\
                    phi=phis, theta=thetas, Npsrs=self.N_p,freqs=self.fT.value)


        #Get Spectra of pulsars
        spectra= []
        for p in psrs:
             sp = hassens.Spectrum(p,freqs=self.fT.value)
             spectra.append(sp)

        self._sensitivitycurve = hassens.DeterSensitivityCurve(spectra)



class Interferometer:
    """
    Class to make an interferometer

    Parameters
    ----------

    name : string
        name of the instrument

    T_obs : float
        the observation time of the PTA in [years]

    load_location : string, optional
        If you want to load an instrument curve from a file, it's the file path
    I_type : string, {'E','A','h'}
        Sets the type of input data.
        'E' is the effective strain spectral density $S_{n}(f)$ ('ENSD'),
        'A' is the amplitude spectral density, $\sqrt{S_{n}(f)}$ ('ASD'),
        'h' is the characteristic strain $h_{n}(f)$ ('h')

    """
    def __init__(self,name,T_obs,**kwargs):
        self.name = name
        self.T_obs = T_obs

        for keys,value in kwargs.items():
            if keys == 'load_location':
                self.load_location = value
            elif keys == 'I_type':
                self.I_type = value

        if hasattr(self,'load_location'):
            Load_Data(self)

    @property
    def T_obs(self):
        return self._T_obs
    @T_obs.setter
    def T_obs(self,value):
        self.var_dict = ['T_obs',value]
        if not isinstance(self._return_value,u.Quantity):
            self._return_value = utils.make_quant(self._return_value,'yr')
        self._T_obs = self._return_value

    @property
    def var_dict(self):
        return self._var_dict
    @var_dict.setter
    def var_dict(self,value):
        utils.Get_Var_Dict(self,value)

    @property
    def fT(self):
        if not hasattr(self,'_fT'):
            if hasattr(self,'_I_data'):
                self._fT = self._I_data[:,0]*u.Hz
            if isinstance(self,SpaceBased):
                self.Set_T_Function_Type()
        return self._fT
    @fT.setter
    def fT(self,value):
        self._fT = value
    @fT.deleter
    def fT(self):
        del self._fT

    @property
    def f_opt(self):
        """The optimal frequency of the instrument ie. the frequecy at the lowest strain"""
        self._f_opt = self.fT[np.argmin(self.h_n_f)]
        return self._f_opt

    @property
    def P_n_f(self):
        """Strain power sensitivity. """
        raise NotImplementedError('Power Spectral Density method must be defined inside SpaceBased or GroundBased classes.')

    @property
    def S_n_f(self):
        """Effective Noise Power Specral Density"""
        if not hasattr(self,'_S_n_f'):
            if hasattr(self,'_I_data'):
                if self._I_Type == 'ASD':
                    S_n_f_sqrt = self._I_data[:,1]
                    self._S_n_f = S_n_f_sqrt**2/u.Hz
                elif self._I_Type == 'ENSD':
                    self._S_n_f = self._I_data[:,1]/u.Hz
                elif self._I_Type == 'h':
                    self._S_n_f = self.h_n_f**2/self.fT
            else:
                raise NotImplementedError('Effective Noise Power Spectral Density method must be defined inside SpaceBased or GroundBased classes.')
        return self._S_n_f
    @S_n_f.deleter
    def S_n_f(self):
        del self._S_n_f

    @property
    def h_n_f(self):
        """Characteristic Strain/effective strain noise amplitude"""
        if not hasattr(self,'_h_n_f'):
            if hasattr(self,'_I_data') and self._I_Type == 'h':
                self._h_n_f = self._I_data[:,1]
            else:
                self._h_n_f = np.sqrt(self.fT*self.S_n_f)
        return self._h_n_f
    @h_n_f.deleter
    def h_n_f(self):
        del self._h_n_f


class GroundBased(Interferometer):
    """
    Class to make a Ground Based Instrument using the Interferometer base class
    """
    def __init__(self,name,T_obs,**kwargs):
        super().__init__(name,T_obs,**kwargs)
        """
        Currently doesn't do anything differently that Instrument object, can be updated if we ever construct at Ground Based PSD...
        """

    @property
    def P_n_f(self):
        """Power Spectral Density. """
        raise NotImplementedError('Uhhh, can only load from a file for right now....')



class SpaceBased(Interferometer):
    """
    Class to make a Space Based interferometer

    Parameters
    ----------
    L : float
        the armlength the of detector in [meters]
    A_acc : float
        the Amplitude of the Acceleration Noise in [meters/second^2]
    f_acc_break_low : float
        the lower break frequency of the acceleration noise in [Hz]
    f_acc_break_high : float
        the higher break frequency of the acceleration noise in [Hz]
    A_IFO : float
        the amplitude of the interferometer

    T_type : string, {'N','A'}
        Picks the transfer function generation method
        'N' uses the numerically approximated method in Robson, Cornish, and Liu, 2019
        'A' uses the analytic fit in Larson, Hiscock, and Hellings, 2000
    Background : Boolean
        Add in a Galactic Binary Confusion Noise
    f_low : float
        Assigned lowest frequency of instrument (default assigns 10^-5Hz)
    f_high : float
        Assigned highest frequency of instrument (default is 1Hz)
    nfreqs : int
        Number of frequencies in logspace the sensitivity is calculated (default is 1e3)

    """
    def __init__(self,name,T_obs,*args,**kwargs):
        super().__init__(name,T_obs,**kwargs)
        self.name = name
        for keys,value in kwargs.items():
            if keys == 'T_type':
                self.T_type = value
            elif keys == 'Background':
                self.Background = value
            elif keys == 'f_low':
                self.f_low = utils.make_quant(value,'Hz')
            elif keys == 'f_high':
                self.f_high = utils.make_quant(value,'Hz')
            elif keys == 'nfreqs':
                self.nfreqs = value

        if not hasattr(self,'nfreqs'):
            self.nfreqs = int(1e3)
        if not hasattr(self,'f_low'):
            self.f_low = 1e-5*u.Hz
        if not hasattr(self,'f_high'):
            self.f_high = 1.0*u.Hz
        if not hasattr(self,'Background'):
            self.Background = False

        if len(args) != 0:
            [L,A_acc,f_acc_break_low,f_acc_break_high,A_IFO,f_IFO_break] = args
            self.L = utils.make_quant(L,'m')
            self.A_acc = utils.make_quant(A_acc,'m/(s*s)')
            self.f_acc_break_low = utils.make_quant(f_acc_break_low,'Hz')
            self.f_acc_break_high = utils.make_quant(f_acc_break_high,'Hz')
            self.A_IFO = utils.make_quant(A_IFO,'m')
            self.f_IFO_break = utils.make_quant(f_IFO_break,'Hz')

        if not hasattr(self,'load_location'):
            if not hasattr(self,'T_type'):
                self.T_type = 'N'
            self.Set_T_Function_Type()

    @property
    def L(self):
        return self._L
    @L.setter
    def L(self,value):
        self.var_dict = ['L',value]
        if not isinstance(self._return_value,u.Quantity):
            self._return_value = utils.make_quant(self._return_value,'m')
        self._L = self._return_value

    @property
    def A_acc(self):
        return self._A_acc
    @A_acc.setter
    def A_acc(self,value):
        self.var_dict = ['A_acc',value]
        if not isinstance(self._return_value,u.Quantity):
            self._return_value = utils.make_quant(self._return_value,'m/s2')
        self._A_acc = self._return_value

    @property
    def f_acc_break_low(self):
        return self._f_acc_break_low
    @f_acc_break_low.setter
    def f_acc_break_low(self,value):
        self.var_dict = ['f_acc_break_low',value]
        if not isinstance(self._return_value,u.Quantity):
            self._return_value = utils.make_quant(self._return_value,'Hz')
        self._f_acc_break_low = self._return_value

    @property
    def f_acc_break_high(self):
        return self._f_acc_break_high
    @f_acc_break_high.setter
    def f_acc_break_high(self,value):
        self.var_dict = ['f_acc_break_high',value]
        if not isinstance(self._return_value,u.Quantity):
            self._return_value = utils.make_quant(self._return_value,'Hz')
        self._f_acc_break_high = self._return_value

    @property
    def A_IFO(self):
        return self._A_IFO
    @A_IFO.setter
    def A_IFO(self,value):
        self.var_dict = ['A_IFO',value]
        if not isinstance(self._return_value,u.Quantity):
            self._return_value = utils.make_quant(self._return_value,'m')
        self._A_IFO = self._return_value

    @property
    def f_IFO_break(self):
        return self._f_IFO_break
    @f_IFO_break.setter
    def f_IFO_break(self,value):
        self.var_dict = ['f_IFO_break',value]
        if not isinstance(self._return_value,u.Quantity):
            self._return_value = utils.make_quant(self._return_value,'Hz')
        self._f_IFO_break = self._return_value

    @property
    def P_n_f(self):
        """Power Spectral Density"""
        if not hasattr(self,'_P_n_f'):
            if not hasattr(self,'_T_Function_Type'):
                self.Set_T_Function_Type()

            P_acc = self.A_acc**2*(1+(self.f_acc_break_low/self.fT)**2)*(1+(self.fT/(self.f_acc_break_high))**4)/(2*np.pi*self.fT)**4 #Acceleration Noise
            P_IMS = self.A_IFO**2*(1+(self.f_IFO_break/self.fT)**4) #Displacement noise of the interferometric TM--to-TM

            f_trans = const.c/2/np.pi/self.L #Transfer frequency
            self._P_n_f = (P_IMS + 2*(1+np.cos(self.fT.value/f_trans.value)**2)*P_acc)/self.L**2/u.Hz
        return self._P_n_f
    @P_n_f.deleter
    def P_n_f(self):
        del self._P_n_f

    @property
    def S_n_f(self):
        """Effective Noise Power Specral Density"""
        if not hasattr(self,'_S_n_f'):
            if hasattr(self,'_I_data'):
                if self._I_Type == 'ASD':
                    S_n_f_sqrt = self._I_data[:,1]
                    self._S_n_f = S_n_f_sqrt**2/u.Hz
                elif self._I_Type == 'ENSD':
                    self._S_n_f = self._I_data[:,1]/u.Hz
                elif self._I_Type == 'h':
                    self._S_n_f = self.h_n_f**2/self.fT
            else:
                S_n_f = self.P_n_f/self.transferfunction**2
                if self.Background:
                    self._S_n_f= S_n_f+self.Add_Background()
                else:
                    self._S_n_f = S_n_f
        return self._S_n_f
    @S_n_f.deleter
    def S_n_f(self):
        del self._S_n_f

    def Load_Transfer_Function(self):
        #Numerical transfer function
        Numerical_Transfer_Function_filedirectory = os.path.join(load_directory,'NumericalTransferFunction/transfer.dat')
        Numerical_Transfer_Function_data = np.loadtxt(Numerical_Transfer_Function_filedirectory)
        self._transferfunctiondata = Numerical_Transfer_Function_data

    def Get_Numeric_Transfer_Function(self):
        if not hasattr(self,'_transferfunctiondata'):
            self.Load_Transfer_Function()

        fc = const.c/(2*self.L)  #light round trip freq
        LISA_Transfer_Function_f = fc*self._transferfunctiondata[:,0]

        idx_f_5 = np.abs(LISA_Transfer_Function_f-self.f_low).argmin()
        idx_f_1 = np.abs(LISA_Transfer_Function_f-self.f_high).argmin()

        #3/10 is normalization 2/5sin(openingangle)
        #Some papers use 3/20, not summing over 2 independent low-freq data channels
        self.transferfunction = np.sqrt(3/10)*self._transferfunctiondata[idx_f_5:idx_f_1,1]
        self.fT = LISA_Transfer_Function_f[idx_f_5:idx_f_1]

    def Get_Analytic_Transfer_Function(self):
        #Response function approximation from Calculation described by Cornish, Robson, Liu 2019
        self.fT = np.logspace(np.log10(self.f_low.value),np.log10(self.f_high.value),self.nfreqs)*u.Hz
        f_L = const.c/2/np.pi/self.L #Transfer frequency
        #3/10 is normalization 2/5sin(openingangle)
        R_f = 3/10/(1+0.6*(self.fT/f_L)**2)
        self.transferfunction = np.sqrt(R_f)

    def Set_T_Function_Type(self):
        if self.T_type == 'n' or self.T_type == 'N':
            self._T_type = 'numeric'
        elif self.T_type == 'a' or self.T_type == 'A':
            self._T_type = 'analytic'
        else:
            print('\nYou can get the transfer function via 2 methods:')
            print(' *To use the numerically approximated method in Robson, Cornish, and Liu, 2019, input "N".')
            print(' *To use the analytic fit in Larson, Hiscock, and Hellings, 2000, input "A".')
            calc_type = input('Please select the calculation type: ')
            self.Set_T_Function_Type(calc_type)
        if hasattr(self,'_T_type'):
            if self._T_type == 'numeric':
                self.Get_Numeric_Transfer_Function()
            if self._T_type == 'analytic':
                self.Get_Analytic_Transfer_Function()


    def Add_Background(self):
        """
        Galactic confusions noise parameters for 6months, 1yr, 2yr, and 4yr
        corresponding to array index 0,1,2,3 respectively
        """
        A = 9e-45
        a = np.array([0.133,0.171,0.165,0.138])
        b = np.array([243,292,299,-221])
        k = np.array([482,1020,611,521])
        g = np.array([917,1680,1340,1680])
        f_k = np.array([0.00258,0.00215,0.00173,0.00113])

        if self.T_obs < 1.*u.yr:
            index = 0
        elif self.T_obs >= 1.*u.yr and self.T_obs < 2.*u.yr:
            index = 1
        elif self.T_obs >= 2.*u.yr and self.T_obs < 4.*u.yr:
            index = 2
        else:
            index = 3
        f = self.fT.value
        return A*np.exp(-(f**a[index])+(b[index]*f*np.sin(k[index]*f)))\
                *(f**(-7/3))*(1 + np.tanh(g[index]*(f_k[index]-f))) #White Dwarf Background Noise

def Load_Data(detector):
    """
    Function to load in a file to initialize any detector.

    Parameters
    ----------
    detector : object
        Instance of a detector class

    """
    if not hasattr(detector,'I_type'):
        print('Is the data:')
        print(' *Effective Noise Spectral Density - "E"')
        print(' *Amplitude Spectral Density- "A"')
        print(' *Effective Strain - "h"')
        detector.I_type = input('Please enter one of the answers in quotations: ')
        Load_Data(detector)

    if detector.I_type == 'E' or detector.I_type == 'e':
        detector._I_Type = 'ENSD'
    elif detector.I_type == 'A' or detector.I_type == 'a':
        detector._I_Type = 'ASD'
    elif detector.I_type == 'h' or detector.I_type == 'H':
        detector._I_Type = 'h'
    else:
        print('Is the data:')
        print(' *Effective Noise Spectral Density - "E"')
        print(' *Amplitude Spectral Density- "A"')
        print(' *Effective Strain - "h"')
        detector.I_type = input('Please enter one of the answers in quotations: ')
        Load_Data(detector)

    detector._I_data = np.loadtxt(detector.load_location)
    detector.fT = detector._I_data[:,0]*u.Hz
