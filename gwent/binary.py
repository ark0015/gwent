import numpy as np
import os
import astropy.constants as const
import astropy.units as u
import scipy.interpolate as interp
from astropy.cosmology import z_at_value
from astropy.cosmology import WMAP9 as cosmo

from waveform import Get_Waveform
import utils

current_path = os.getcwd()
splt_path = current_path.split("/")
top_path_idx = splt_path.index('gwent')
top_directory = "/".join(splt_path[0:top_path_idx+1])

class BinaryBlackHole:
    def __init__(self,M,q,z,**kwargs):
        """args order: M,q,chi1,chi2,z,inc
            kwargs: f_low=1e-5,nfreqs=int(1e3)"""
        self.M = M
        self.q = q
        self.z = z

        for keys,value in kwargs.items():
            if keys == 'load_location':
                self.load_location = value

        if hasattr(self,'load_location'):
            self.Load_Data()

    @property
    def M(self):
        self._M = utils.make_quant(self._M,'M_sun')
        return self._M
    @M.setter
    def M(self,value):
        self.var_dict = ['M',value]
        self._M = self._return_value

    @property
    def q(self):
        return self._q
    @q.setter
    def q(self,value):
        self.var_dict = ['q',value]
        self._q = self._return_value

    @property
    def z(self):
        return self._z
    @z.setter
    def z(self,value):
        self.var_dict = ['z',value]
        self._z = self._return_value

    @property
    def h_f(self):
        if not hasattr(self,'_h_f'):
            raise NotImplementedError('The strain must be defined inside BBHFrequencyDomain or BBHTimeDomain classes.')
        return self._h_f
    @h_f.setter
    def h_f(self,value):
        self._h_f = value

    @property
    def f(self):
        if not hasattr(self,'_f'):
            raise NotImplementedError('Interferometer frequency must be defined inside SpaceBased or GroundBased classes.')
        return self._f
    @f.setter
    def f(self,value):
        self._f = value

    @property
    def var_dict(self):
        return self._var_dict
    @var_dict.setter
    def var_dict(self,value):
        utils.Get_Var_Dict(self,value)
    
    def Load_Data(self):
        if hasattr(self,'load_location'):
            if os.path.exists(self.load_location):
                self._load_data = np.loadtxt(self.load_location)
            else:
                raise IOError('File %s does not exist, please assign load_location a correct filepath.' %self.load_location)
        else:
            raise ValueError('load_location is not assigned, please set with name_of_BBH.load_location="path/to/file".')

class BBHFrequencyDomain(BinaryBlackHole):
    def __init__(self,*args,**kwargs):
        """args order: M,q,chi1,chi2,z,inc
            kwargs: f_low=1e-5,nfreqs=int(1e3)"""

        [M,q,chi1,chi2,z,inc] = args
        super().__init__(M,q,z,**kwargs)

        self.chi1 = chi1
        self.chi2 = chi2
        self.inc = inc

        for keys,value in kwargs.items():
            if keys == 'f_low':
                self.f_low = value
            elif keys == 'f_high':
                self.f_high = value
            elif keys == 'nfreqs':
                self.nfreqs = value
            elif keys == 'instrument':
                self.instrument = value
                self.Check_Freq_Evol()
        if not hasattr(self,'nfreqs'):
            self.nfreqs = int(1e3)
        if not hasattr(self,'f_low'):
            self.f_low = 1e-5*u.Hz

        self.Get_Fitcoeffs()

    @property
    def chi1(self):
        return self._chi1
    @chi1.setter
    def chi1(self,value):
        self.var_dict = ['chi1',value]
        self._chi1 = self._return_value

    @property
    def chi2(self):
        return self._chi2
    @chi2.setter
    def chi2(self,value):
        self.var_dict = ['chi2',value]
        self._chi2 = self._return_value

    @property
    def inc(self):
        return self._inc
    @inc.setter
    def inc(self,value):
        self.var_dict = ['inc',value]
        self._inc = self._return_value

    @property
    def instrument(self):
        return self._instrument
    @instrument.setter
    def instrument(self,value):
        self._instrument = value

    @property
    def h_gw(self):
        if not hasattr(self,'_h_gw'):
            if not hasattr(self,'f_init'):
                if hasattr(self,'_instrument'):
                    self.Check_Freq_Evol()
                else:
                    raise ValueError('No instrument assigned, please fix it. '\
                        'Try: "source.instrument = instrument".')
                self._h_gw = Get_Mono_Strain(self,self.instrument.f_opt).to('')
            else:
                self._h_gw = Get_Mono_Strain(self,self.f_init).to('')
        return self._h_gw
    @h_gw.setter
    def h_gw(self,value):
        self._h_gw = value
    @h_gw.deleter
    def h_gw(self):
        del self._h_gw

    @property
    def h_f(self):
        if not hasattr(self,'_h_f'):
            if not (hasattr(self,'_phenomD_f') and hasattr(self,'_phenomD_h')):
                self.Get_PhenomD_Strain()
            [_,self._h_f] = Strain_Conv(self,self._phenomD_f,self._phenomD_h)
        return self._h_f
    @h_f.deleter
    def h_f(self):
        del self._h_f

    @property
    def f(self):
        if not hasattr(self,'_f'):
            if not (hasattr(self,'_phenomD_f') and hasattr(self,'_phenomD_h')):
                self.Get_PhenomD_Strain()
            [self._f,_] = Strain_Conv(self,self._phenomD_f,self._phenomD_h)
        return self._f
    @f.deleter
    def f(self):
        del self._f

    def Get_Fitcoeffs(self):
        #load QNM fitting files for speed later
        fit_coeffs_filedirectory = top_directory + '/docs/LoadFiles/PhenomDFiles/fitcoeffsWEB.dat'
        self._fitcoeffs = np.loadtxt(fit_coeffs_filedirectory)

    def Get_PhenomD_Strain(self):
        if not hasattr(self,'_fitcoeffs'):
            self.Get_Fitcoeffs()
        [self._phenomD_f,self._phenomD_h] = Get_Waveform(self)

    def Get_Time_From_Merger(self,f_obs):
        """Takes in an initally observed frequency, outputs the binary's time
            from merger.
        """
        m_conv = const.G/const.c**3 #Converts M = [M] to M = [sec]
        eta = self.q/(1+self.q)**2

        M_time = self.M.to('kg')*m_conv
        M_chirp = eta**(3/5)*M_time

        f_obs_source = f_obs*(1+self.z)
        return 5*(M_chirp)**(-5/3)*(8*np.pi*f_obs_source)**(-8/3)

    def Get_Source_Freq(self,tau):
        """Takes in a time from merger (tau) and calculates the binary's
            GW frequency at that time. Assumes tau is in the source frame
        """
        m_conv = const.G/const.c**3 #Converts M = [M] to M = [sec]
        eta = self.q/(1+self.q)**2

        M_time = self.M.to('kg')*m_conv
        M_chirp = eta**(3/5)*M_time

        return 1./8./np.pi/M_chirp*(5*M_chirp/tau)**(3./8.)

    def Check_Freq_Evol(self):
        #####################################
        #If the initial observed time from merger is less than the time observed
        #(ie t_init-T_obs < 0 => f_evolve is complex),
        #the BBH will or has already merged during the observation

        #If the initial observed time from merger is greater than the time observed
        #(ie t_init-T_obs > 0 => f_evolve is real),
        #And if the frequency of the binary does evolve over more than one bin,
        #(ie f_T_obs-f_init < 1/T_obs), it is monochromatic, so we set the frequency
        #to the optimal frequency of the detector

        #Otherwise it is chirping and evolves over the observation and we
        #set the starting frequency we observe it at to f(Tobs), which is the 
        #frequency at an observation time before merger
        #####################################
        m_conv = const.G/const.c**3 #Converts M = [M] to M = [sec]
        eta = self.q/(1+self.q)**2

        M_time = self.M.to('kg')*m_conv
        M_chirp_source = eta**(3/5)*M_time

        T_obs = utils.make_quant(self.instrument.T_obs,'s')
        T_obs_source = T_obs/(1+self.z)
        #print('T_obs_source: ',T_obs_source.to('yr'))

        """
        #Assumes t_init is in source frame, can either be randomly drawn
        t_init_source = np.random.uniform(0,100)*u.yr
        """
        #Assumes f_init is the optimal frequency in the instrument frame to get t_init_source
        self.f_init = self.instrument.f_opt
        t_init_source = self.Get_Time_From_Merger(self.f_init)

        #f(T_obs), the frequency of the source at T_obs before merger
        f_T_obs_source = self.Get_Source_Freq(T_obs_source)
        #f(T_obs) in the instrument frame
        self.f_T_obs = f_T_obs_source/(1+self.z)

        #t_init_source = make_quant(t_init_source,'s')
        #print('t_init_source: ',t_init_source.to('yr'))

        #f_init_source = self.Get_Source_Freq(t_init_source)
        #print('f_init_source: ',f_init_source)
        
        #self.f_init = f_init_source/(1+self.z)
        #print('f_init_inst: ',self.f_init)
        
        #f_after_T_obs_source = self.Get_Source_Freq((t_init_source-T_obs_source))
        #print('f_end_source: ',f_after_T_obs_source)
        
        #self.f_T_obs = f_after_T_obs_source/(1+self.z)
        #print('f_T_obs_inst: ',self.f_T_obs)
        
        #delf_obs_source_exact = f_after_T_obs_source-f_init_source
        #print('delf_source: ',delf_obs_source_exact)
        
        #from eqn 41 from Hazboun,Romano, and Smith (2019) https://arxiv.org/abs/1907.04341
        #Uses binomial expansion of f_T_obs_inst - f_init_inst
        #Will not ever be imaginary, so probably better to use
        delf_obs_source_approx = 1./8./np.pi/M_chirp_source*(5*M_chirp_source/t_init_source)**(3./8.)*(3*T_obs_source/8/t_init_source)
        #print('delf_Jeff: ',delf_obs_source_approx)
        
        delf_obs =  delf_obs_source_approx/(1+self.z)
        #print('delf_obs: ',delf_obs)

        """
        #Old way I was doing this....
        M_redshifted_time = self.M.to('kg')*(1+self.z)*m_conv
        M_chirp = eta**(3/5)*M_redshifted_time
        t_init = 5*(M_chirp)**(-5/3)*(8*np.pi*self.instrument.f_opt)**(-8/3)
        #print('t_init: ', t_init.to('yr'))
        #f(t) from eqn 40
        f_evolve = 1./8./np.pi/M_chirp*(5*M_chirp/(t_init-T_obs))**(3./8.)
        f_T_obs = 1./8./np.pi/M_chirp*(5*M_chirp/T_obs)**(3./8.)
        print(f_T_obs)
        #from eqn 41 from Hazboun,Romano, and Smith (2019) https://arxiv.org/abs/1907.04341
        delf = 1./8./np.pi/M_chirp*(5*M_chirp/t_init)**(3./8.)*(3*T_obs/8/t_init)
        print('delf old: ',delf)
        print('')
        """
        
        if delf_obs < (1/T_obs):
            self.ismono = True
        else:
            self.ismono = False

    


class BBHTimeDomain(BinaryBlackHole):
    def __init__(self,*args,**kwargs):
        """args order: M,q,z"""
        super().__init__(*args,**kwargs)
        self.Get_hf_from_hcross_hplus()

    @property
    def t(self):
        if not hasattr(self,'_t'):
            self._t = self._load_data[:,0]
        self._t = utils.make_quant(self._t,'s')
        return self._t

    @property
    def h_plus_t(self):
        if not hasattr(self,'_h_plus_t'):
            self._h_plus_t = self._load_data[:,1]
        return self._h_plus_t

    @property
    def h_cross_t(self):
        if not hasattr(self,'_h_cross_t'):
            self._h_cross_t = self._load_data[:,1]
        return self._h_cross_t

    @property
    def h_f(self):
        if not hasattr(self,'_h_f'):
            [natural_f,natural_h] = self.Get_hf_from_hcross_hplus()
            [_,self._h_f] = Strain_Conv(self,natural_f,natural_h)
        return self._h_f
    @h_f.deleter
    def h_f(self):
        del self._h_f

    @property
    def f(self):
        if not hasattr(self,'_f'):
            [natural_f,natural_h] = self.Get_hf_from_hcross_hplus()
            [self._f,_] = Strain_Conv(self,natural_f,natural_h)
        return self._f
    @f.deleter
    def f(self):
        del self._f


    def Get_hf_from_hcross_hplus(self,interp_res='coarse',windowing='left'):
        """Converts dimensionless, time domain strain to frequency space

        Parameters
        ----------
        interp_res : {'coarse','fine'}, optional
            'coarse' 

        """

        #Interpolate time to evenly sampled data, can be fine or coarse
        diff_t = np.diff(self.t.value)
        if interp_res == 'fine':
            dt = min(diff_t)
        elif interp_res == 'coarse':
            dt = max(diff_t)

        interp_t = np.arange(self.t[0].value,self.t[-1].value,dt)
        #interpolate strain to evenly sampled data for FFT
        h_cross_t = interp.interp1d(self.t,self.h_cross_t,kind='cubic')
        h_plus_t = interp.interp1d(self.t,self.h_plus_t,kind='cubic')
        interp_h_cross_t = h_cross_t(interp_t)
        interp_h_plus_t = h_plus_t(interp_t)

        #Filter/Window
        hann_window = np.hanning(len(interp_t)) #Two sided
        if windowing == 'left':
            #########################
            """Applies window to first (left) half"""
            first_half = hann_window[:int(len(interp_t)/2)] # Only need tapering on first half of waveform
            second_half = np.ones(len(interp_t)-len(first_half)) #no windowing on second half of waveform
            #########################
            window = np.append(first_half,second_half) # Only apply window to first half of waveform
        elif windowing == 'right':
            #########################
            """Applies window to second (right) half"""
            second_half = hann_window[int(len(interp_t)/2):] # Only need tapering on second half of waveform
            first_half = np.ones(len(interp_t)-len(second_half)) #no windowing on first half of waveform
            #########################
            window = np.append(first_half,second_half)
        elif windowing == 'all':
            window = hann_window
        #Window!     
        win_h_cross_t = np.multiply(interp_h_cross_t,window)
        win_h_plus_t = np.multiply(interp_h_plus_t,window)

        #FFT the two polarizations
        h_cross_f = np.fft.fft(win_h_cross_t)
        h_plus_f = np.fft.fft(win_h_plus_t)
        freqs = np.fft.fftfreq(len(interp_t),d=dt)

        #cut = np.abs(freqs).argmax() #Cut off the negative frequencies
        f_cut_low = 3e-3 #Low Cutoff frequency
        f_cut_high = 1.5e-1 #High Cutoff frequency 
        cut_low = np.abs(freqs-f_cut_low).argmin() #Cut off frequencies lower than a frequency
        cut_high = np.abs(freqs-f_cut_high).argmin() #Cut off frequencies higher than a frequency
        #cut=int(len(freqs)*0.9) #Cut off percentage of frequencies
        h_cross_f = h_cross_f[cut_low:cut_high]
        h_plus_f = h_plus_f[cut_low:cut_high]
        natural_f = freqs[cut_low:cut_high]
        
        #Combine them for raw spectral power
        natural_h_f = np.sqrt((np.abs(h_cross_f))**2 + (np.abs(h_plus_f))**2)
        return [natural_f,natural_h_f]



def Strain_Conv(source,natural_f,natural_h):
    """Converts frequency and strain in natural units (G=c=1) to Hertz and dimensionless, respectively.

    Parameters
    ----------
    source
        Instance of gravitational wave source class
    natural_f : array [Mf]
        the frequency of the source in natural units (G=c=1)
    natural_h : array [Mf]
        the strain of the source in natural units (G=c=1)


    """
    DL = cosmo.luminosity_distance(source.z)
    DL = DL.to('m')

    m_conv = const.G/const.c**3 #Converts M = [M] to M = [sec]
    M_redshifted_time = source.M.to('kg')*(1+source.z)*m_conv
    
    #frequency and strain of source in detector frame
    freq_conv = 1/M_redshifted_time
    #Normalized factor to match Stationary phase approx at low frequencies?
    #Changed from sqrt(5/16/pi)
    strain_conv = np.sqrt(1/4/np.pi)*(const.c/DL)*M_redshifted_time**2
    
    f = natural_f*freq_conv
    h_f = natural_h*strain_conv
    return [f,h_f]

def Get_Char_Strain(source):
    """Converts source strain to characteristic strain

    Parameters
    ----------
    source
        Instance of gravitational wave source class

    """
    h_char = np.sqrt(4*source.f**2*source.h_f**2)
    return h_char

def Get_Mono_Strain(source,f_gw,strain_const='Averaged'):
    """Calculates the strain from a binary black hole.

    Parameters
    ----------
    f_gw : float
        The source frequency of the gravitational wave.
    strain_const : {'Averaged','UseInc','Optimal'}, optional
        'Averaged' gives the sky and inclination averaged strain from Robson et al. 2019 (eqn 27) https://arxiv.org/pdf/1803.01944.pdf
        'UseInc' uses the source inclination value from Rosado, Sesana, and Gair (2015) https://arxiv.org/abs/1503.04803, will be innacurate because nothing else accounts for inclination
        'Optimal' gives the optimally oriented, face-on, inclination (ie. inc=0) value

    Returns
    -------
    float
        the strain of a monochromatic source in the dector frame

    """
    f_gw = utils.make_quant(f_gw,'Hz')
    if isinstance(strain_const,str):
        DL = cosmo.luminosity_distance(source.z)
        DL = DL.to('m')

        #Converts M = [M] to M = [sec]
        m_conv = const.G/const.c**3 

        eta = source.q/(1+source.q)**2
        M_redshifted_time = source.M.to('kg')*(1+source.z)*m_conv
        M_chirp = eta**(3/5)*M_redshifted_time

        if strain_const == 'UseInc':
            a = 1+np.cos(source.inc)**2
            b = -2*np.cos(source.inc)
            const_val = 2*np.sqrt(.5*(a**2+b**2))
        elif strain_const == 'Optimal':
            const_val = 4.
        elif strain_const == 'Averaged':
            const_val = 8/np.sqrt(5)
        else:
            raise ValueError('Can only use "UseInc" or "Averaged" monochromatic strain calculation.')

        return const_val*(const.c/DL)*(np.pi*f_gw)**(2./3.)*M_chirp**(5./3.)
    else:
        raise ValueError('Can only use "UseInc" or "Averaged" monochromatic strain calculation.')