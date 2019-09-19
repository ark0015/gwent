import numpy as np

def Get_Waveform(source,pct_of_peak=0.01):
    """Uses Mass Ratio (q <= 18), aligned spins (abs(a/m)~0.85 or when q=1 abs(a/m)<0.98),
    fitting coefficients for QNM type, and sampling rate
    Returns the frequency, the Phenom amplitude of the inspiral-merger-ringdown
    Uses methods found in <https://arxiv.org/abs/1508.07253> and <https://arxiv.org/abs/1508.07250>

    Parameters
    ----------
    source : object
            source object from StrainandNoise, contains all source parameters
    pct_of_peak : float, optional
            the percentange of the strain at merger that dictates the maximum frequency the waveform is calculated at in geometrized units (G=c=1)

    Returns
    -------
    Mf : numpy array of floats
        the waveform frequencies in geometrized units (G=c=1)
    fullwaveform : numpy array of floats
        the waveform strain in geometrized units (G=c=1)

    """
    f_low = source.f_low
    N = source.nfreqs
    q = source.q
    x1 = source.chi1
    x2 = source.chi2
    fitcoeffs = source._fitcoeffs

    #M = m1+m2 #Total Mass
    #q = m2/m1 #Mass Ratio: Paper tested up to 18
    #eta = m1*m2/M**2 reduced mass: Paper tested up to 0.05 (q=18)
    eta = q/(q+1)**2
    x_PN = chi_PN(eta,x1,x2) #PN reduced spin parameter
    a_f = a_final(x1,x2,q,eta) #dimensionless spin

    ##################
    #Finds f_ringdown and f_damp from fit taken from <https://arxiv.org/abs/gr-qc/0512160>
    n = 0      #QNM indices
    l = 2
    m = 2
    numn = 3   #number of n's included in the table

    index = (l-2)*(2*l+1)*numn + (l-m)*numn + n
    f_fit = fitcoeffs[index][3:6]
    q_fit = fitcoeffs[index][6:9]

    omega_RD = f_fit[0]+f_fit[1]*(1-a_f)**f_fit[2]          #M omega_{lmn}
    tau = 2*(q_fit[0]+q_fit[1]*(1-a_f)**q_fit[2])/omega_RD  #tau_{lmn}/M = 2 Q_{lmn}/(M omega_{lmn})
    ########################
    f_RD = omega_RD/2/np.pi
    f_damp = 1/tau/2/np.pi

    Gamma1 = Lambda(eta,x_PN,4)
    Gamma2 = Lambda(eta,x_PN,5)
    Gamma3 = Lambda(eta,x_PN,6)

    f_peak = Calc_f_peak(f_RD,f_damp,[Gamma1,Gamma2,Gamma3])

    f1 = 0.014
    f3 = f_peak
    f2 = (f1+f3)/2

    cutoffFreq = Find_Cutoff_Freq(f_RD,f_damp,[Gamma1,Gamma2,Gamma3],pct_of_peak=pct_of_peak)

    #If lowest frequency is greater than cutoffFreq, then raise error.
    if f_low >= cutoffFreq:
        raise ValueError('Lower frequency bound (ie. f_low) must be lower than that of the merger ringdown.')
    
    Mf = np.logspace(np.log10(f_low),np.log10(cutoffFreq),N)

    v1 = A_insp(f1,eta,x1,x2,x_PN)
    v2 = Lambda(eta,x_PN,3)
    v3 = A_MR(f3,f_RD,f_damp,[Gamma1,Gamma2,Gamma3])
    fund1 = DA_insp(f1,eta,x1,x2,x_PN)
    fund3 = DA_MR(f3,f_RD,f_damp,[Gamma1,Gamma2,Gamma3])

    #############################
    #Calculate Solutions to eqn 21 in intermediate region
    Del_solns = A_intermediate(f1,f2,f3,v1,v2,v3,fund1,fund3) # Solutions to eqn 21

    ##############################
    #Calculate all sections of waveform and Paste together
    indxf1 = np.argmin(np.abs(Mf-f1))
    indxfpeak = np.argmin(np.abs(Mf-f_peak))

    tmpinspiral = A_norm(Mf[0:indxf1+1],eta)*A_insp(Mf[0:indxf1+1],eta,x1,x2,x_PN)
    tmpintermediate = A_norm(Mf[indxf1+1:indxfpeak],eta)*A_int(Mf[indxf1+1:indxfpeak],Del_solns)
    tmpmergerringdown = A_norm(Mf[indxfpeak:],eta)*A_MR(Mf[indxfpeak:],f_RD,f_damp,[Gamma1,Gamma2,Gamma3])
    fullwaveform = np.hstack((tmpinspiral,tmpintermediate,tmpmergerringdown))

    return [Mf,fullwaveform]

def A_norm(freqs,eta):
    """Calculates the constant scaling factor A_0

    Parameters
    ----------
    freqs : array
        The frequencies in Natural units (Mf, G=c=1) of the waveform
    eta : float
        The reduced mass ratio

    """
    const = np.sqrt(2*eta/3/np.pi**(1/3))
    return const*freqs**-(7/6)


def A_insp(freqs,eta,x1,x2,X_PN):
    """Calculates the Inspiral Amplitude

    Parameters
    ----------
    freqs : array
        The frequencies in Natural units (Mf, G=c=1) of the waveform
    eta : float
        The reduced mass ratio
    x1 : float
        The dimensionless spin parameter abs(a/m) for black hole m1.
    x2 : float
        The dimensionless spin parameter abs(a/m) for black hole m2.
    x_PN : float
        The PN reduced spin parameter

    """
    A_PN = 0.0
    A_higher = 0.0
    for i in range(7):
        A_PN = A_PN + PN_coeffs(eta,x1,x2,i)*(np.pi*freqs)**(i/3)
        if i >= 1 and i <= 3:
            A_higher = A_higher + Lambda(eta,X_PN,i-1)*freqs**((6+i)/3)
    return (A_PN + A_higher)


def DA_insp(freqs,eta,x1,x2,X_PN):
    """Calculates Derivative of the inspiral amplitude.

    Parameters
    ----------
    freqs : array
        The frequencies in Natural units (Mf, G=c=1) of the waveform
    eta : float
        The reduced mass ratio
    x1 : float
        The dimensionless spin parameter abs(a/m) for black hole m1.
    x2 : float
        The dimensionless spin parameter abs(a/m) for black hole m2.
    x_PN : float
        The PN reduced spin parameter

    """
    DA_PN = 0.0
    DA_higher = 0.0
    for i in range(7):
        PN_const = np.pi**(i/3)*(i/3)*PN_coeffs(eta,x1,x2,i)
        DA_PN = DA_PN + PN_const*(freqs)**((i-3)/3)
        if i >= 1 and i <= 3:
            higher_const = ((6+i)/3)*Lambda(eta,X_PN,i-1)
            DA_higher = DA_higher + higher_const*freqs**((i+3)/3)
        
    return DA_PN + DA_higher


def A_MR(freqs,f_RD,f_damp,Gammas):
    """Calculates the Normalized Merger-Ringdown Amplitude
    
    Parameters
    ----------
    freqs : array
        The frequencies in Natural units (Mf, G=c=1) of the waveform
    f_RD : float
        Frequency of the Ringdown transition
    f_damp : float
        Damping frequency
    Gammas : array-like
        Normalizes lorentzian to correct shape

    """
    varf = freqs-f_RD
    fg_d = Gammas[2]*f_damp 
    return (Gammas[0]*fg_d)/(varf**2+fg_d**2)*np.exp(-(Gammas[1]/fg_d)*varf)


def DA_MR(freqs,f_RD,f_damp,Gammas):
    """Calculates Derivative of the Merger-Ringdown Amplitude
    
    Parameters
    ----------
    freqs : array
        The frequencies in Natural units (Mf, G=c=1) of the waveform
    f_RD : float
        Frequency of the Ringdown transition
    f_damp : float
        Damping frequency
    Gammas : array-like
        Normalizes lorentzian to correct shape

    """
    varf = freqs-f_RD
    fg_d = Gammas[2]*f_damp
    A_MR_0 = A_MR(freqs,f_RD,f_damp,Gammas)
    return -A_MR_0*(2*varf/(varf**2+fg_d**2)+Gammas[1]/fg_d)


def A_intermediate(f1,f2,f3,v1,v2,v3,d1,d3):
    """Solves system of equations for intermediate amplitude matching"""
    Mat = np.array([[1., f1, f1**2, f1**3, f1**4],[1., f2, f2**2, f2**3, f2**4],[1., f3, f3**2, f3**3, f3**4], \
            [0., 1., 2*f1, 3*f1**2, 4*f1**3],[0., 1., 2*f3, 3*f3**2, 4*f3**3]],dtype='float')
    a = np.array([v1,v2,v3,d1,d3],dtype='float')
    return np.linalg.solve(Mat,a)


def A_int(freqs,delt):
    """Calculates the Intermediate Amplitude

    Parameters
    ----------
    freqs : array
        The frequencies in Natural units (Mf, G=c=1) of the waveform
    delt : array
        Coefficient solutions to match the inspiral to the merger-ringdown portion of the waveform

    """
    return (delt[0]+delt[1]*freqs+delt[2]*freqs**2+delt[3]*freqs**3+delt[4]*freqs**4)


def Lambda(eta,x_PN,lmbda):
    """Gets the Lambdas from Eqn 31 in <https://arxiv.org/abs/1508.07253>

    Parameters
    ----------
    eta : float
        The reduced mass ratio
    x_PN : float
        The PN reduced spin parameter
    lmbda : int
        Iterator for different Lambda variables using the zeta function

    """
    xi = x_PN-1
    xi2 = xi*xi
    xi3 = xi2*xi
    eta2 = eta*eta
    if lmbda == 0:  #rho1
        coeffs = zeta(0)
    elif lmbda == 1: #rho2
        coeffs = zeta(1)
    elif lmbda == 2: #rho3
        coeffs = zeta(2)
    elif lmbda == 3: #v2
        coeffs = zeta(3)
    elif lmbda == 4: #gamma1
        coeffs = zeta(4)
    elif lmbda == 5: #gamma2
        coeffs = zeta(5)
    elif lmbda == 6: #gamma3
        coeffs = zeta(6)
    
    return coeffs[0] + coeffs[1]*eta + \
        (coeffs[2] + coeffs[3]*eta + coeffs[4]*eta2)*xi + \
        (coeffs[5] + coeffs[6]*eta + coeffs[7]*eta2)*xi2 + \
        (coeffs[8] + coeffs[9]*eta + coeffs[10]*eta2)*xi3


def zeta(k):
    """Coefficients in table 5 of <https://arxiv.org/abs/1508.07253>"""
    if k == 0: #rho 1
        coeffs = [3931.9, -17395.8, 3132.38, 343966.0, -1.21626e6, -70698.0, 1.38391e6, -3.96628e6, -60017.5, 803515.0, -2.09171e6]
    elif k == 1: #rho 2
        coeffs = [-40105.5, 112253.0, 23561.7, -3.47618e6, 1.13759e7, 754313.0, -1.30848e7, 3.64446e7, 596227.0, -7.42779e6, 1.8929e7]
    elif k == 2: #rho 3
        coeffs = [83208.4, -191238.0, -210916.0, 8.71798e6, -2.69149e7, -1.98898e6, 3.0888e7, -8.39087e7, -1.4535e6, 1.70635e7, -4.27487e7]
    elif k == 3: #v 2
        coeffs = [0.814984, 2.57476, 1.16102, -2.36278, 6.77104, 0.757078, -2.72569, 7.11404, 0.176693, -0.797869, 2.11624]
    elif k == 4: #gamma 1
        coeffs = [0.0069274, 0.0302047, 0.00630802, -0.120741, 0.262716, 0.00341518, -0.107793, 0.27099, 0.000737419, -0.0274962, 0.0733151]
    elif k == 5: #gamma 2
        coeffs = [1.01034, 0.000899312, 0.283949, -4.04975, 13.2078, 0.103963, -7.02506, 24.7849, 0.030932, -2.6924, 9.60937]
    elif k == 6: #gamma 3
        coeffs = [1.30816, -0.00553773, -0.0678292, -0.668983, 3.40315, -0.0529658, -0.992379, 4.82068, -0.00613414, -0.384293, 1.75618]
    return coeffs


def PN_coeffs(eta,x1,x2,i):
    """Gets the PN Amplitude coefficients

    Parameters
    ----------
    eta : float
        The reduced mass ratio
    x1 : float
        The dimensionless spin parameter abs(a/m) for black hole m1.
    x2 : float
        The dimensionless spin parameter abs(a/m) for black hole m2.
    q : float
        The mass ratio m1/m2, m1<=m2
    i : int
        iterator to dictate which PN Amplitude to use

    Notes
    -----
    Coefficients in appendix B (eqns B14-B20) of <https://arxiv.org/abs/1508.07253>

    """
    delta = np.sqrt(1.0-4.0*eta)
    chi_s = (x1+x2)/2.0
    chi_a = (x1-x2)/2.0
    if i == 0:
        A_i = 1
    elif i == 1:
        A_i = 0
    elif i == 2:
        A_i = -323/224 + (451/168)*eta
    elif i == 3:
        A_i = (27/8)*delta*chi_a + (27/8 -(11/6)*eta)*chi_s
    elif i == 4:
        A_i = -27312085/8128512 -(1975055/338688)*eta + (105271/24192)*eta**2 + \
            (-81/32+8*eta)*chi_a**2 - 81/16*delta*chi_a*chi_s + (-81/32+17/8*eta)*chi_s**2
    elif i == 5:
        A_i = -85*np.pi/64 + 85*np.pi/16*eta + (285197/16128-1579/4032*eta)*delta*chi_a + \
            (285197/16128 - 15317/672*eta - 2227/1008*eta**2)*chi_s
    elif i == 6:
        A_i = -177520268561/8583708672 + (545384828789/5007163392 - 205*np.pi**2/48)*eta - \
            3248849057/178827264*eta**2 + 34473079/6386688*eta**3 + (1614569/64512-1873643/16128*eta+2167/42*eta**2)*chi_a**2 + \
            (31*np.pi/12 - 7*np.pi/3*eta)*chi_s + (1614569/64512-61391/1344*eta+57451/4032*eta**2)*chi_s**2 + \
            delta*chi_a*(31*np.pi/12+(1614569/32256-165961/2688*eta)*chi_s)
    return A_i


def Calc_f_peak(f_RD,f_damp,Gammas):
    """Calculates the frequency at the peak of the merger

    Parameters
    ----------
    f_RD : float
        Frequency of the Ringdown transition
    f_damp : float
        Damping frequency
    Gammas : array-like
        Normalizes lorentzian to correct shape
    
    Notes
    -----
    There is a problem with this expression from the paper becoming imaginary if gamma2 >= 1 
    so if gamma2 >= 1 then set the square root term to zero.

    """
    if Gammas[1] <= 1:
        f_max = np.abs(f_RD+f_damp*Gammas[2]*(np.sqrt(1-Gammas[1]**2)-1)/Gammas[1])
    else:
        f_max = np.abs(f_RD+(f_damp*Gammas[2]*-1)/Gammas[1])
    return f_max


def Find_Cutoff_Freq(f_RD,f_damp,Gammas,pct_of_peak=0.0001):
    """Cutoff signal when the amplitude is a factor of 10 below the value at f_RD

    Parameters
    ----------
    f_RD : float
        Frequency of the Ringdown transition
    f_damp : float
        Damping frequency
    Gammas : array-like
        Normalizes lorentzian to correct shape

    pct_of_peak : float, optional
        the percentange of the strain at merger that dictates the maximum 
        frequency the waveform is calculated at in geometrized units (G=c=1) 

    """
    tempfreqs = np.logspace(np.log10(f_RD),np.log10(10*f_RD),100)
    cutoffAmp = pct_of_peak*A_MR(f_RD,f_RD,f_damp,[Gammas[0],Gammas[1],Gammas[2]])
    merger_ringdown_Amp = A_MR(tempfreqs,f_RD,f_damp,[Gammas[0],Gammas[1],Gammas[2]])
    cutoffindex = np.argmin(np.abs(cutoffAmp-merger_ringdown_Amp))
    return tempfreqs[cutoffindex]

def a_final(x1,x2,q,eta):
    """The Final spin of the binary remnant black hole

    Parameters
    ----------
    x1 : float
        The dimensionless spin parameter abs(a/m) for black hole m1.
    x2 : float
        The dimensionless spin parameter abs(a/m) for black hole m2.
    q : float
        The mass ratio m1/m2, m1<=m2
    eta : float
        The reduced mass ratio

    Notes
    -----
    Uses eq. 3 in <https://arxiv.org/abs/0904.2577>, changed to match our q convention
    a=J/M**2 where J = x1*m1**2 + x2*m2**2

    """
    a = (q**2*x1+x2)/(q**2+1)
    s4 = -0.1229
    s5 = 0.4537
    t0 = -2.8904
    t2 = -3.5171
    t3 = 2.5763
    return a + s4*a**2*eta + s5*a*eta**2 + t0*a*eta + 2*np.sqrt(3)*eta + t2*eta**2 + t3*eta**3


def chi_PN(eta,x1,x2):
    """Calculates the PN reduced spin parameter

    Parameters
    ----------
    eta : float
        The reduced mass ratio
    x1 : float
        The dimensionless spin parameter abs(a/m) for black hole m1.
    x2 : float
        The dimensionless spin parameter abs(a/m) for black hole m2.

    Notes
    -----
    See Eq 5.9 in <https://arxiv.org/abs/1107.1267v2>

    """
    delta = np.sqrt(1.0-4.0*eta)
    chi_s = (x1+x2)/2.0
    chi_a = (x1-x2)/2.0
    return chi_s*(1.0-eta*76.0/113.0) + delta*chi_a