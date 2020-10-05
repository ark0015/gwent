"""Physical constants

"""
import scipy.constants
import scipy.special


c = scipy.constants.c
G = scipy.constants.G
g = scipy.constants.g
E0 = scipy.constants.epsilon_0
hbar = scipy.constants.hbar
kB = scipy.constants.k
yr = scipy.constants.year
AU = scipy.constants.astronomical_unit
parsec = scipy.constants.parsec
Mpc = scipy.constants.parsec * 1e6

# FIXME: use astropy for the following astrophysical/cosmological
# constants
R_earth = 6.3781e6
SolarMassParameter = 1.3271244e20
MSol = 1.3271244e20 / scipy.constants.G
# http://physics.nist.gov/cuu/Constants/
H0 = 67110
# http://arxiv.org/pdf/1303.5076v3.pdf
omegaM = 0.3175,
omegaLambda = 1 - 0.3175

# bessel function zeros
BESSEL_ZEROS = scipy.special.jn_zeros(1, 300)
J0M = scipy.special.jn(0, BESSEL_ZEROS)
