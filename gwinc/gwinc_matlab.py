import os
import copy
import tempfile
import logging

import scipy.io

from . import const
from . import suspension
from .struct import Struct

##################################################

MATLAB_ENGINE = None

class Matlab:
    def __init__(self, gwincpath=None):
        """Start a MATLAB engine for GWINC processing

         engine is provided, the GWINC path is added
        to it's path.

        """
        global MATLAB_ENGINE

        if MATLAB_ENGINE:
            return

        import matlab.engine

        if not gwincpath:
            gwincpath = os.path.expanduser(os.getenv('GWINCPATH', 'gwinc'))
        if not os.path.exists(os.path.join(gwincpath, 'gwinc.m')):
            raise IOError("Invalid MATLAB GWINC path: '{}'".format(gwincpath))

        logging.info("starting MATLAB engine...")
        MATLAB_ENGINE = matlab.engine.start_matlab()

        logging.info("gwinc path: {}".format(gwincpath))
        MATLAB_ENGINE.addpath(gwincpath)


    @property
    def eng(self):
        return MATLAB_ENGINE

    @property
    def workspace(self):
        return MATLAB_ENGINE.workspace

    def addpath(self, *args):
        return MATLAB_ENGINE.addpath(*args)

    def eval(self, *args, **kwargs):
        return MATLAB_ENGINE.eval(*args, **kwargs)


    def load_array(self, var, array):
        """Load numpy array into workspace as vector.

        `var` is name of workspace variable, and `array` is numpy
        ndarray.

        """
        # this stupidity because you can't just give the engine a np.ndarray
        MATLAB_ENGINE.workspace[var] = array.tolist()
        MATLAB_ENGINE.eval('{0} = cell2mat({0});'.format(var), nargout=0)


    def load_struct(self, var, struct):
        """Load pygwinc.Struct array into workspace as vector.

        `var` is name of workspace variable, and `struct` is
        pygwinc.Struct.

        """
        # similar stupidity prevents this (this time recarrays in the dict):
        #matlab.workspace['ifo'] = ifo.to_dict(array=True)
        with tempfile.NamedTemporaryFile(suffix='.mat') as f:
            scipy.io.savemat(f, struct.to_dict(array=True))
            MATLAB_ENGINE.eval("{} = load('{}');".format(var, f.name), nargout=0)


    def extract(self, *wvars):
        """Extract workspace variables from engine.

        Returns dict with wvars as keys.

        """
        assert len(wvars) > 0
        with tempfile.NamedTemporaryFile(suffix='.mat') as f:
            MATLAB_ENGINE.save(f.name, *wvars, nargout=0)
            data = scipy.io.loadmat(f, squeeze_me=True, struct_as_record=False)
        if len(wvars) == 1:
            return data[wvars[0]]
        else:
            return data

##################################################

NOISE_NAME_MAP = {
    'Quantum': 'Quantum Vacuum',
    'Newtonian': 'Newtonian Gravity',
    'CoatBrown': 'Coating Brownian',
    'CoatTO': 'Coating Thermo-Optic',
    'SubBrown': 'Substrate Brownian',
    'SubTE': 'Substrate Thermo-Elastic',
    'SuspThermal': 'Suspension Thermal',
    'ResGas': 'Excess Gas',
}


def ifo_matlab_transform(ifo):
    """Prep the ifo structure for use with MATLAB gwinc

    * add "Constants" sub-Struct
    * copy Temp to Constants
    * change Suspension.FiberType string into number

    """
    # add constants
    CONSTS = {k:v for k, v in const.__dict__ if not k.startswith('__')}
    ifo.Constants = Struct.from_dict(CONSTS)

    # copy tempurature into Constants
    ifo.Constants.Temp = ifo.Infrastructure.Temp

    # translate FiberType string to int
    ifo.Suspension.FiberType = suspension.FIBER_TYPES.index(ifo.Suspension.FiberType)

    return ifo


def _rename_noises(d):
    nd = {}
    for k,v in d.items():
        try:
            nk = NOISE_NAME_MAP[k]
        except KeyError:
            nk = k
        if isinstance(v, dict):
            nd[nk] = _rename_noises(v)
        else:
            nd[nk] = v
    return nd


def gwinc_matlab(f, ifoin, plot=False):
    """Execute gwinc in MATLAB with the specified ifo model.

    This uses the python matlab.engine (see Matlab class) to calculate
    noises with MATLAB gwinc (gwinc directory specified by GWINCPATH
    environment variable).

    Returns `score` (dict), `noises` (dict), and `ifo` (Struct) as
    returned from MATLAB.

    If `plot` is True will cause MATLAB to produce it's own plot for
    the noise budget.

    """
    ifo = copy.deepcopy(ifoin)

    ifo_matlab_transform(ifo)

    matlab = Matlab()

    matlab.load_array('f', f)
    matlab.load_struct('ifo', ifo)

    if plot:
        plot_flag = '3'
    else:
        plot_flag = '0'

    cmd = "[score, noises, ifo] = gwinc(f, [], ifo, SourceModel, {});".format(plot_flag)
    matlab.eval(cmd, nargout=0)

    data = matlab.extract('score', 'noises', 'ifo')

    score = data['score']
    mnoises = Struct.from_matstruct(data['noises']).to_dict()
    ##### blow out 'MirrorThermal' sub-dict
    for n,d in mnoises['MirrorThermal'].items():
        if n == 'Total':
            continue
        mnoises[n] = d
    del mnoises['MirrorThermal']
    #####
    noises = _rename_noises(mnoises)
    ifo = Struct.from_matstruct(data['ifo'])

    return score, noises, ifo
