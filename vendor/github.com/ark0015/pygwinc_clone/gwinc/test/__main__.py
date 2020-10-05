import os
import sys
import signal
import pickle
import subprocess
import hashlib
import numpy as np
import matplotlib.pyplot as plt
import argparse

import logging
logging.basicConfig(format='%(message)s',
                    level=os.getenv('LOG_LEVEL', logging.INFO))

from .. import load_ifo
from ..gwinc import gwinc
from ..gwinc_matlab import gwinc_matlab
try:
    import inspiral_range
except ImportError:
    inspiral_range = None


FLO = 5
FHI = 6000
NPOINTS = 3000


def path_hash(path):
    """Calculate SHA1 hash of path, either directory or file"""
    if not path or not os.path.exists(path):
        return
    path = os.path.expanduser(path)
    if os.path.isdir(path):
        d = path
        f = '.'
    else:
        d = os.path.dirname(path)
        f = os.path.basename(path)
    CWD = os.getcwd()
    os.chdir(d)
    cmd = 'find {} -type f ! -wholename "*/.*" -print0 | sort -z | xargs -0 sha1sum | sha1sum'.format(f)
    sha1sum_out = subprocess.check_output(cmd, shell=True)
    sha1sum = sha1sum_out.split()[0]
    os.chdir(CWD)
    return sha1sum.decode()

##################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plot', '-p', action='store_true', help='plot differences')
    parser.add_argument('--save', '-s', help='save plot to file')
    parser.add_argument('--recalc', '-r', action='store_true', help='recalculate all traces')
    parser.add_argument('--tolerance', '-t', help='fractional tolerance', type=float, default=1e-6)
    parser.add_argument('--skip', '-k', action='append', help='traces to skip comparing')
    parser.add_argument('IFO', help='IFO name or description file')
    args = parser.parse_args()

    logging.info("loading IFO '{}'...".format(args.IFO))
    Budget, ifo, freq, plot_style = load_ifo(args.IFO)

    freq = np.logspace(np.log10(FLO), np.log10(FHI), NPOINTS)

    ##############################
    # matgwinc processing

    mdata_pkl = os.path.join(os.path.dirname(__file__), '{}.pkl'.format(args.IFO))

    ifo_hash = hashlib.sha1(ifo.to_txt().encode()).hexdigest()
    gwinc_hash = path_hash(os.getenv('GWINCPATH'))
    if not gwinc_hash:
        logging.warning("GWINCPATH not specified or does not exist; skipping check for changes to matgwinc code.")

    mrecalc = args.recalc

    if os.path.exists(mdata_pkl) and not mrecalc:
        mrecalc = False
        logging.info("loading matgwinc data {}...".format(mdata_pkl))
        with open(mdata_pkl, 'rb') as f:
            if sys.version_info.major > 2:
                mdata = pickle.load(f, encoding='latin1')
            else:
                mdata = pickle.load(f)

        if mdata['ifo_hash'] != ifo_hash:
            logging.info("ifo hash has changed: {}".format(ifo_hash))
            mrecalc = True
        if gwinc_hash and mdata['gwinc_hash'] != gwinc_hash:
            logging.info("matgwinc hash has changed: {}".format(gwinc_hash))
            mrecalc = True

    if mrecalc:
        logging.info("calculating matgwinc noises...")
        try:
            mscore, mnoises, mifo = gwinc_matlab(freq, ifo)
        except (ImportError, OSError):
            sys.exit("MATLAB engine not available.")
        mdata = dict(score=mscore, noises=mnoises, ifo=mifo, ifo_hash=ifo_hash, gwinc_hash=gwinc_hash)
        with open(mdata_pkl, 'wb') as f:
            pickle.dump(mdata, f)

    mnoises = mdata['noises']

    ##############################
    # pygwinc processing

    logging.info("calculating pygwinc noises...")
    score, noises, ifo = gwinc(freq, ifo)

    ##############################
    # calc inspiral ranges

    if inspiral_range:
        logging.info("calculating inspiral ranges...")

        range_func = inspiral_range.range
        H = inspiral_range.waveform.CBCWaveform(freq)

        mfom = range_func(freq, mnoises['Total'], H=H)
        _, mnoises['int73'] = inspiral_range.int73(freq, mnoises['Total'])
        logging.info("matgwinc range: {:.2f} Mpc".format(mfom))

        fom = range_func(freq, noises['Total'], H=H)
        _, noises['int73'] = inspiral_range.int73(freq, noises['Total'])
        logging.info("pygwinc range: {:.2f} Mpc".format(fom))

        fom_title = """inspiral {func} {m1}/{m2} Msol:
matgwinc: {mfom:.2f} Mpc
pygwinc: {fom:.2f} Mpc""".format(
            func=range_func.__name__,
            m1=H.params['m1'],
            m2=H.params['m2'],
            mfom=mfom,
            fom=fom,
        )

    else:
        fom_title = ''

    ##############################
    # find differences

    skip = args.skip
    fractional_tolerance = args.tolerance

    diffs = {}
    for name, noise in noises.items():
        if name in ['Freq']:
            continue
        if skip and name in skip:
            logging.warning("SKIPPING TEST: '{}'".format(name))
            continue

        try:
            mnoise = mnoises[name]
        except KeyError:
            continue
        # logging.info("compare: {}".format(name))

        mn = mnoise
        pn = noise
        # _, mn = inspiral_range.int73(freq, mnoise)
        # _, pn = inspiral_range.int73(freq, noise)

        diff = np.sqrt(mn) - np.sqrt(pn)
        frac = abs(diff / np.sqrt(pn))

        if max(frac) < fractional_tolerance:
            continue

        logging.warning("EXCESSIVE DIFFERENCE: {:{w}} {:6.1f} ppm".format(
            name, max(frac)*1e6, w=max([len(n) for n in noises])))
        # logging.warning("  max: {:e}, min: {:e}".format(max(frac), min(frac)))

        diffs[name] = (mn, pn, frac)

    ##############################
    # plot

    if args.plot:
        spec = (len(diffs)+1, 2)
        sharex = None
        for i, name in enumerate(diffs):
            mn, pn, frac = diffs[name]

            axl = plt.subplot2grid(spec, (i, 0), sharex=sharex)
            axl.loglog(freq, np.sqrt(pn), label='pygwinc')
            axl.loglog(freq, np.sqrt(mn), label='matgwinc')
            axl.grid()
            axl.legend(loc='upper right')
            axl.set_ylabel(name)
            if i == 0:
                sharex = axl

            axr = plt.subplot2grid(spec, (i, 1), sharex=sharex)
            axr.loglog(freq, frac)
            axr.grid()
            axr.axhline(y=max(frac), color='r', linestyle='--')
            axr.text(max(freq)+4000, max(frac), '{:.1f} ppm'.format(max(frac)*1e6),
                     horizontalalignment='left', verticalalignment='center',
                     color='red')

        if diffs:
            axl.set_xlabel("frequency [Hz]")
            axr.set_xlabel("frequency [Hz]")
    
            plt.suptitle("""{} mat/py gwinc noise comparison
noises that differ by more than {} ppm [(mat-py)/py]
{}""".format(args.IFO, fractional_tolerance*1e6, fom_title))

            if args.save:
                plt.gcf().set_size_inches(11, (len(diffs)+1)*4)
                plt.savefig(args.save)
            else:
                plt.show()

        else:
            logging.warning("All tests passed, so no plot was generated")

    ##############################

    if len(diffs) > 0:
        return 1
    return 0

##################################################

if __name__ == '__main__':
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    sys.exit(main())
