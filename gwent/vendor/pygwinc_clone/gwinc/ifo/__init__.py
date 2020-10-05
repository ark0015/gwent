import os
import logging

from ..struct import load_struct, STRUCT_EXT
from ..util import load_module


PLOT_STYLE = dict(
    ylabel=u"Strain [1/\u221AHz]",
)


def available_ifos():
    """List available included pre-defined IFOs

    """
    ifos = []
    root = os.path.dirname(__file__)
    for f in os.listdir(root):
        if os.path.isdir(os.path.join(root, f)) and f[0] != '_':
            ifos.append(f)
    return sorted(ifos)


def load_ifo(name_or_path):
    """Load GWINC IFO Budget by name or from file.

    Named IFOs should correspond to one of the IFOs available in the
    gwinc package (see gwinc.available_ifos()).  If a path is provided
    it should either be a budget package (directory) or module (ending
    in .py), or an IFO struct (see gwinc.load_struct()).  In the
    latter case the base aLIGO budget definition will be used.

    Returns primary Budget class, ifo structure, frequency array, and
    plot style dictionary, with the last three being None if they are
    not defined in the budget.

    """
    ifo = None

    if os.path.exists(name_or_path):
        path = name_or_path.rstrip('/')
        bname, ext = os.path.splitext(os.path.basename(path))

        if ext in STRUCT_EXT:
            logging.info("loading struct {}...".format(path))
            ifo = load_struct(path)
            bname = 'aLIGO'
            modname = 'gwinc.ifo.aLIGO'
            logging.info("loading budget {}...".format(modname))

        else:
            modname = path
            logging.info("loading module path {}...".format(modname))

    else:
        if name_or_path not in available_ifos():
            raise RuntimeError("Unknonw IFO '{}' (available IFOs: {}).".format(
                name_or_path,
                available_ifos(),
            ))
        bname = name_or_path
        modname = 'gwinc.ifo.'+name_or_path
        logging.info("loading module {}...".format(modname))

    mod, modpath = load_module(modname)

    Budget = getattr(mod, bname)
    ifo = getattr(mod, 'IFO', ifo)
    ifopath = os.path.join(modpath, 'ifo.yaml')
    if not ifo and ifopath:
        ifo = load_struct(ifopath)
    freq = getattr(mod, 'FREQ', None)
    plot_style = getattr(mod, 'PLOT_STYLE', PLOT_STYLE)

    return Budget, ifo, freq, plot_style
