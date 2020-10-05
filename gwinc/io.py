import datetime
import h5py


SCHEMA = 'GWINC noise budget'
SCHEMA_VERSION = 1


def _write_trace_recursive(grp, traces):
    for name, trace in traces.items():
        if isinstance(trace, dict):
            tgrp = grp.create_group(name)
            _write_trace_recursive(tgrp, trace)
        else:
            data, style = trace
            dset = grp.create_dataset(name, data=data)
            for key, val in style.items():
                dset.attrs[key] = val


def save_hdf5(path, freq, traces, **kwargs):
    """Save GWINC traces dict to HDF5 file.

    See HDF5_SCHEMA.

    """
    with h5py.File(path, 'w') as f:
        f.attrs['SCHEMA'] = SCHEMA
        f.attrs['SCHEMA_VERSION'] = SCHEMA_VERSION
        # FIXME: add budget code hash or something
        f.attrs['date'] = datetime.datetime.now().isoformat()
        for key, val in kwargs.items():
            f.attrs[key] = val
        f.create_dataset('Freq', data=freq)
        tgrp = f.create_group('traces')
        _write_trace_recursive(tgrp, traces)


def _read_trace_recursive(element):
    trace = {}
    for name, item in element.items():
        if isinstance(item, h5py.Group):
            trace[name] = _read_trace_recursive(item)
        else:
            trace[name] = item.value, dict(item.attrs.items())
    return trace


def load_hdf5(path):
    """Load traces from HDF5 file.

    Returns a recursive traces dictionary.  See HDF5_SCHEMA.

    """
    with h5py.File(path, 'r') as f:
        # FIXME: check SCHEMA name/version
        freq = f['Freq'].value
        traces = _read_trace_recursive(f['/traces'])
        attrs = dict(f.attrs.items())
        return freq, traces, attrs
