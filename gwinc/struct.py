import os
import re
import io
import yaml
import numpy as np
from scipy.io import loadmat
from scipy.io.matlab.mio5_params import mat_struct


# HACK: fix loading number in scientific notation
#
# https://stackoverflow.com/questions/30458977/yaml-loads-5e-6-as-string-and-not-a-number
#
# An apparent bug in python-yaml prevents it from regognizing
# scientific notation as a float.  The following is a modified version
# of the parser that recognize scientific notation appropriately.
yaml_loader = yaml.SafeLoader
yaml_loader.add_implicit_resolver(
    'tag:yaml.org,2002:float',
    re.compile('''^(?:
     [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
    |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
    |\\.[0-9_]+(?:[eE][-+][0-9]+)?
    |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
    |[-+]?\\.(?:inf|Inf|INF)
    |\\.(?:nan|NaN|NAN))$''', re.X),
    list('-+0123456789.'))


def dictlist2recarray(l):
    def dtype(v):
        if isinstance(v, int):
            return float
        else:
            return type(v)
    # get dtypes from first element dict
    dtypes = [(k, dtype(v)) for k,v in l[0].items()]
    values = [tuple(el.values()) for el in l]
    out = np.array(values, dtype=dtypes)
    return out.view(np.recarray)


class Struct(object):
    """Matlab struct-like object

    This is a simple implementation of a MATLAB struct-like object
    that stores values as attributes of a simple class: and allows
    assigning to attributes recursively, e.g.:

    >>> s = Struct()
    >>> s.a = 4
    >>> s.b = Struct()
    >>> s.b.c = 8

    Various classmethods allow creating one of these objects from YAML
    file, a nested dict, or a MATLAB struct object.

    """

    # FIXME: This would be a way to allow setting nested struct
    # attributes, e.g.:
    #
    # >>> s = Struct()
    # >>> s.a.b.c = 4
    #
    # Usage of __getattr__ like this is dangerous and creates
    # non-intuitive behavior (i.e. an empty struct is returned when
    # accessing attributes that don't exist).  Is there a way to
    # accomplish this without that adverse side affect?
    #
    # def __getattr__(self, name):
    #     if name not in self.__dict__:
    #         self.__dict__[name] = Struct()
    #     return self.__dict__[name]

    ##########

    def __init__(self, **kwargs):
        """Arguments can pre-fill the structure

        """
        self.__dict__.update(kwargs)

    def __getitem__(self, key):
        """Get a (possibly nested) value from the struct.

        """
        if '.' in key:
            k, r = key.split('.', 1)
            # FIXME: this is inelegant.  better done with regexp?
            if len(k.split('[')) > 1:
                kl, i = k.split('[')
                i = int(i.strip(']'))
                return self.__dict__[kl][i][r]
            return self.__dict__[k][r]
        else:
            return self.__dict__[key]

    def get(self, key, default):
        """Get a (possibly nested) value from the struct, or default.

        """
        try:
            return self[key]
        except KeyError:
            return default

    def __setitem__(self, key, value):
        if '.' in key:
            k, r = key.split('.', 1)
            self.__dict__[k][r] = value
        else:
            self.__dict__[key] = value

    def setdefault(self, key, default):
        return self.__dict__.setdefault(key, default)

    def items(self):
        return self.__dict__.items()

    def keys(self):
        return self.__dict__.keys()

    def values(self):
        return self.__dict__.values()

    def __contains__(self, key):
        return key in self.__dict__


    def to_dict(self, array=False):
        """Return nested dictionary representation of Struct.

        If `array` is True any lists encountered will be turned into
        numpy arrays, and lists of Structs will be turned into record
        arrays.  This is needed to convert to structure arrays in
        matlab.

        """
        d = {}
        for k,v in self.__dict__.items():
            if isinstance(v, type(self)):
                d[k] = v.to_dict(array=array)
            else:
                if isinstance(v, list):
                    try:
                        # this should fail if the elements of v are
                        # not Struct
                        # FIXME: need cleaner way to do this
                        v = [i.to_dict(array=array) for i in v]
                        if array:
                            v = dictlist2recarray(v)
                    except AttributeError:
                        if array:
                            v = np.array(v)
                elif isinstance(v, int):
                    v = float(v)
                d[k] = v
        return d

    def to_yaml(self, path=None):
        """Return YAML representation of Struct.

        Write YAML to `path` if specified.

        """
        y = yaml.dump(self.to_dict(), default_flow_style=False)
        if path:
            with open(path, 'w') as f:
                f.write(y)
        else:
            return y

    # def __repr__(self):
    #     return self.to_yaml().strip('\n')

    def __str__(self):
        return '<GWINC Struct: {}>'.format(list(self.__dict__.keys()))

    def __iter__(self):
        return iter(self.__dict__)

    def walk(self):
        """Iterate over all leaves in the struct tree.

        """
        for k,v in self.__dict__.items():
            if isinstance(v, type(self)):
                for sk,sv in v.walk():
                    yield k+'.'+sk, sv
            else:
                try:
                    for i,vv in enumerate(v):
                        for sk,sv in vv.walk():
                            yield '{}[{}].{}'.format(k,i,sk), sv
                except (AttributeError, TypeError):
                    yield k, v


    def diff(self, other):
        """Return tuple of differences between target IFO.

        Returns list of (key, value, other_value) tuples.  Value is
        None if key not present.

        """
        diffs = []
        for k, ov in other.walk():
            v = self.get(k, None)
            if ov != v and ov is not v:
                diffs.append((k, v, ov))
        for k, v in self.walk():
            ov = other.get(k, None)
            if ov is None:
                diffs.append((k, v, ov))
        return diffs


    def to_txt(self, path=None, fmt='0.6e', delimiter=': ', end=''):
        """Return text represenation of Struct, one element per line.

        Struct keys use '.' to indicate hierarchy.  The `fmt` keyword
        controls the formatting of numeric values.  MATLAB code can be
        generated with the following parameters:

        >>> ifo.to_txt(delimiter=' = ', end=';')

        Write text to `path` if specified.

        """
        txt = io.StringIO()
        for k, v in sorted(self.walk()):
            if isinstance(v, (int, float, complex)):
                base = fmt
            elif isinstance(v, (list, np.ndarray)):
                if isinstance(v, list):
                    v = np.array(v)
                v = np.array2string(v, separator='', max_line_width=np.Inf, formatter={'all':lambda x: "{:0.6e} ".format(x)})
                base = 's'
            else:
                base = 's'
            txt.write(u'{key}{delimiter}{value:{base}}{end}\n'.format(
                key=k, value=v, base=base,
                delimiter=delimiter,
                end=end,
            ))
        if path:
            with open(path, 'w') as f:
                f.write(txt.getvalue())
        else:
            return txt.getvalue()


    @classmethod
    def from_dict(cls, d):
        """Create Struct from nested dict.

        """
        c = cls()
        for k,v in d.items():
            if type(v) == dict:
                c.__dict__[k] = Struct.from_dict(v)
            else:
                try:
                    c.__dict__[k] = list(map(Struct.from_dict, v))
                except (AttributeError, TypeError):
                    c.__dict__[k] = v
        return c


    @classmethod
    def from_yaml(cls, y):
        """Create Struct from YAML string.

        """
        d = yaml.load(y)
        return cls.from_dict(d)


    @classmethod
    def from_matstruct(cls, s):
        """Create Struct from scipy.io.matlab mat_struct object.

        """
        c = cls()
        try:
            s = s['ifo']
        except:
            pass
        for k,v in s.__dict__.items():
            if k in ['_fieldnames']:
                # skip these fields
                pass
            elif type(v) is mat_struct:
                c.__dict__[k] = Struct.from_matstruct(v)
            else:
                # handle lists of Structs
                try:
                    c.__dict__[k] = list(map(Struct.from_matstruct, v))
                except:
                    c.__dict__[k] = v
                    # try:
                    #     c.__dict__[k] = float(v)
                    # except:
                    #     c.__dict__[k] = v
        return c


    @classmethod
    def from_file(cls, path):
        """Load Struct from .yaml or MATLAB .mat file.

        File type will be determined by extension.

        """
        (root, ext) = os.path.splitext(path)

        with open(path, 'r') as f:
            if ext in ['.yaml', '.yml']:
                d = yaml.load(f, Loader=yaml_loader)
                return cls.from_dict(d)
            elif ext == '.mat':
                s = loadmat(f, squeeze_me=True, struct_as_record=False)
                return cls.from_matstruct(s)
            else:
                raise IOError("Unknown file type: {}".format(ext))


def load_struct(path):
    """Load struct from YAML or MATLAB file.

    Files may be either .yaml, .mat or .m.  For .m files, the file is
    expected to include either an object or function that corresponds
    to the basename of the file.  The MATLAB engine will be invoked to
    execute the .m code and extract the resultant IFO data.

    """
    root, ext = os.path.splitext(path)

    if ext == '.m':
        from ..gwinc_matlab import Matlab
        matlab = Matlab()
        matlab.addpath(os.path.dirname(path))
        func_name = os.path.basename(root)
        matlab.eval("ifo = {};".format(func_name), nargout=0)
        ifo = matlab.extract('ifo')
        return Struct.from_matstruct(ifo)

    else:
        return Struct.from_file(path)


# accepted extension types for struct files
STRUCT_EXT = ['.yaml', '.yml', '.mat', '.m']
