import os
import importlib


def lpath(file0, file1):
    """Return path of file1 when expressed relative to file0.

    For instance, if file0 is "/path/to/file0" and file1 is
    "../for/file1" then what is returned is "/path/for/file1".

    This is useful for resolving paths within packages with e.g.:

      rpath = lpath(__file__, '../resource.yaml')

    """
    return os.path.abspath(os.path.join(os.path.dirname(file0), file1))


def load_module(name_or_path):
    """Load module from name or path.

    Return loaded module and module path.

    """
    if os.path.exists(name_or_path):
        path = name_or_path.rstrip('/')
        modname = os.path.splitext(os.path.basename(path))[0]
        if os.path.isdir(path):
            path = os.path.join(path, '__init__.py')
        spec = importlib.util.spec_from_file_location(modname, path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
    else:
        mod = importlib.import_module(name_or_path)
    try:
        path = mod.__path__[0]
    except AttributeError:
        path = mod.__file__
    return mod, path
