import astropy.units as u
import numpy as np


def make_quant(param, default_unit):
    """Convenience function to intialize a parameter as an astropy quantity.

    Parameters
    ----------
    param: float, or Astropy Quantity
        Parameter to initialize
    default_unit: str
        Astropy unit string, sets as default for param.

    Returns
    -------
        an astropy quantity

    Examples
    --------
        self.f0 = make_quant(f0,'MHz')

    Notes
    -----
    Taken from <https://github.com/Hazboun6/hasasia/blob/master/hasasia/sensitivity.py#L834>

    """
    default_unit = u.core.Unit(default_unit)
    if hasattr(param, "unit"):
        try:
            quantity = param.to(default_unit)
        except u.UnitConversionError:
            raise ValueError(
                "Quantity {0} with incompatible unit {1}".format(param, default_unit)
            )
    else:
        quantity = param * default_unit

    return quantity


def Get_Var_Dict(obj, value):
    """Updates and initializes variable dictionaries used to keep track of
    current values and variable minima and maxima.

    Parameters
    ----------
    obj: object
        Instance of class with parameter variables
    value: array-like
        value(s) that are assigned into dictionary

    Notes
    -----
    value contains the variable name in the first index
    the next is the current value of the variable
    the last two are optional and contain the variable min and max

    Examples
    --------
    ``obj.var_dict = ['M',value]``
        where obj is in this case an instance of a BinaryBlackHole

    """
    if not hasattr(obj, "var_dict"):
        obj._var_dict = {}
    if isinstance(value, list):
        if len(value) == 2 and isinstance(value[0], str):
            var_name = value[0]
            vals = value[1]
            if isinstance(vals, u.Quantity):
                no_unit_vals = vals.value
            else:
                no_unit_vals = vals

            if isinstance(no_unit_vals, list) and len(no_unit_vals) == 3:
                if hasattr(obj, "n_p"):
                    if obj.n_p == len(no_unit_vals):
                        raise ValueError(LenError_1())
                if (
                    isinstance(vals[0], (float, int, u.Quantity))
                    and isinstance(vals[1], (float, int, u.Quantity))
                    and isinstance(vals[2], (float, int, u.Quantity))
                ):
                    obj._return_value = vals[0]
                    obj._var_dict[var_name] = {
                        "val": vals[0],
                        "min": vals[1],
                        "max": vals[2],
                        "sampled": False,
                    }
                else:
                    raise ValueError(DictError_3())
            elif (
                isinstance(no_unit_vals, (list, np.ndarray)) and len(no_unit_vals) != 3
            ):
                if var_name in obj._var_dict.keys():
                    obj._var_dict[var_name]["val"] = vals
                else:
                    if len(no_unit_vals) == 2:
                        obj.var_dict[var_name] = {
                            "val": vals,
                            "min": None,
                            "max": None,
                            "sampled": True,
                        }
                    else:
                        if hasattr(obj, "n_p"):
                            if obj.n_p == len(no_unit_vals):
                                obj.var_dict[var_name] = {
                                    "val": vals,
                                    "min": None,
                                    "max": None,
                                    "sampled": False,
                                }
                            else:
                                raise ValueError(LenError_2())
                        else:
                            raise ValueError(LenError_2())

                obj._return_value = vals
            elif isinstance(no_unit_vals, (float, int, np.int64)):
                if var_name in obj._var_dict.keys():
                    obj._var_dict[var_name]["val"] = vals
                else:
                    obj.var_dict[var_name] = {
                        "val": vals,
                        "min": None,
                        "max": None,
                        "sampled": False,
                    }
                obj._return_value = vals
            else:
                raise ValueError(DictError_2())
        else:
            raise ValueError(DictError_Full())
    else:
        raise ValueError(DictError_Full())


def DictError_Full():
    return 'Must assign either: \n\
    - A name and value in a list (ie. ["name",val]), or \n\
    - A name, a value, a minimum value, and maximum value in a list (ie. ["name",val,min,max]), \n\
    where where name is a string, and val,min,and max are either floats, ints, or an astropy Quantity.'


def DictError_3():
    return 'Must assign a name, a value, a minimum value, and maximum value in a list (ie. ["name",val,min,max]), \n\
    where name is a string, and val, min, and max are either floats, ints, or astropy Quantities.'


def DictError_2():
    return 'Must assign a name and value in a list (ie. ["name",val]) \n\
    where name is a string, and val is either a float, an int, or an astropy Quantity.'


def LenError_1():
    return "Could not tell if values are pulsar values or min/maxes. Try using more than 3 pulsars."


def LenError_2():
    return "To assign an array of values, it must be either [min,max], or an array of individual pulsar parameters of length n_p."
