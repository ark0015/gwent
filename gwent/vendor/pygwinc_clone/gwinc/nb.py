import os
import logging
import itertools
import importlib
import importlib.util
import collections
import numpy as np
import scipy.interpolate


def quadsum(data):
    """Calculate quadrature sum of list of data arrays.

    Provided data are assumed to be power-referred, so this is a
    simple point-by-point sum.

    NaNs in sum elements do not contribute to sum.

    """
    return np.nansum(data, 0)



class BudgetItem:
    """GWINC BudgetItem class

    """
    def load(self):
        """Overload method for initial loading of static data.

        """
        return None

    def update(self, **kwargs):
        """Overload method for updating data needed to calculate final PSD.

        By default any keyword arguments provided are written directly
        as attribute variables (as with __init__).

        """
        for key, val in kwargs.items():
            setattr(self, key, val)

    def calc(self):
        """Overload method for calculation of final PSD.

        Should return an array of power-referenced values evaluated at
        all evaluation frequencies (self.freq).

        """
        return None

    ##########

    def __init__(self, freq, **kwargs):
        """Initialize budget item.

        Primary argument is the evaluation frequency array.  Any
        keyword arguments provided are simple written as attribute
        variables in the initialized object.

        """
        self.__freq = freq
        for key, val in kwargs.items():
            setattr(self, key, val)

    @property
    def name(self):
        """"Name of this BudgetItem class."""
        return self.__class__.__name__

    def __str__(self):
        # FIXME: provide info on internal state (load/update/calc/etc.)
        return '<{} {}>'.format(
            ' '.join([c.__name__ for c in self.__class__.__bases__]),
            self.name,
        )

    @property
    def freq(self):
        """Evaluation frequency array supplied at initialization."""
        return self.__freq

    def interpolate(self, freq, data):
        """Interpolate data to the evaluation frequencies.

        """
        func = scipy.interpolate.interp1d(
            freq, data,
            kind='nearest',
            copy=False,
            assume_sorted=True,
            bounds_error=False,
            fill_value=np.nan,
        )
        return func(self.freq)


class Calibration(BudgetItem):
    """GWINC Calibration class

    BudgetItem that represents a calibration transfer function for a
    Noise.  The calc() method should return a transfer function
    amplitude array evaluated at the evaluation frequencies supplied
    at initialization and available in the `freq` array attribute
    (self.freq).

    """
    def __call__(self, data):
        """Calibrate input data.

        Returns calibrated version of input data array,
        e.g. point-by-point product of data and calibration arrays.

        """
        cal = self.calc()
        assert data.shape == cal.shape, \
            "data shape does not match calibration ({} != {})".format(data.shape, cal.shape)
        return data * cal


class Noise(BudgetItem):
    """GWINC Noise class

    BudgetItem that represents a PSD noise calculation.  The calc()
    method should return the noise PSD spectrum array evaluated at the
    evaluation frequencies supplied at initialization and available in
    the `freq` array attribute (self.freq).

    """

    style = {}
    """Trace plot style dictionary"""

    def calc_trace(self, calibration=None, calc=True):
        """Returns noise (PSD, style) tuple.

        If `calibration` is not None it is assumed to be a
        len(self.freq) array that will be multiplied to the output
        PSD.

        If calc=False, the noise will not be calculated and the PSD
        will be None.  This is useful for just getting style the
        style.

        """
        if calc:
            data = self.calc()
            if calibration is not None:
                data *= calibration
        else:
            data = None
        return data, self.style


class Budget(Noise):
    """GWINC Budget class

    This is a Noise that represents the budget of multiple sub noises.

    The `noises` attribute of this class should be list constituent
    Noise classes.  Each element can be either a single Noise class,
    or a tuple of (Noise, Calibration) classes, e.g.:

    noises = [
        Thermal,
        (Shot, Sensing),
    ]

    When this object is initialized, all sub noises and calibrations
    are initialized.  Pre-defined load() and update() methods call the
    load() and update() methods of all sub noises and calibrations.
    When calc() is called, the PSD is calculated for all sub noises,
    the relevant calibration is evaluated and applied, and the
    quadrature sum of all calibrated consituent noises is returned.

    Additionally a `references` attribute may be definied, similar to
    the `noises` attribute described above except that the specified
    noises do not contribute to the overall budget total.

    NOTE: an `ifo` attribute is always passed as an initialization
    argument to sub noises.

    """

    noises = []
    """List of constituent noise classes"""

    references = []
    """List of reference nosie classes"""

    def __init__(self, *args, noises=None, **kwargs):
        """Initialize Budget object.

        See BudgetItem for base initialization arguments.

        If a `noises` keyword argument is provided it should be an
        iterable of noise names (constituent or reference) which will
        be used to filter the noises initialized in this budget.

        """
        super().__init__(*args, **kwargs)
        # store args and kwargs for later use
        self.args = args
        self.kwargs = kwargs
        # FIXME: special casing the IFO here, in case it's defined as
        # a class attribute rather than passed at initialization.  we
        # do this because we're not defining a standard way to extract
        # IFO variables that get passed around in a reasonable way.
        # how can we clarify this?
        if 'ifo' not in kwargs and getattr(self, 'ifo', None):
            self.kwargs['ifo'] = getattr(self, 'ifo', None)
        # all noise objects keyed by name
        self._noise_objs = collections.OrderedDict()
        # all cal objects keyed by name
        self._cal_objs = {}
        # noise to calibration mapping
        self._noise_cal = {}
        # set of all constituent budget noise names
        self._budget_noises = set()
        # initialize all noise objects
        for nc in self.noises:
            name, noise_obj, cal = self.__init_noise(nc)
            if noises and name not in noises:
                continue
            self.__add_noise(noise_obj, cal)
            self._budget_noises.add(name)
        for nc in self.references:
            name, noise_obj, cal = self.__init_noise(nc)
            if noises and name not in noises:
                continue
            self.__add_noise(noise_obj, cal)
        # error if requested noise is not present
        if noises:
            sset = set(noises)
            nset = set([name for name in self._noise_objs.keys()])
            if not sset <= nset:
                raise AttributeError("unknown noise terms: {}".format(' '.join(sset-nset)))

    def __init_noise(self, nc):
        cal = None
        if isinstance(nc, (list, tuple)):
            noise, cal = nc[:2]
        else:
            noise = nc
        noise_obj = noise(
            *self.args,
            **self.kwargs
        )
        return noise_obj.name, noise_obj, cal

    def __add_noise(self, noise_obj, cal):
        logging.debug("init {}".format(noise_obj))
        # instantiate the noise object
        name = noise_obj.name
        self._noise_objs[name] = noise_obj
        if cal is not None:
            # if a cal object is specified, instantiate and store
            cal_obj = cal(*self.args, **self.kwargs)
            if cal_obj.name not in self._cal_objs:
                self._cal_objs[cal_obj.name] = cal_obj
            self._noise_cal[name] = cal_obj.name

    def __getitem__(self, name):
        try:
            return self._noise_objs[name]
        except KeyError:
            try:
                return self._cal_objs[name]
            except KeyError:
                raise KeyError("unknown noise or cal name '{}".format(name))

    def keys(self):
        """Iterate over budget noise names."""
        return iter(self._noise_objs.keys())

    def values(self):
        """Iterate over budget noise objects."""
        return iter(self._noise_objs.values())

    def items(self):
        """Iterate over budget noise (name, object) tuples."""
        return iter(self._noise_objs.items())

    def __iter__(self):
        return iter(self.keys())

    def walk(self):
        """Walk recursively through every BudgetItem in the budget.

        This includes Noise, Calibration and Budget objects, as well
        as any decendents of Budget objects.

        For each leaf item yields a tuple of all ancestor objects,
        e.g.:

          (self)
          (self, BudgetItem)
          (self, ChildBudget1)
          (self, ChildBudget1, BudgetItem)
          ...

        """
        yield (self,)
        for item in itertools.chain(
                self._cal_objs.values(),
                self._noise_objs.values()):
            if isinstance(item, Budget):
                for i in item.walk():
                    yield (self,) + i
            else:
                yield (self, item)

    def load(self):
        """Load all noise and cal objects."""
        for name, item in itertools.chain(
                self._cal_objs.items(),
                self._noise_objs.items()):
            logging.debug("load {}".format(item))
            item.load()

    def update(self, **kwargs):
        """Update all noise and cal objects with supplied kwargs."""
        for name, item in itertools.chain(
                self._cal_objs.items(),
                self._noise_objs.items()):
            logging.debug("update {}".format(item))
            item.update(**kwargs)

    def cal_for_noise(self, name):
        """Return the calibration object for named noise."""
        try:
            return self._cal_objs[self._noise_cal[name]]
        except KeyError:
            return None

    def calc_noise(self, name):
        """Return calibrated individual noise term.

        The noise PSD and calibration transfer functions are
        calculated, and the calibrated noise array is returned.

        """
        noise = self._noise_objs[name]
        nd = noise.calc()
        cal = self.cal_for_noise(name)
        if cal:
            cd = cal.calc()
            return nd * cd
        else:
            return nd

    def calc(self):
        """Calculate sum of all noises.

        """
        data = [self.calc_noise(name) for name in self._noise_objs.keys() if name in self._budget_noises]
        return quadsum(data)

    def calc_trace(self, calibration=None, calc=True):
        """Returns a dictionary of noises traces, keyed by noise names.

        Values are (data, style) trace tuples (see Noise.calc_trace).
        The key of the budget sum total is 'Total'.  The values of sub
        budgets are themselves dictionaries returned from
        calc_trace() of the sub budget.

        If `calibration` is not None it is assumed to be a
        len(self.freq) array that will be multiplied to the output
        PSD of the budget and all sub noises.

        If calc=False, the noise will not be calculated and the PSD
        will be None.  This is useful for just getting style the
        style.

        """
        # start by creating an empty OrderedDict used for outputing trace data
        # or style info with the following order:
        #   references
        #   total
        #   constituents
        d = collections.OrderedDict()
        # allocate references
        for name, noise in self._noise_objs.items():
            if name in self._budget_noises:
                continue
            d[name] = noise.calc_trace(calc=False)
        # allocate total
        if self._budget_noises:
            d['Total'] = None, self.style
        # allocate constituent
        for name, noise in self._noise_objs.items():
            if name not in self._budget_noises:
                continue
            d[name] = noise.calc_trace(calc=False)
        # if we're not calc'ing, just return the dict now
        if not calc:
            return d
        # calc all noises
        for name, noise in self._noise_objs.items():
            # extract/calc the budget-level calibration for this noise
            cal_obj = self.cal_for_noise(name)
            if cal_obj:
                cal = cal_obj.calc()
            else:
                cal = np.ones_like(self.freq)
            # then multiply by the supplied calibration
            if calibration is not None:
                cal *= calibration
            d[name] = noise.calc_trace(calibration=cal, calc=True)
        # calc budget total
        constituent_data = []
        for name in self._budget_noises:
            if isinstance(d[name], dict):
                data = d[name]['Total'][0]
            else:
                data = d[name][0]
            constituent_data.append(data)
        d['Total'] = quadsum(constituent_data), self.style
        return d
