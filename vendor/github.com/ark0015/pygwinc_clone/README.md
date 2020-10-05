[![pipeline status](https://git.ligo.org/gwinc/pygwinc/badges/master/pipeline.svg)](https://git.ligo.org/gwinc/pygwinc/commits/master)

# Python port of GW Interferometer Noise Calculator

![gwinc](https://gwinc.docs.ligo.org/pygwinc/aLIGO.png)

This is a collection of mostly analytic noise calculations (e.g. quantum, thermal)


## basic usage

`pygwinc` creates noise budgets based on detector descriptions
provided in either .yml or .mat files (see below).  Once the detector
description is loaded, the noise budget can be calculated and plotted:
```python
>>> import gwinc
>>> import numpy as np
>>> freq = np.logspace(1, 3, 1000)
>>> Budget, ifo, freq_, plot_style = gwinc.load_ifo('aLIGO')
>>> ifo = gwinc.precompIFO(freq, ifo)
>>> traces = Budget(freq, ifo=ifo).calc_trace()
>>> fig = gwinc.plot_noise(freq, traces, **plot_style)
>>> fig.show()
```


## command line interface

You can make gwinc plots directly from the command line by executing
the package directly:
```shell
$ python3 -m gwinc aLIGO
```


## detector description files

`pygwinc` can load budget descriptions in different formats: python
package/module, .yaml YAML file, and MATLAB gwinc .mat or .m files.

`pygwinc` includes budgets for various canonical detectors:

* [aLIGO](https://git.ligo.org/gwinc/pygwinc/blob/master/gwinc/ifo/aLIGO)
* [A+](https://git.ligo.org/gwinc/pygwinc/blob/master/gwinc/ifo/Aplus)
* [Voyager](https://git.ligo.org/gwinc/pygwinc/blob/master/gwinc/ifo/Voyager)
* [Cosmic Explorer](https://git.ligo.org/gwinc/pygwinc/blob/master/gwinc/ifo/CE)


## noise budgets


GWINC provides an `nb` package for defining arbitrary noise budgets:

```python
import numpy as np
from gwinc import nb
from gwinc import noise


class ExcessGas(nb.Noise):
    """Excess gas"""
    style = dict(
        label='Excess Gas',
        color='#ad900d',
        linestyle='--',
    )

    def calc(self):
        return noise.residualgas.gas(self.freq, self.ifo)


class MeasuredNoise(nb.Noise):
    """My measured noise"""
    style = dict(
        label='Measured Noise',
        color='#838209',
        linestyle='-',
    )

    def load(self):
        psd, freq = np.loadtxt('/path/to/measured/psd.txt')
        self.data = self.interpolate(f, psd)

    def calc(self):
        return self.data


class MyBudget(nb.Budget):
    noises = [
        ExcessGas,
        MeasuredNoise,
    ]
```


## comparison with MATLAB gwinc

`pygwinc` includes the ability use MATLAB gwinc directly via the
MATLAB python interface (see the CLI '--matlab' option above).  This
also allows for easy direct comparison between the pygwinc and
matgwinc noise budgets.

If you have a local checkout of matgwinc (at e.g. /path/to/gwinc) and
a local installation of MATLAB and it's python interface (at
e.g. /opt/matlab/python/lib/python3.6/site-packages) you can run the
comparison as so:
```shell
$ export GWINCPATH=/path/to/matgwinc
$ export PYTHONPATH=/opt/matlab/python/lib/python3.6/site-packages
$ python3 -m gwinc.test -p aLIGO
```
This will produce a summary page of the various noise spectra that
differ between matgwinc and pygwinc.

Latest comparison plots from continuous integration:

* [aLIGO comparison](https://gwinc.docs.ligo.org/pygwinc/aLIGO_test.png)
* [A+ comparison](https://gwinc.docs.ligo.org/pygwinc/A+_test.png)
