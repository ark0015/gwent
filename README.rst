=====
gwent
=====


.. image:: https://img.shields.io/pypi/v/gwent.svg
        :target: https://pypi.python.org/pypi/gwent

.. image:: https://img.shields.io/travis/ark0015/gwent.svg
        :target: https://travis-ci.org/ark0015/gwent

.. image:: https://readthedocs.org/projects/gwent/badge/?version=latest
        :target: https://gwent.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status


Gravitational Wave dEtector desigN Toolkit.

Generates strain sensitivity curves and Waterfall plots for various gravitational wave detector designs.

.. image:: https://raw.githubusercontent.com/ark0015/gwent/master/data/full_waterfall_plots_lb.png
        :align: center
        :alt: gwent Waterfall Plots

* Free software: MIT license
* Documentation: https://gwent.readthedocs.io.


Features
--------
Calculates the sensitivity curves for various designs of pulsar timing arrays, space-based detectors, and ground-based detectors.
This includes:

* NANOGrav
* SKA
* LISA
* aLIGO
* Voyager
* and more!

Calculates the strain from coalescing black hole binaries. It contains functionality for different source descriptions:

* Slowly-evolving sources, ie. BHBs early in their inspiral where they appear to not change in frequency.
* Rapidly-evolving sources, ie. BHBs in the final stages of coalescence. 

  * Uses a fully Pythonic implementation of the phenomenological model ``IMRPhenomD`` to accurately represent the inspiral, merger, and ringdown of the BHB.

Calculates the matched-filtered signal-to-noise ratio (SNR) to help assess the detectability of any BHB source configuration by any represented gravitational wave detector.

* Includes robust plotting methods to represent these SNRs.


Getting Started
---------------
``gwent`` is available on the Python Package Inventory, so the preferred method to install ``gwent`` is to install it with ``pip``, as it will always install the most recent stable release.

.. code-block:: console

    $ pip install gwent

README Figure and Data
----------------------
If you are looking for quick data, we conveniently place the figure above in the `data <https://github.com/ark0015/gwent/tree/master/data>`_ folder on the Github repo. There you can also find the raw data used for this figure in ``.npz`` format. To load this data, simply use ``np.load(filename)``, and the data can be accessed by the kwargs ``'mass'``, ``'redshift'``, and ``'snr'``. E.g., 

.. code-block:: python

    import numpy as np
    import gwent
    from gwent.snrplot import Plot_SNR
    loaded_file = np.load(filename)
    Plot_SNR('M',load_file['mass'],'z',load_file['redshift'],load_file['snr'])
    
Publication
-----------
This work and methodology is available on arXiv_. If you use ``gwent``, please cite this work using the following:

.._arXiv: https://arxiv.org/abs/2010.02135

.. code-block:: tex

    @ARTICLE{2020arXiv201002135K,
           author = {{Kaiser}, Andrew R. and {McWilliams}, Sean T.},
            title = "{Sensitivity of present and future black-hole binary observations across the gravitational wave spectrum}",
          journal = {arXiv e-prints},
         keywords = {General Relativity and Quantum Cosmology, Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Instrumentation and Methods for Astrophysics},
             year = 2020,
            month = oct,
              eid = {arXiv:2010.02135},
            pages = {arXiv:2010.02135},
    archivePrefix = {arXiv},
           eprint = {2010.02135},
     primaryClass = {gr-qc},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2020arXiv201002135K},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

Credits
-------
We utilize and include within the package a specific commit of ``pygwinc`` found at https://git.ligo.org/gwinc/pygwinc to create many of the ground-based gravitational wave detector sensitivity curves. At the time of creation, there is no ``pygwinc`` availability on PyPI, so we explicitly include the necessary portions of the code within.

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
