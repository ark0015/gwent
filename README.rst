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

.. image:: https://raw.githubusercontent.com/ark0015/gwent/master/docs/calcSNR_tutorial_files/full_waterfall_plots_lb.png
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

	* Uses a fully Pythonic implementation of the phenomenological model `IMRPhenomD` to accurately represent the inspiral, merger, and ringdown of the BHB.

Calculates the matched-filtered signal-to-noise ratio (SNR) to help assess the detectability of any BHB source configuration by any represented gravitational wave detector.

* Includes robust plotting methods to represent these SNRs.


Getting Started
---------------
`gwent` is available on the Python Package Inventory, so the preferred method to install `gwent` is to install it with `pip`, as it will always install the most recent stable release.

.. code-block:: console

    $ pip install gwent


To install `pygwinc`, a GitLab hosted package necessary to fully utilize `gwent`, run this command in your terminal:

.. code-block:: console

    $ pip install git+https://git.ligo.org/gwinc/pygwinc.git@65396ee42e851ab7189618cabe1c12081b5d982e#egg=pygwinc

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
