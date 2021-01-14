---
title: 'Gwent: A Python package for assessing gravitational wave detector designs'
tags:
  - Python
  - astronomy
  - gravitational waves
  - gravitational wave detectors
authors:
  - name: Andrew R. Kaiser^[Corresponding author]
    orcid: 0000-0002-3654-980X
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Sean T. McWilliams
    orcid: 0000-0003-2397-8290
    affiliation: "1, 2"
  - name: Jeffrey S. Hazboun
    orcid: 0000-0003-2742-3321
    affiliation: 3
affiliations:
 - name: West Virginia University
   index: 1
 - name: Center for Gravitational Waves and Cosmology
   index: 2
 - name: University of Washington Bothell
   index: 3
date: 23 September 2020
bibliography: paper.bib
---

# Summary

Black-holes are known to span at least nine orders of magnitude in mass: from the stellar-mass objects born from the most massive stars, to supermassive black-holes at the center of galaxies. 
Regardless of the mass scale, all of these objects are expected to form binaries and eventually emit observable gravitational waves (GWs). 
The observation of GWs from black-hole binaries (BHBs) allows us to probe the Universe on a wide range of scales.
Because of the broad GW spectrum from the BHBs across mass scales, there is no single instrument capable of observing the entire scope.
Stellar BHBs up to hundreds of solar masses can be observed from the ground, whereas more massive objects can only be observed from space. 
Massive binaries that are millions of solar masses require space-based detectors, whereas supermassive binaries that are billions of solar masses require detectors larger than our solar system, and are currently being pursued with pulsar timing arrays (PTAs).

``gwent`` is a Python [@Python] package for modeling the sensitivities of current and future generations of GW detectors across the entire GW spectrum of coalescing BHBs. 
It provides methods based on the formalism in [@Kaiser:2020] to generate sensitivity curves for PTAs using a novel realistic PTA sensitivity curve generator [@Hazboun:2019; @Hasasia:2019], space-based interferometers using adaptive models that can represent a wide range of proposed detector designs [@AmaroSeoane:2017], and ground-based interferometers using realistic models that can reproduce current [@Abbott:2016], second, and third generation designs [@Hild:2011], as well as novel variations of the essential design parameters. 
To model the signal from BHBs at any mass scale, ``gwent`` uses GW waveforms capable of modeling each piece of a binary's evolution for sources with varying mass ratios and spins [@Khan:2016; @Husa:2016]. 
It uses standard Python packages, such as ``NumPy`` [@Numpy], ``SciPy`` [@Scipy], ``Matplotlib`` [@Matplotlib] and ``Astropy`` [@Astropy] to construct these GW signals and detector models. 
Additionally, ``gwent`` utilizes and explicitly includes within the package a specific commit of ``pygwinc`` [@Pygwinc] to create many of the ground-based GW detector sensitivity curves as at the time of creation, there is no ``pygwinc`` availability on ``PyPI``.

Using this adaptable framework, we produce signal-to-noise ratios (SNRs) for the combination of any modeled parameter, associated with either the detector or the source. 
By allowing variation across each detector and source parameter, we can pinpoint the most important factors to determining the optimal performance for particular instrument designs. 
The adaptability of our detector and signal models can easily be extended to new detector designs and other models of GW signals.

This software is designed to be used by astronomers to assess a BHB's detectability in any proposed detector. 
In an effort to make the package useful to a broader community, ``gwent`` includes many standard detectors. 
To be useful to a more advanced audience, ``gwent`` includes the ability to use cutting edge GW source waveform models available in ``lalsuite`` [@lalsuite], and user-provided GW detectors to assess their feasibility. 
It has already been used in a number of scientific publications [@Chen:2020; @Campeti:2020] and a GW detector design proposal [@AmaroSeoane:2020].

# Acknowledgments

We thank Maura McLaughlin for helpful feedback, support, and guidance throughout the work.
This work was supported by NSF PFC Grant PHY-1430284 and NSF CAREER Grant PHY-1945130.

# References