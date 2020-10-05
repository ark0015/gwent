#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from setuptools import find_packages, setup


version = '0.1'


setup_args = dict(
    name             = 'pygwinc',
    version          = version,
    url              = 'https://git.ligo.org/gwinc/pygwinc',
    author           = 'LIGO Laboratory',
    author_email     = 'jrollins@ligo.caltech.edu ',
    description      = "Gravitation Wave Interferometer Noise Calculator",
    license          = 'Copyright 2018 LIGO Laboratory',
    keywords         = 'Noise, LIGO, Gravitational Wave,',
    classifiers = [
        'Topic :: Scientific/Engineering',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],

    install_requires = [
        'h5py',
        'matplotlib',
        'numpy',
        'pyyaml',
        'scipy',
    ],

    packages = find_packages(
        exclude = ['docs',],
    ),

    entry_points={
        'console_scripts': [
            'gwinc = gwinc.__main__:main',
        ],
    },

    include_package_data = True,
    zip_safe = False,
)

if __name__ == "__main__":
    setup(**setup_args)
