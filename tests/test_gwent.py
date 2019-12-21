#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `gwent` package."""


import pytest

from gwent import binary
from gwent import detector
from gwent import snr
from gwent import snrplot
from gwent import utils
from gwent import waveform


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string
