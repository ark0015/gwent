#!/usr/bin/env python
# -*- coding: utf-8 -*-

import astropy.units as u
import numpy as np
import pytest

from gwent import detector, utils

unit_1 = 1.0 * u.Hz


def test_conversion_error():
    with pytest.raises(ValueError):
        utils.make_quant(unit_1, "s")


# LISA Proposal 2
# Values from Robson, Cornish, and Liu 2019 https://arxiv.org/abs/1803.01944 using the Transfer Function Approximation within.
L = 2.5 * u.Gm  # armlength in Gm
L = L.to("m")
LISA_T_obs = 4.0 * u.yr
f_acc_break_low = 0.4 * u.mHz.to("Hz") * u.Hz
f_acc_break_high = 8.0 * u.mHz.to("Hz") * u.Hz
f_IMS_break = 2.0 * u.mHz.to("Hz") * u.Hz
A_acc = 3e-15 * u.m / u.s / u.s
A_IMS = 1.5e-11 * u.m
Background = True


@pytest.fixture
def LISA_prop2():
    LISA_prop2 = detector.SpaceBased(
        "LISA Approximate",
        LISA_T_obs,
        L,
        A_acc,
        f_acc_break_low,
        f_acc_break_high,
        A_IMS,
        f_IMS_break,
        Background=Background,
        Background_model=1,
        T_type="A",
    )
    return LISA_prop2


def test_DictError_2(LISA_prop2):
    unit_2 = ["L", np.array([1, 2, 3])]
    with pytest.raises(ValueError):
        utils.Get_Var_Dict(LISA_prop2, unit_2)


def test_DictError_3(LISA_prop2):
    unit_3 = ["L", [1, [2, 3], 3]]
    with pytest.raises(ValueError):
        utils.Get_Var_Dict(LISA_prop2, unit_3)


def test_DictError_Full(LISA_prop2):
    unit_4 = np.array(["L", 1])
    with pytest.raises(ValueError):
        utils.Get_Var_Dict(LISA_prop2, unit_4)
