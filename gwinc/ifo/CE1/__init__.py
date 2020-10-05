from gwinc.ifo.noises import *

class Newtonian(nb.Budget):
    """Newtonian Gravity

    """

    name = 'Newtonian'

    style = dict(
        label='Newtonian',
        color='#15b01a',
    )

    noises = [
        NewtonianRayleigh,
        NewtonianBody,
        NewtonianInfrasound,
    ]

class Coating(nb.Budget):
    """Coating Thermal

    """

    name = 'Coating'

    style = dict(
        label='Coating Thermal',
        color='#fe0002',
    )

    noises = [
        CoatingBrownian,
        CoatingThermoOptic,
    ]

class Substrate(nb.Budget):
    """Substrate Thermal

    """

    name = 'Substrate'

    style = dict(
        label='Substrate Thermal',
        color='#fb7d07',
        linestyle='--',
    )

    noises = [
        SubstrateBrownian,
        SubstrateThermoElastic,
    ]

class CE1(nb.Budget):

    name = 'Cosmic Explorer 1'

    noises = [
        QuantumVacuum,
        Seismic,
        Newtonian,
        SuspensionThermal,
        Coating,
        Substrate,
        ExcessGas,
    ]
