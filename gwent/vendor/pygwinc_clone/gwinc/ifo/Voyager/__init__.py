from gwinc.ifo.noises import *


class Voyager(nb.Budget):

    name = 'Voyager'

    noises = [
        QuantumVacuum,
        Seismic,
        Newtonian,
        SuspensionThermal,
        CoatingBrownian,
        CoatingThermoOptic,
        ITMThermoRefractive,
        ITMCarrierDensity,
        SubstrateBrownian,
        SubstrateThermoElastic,
        ExcessGas,
    ]
