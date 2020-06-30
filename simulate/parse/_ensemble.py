"""
_ensemble.py: Parses ensemble options.

Copyright (c) 2020 Charles Li // UCSB, Department of Chemical Engineering

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from ast import literal_eval
import os

from simtk.openmm.app import Simulation
from simtk.openmm import VerletIntegrator, LangevinIntegrator
from simtk.unit import amu, nanometer, picosecond
from openmmtools.integrators import VelocityVerletIntegrator, NoseHooverChainVelocityVerletIntegrator
import numpy as np

from ._options import _Options
from ._ensemble_integrator import *
from ._ensemble_barostat import *
from ._ensemble_reporter import *
from ._ensemble_minimizeenergy import *
from ._ensemble_average import *
from ._ensemble_temperatureramp import *


class _EnsembleOptions(_Options):

    _SECTION_NAME = "_Ensemble"

    # =========================================================================

    def __init__(self, simulations_options):
        super(_EnsembleOptions, self).__init__()
        self._create_integrator_options()
        self._create_reporter_options()
        self.simulations_options = simulations_options
        self.integrator_options = None
        self.integrator = None
        self.reporters = []
        self.minimize_energy_options = None
        self.steps = 0
        self.saveState = None
        self.loadState = None

    def _create_options(self):
        super(_EnsembleOptions, self)._create_options()
        self._OPTIONS['integrator'] = self._parse_integrator
        self._OPTIONS['steps'] = self._parse_steps
        self._OPTIONS['saveState'] = self._parse_save_state
        self._OPTIONS['loadState'] = self._parse_load_state

    def _create_sections(self):
        super(_EnsembleOptions, self)._create_sections()
        self._SECTIONS['reporters'] = self._parse_reporters
        self._SECTIONS['minimizeEnergy'] = self._parse_minimize_energy

    def _create_integrator_options(self):
        self._INTEGRATOR_OPTIONS = {'VerletIntegrator': VerletIntegratorOptions,
                                    'VelocityVerletIntegrator': VelocityVerletIntegratorOptions}

    def _create_reporter_options(self):
        self._REPORTER_OPTIONS = {'DCDReporter': DCDReporterOptions,
                                  'PDBReporter': PDBReporterOptions,
                                  'StateDataReporter': StateDataReporterOptions,
                                  'EnergyReporter': EnergyReporterOptions,
                                  'PotentialEnergyReporter': PotentialEnergyReporterOptions,
                                  'KineticEnergyReporter': KineticEnergyReporterOptions,
                                  'RadiusOfGyrationReporter': RadiusOfGyrationReporterOptions,
                                  'EndToEndDistanceReporter': EndToEndDistanceReporterOptions,
                                  'CheckpointReporter': CheckpointReporterOptions}

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.integrator is None:
            self._incomplete_error('integrator')

    # =========================================================================

    def _parse_integrator(self, *args):
        integrator_name = args[0]
        line_deque = args[1]
        integrator_options = self._INTEGRATOR_OPTIONS[integrator_name]()
        integrator_options.parse(line_deque.popleft())
        self.integrator_options = integrator_options
        self.integrator = integrator_options.integrator()

    def _parse_reporters(self, *args):
        line_deque = args[1].popleft()
        while len(line_deque) > 0:
            reporter_name = line_deque.popleft()
            reporter_options = self._REPORTER_OPTIONS[reporter_name](self)
            reporter_options.parse(line_deque.popleft())
            self.reporters.append(reporter_options.reporter())

    def _parse_minimize_energy(self, *args):
        line_deque = args[1].popleft()
        minimize_energy_options = MinimizeEnergyOptions(self)
        minimize_energy_options.parse(line_deque)
        self.minimize_energy_options = minimize_energy_options

    def _parse_steps(self, *args):
        self.steps = literal_eval(args[0])

    def _parse_save_state(self, *args):
        self.saveState = self._create_filepath(args[0])

    def _parse_load_state(self, *args):
        self.loadState = self._create_filepath(args[0])

    # =========================================================================

    # Helper methods for parsing options

    def _create_filepath(self, filepath):
        directory = self.simulations_options.input_options.directory
        return os.path.join(directory, filepath)

    # =========================================================================

    def create_simulation(self, topology, system):
        simulation = Simulation(topology, system, self.integrator)
        for reporter in self.reporters:
            simulation.reporters.append(reporter)
        return simulation

    def create_integrator(self):
        return self.integrator_options.integrator()


class NVEOptions(_EnsembleOptions):

    _SECTION_NAME = 'NVE'

    # =========================================================================

    def __init__(self, simulations_options):
        super(NVEOptions, self).__init__(simulations_options)
        self.average_options = None

    def _create_sections(self):
        super(NVEOptions, self)._create_sections()
        self._SECTIONS['average'] = self._parse_average

    # =========================================================================

    def _check_for_incomplete_input(self):
        super(NVEOptions, self)._check_for_incomplete_input()

    # =========================================================================

    def _parse_average(self, *args):
        line_deque = args[1].popleft()
        average_options = AverageOptions()
        average_options.parse(line_deque)
        self.average_options = average_options


class NVTOptions(NVEOptions):

    _SECTION_NAME = 'NVT'

    def __init__(self, simulations_options):
        super(NVTOptions, self).__init__(simulations_options)
        self.thermostat = None
        self.temperature_ramp_options = None

    def _create_options(self):
        super(NVTOptions, self)._create_options()
        self._OPTIONS['thermostat'] = self._parse_thermostat

    def _create_sections(self):
        super(NVTOptions, self)._create_sections()
        self._SECTIONS['temperatureRamp'] = self._parse_temperature_ramp

    def _create_integrator_options(self):
        super(NVTOptions, self)._create_integrator_options()
        self._INTEGRATOR_OPTIONS['LangevinIntegrator'] = LangevinIntegratorOptions

    def _create_thermostat_options(self):
        self._THERMOSTAT_OPTIONS = {}

    # =========================================================================

    INTEGRATORS_NO_THERMOSTAT = (VerletIntegrator,
                                 VelocityVerletIntegrator)
    INTEGRATORS_WITH_THERMOSTAT = (LangevinIntegrator,
                                   NoseHooverChainVelocityVerletIntegrator)
    
    def _check_for_incomplete_input(self):
        super(NVTOptions, self)._check_for_incomplete_input()
        if isinstance(self.integrator, self.INTEGRATORS_NO_THERMOSTAT):
            if self.thermostat is None:
                raise ValueError("Thermostat must be used with the specified "
                                 "integrator in the {} ensemble".format(self._SECTION_NAME))
        if isinstance(self.integrator, self.INTEGRATORS_WITH_THERMOSTAT):
            if self.thermostat is not None:
                raise ValueError("A thermostat does not need to be used with specified "
                                 "integrator in the {} ensemble".format(self._SECTION_NAME))
        if self.temperature_ramp_options is not None and self.average_options is not None:
            raise ValueError("The average section cannot be used with the temperature ramp section.")

    # =========================================================================

    def _parse_thermostat(self, *args):
        thermostat_name = args[0]
        line_deque = args[1]
        thermostat_options = self._THERMOSTAT_OPTIONS[thermostat_name]()
        thermostat_options.parse(line_deque.popleft())
        self.thermostat = thermostat_options.thermostat()

    def _parse_temperature_ramp(self, *args):
        line_deque = args[1].popleft()
        temperature_ramp_options = TemperatureRampOptions()
        temperature_ramp_options.parse(line_deque)
        self.temperature_ramp_options = temperature_ramp_options

    # =========================================================================

    def create_simulation(self, topology, system):
        if self.thermostat is not None:
            system.addForce(self.thermostat)
        return super(NVTOptions, self).create_simulation(topology, system)


class NPTOptions(NVTOptions):

    _SECTION_NAME = 'NPT'

    def __init__(self, simulations_options):
        super(NPTOptions, self).__init__(simulations_options)
        self.barostat = None

    def _create_options(self):
        super(NPTOptions, self)._create_options()
        self._OPTIONS['barostat'] = self._parse_barostat

    # =========================================================================

    def _check_for_incomplete_input(self):
        super(NPTOptions, self)._check_for_incomplete_input()
        if self.barostat is None:
            self._incomplete_error('barostat')

    # =========================================================================

    BAROSTAT_OPTIONS = {'MonteCarloBarostat': MonteCarloBarostatOptions}

    def _parse_barostat(self, *args):
        barostat_name = args[0]
        line_deque = args[1]
        barostat_options = self.BAROSTAT_OPTIONS[barostat_name]()
        barostat_options.parse(line_deque.popleft())
        self.barostat = barostat_options.barostat()

    # =========================================================================
    
    def create_simulation(self, topology, system):
        system.addForce(self.barostat)
        self.barostat.setForceGroup(system.getNumForces() - 1)
        return super(NPTOptions, self).create_simulation(topology, system)


class RNEMDOptions(_EnsembleOptions):

    _SECTION_NAME = 'RNEMD'

    def __init__(self, simulations_options):
        super(RNEMDOptions, self).__init__(simulations_options)
        self.thermostat = None
        self.numSlabs = None
        self.swapFrequency = None

    def _create_options(self):
        super(RNEMDOptions, self)._create_options()
        self._OPTIONS['thermostat'] = self._parse_thermostat
        self._OPTIONS['numSlabs'] = self._parse_num_slabs
        self._OPTIONS['swapFrequency'] = self._parse_swap_frequency

    def _create_integrator_options(self):
        super(RNEMDOptions, self)._create_integrator_options()
        self._INTEGRATOR_OPTIONS['LangevinIntegrator'] = LangevinIntegratorOptions

    def _create_reporter_options(self):
        super(RNEMDOptions, self)._create_reporter_options()
        self._REPORTER_OPTIONS['RNEMDReporter'] = RNEMDReporterOptions
        self._REPORTER_OPTIONS['RNEMDVelocityReporter'] = RNEMDVelocityReporterOptions

    # =========================================================================

    def _check_for_incomplete_input(self):
        super(RNEMDOptions, self)._check_for_incomplete_input()
        if self.numSlabs is None:
            self._incomplete_error('numSlabs')
        if self.swapFrequency is None:
            self._incomplete_error('swapFrequency')

    # =========================================================================

    def _parse_thermostat(self, *args):
        pass

    def _parse_num_slabs(self, *args):
        self.numSlabs = literal_eval(args[0])

    def _parse_swap_frequency(self, *args):
        self.swapFrequency = literal_eval(args[0])

    # =========================================================================

    def create_simulation(self, topology, system):
        if self.thermostat is not None:
            system.addForce(self.thermostat)
        simulation = super(RNEMDOptions, self).create_simulation(topology, system)
        simulation.context.totalMomentumExchanged = np.array([0.0] * 3) * amu*nanometer/picosecond
        simulation.context.xMomentumExchanged = 0.0*amu*nanometer/picosecond
        simulation.context.yMomentumExchanged = 0.0*amu*nanometer/picosecond
        simulation.context.zMomentumExchanged = 0.0*amu*nanometer/picosecond
        return simulation
