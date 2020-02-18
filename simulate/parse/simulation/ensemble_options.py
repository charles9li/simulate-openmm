from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from copy import copy

from simtk.openmm.app import Simulation
from simtk.openmm import VerletIntegrator, LangevinIntegrator
from simtk.unit import amu, nanometer, picosecond
from openmmtools.integrators import VelocityVerletIntegrator, NoseHooverChainVelocityVerletIntegrator

from simulate.parse._options import _Options
from simulate.parse.simulation.integrator_options import *
from simulate.parse.simulation.barostat_options import *
from simulate.parse.simulation.reporter_options import *
from simulate.parse.simulation.minimize_energy_options import *
from simulate.parse.simulation.average_options import *


class _EnsembleOptions(_Options):

    def __init__(self):
        super(_EnsembleOptions, self).__init__()
        self.integrator = None
        self.reporters = []
        self.minimize_energy_options = None
        self.steps = 0

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.integrator is None:
            self._incomplete_error('integrator')

    # =========================================================================

    INTEGRATOR_OPTIONS = {'VerletIntegrator': VerletIntegratorOptions,
                          'VelocityVerletIntegrator': VelocityVerletIntegratorOptions}
    REPORTER_OPTIONS = {'DCDReporter': DCDReporterOptions,
                        'StateDataReporter': StateDataReporterOptions}

    def _parse_integrator(self, *args):
        integrator_name = args[0]
        line_deque = args[1]
        integrator_options = self.INTEGRATOR_OPTIONS[integrator_name]()
        integrator_options.parse(line_deque.popleft())
        self.integrator = integrator_options.integrator()

    def _parse_reporters(self, *args):
        line_deque = args[1].popleft()
        while len(line_deque) > 0:
            reporter_name = line_deque.popleft()
            reporter_options = self.REPORTER_OPTIONS[reporter_name]()
            reporter_options.parse(line_deque.popleft())
            self.reporters.append(reporter_options.reporter())

    def _parse_minimize_energy(self, *args):
        line_deque = args[1].popleft()
        minimize_energy_options = MinimizeEnergyOptions()
        minimize_energy_options.parse(line_deque)
        self.minimize_energy_options = minimize_energy_options

    def _parse_steps(self, *args):
        self.steps = literal_eval(args[0])

    OPTIONS = {'integrator': _parse_integrator,
               'steps': _parse_steps}

    SECTIONS = {'reporters': _parse_reporters,
                'minimizeEnergy': _parse_minimize_energy}

    # =========================================================================

    def simulation(self, topology, system):
        simulation = Simulation(topology, system, self.integrator)
        for reporter in self.reporters:
            simulation.reporters.append(reporter)
        return simulation


class NVEOptions(_EnsembleOptions):

    SECTION_NAME = 'NVE'

    def __init__(self):
        super(NVEOptions, self).__init__()
        self.average_options = None

    # =========================================================================

    def _check_for_incomplete_input(self):
        super(NVEOptions, self)._check_for_incomplete_input()

    # =========================================================================

    def _parse_average(self, *args):
        line_deque = args[1].popleft()
        average_options = AverageOptions()
        average_options.parse(line_deque)
        self.average_options = average_options

    SECTIONS = copy(_EnsembleOptions.SECTIONS)
    SECTIONS['average'] = _parse_average


class NVTOptions(NVEOptions):

    SECTION_NAME = 'NVT'

    def __init__(self):
        super(NVTOptions, self).__init__()
        self.thermostat = None

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
                                 "integrator in the {} ensemble".format(self.SECTION_NAME))
        if isinstance(self.integrator, self.INTEGRATORS_WITH_THERMOSTAT):
            if self.thermostat is not None:
                raise ValueError("A thermostat does not need to be used with specified "
                                 "integrator in the {} ensemble".format(self.SECTION_NAME))

    # =========================================================================

    INTEGRATOR_OPTIONS = copy(NVEOptions.INTEGRATOR_OPTIONS)
    INTEGRATOR_OPTIONS['LangevinIntegrator'] = LangevinIntegratorOptions

    THERMOSTAT_OPTIONS = {}

    def _parse_thermostat(self, *args):
        thermostat_name = args[0]
        line_deque = args[1]
        thermostat_options = self.THERMOSTAT_OPTIONS[thermostat_name]()
        thermostat_options.parse(line_deque.popleft())
        self.thermostat = thermostat_options.thermostat()

    OPTIONS = copy(NVEOptions.OPTIONS)
    OPTIONS['thermostat'] = _parse_thermostat

    # =========================================================================

    def simulation(self, topology, system):
        if self.thermostat is not None:
            system.addForce(self.thermostat)
        return super(NVTOptions, self).simulation(topology, system)


class NPTOptions(NVTOptions):

    SECTION_NAME = 'NPT'

    def __init__(self):
        super(NPTOptions, self).__init__()
        self.barostat = None

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

    OPTIONS = copy(NVTOptions.OPTIONS)
    OPTIONS['barostat'] = _parse_barostat

    # =========================================================================
    
    def simulation(self, topology, system):
        system.addForce(self.barostat)
        return super(NPTOptions, self).simulation(topology, system)


class RNEMDOptions(NVTOptions):

    SECTION_NAME = 'RNEMD'

    def __init__(self):
        super(RNEMDOptions, self).__init__()
        self.numSlabs = None
        self.swapFrequency = None

    # =========================================================================

    def _check_for_incomplete_input(self):
        super(super(RNEMDOptions, self), self)._check_for_incomplete_input()
        if self.numSlabs is None:
            self._incomplete_error('numSlabs')
        if self.swapFrequency is None:
            self._incomplete_error('swapFrequency')

    # =========================================================================

    # TODO: add RNEMDReporter option when done
    REPORTER_OPTIONS = copy(NVEOptions.REPORTER_OPTIONS)

    def _parse_num_slabs(self, *args):
        self.numSlabs = literal_eval(args[0])

    def _parse_swap_frequency(self, *args):
        self.swapFrequency = literal_eval(args[0])

    OPTIONS = copy(_EnsembleOptions.OPTIONS)
    OPTIONS['numSlabs'] = _parse_num_slabs
    OPTIONS['swapFrequency'] = _parse_swap_frequency

    # =========================================================================

    def simulation(self, topology, system):
        if self.thermostat is not None:
            system.addForce(self.thermostat)
        simulation = super(super(RNEMDOptions, self), self).simulation(topology, system)
        simulation.context.totalMomentumExchanged = 0.0*amu*nanometer/picosecond
        return simulation
