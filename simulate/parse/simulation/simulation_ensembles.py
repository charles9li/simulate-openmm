from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from simulate.parse._options import _Options
from .position_options import FileOptions
from .velocity_options import SetVelocitiesToTemperatureOptions
from .ensemble_options import NVEOptions, NVTOptions, NPTOptions, RNEMDOptions


class SimulationEnsembles(_Options):

    SECTION_NAME = "SimulationEnsembles"

    # =========================================================================

    def __init__(self):
        super(SimulationEnsembles, self).__init__()
        self.position_options = None
        self.velocity_options = None
        self.ensembles = []

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.position_options is None:
            self._incomplete_error('initialPositions')

    # =========================================================================

    POSITION_METHODS = {'File': FileOptions}

    VELOCITY_METHODS = {'SetVelocitiesToTemperature': SetVelocitiesToTemperatureOptions}

    ENSEMBLE_METHODS = {'NVE': NVEOptions,
                        'NVT': NVTOptions,
                        'NPT': NPTOptions,
                        'RNEMD': RNEMDOptions}

    def _parse_initial_positions(self, *args):
        option_value = args[0]
        line_deque = args[1]
        position_options = self.POSITION_METHODS[option_value]()
        position_options.parse(line_deque.popleft())
        self.position_options = position_options

    def _parse_initial_velocities(self, *args):
        option_value = args[0]
        line_deque = args[1]
        velocity_options = self.VELOCITY_METHODS[option_value]()
        velocity_options.parse(line_deque.popleft())
        self.velocity_options = velocity_options

    def _parse_ensemble(self, *args):
        ensemble_name = args[0]
        line_deque = args[1]
        ensemble_options = self.ENSEMBLE_METHODS[ensemble_name]()
        ensemble_options.parse(line_deque.popleft())
        self.ensembles.append(ensemble_options)

    OPTIONS = {'initialPositions': _parse_initial_positions,
               'initialVelocities': _parse_initial_velocities}

    SECTIONS = {'NVE': _parse_ensemble,
                'NVT': _parse_ensemble,
                'NPT': _parse_ensemble,
                'RNEMD': _parse_ensemble}

    # =========================================================================

    def set_positions(self, simulation):
        self.position_options.set_positions(simulation)

    def set_velocities(self, simulation):
        if self.velocity_options is not None:
            self.velocity_options.set_velocities(simulation)
