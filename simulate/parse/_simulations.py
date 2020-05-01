"""
_simulations.py: Parses options for the simulation(s).

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

from ._options import _Options
from ._simulations_positions import FileOptions, SubrandomParticlePositions, DodecaneAcrylatePositionOptions
from ._simulations_velocity import SetVelocitiesToTemperatureOptions
from ._ensemble import NVEOptions, NVTOptions, NPTOptions, RNEMDOptions


class SimulationsOptions(_Options):

    _SECTION_NAME = "SimulationEnsembles"

    # =========================================================================

    def __init__(self, input_options):
        super(SimulationsOptions, self).__init__()
        self.input_options = input_options
        self.position_options = None
        self.velocity_options = None
        self.ensembles = []

    def _create_options(self):
        super(SimulationsOptions, self)._create_options()
        self._OPTIONS['initialPositions'] = self._parse_initial_positions
        self._OPTIONS['initialVelocities'] = self._parse_initial_velocities

    def _create_sections(self):
        super(SimulationsOptions, self)._create_sections()
        self._SECTIONS['NVE'] = self._parse_ensemble
        self._SECTIONS['NVT'] = self._parse_ensemble
        self._SECTIONS['NPT'] = self._parse_ensemble
        self._SECTIONS['RNEMD'] = self._parse_ensemble

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.position_options is None:
            self._incomplete_error('initialPositions')

    # =========================================================================

    POSITION_OPTIONS = {'File': FileOptions,
                        'SubrandomParticlePositions': SubrandomParticlePositions,
                        'DodecaneAcrylatePositions': DodecaneAcrylatePositionOptions}

    VELOCITY_OPTIONS = {'SetVelocitiesToTemperature': SetVelocitiesToTemperatureOptions}

    ENSEMBLE_OPTIONS = {'NVE': NVEOptions,
                        'NVT': NVTOptions,
                        'NPT': NPTOptions,
                        'RNEMD': RNEMDOptions}

    def _parse_initial_positions(self, *args):
        option_value = args[0]
        line_deque = args[1]
        position_options = self.POSITION_OPTIONS[option_value](self)
        position_options.parse(line_deque.popleft())
        self.position_options = position_options

    def _parse_initial_velocities(self, *args):
        option_value = args[0]
        line_deque = args[1]
        velocity_options = self.VELOCITY_OPTIONS[option_value]()
        velocity_options.parse(line_deque.popleft())
        self.velocity_options = velocity_options

    def _parse_ensemble(self, *args):
        ensemble_name = args[0]
        line_deque = args[1]
        ensemble_options = self.ENSEMBLE_OPTIONS[ensemble_name](self)
        ensemble_options.parse(line_deque.popleft())
        self.ensembles.append(ensemble_options)

    # =========================================================================

    def set_positions(self, simulation, *args):
        self.position_options.set_positions(simulation, *args)

    def set_velocities(self, simulation):
        if self.velocity_options is not None:
            self.velocity_options.set_velocities(simulation)
