"""
_simulations_velocity.py: Parses velocity initialization options for simulations.

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

from simtk.unit import kelvin

from ._options import _Options


class _VelocityOptions(_Options):

    def __init__(self):
        super(_VelocityOptions, self).__init__()

    # =========================================================================

    def set_velocities(self, simulation):
        pass


class SetVelocitiesToTemperatureOptions(_VelocityOptions):

    _SECTION_NAME = 'SetVelocitiesToTemperature'

    # =========================================================================

    def __init__(self):
        super(SetVelocitiesToTemperatureOptions, self).__init__()
        self.temperature = None
        self.randomSeed = None

    def _create_options(self):
        super(SetVelocitiesToTemperatureOptions, self)._create_options()
        self._OPTIONS['temperature'] = self._parse_temperature

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.temperature is None:
            self._incomplete_error('temperature')

    # =========================================================================

    def _parse_temperature(self, *args):
        self.temperature = literal_eval(args[0])*kelvin

    # =========================================================================

    def set_velocities(self, simulation):
        if self.randomSeed is None:
            simulation.context.setVelocitiesToTemperature(self.temperature)
        else:
            simulation.context.setVelocitiesToTemperature(self.temperature, self.randomSeed)
