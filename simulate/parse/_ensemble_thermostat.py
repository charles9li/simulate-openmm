"""
_ensemble_thermostat.py: Parses thermostat options for an ensemble;

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

from simtk.unit import kelvin, picosecond

from ._options import _Options

__all__ = ['AndersenThermostatOptions']


class AndersenThermostatOptions(_Options):

    _SECTION_NAME = "AndersenThermostat"

    # =========================================================================

    def __init__(self):
        super(AndersenThermostatOptions, self).__init__()
        self.defaultTemperature = None
        self.defaultCollisionFrequency = 1.0/picosecond

    def _create_options(self):
        super(AndersenThermostatOptions, self)._create_options()
        self._OPTIONS['defaultTemperature'] = self._parse_default_temperature
        self._OPTIONS['defaultCollisionFrequency'] = self._parse_default_collision_frequency

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.defaultTemperature is None:
            self._incomplete_error('defaultTemperature')
        if self.defaultCollisionFrequency is None:
            self._incomplete_error('defaultCollisionFrequency')

    # =========================================================================

    def _parse_default_temperature(self, *args):
        self.defaultTemperature = literal_eval(args[0])*kelvin

    def _parse_default_collision_frequency(self, *args):
        self.defaultCollisionFrequency = literal_eval(args[0])/picosecond

    # =========================================================================

    def thermostat(self):
        from simtk.openmm import AndersenThermostat
        return AndersenThermostat(self.defaultTemperature, self.defaultCollisionFrequency)
