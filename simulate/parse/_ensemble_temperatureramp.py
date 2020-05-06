"""
_ensemble_temperatureramp.py: Parses options for temperature ramping.

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

__all__ = ['TemperatureRampOptions']


class TemperatureRampOptions(_Options):

    _SECTION_NAME = "TemperatureRamp"

    # =========================================================================

    def __init__(self):
        super(TemperatureRampOptions, self).__init__()
        self.startTemperature = None
        self.endTemperature = None
        self.stepTemperature = None
        self.stepsPerInterval = None

    def _create_options(self):
        super(TemperatureRampOptions, self)._create_options()
        self._OPTIONS['startTemperature'] = self._parse_start_temperature
        self._OPTIONS['endTemperature'] = self._parse_end_temperature
        self._OPTIONS['stepTemperature'] = self._parse_step_temperature
        self._OPTIONS['stepsPerInterval'] = self._parse_steps_per_interval

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.startTemperature is None:
            self._incomplete_error('startTemperature')
        if self.endTemperature is None:
            self._incomplete_error('endTemperature')
        if self.stepTemperature is None:
            self._incomplete_error('stepTemperature')
        if self.stepsPerInterval is None:
            self._incomplete_error('stepsPerInterval')
        if (self.endTemperature - self.startTemperature)*self.stepTemperature < 0.0*kelvin**2:
            raise ValueError("Final temperature will never be reached with the temperature step.")

    # =========================================================================

    def _parse_start_temperature(self, *args):
        self.startTemperature = literal_eval(args[0])*kelvin

    def _parse_end_temperature(self, *args):
        self.endTemperature = literal_eval(args[0])*kelvin

    def _parse_step_temperature(self, *args):
        self.stepTemperature = literal_eval(args[0])*kelvin

    def _parse_steps_per_interval(self, *args):
        self.stepsPerInterval = literal_eval(args[0])
