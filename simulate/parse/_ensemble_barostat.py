"""
_ensemble_barostat.py: Parses barostat options for an ensemble;

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

from simtk.unit import bar, kelvin

from ._options import _Options

__all__ = ['MonteCarloBarostatOptions']


class MonteCarloBarostatOptions(_Options):

    _SECTION_NAME = "MonteCarloBarostat"

    # =========================================================================

    def __init__(self):
        super(MonteCarloBarostatOptions, self).__init__()
        self.defaultTemperature = None
        self.defaultPressure = None
        self.frequency = 25

    def _create_options(self):
        super(MonteCarloBarostatOptions, self)._create_options()
        self._OPTIONS['defaultTemperature'] = self._parse_default_temperature
        self._OPTIONS['defaultPressure'] = self._parse_default_pressure
        self._OPTIONS['frequency'] = self._parse_frequency

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.defaultTemperature is None:
            self._incomplete_error('defaultTemperature')
        if self.defaultPressure is None:
            self._incomplete_error('defaultPressure')

    # =========================================================================

    def _parse_default_temperature(self, *args):
        self.defaultTemperature = literal_eval(args[0])*kelvin

    def _parse_default_pressure(self, *args):
        self.defaultPressure = literal_eval(args[0])*bar

    def _parse_frequency(self, *args):
        self.frequency = literal_eval(args[0])

    # =========================================================================

    def barostat(self):
        from simtk.openmm import MonteCarloBarostat
        return MonteCarloBarostat(self.defaultPressure, self.defaultTemperature, self.frequency)


class MonteCarloAnisotropicBarostatOptions(MonteCarloBarostatOptions):

    _SECTION_NAME = "MonteCarloAnisotropicBarostat"

    # =========================================================================

    def __init__(self):
        super(MonteCarloAnisotropicBarostatOptions, self).__init__()
        self.scaleX = True
        self.scaleY = True
        self.scaleZ = True

    def _create_options(self):
        super(MonteCarloAnisotropicBarostatOptions, self)._create_options()
        self._OPTIONS['scaleX'] = self._parse_scale_x
        self._OPTIONS['scaleY'] = self._parse_scale_y
        self._OPTIONS['scaleZ'] = self._parse_scale_z

    # =========================================================================

    def _parse_scale_x(self, *args):
        self.scaleX = literal_eval(args[0])

    def _parse_scale_y(self, *args):
        self.scaleY = literal_eval(args[0])

    def _parse_scale_z(self, *args):
        self.scaleZ = literal_eval(args[0])

    # =========================================================================

    def barostat(self):
        from simtk.openmm import MonteCarloAnisotropicBarostat, Vec3
        pressure_vec = Vec3(self.defaultPressure, self.defaultPressure, self.defaultPressure)
        return MonteCarloAnisotropicBarostat(pressure_vec, self.defaultPressure,
                                             self.scaleX, self.scaleY, self.scaleZ, self.frequency)
