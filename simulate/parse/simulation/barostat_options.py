from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from ast import literal_eval

from simtk.unit import bar, kelvin

from simulate.parse._options import _Options


class MonteCarloBarostatOptions(_Options):

    SECTION_NAME = "MonteCarloBarostat"

    # =========================================================================

    def __init__(self):
        super(MonteCarloBarostatOptions, self).__init__()
        self.defaultTemperature = None
        self.defaultPressure = None
        self.frequency = 25

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.defaultPressure is None:
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

    OPTIONS = {'defaultTemperature': _parse_default_temperature,
               'defaultPressure': _parse_default_pressure,
               'frequency': _parse_frequency}

    # =========================================================================

    def barostat(self):
        from simtk.openmm import MonteCarloBarostat
        return MonteCarloBarostat(self.defaultPressure, self.defaultTemperature, self.frequency)
