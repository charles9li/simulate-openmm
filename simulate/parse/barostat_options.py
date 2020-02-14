from ast import literal_eval

from simtk.openmm import MonteCarloBarostat
from simtk.unit import bar, kelvin

from ._options import _Options


class MonteCarloBarostatOptions(_Options):

    SECTION_NAME = "MonteCarloBarostat"

    # =========================================================================

    def __init__(self):
        super(MonteCarloBarostatOptions, self).__init__()
        self.defaultTemperature = None
        self.defaultPressure = None
        self.frequency = 25

    def _no_option_specified_exception(self, option_name):
        raise ValueError("No {} specified for {}".format(option_name, self.SECTION_NAME))

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
        if self.defaultTemperature is None:
            self._no_option_specified_exception("temperature")
        if self.defaultPressure is None:
            self._no_option_specified_exception("pressure")
        return MonteCarloBarostat(self.defaultPressure, self.defaultTemperature, self.frequency)
