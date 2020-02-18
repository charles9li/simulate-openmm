from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from ast import literal_eval

from simtk.unit import kelvin

from simulate.parse._options import _Options


class _VelocityOptions(_Options):

    def __init__(self):
        super(_VelocityOptions, self).__init__()

    def set_velocities(self, simulation):
        pass


class SetVelocitiesToTemperatureOptions(_VelocityOptions):

    SECTION_NAME = 'SetVelocitiesToTemperature'

    def __init__(self):
        super(SetVelocitiesToTemperatureOptions, self).__init__()
        self.temperature = None
        self.randomSeed = None

    def _check_for_incomplete_input(self):
        if self.temperature is None:
            self._incomplete_error('temperature')

    def _parse_temperature(self, *args):
        self.temperature = literal_eval(args[0])*kelvin

    OPTIONS = {'temperature': _parse_temperature}

    def set_velocities(self, simulation):
        if self.randomSeed is None:
            simulation.context.setVelocitiesToTemperature(self.temperature)
        else:
            simulation.context.setVelocitiesToTemperature(self.temperature, self.randomSeed)
