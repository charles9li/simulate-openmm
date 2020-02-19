from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from ast import literal_eval

from simulate.parse._options import _Options


class AverageOptions(_Options):

    SECTION_NAME = "AverageOptions"
    
    def __init__(self):
        super(AverageOptions, self).__init__()
        self.volume = False
        self.energy = False

    def _parse_volume(self, *args):
        self.volume = literal_eval(args[0])

    def _parse_energy(self, *args):
        self.energy = literal_eval(args[0])

    OPTIONS = {'volume': _parse_volume,
               'energy': _parse_energy}
