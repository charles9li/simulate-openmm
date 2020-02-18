from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from ast import literal_eval

from simtk.unit import kilojoule_per_mole

from simulate.parse._options import _Options


class MinimizeEnergyOptions(_Options):
    
    def __init__(self):
        super(MinimizeEnergyOptions, self).__init__()
        self.tolerance = 10*kilojoule_per_mole
        self.maxIterations = 0

    def _parse_tolerance(self, *args):
        self.tolerance = literal_eval(args[0])*kilojoule_per_mole

    def _parse_max_iterations(self, *args):
        self.maxIterations = literal_eval(args[0])

    OPTIONS = {'tolerance': _parse_tolerance,
               'maxIterations': _parse_max_iterations}
