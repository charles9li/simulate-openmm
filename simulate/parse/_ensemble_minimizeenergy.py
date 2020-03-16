"""
_ensemble_minimizeenergy.py: Parses energy minimization options for an ensemble.

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

from simtk.unit import kilojoule_per_mole

from ._options import _Options


class MinimizeEnergyOptions(_Options):

    _SECTION_NAME = "MinimizeEnergyOptions"

    # =========================================================================
    
    def __init__(self):
        super(MinimizeEnergyOptions, self).__init__()
        self.tolerance = 10*kilojoule_per_mole
        self.maxIterations = 0
        self.PDBFile = None

    def _create_options(self):
        super(MinimizeEnergyOptions, self)._create_options()
        self._OPTIONS['tolerance'] = self._parse_tolerance
        self._OPTIONS['maxIterations'] = self._parse_max_iterations
        self._OPTIONS['PDBFile'] = self._parse_PDB_file

    # =========================================================================

    def _parse_tolerance(self, *args):
        self.tolerance = literal_eval(args[0])*kilojoule_per_mole

    def _parse_max_iterations(self, *args):
        self.maxIterations = literal_eval(args[0])

    def _parse_PDB_file(self, *args):
        self.PDBFile = args[0]