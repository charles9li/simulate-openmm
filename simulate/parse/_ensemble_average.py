"""
_ensemble_average.py: Parses averaging options for an ensemble.

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

from ._options import _Options


class AverageOptions(_Options):

    _SECTION_NAME = "AverageOptions"

    # =========================================================================
    
    def __init__(self):
        super(AverageOptions, self).__init__()
        self.volume = False
        self.energy = False

    # =========================================================================

    def _parse_volume(self, *args):
        self.volume = literal_eval(args[0])

    def _parse_energy(self, *args):
        self.energy = literal_eval(args[0])

    _OPTIONS = {'volume': _parse_volume,
                'energy': _parse_energy}
