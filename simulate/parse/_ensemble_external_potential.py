"""
_ensemble_external_potential.py: Parses external potential options.

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

from simtk.unit import nanometer
from simtk.openmm import CustomCentroidBondForce
import numpy as np

from ._options import _Options


class SinusoidalOptions(_Options):

    _SECTION_NAME = "Sinusoidal"

    # =========================================================================

    def __init__(self):
        super(SinusoidalOptions, self).__init__()
        self.num_period = 1
        self.period = None
        self.amplitude = 1.0
        self.axis = 'x'
        self.residue = None
        self.indices = None

    def _create_options(self):
        super(SinusoidalOptions, self)._create_options()
        self._OPTIONS['numPeriod'] = self._parse_num_period
        self._OPTIONS['period'] = self._parse_period
        self._OPTIONS['amplitude'] = self._parse_amplitude
        self._OPTIONS['axis'] = self._parse_axis
        self._OPTIONS['residue'] = self._parse_residue
        self._OPTIONS['indices'] = self._parse_indices

        # =========================================================================

    def _check_for_incomplete_input(self):
        if self.residue is None:
            self._incomplete_error('residue')

    # =========================================================================

    def _parse_num_period(self, *args):
        self.num_period = literal_eval(args[0])

    def _parse_period(self, *args):
        self.period = literal_eval(args[0])

    def _parse_amplitude(self, *args):
        self.amplitude = literal_eval(args[0])

    def _parse_axis(self, *args):
        self.axis = args[0]

    def _parse_residue(self, *args):
        self.residue = args[0]

    def _parse_indices(self, *args):
        self.indices = [literal_eval(i) for i in args[0].split()]

    # =========================================================================

    def add_potential_to_system(self, topology, system):
        if self.period is None:
            a, _, _ = system.getDefaultPeriodicBoxVectors()
            axis_index_map = {'x': 0, 'y': 1, 'z': 2}
            period = a[axis_index_map[self.axis]].value_in_unit(nanometer)
        else:
            period = self.period
        energy_function = "{amp}*sin(2*{pi}*{nperiod}*{axis}1/{L})".format(amp=self.amplitude,
                                                                           pi=np.pi,
                                                                           nperiod=self.num_period,
                                                                           axis=self.axis,
                                                                           L=period)
        force = CustomCentroidBondForce(1, energy_function)
        for residue in topology.residues():
            if residue.name == self.residue:
                indices = [atom.index for atom in residue.atoms()]
                particles = []
                if self.indices is None:
                    particles = indices
                else:
                    for i in self.indices:
                        particles.append(indices[i])
                group = force.addGroup(particles)
                force.addBond([group])
        system.addForce(force)
        force.setForceGroup(system.getNumForces() - 1)
        return system
