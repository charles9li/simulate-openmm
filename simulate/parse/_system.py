"""
_system.py: Parses options for the system.

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

from simtk.openmm import app, NonbondedForce
from simtk.unit import nanometer

from ._options import _Options
from ._system_topology import GromacsTopologyOptions, DodecaneAcrylateTopologyOptions


class SystemOptions(_Options):

    _SECTION_NAME = 'system'

    # =========================================================================

    def __init__(self):
        super(SystemOptions, self).__init__()
        self.topology_options = None
        self.nonbondedMethod = app.NoCutoff
        self.nonbondedCutoff = 0.9*nanometer
        self.constraints = None
        self.rigidWater = True
        self.implicitSolvent = None
        self.soluteDielectric = 1.0
        self.solventDielectric = 78.5
        self.ewaldErrorTolerance = 0.0005
        self.removeCMMotion = True
        self.hydrogenMass = None
        self.useDispersionCorrection = True

    def _create_options(self):
        super(SystemOptions, self)._create_options()
        self._OPTIONS['topology'] = self._parse_topology
        self._OPTIONS['nonbondedMethod'] = self._parse_nonbonded_method
        self._OPTIONS['nonbondedCutoff'] = self._parse_nonbonded_cutoff
        self._OPTIONS['constraints'] = self._parse_constraints
        self._OPTIONS['ewaldErrorTolerance'] = self._parse_ewald_error_tolerance
        self._OPTIONS['useDispersionCorrection'] = self._parse_use_dispersion_correction

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.topology_options is None:
            self._incomplete_error('topology')

    # =========================================================================

    NONBONDED_METHODS = {'NoCutoff': app.NoCutoff,
                         'CutoffPeriodic': app.CutoffPeriodic,
                         'CutoffNonPeriodic': app.CutoffNonPeriodic,
                         'Ewald': app.Ewald,
                         'PME': app.PME,
                         'LJPME': app.LJPME}
    TOPOLOGY_OPTIONS = {'GromacsTopology': GromacsTopologyOptions,
                        'DodecaneAcrylateTopology': DodecaneAcrylateTopologyOptions}
    CONSTRAINTS = {'None': None,
                   'HBonds': app.HBonds,
                   'AllBonds': app.AllBonds,
                   'HAngles': app.HAngles}

    def _parse_topology(self, *args):
        option_value = args[0]
        line_deque = args[1]
        topology_options = self.TOPOLOGY_OPTIONS[option_value]()
        topology_options.parse(line_deque.popleft())
        self.topology_options = topology_options

    def _parse_nonbonded_method(self, *args):
        option_value = args[0]
        try:
            self.nonbondedMethod = self.NONBONDED_METHODS[option_value]
        except KeyError:
            raise ValueError("{} is not a valid option for nonbondedMethod.".format(option_value))

    def _parse_nonbonded_cutoff(self, *args):
        self.nonbondedCutoff = literal_eval(args[0])*nanometer

    def _parse_constraints(self, *args):
        self.constraints = self.CONSTRAINTS[args[0]]

    def _parse_ewald_error_tolerance(self, *args):
        self.ewaldErrorTolerance = literal_eval(args[0])

    def _parse_use_dispersion_correction(self, *args):
        self.useDispersionCorrection = literal_eval(args[0])

    # =========================================================================

    def topology(self):
        return self.topology_options.topology()

    def create_system(self):

        # Create system
        system = self.topology_options.create_system(nonbondedMethod=self.nonbondedMethod,
                                                     nonbondedCutoff=self.nonbondedCutoff,
                                                     constraints=self.constraints,
                                                     rigidWater=self.rigidWater,
                                                     implicitSolvent=self.implicitSolvent,
                                                     soluteDielectric=self.soluteDielectric,
                                                     solventDielectric=self.solventDielectric,
                                                     ewaldErrorTolerance=self.ewaldErrorTolerance,
                                                     removeCMMotion=self.removeCMMotion,
                                                     hydrogenMass=self.hydrogenMass)

        # Set dispersion correction
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                force.setUseDispersionCorrection(self.useDispersionCorrection)

        # Set force groups
        for i in range(system.getNumForces()):
            system.getForce(i).setForceGroup(i)

        return system
