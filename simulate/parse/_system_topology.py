"""
_system_topology.py: Parses topology information for a system.

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

import os
from ast import literal_eval

import numpy as np
from simtk.openmm.app import ForceField, NoCutoff, Topology
from simtk.unit import nanometer
from parmed import gromacs

from ._options import _Options
from ._system_topology_chain import ChainOptions


class _TopologyOptions(_Options):

    # =========================================================================

    _SECTION_NAME = '_TopologyOptions'

    # =========================================================================

    def __init__(self):
        super(_TopologyOptions, self).__init__()
        self._topology = None

    # =========================================================================

    def topology(self):
        return self._topology

    def create_system(self, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*nanometer,
                      constraints=None, rigidWater=True, implicitSolvent=None,
                      soluteDielectric=1.0, solventDielectric=78.5,
                      ewaldErrorTolerance=0.0005, removeCMMotion=True,
                      hydrogenMass=None):
        pass


class AmberTopologyOptions(_TopologyOptions):

    _SECTION_NAME = 'AmberTopologyOptions'

    def __init__(self):
        super(AmberTopologyOptions, self).__init__()
        raise NotImplementedError("'{}' is not supported yet.".format(self._SECTION_NAME))


class GromacsTopologyOptions(_TopologyOptions):

    _SECTION_NAME = 'GromacsTopologyOptions'

    # =========================================================================

    def __init__(self):
        super(GromacsTopologyOptions, self).__init__()
        self.topFilename = None
        self.groFilename = None
        self._gromacs_topology = None

    def _create_options(self):
        super(GromacsTopologyOptions, self)._create_options()
        self._OPTIONS['topFilename'] = self._parse_top_filename
        self._OPTIONS['groFilename'] = self._parse_gro_filename

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.topFilename is None:
            self._incomplete_error('topFilename')
        if self.groFilename is None:
            self._incomplete_error('groFilename')

    # =========================================================================

    def _parse_top_filename(self, *args):
        self.topFilename = args[0]

    def _parse_gro_filename(self, *args):
        self.groFilename = args[0]

    # =========================================================================

    def topology(self):
        self._create_gromacs_topology()
        return self._gromacs_topology.topology

    def _create_gromacs_topology(self):
        if self._gromacs_topology is None:
            gro = gromacs.GromacsGroFile.parse(self.groFilename)
            self._gromacs_topology = gromacs.GromacsTopologyFile(self.topFilename)
            self._gromacs_topology.box = gro.box

    def create_system(self, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*nanometer,
                      constraints=None, rigidWater=True, implicitSolvent=None,
                      soluteDielectric=1.0, solventDielectric=78.5,
                      ewaldErrorTolerance=0.0005, removeCMMotion=True,
                      hydrogenMass=None):
        self._create_gromacs_topology()
        return self._gromacs_topology.createSystem(nonbondedMethod=nonbondedMethod,
                                                   nonbondedCutoff=nonbondedCutoff,
                                                   constraints=constraints, rigidWater=rigidWater, implicitSolvent=implicitSolvent,
                                                   soluteDielectric=soluteDielectric, solventDielectric=solventDielectric,
                                                   ewaldErrorTolerance=ewaldErrorTolerance, removeCMMotion=removeCMMotion,
                                                   hydrogenMass=hydrogenMass)
            
            
class DodecaneAcrylateTopologyOptions(_TopologyOptions):

    # =========================================================================

    _SECTION_NAME = "DodecaneAcrylateTopologyOptions"

    # =========================================================================

    # Paths to forcefield files

    _data_directory = os.path.join(os.path.dirname(__file__), 'data')

    TRAPPEUA_FF_PATH = os.path.join(_data_directory, "trappeua-acrylates.xml")

    # =========================================================================
    
    def __init__(self):
        super(DodecaneAcrylateTopologyOptions, self).__init__()
        self.force_field = ForceField(self.TRAPPEUA_FF_PATH)
        self.numDodecane = 0
        self.box_vectors = None
        self.chains = []
        self.id_to_sequence = {}

    def _create_options(self):
        super(DodecaneAcrylateTopologyOptions, self)._create_options()
        self._OPTIONS['numDodecane'] = self._parse_num_dodecane
        self._OPTIONS['box'] = self._parse_box

    def _create_sections(self):
        super(DodecaneAcrylateTopologyOptions, self)._create_sections()
        self._SECTIONS['chain'] = self._parse_chain

    # =========================================================================

    def _parse_num_dodecane(self, *args):
        self.numDodecane = literal_eval(args[0])

    def _parse_box(self, *args):
        a, b, c = args[0].split(' ')
        self.box_vectors = np.array([[literal_eval(a), 0.0, 0.0],
                                     [0.0, literal_eval(b), 0.0],
                                     [0.0, 0.0, literal_eval(c)]])*nanometer

    def _parse_chain(self, *args):
        line_deque = args[1]
        chain_options = ChainOptions()
        chain_options.parse(line_deque.popleft())
        self.chains.append(chain_options)

    # =========================================================================

    def topology(self):
        self._create_dodecane_acrylate_topology()
        return self._topology

    def _create_dodecane_acrylate_topology(self):
        if self._topology is None:
            topology = Topology()
            if self.box_vectors is not None:
                topology.setPeriodicBoxVectors(self.box_vectors)
            for chain_option in self.chains:
                id_to_sequence = chain_option.add_chain_to_topology(topology)
                self.id_to_sequence.update(id_to_sequence)
            self._topology = topology

    def create_system(self, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*nanometer,
                      constraints=None, rigidWater=True, implicitSolvent=None,
                      soluteDielectric=1.0, solventDielectric=78.5,
                      ewaldErrorTolerance=0.0005, removeCMMotion=True,
                      hydrogenMass=None):
        self._create_dodecane_acrylate_topology()
        return self.force_field.createSystem(self._topology, nonbondedMethod=nonbondedMethod,
                                             nonbondedCutoff=nonbondedCutoff,
                                             constraints=constraints, rigidWater=rigidWater,
                                             implicitSolvent=implicitSolvent, soluteDielectric=soluteDielectric,
                                             solventDielectric=solventDielectric, ewaldErrorTolerance=ewaldErrorTolerance,
                                             removeCMMotion=removeCMMotion, hydrogenMass=hydrogenMass)
