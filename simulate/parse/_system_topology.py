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
from simtk.openmm.app import Element, ForceField, NoCutoff, Topology
from simtk.unit import nanometer

from ._options import _Options
from ._system_topology_chain import ChainOptions, BranchedChainOptions


class _TopologyOptions(_Options):

    # =========================================================================

    _SECTION_NAME = '_Topology'

    # =========================================================================

    def __init__(self, system_options):
        super(_TopologyOptions, self).__init__()
        self.system_options = system_options
        self._topology = None

    # =========================================================================

    # Helper functions for parsing input

    def _create_filepath(self, filepath):
        directory = self.system_options.input_options.directory
        return os.path.join(directory, filepath)

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

    _SECTION_NAME = 'AmberTopology'

    def __init__(self, system_options):
        super(AmberTopologyOptions, self).__init__(system_options)
        raise NotImplementedError("'{}' is not supported yet.".format(self._SECTION_NAME))


class GromacsTopologyOptions(_TopologyOptions):

    _SECTION_NAME = 'GromacsTopology'

    # =========================================================================

    def __init__(self, system_options):
        super(GromacsTopologyOptions, self).__init__(system_options)
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
        self.topFilename = self._create_filepath(args[0])

    def _parse_gro_filename(self, *args):
        self.groFilename = self._create_filepath(args[0])

    # =========================================================================

    def topology(self):
        self._create_gromacs_topology()
        return self._gromacs_topology.topology

    def _create_gromacs_topology(self):
        from parmed import gromacs
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

    _SECTION_NAME = "DodecaneAcrylateTopology"

    # =========================================================================

    # Paths to forcefield files

    data_directory = os.path.join(os.path.dirname(__file__), 'data')

    TRAPPEUA_FF_PATH = os.path.join(data_directory, "trappeua-acrylates.xml")
    OPLS_AA_PATH = os.path.join(data_directory, "opls_aa.xml")

    # =========================================================================
    
    def __init__(self, system_options):
        super(DodecaneAcrylateTopologyOptions, self).__init__(system_options)
        self.forceField_str = "TraPPE-UA"
        self.forceField = ForceField(self.TRAPPEUA_FF_PATH)
        self.numDodecane = 0
        self.numSqualane = 0
        self.dodecaneInstructions = None
        self.box_vectors = None
        self.chains = []
        self.branched_chains = []
        self.id_to_sequence = {}

    def _create_options(self):
        super(DodecaneAcrylateTopologyOptions, self)._create_options()
        self._OPTIONS['forceField'] = self._parse_force_field
        self._OPTIONS['numDodecane'] = self._parse_num_dodecane
        self._OPTIONS['numSqualane'] = self._parse_num_squalane
        self._OPTIONS['dodecaneInstructions'] = self._parse_dodecane_instructions
        self._OPTIONS['box'] = self._parse_box

    def _create_sections(self):
        super(DodecaneAcrylateTopologyOptions, self)._create_sections()
        self._SECTIONS['chain'] = self._parse_chain
        self._SECTIONS['BranchedChain'] = self._parse_branched_chain

    # =========================================================================

    def _parse_force_field(self, *args):
        if args[0] == 'TraPPE-UA':
            self.forceField = ForceField(self.TRAPPEUA_FF_PATH)
        elif args[0] == 'OPLS-AA':
            self.forceField = ForceField(self.OPLS_AA_PATH)
        else:
            raise ValueError("Invalid force field.")
        self.forceField_str = args[0]

    def _parse_num_dodecane(self, *args):
        self.numDodecane = literal_eval(args[0])

    def _parse_num_squalane(self, *args):
        self.numSqualane = literal_eval(args[0])

    def _parse_dodecane_instructions(self, *args):
        self.dodecaneInstructions = [instruction.strip() for instruction in args[0].split('/')]

    def _parse_box(self, *args):
        a, b, c = args[0].split(' ')
        self.box_vectors = np.array([[literal_eval(a), 0.0, 0.0],
                                     [0.0, literal_eval(b), 0.0],
                                     [0.0, 0.0, literal_eval(c)]])*nanometer

    def _parse_chain(self, *args):
        line_deque = args[1]
        chain_options = ChainOptions(self)
        chain_options.parse(line_deque.popleft())
        self.chains.append(chain_options)

    def _parse_branched_chain(self, *args):
        line_deque = args[1]
        branched_chain_options = BranchedChainOptions(self)
        branched_chain_options.parse(line_deque.popleft())
        self.branched_chains.append(branched_chain_options)

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
            for branched_chain_option in self.branched_chains:
                branched_chain_option.add_chain_to_topology(topology)
            for _ in range(self.numDodecane):
                dodecane_id = self._add_dodecane_to_topology(topology)
                self.id_to_sequence[dodecane_id] = "C12"
            for _ in range(self.numSqualane):
                squalane_id = self._add_squalane_to_topology(topology)
                self.id_to_sequence[squalane_id] = "squalane"
            self._topology = topology

    def _add_dodecane_to_topology(self, topology):

        # Carbon element
        carbon_element = Element.getBySymbol('C')
        hydrogen_element = Element.getBySymbol('H')

        chain = topology.addChain("{}-C12".format(topology.getNumChains() + 1))
        residue = topology.addResidue("C12", chain)
        prev_atom = topology.addAtom("C", carbon_element, residue)
        if self.forceField_str == "TraPPE-UA":
            for i in range(11):
                curr_atom = topology.addAtom("C{}".format(i + 1), carbon_element, residue)
                topology.addBond(prev_atom, curr_atom)
                prev_atom = curr_atom
        else:
            H_counter = 0
            for _ in range(3):
                H = topology.addAtom("H{}".format(H_counter), hydrogen_element, residue)
                topology.addBond(H, prev_atom)
            for i in range(11):
                curr_atom = topology.addAtom("C{}".format(i + 1), carbon_element, residue)
                topology.addBond(prev_atom, curr_atom)
                for _ in range(2):
                    H = topology.addAtom("H{}".format(H_counter), hydrogen_element, residue)
                    topology.addBond(H, curr_atom)
                    prev_atom = curr_atom
            H = topology.addAtom("H{}".format(H_counter), hydrogen_element, residue)
            topology.addBond(H, prev_atom)

        return chain.id

    def _add_squalane_to_topology(self, topology):

        # Carbon element
        carbon_element = Element.getBySymbol('C')
        hydrogen_element = Element.getBySymbol('H')

        chain = topology.addChain("{}-squalane".format(topology.getNumChains() + 1))
        residue = topology.addResidue("squalane", chain)
        prev_atom = None
        if self.forceField_str == "TraPPE-UA":
            atom_index = 0
            for _ in range(3):
                C1 = topology.addAtom("C{}".format(atom_index), carbon_element, residue)
                atom_index += 1
                if prev_atom is not None:
                    topology.addBond(prev_atom, C1)
                C2 = topology.addAtom("C{}".format(atom_index), carbon_element, residue)
                atom_index += 1
                topology.addBond(C1, C2)
                C3 = topology.addAtom("C{}".format(atom_index), carbon_element, residue)
                atom_index += 1
                topology.addBond(C2, C3)
                C4 = topology.addAtom("C{}".format(atom_index), carbon_element, residue)
                atom_index += 1
                topology.addBond(C2, C4)
                C5 = topology.addAtom("C{}".format(atom_index), carbon_element, residue)
                atom_index += 1
                topology.addBond(C4, C5)
                prev_atom = C5
            for _ in range(3):
                C1 = topology.addAtom("C{}".format(atom_index), carbon_element, residue)
                atom_index += 1
                topology.addBond(prev_atom, C1)
                C2 = topology.addAtom("C{}".format(atom_index), carbon_element, residue)
                atom_index += 1
                topology.addBond(C1, C2)
                C3 = topology.addAtom("C{}".format(atom_index), carbon_element, residue)
                atom_index += 1
                topology.addBond(C2, C3)
                C4 = topology.addAtom("C{}".format(atom_index), carbon_element, residue)
                atom_index += 1
                topology.addBond(C3, C4)
                C5 = topology.addAtom("C{}".format(atom_index), carbon_element, residue)
                atom_index += 1
                topology.addBond(C3, C5)
                prev_atom = C5
        else:
            raise NotImplementedError("OPLS-AA not implemented for squalane")

        return chain.id

    def create_system(self, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*nanometer,
                      constraints=None, rigidWater=True, implicitSolvent=None,
                      soluteDielectric=1.0, solventDielectric=78.5,
                      ewaldErrorTolerance=0.0005, removeCMMotion=True,
                      hydrogenMass=None):
        self._create_dodecane_acrylate_topology()
        return self.forceField.createSystem(self._topology, nonbondedMethod=nonbondedMethod,
                                            nonbondedCutoff=nonbondedCutoff,
                                            constraints=constraints, rigidWater=rigidWater,
                                            implicitSolvent=implicitSolvent, soluteDielectric=soluteDielectric,
                                            solventDielectric=solventDielectric, ewaldErrorTolerance=ewaldErrorTolerance,
                                            removeCMMotion=removeCMMotion, hydrogenMass=hydrogenMass)
