"""
_system_topology_chain.py: Parses chain options for the DodecaneAcrylateTopology
section.

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
import os

import numpy as np
from simtk.openmm.app import Element, PDBFile, Topology

from ._options import _Options


class _StructureParser(object):

    def __init__(self, structure_string):
        self._structure_string = "".join(structure_string.split())
        self._index = 0
        tokenized_string = self._parse_structure_string(0)
        self._chain_sequence = []
        self._parse_tokenized_string(tokenized_string)

    def __str__(self):
        return self._structure_string

    def get_chain_sequence(self):
        return self._chain_sequence

    def _parse_structure_string(self, level):
        tokenized_string = []
        while self._index < len(self._structure_string):
            char = self._structure_string[self._index]
            if char.isdigit():
                tokenized_string.append(self._parse_multiplier())
            elif char == 'm' or char == 'A':
                tokenized_string.append(self._parse_monomer())
            elif char == '+':
                self._index += 1
            elif char == '(':
                self._index += 1
                tokenized_string.append(self._parse_structure_string(1))
            elif char == ')':
                level -= 1
                break
            else:
                raise ValueError("Invalid structure.")
        if level != 0:
            raise ValueError("Mismatched parentheses.")
        return tokenized_string

    def _parse_multiplier(self):
        multiplier_str = ""
        while self._index < len(self._structure_string):
            char = self._structure_string[self._index]
            if char.isdigit():
                multiplier_str += char
                self._index += 1
            elif char == '*':
                multiplier = literal_eval(multiplier_str)
                self._index += 1
                return multiplier
            else:
                raise ValueError("Invalid structure.")

    def _parse_monomer(self):

        # Parse monomer type
        monomer_str = ""
        if self._structure_string[self._index] == 'A':
            monomer_str += self._structure_string[self._index]
            self._index += 1
        elif self._structure_string[self._index:self._index + 2] == 'mA':
            monomer_str += self._structure_string[self._index:self._index + 2]
            self._index += 2
        else:
            raise ValueError("Invalid structure.")

        # Parse end chain length
        while self._index < len(self._structure_string):
            char = self._structure_string[self._index]
            if char.isdigit():
                monomer_str += char
                self._index += 1
            else:
                break
        return monomer_str

    def _parse_tokenized_string(self, tokenized_string):
        if isinstance(tokenized_string, str):
            self._chain_sequence.append(tokenized_string)
        elif isinstance(tokenized_string, list):
            i = 0
            multiplier = 1
            while i < len(tokenized_string):
                token = tokenized_string[i]
                i += 1
                if isinstance(token, int):
                    multiplier *= token
                elif isinstance(token, str) or isinstance(token, str):
                    for _ in range(multiplier):
                        self._parse_tokenized_string(token)
                    multiplier = 1
                else:
                    raise ValueError("Invalid sequence.")


class ChainOptions(_Options):

    _SECTION_NAME = "Chain"

    # =========================================================================

    BORON = Element.getBySymbol('B')
    CARBON = Element.getBySymbol('C')
    NITROGEN = Element.getBySymbol('N')
    OXYGEN = Element.getBySymbol('O')

    # =========================================================================

    def __init__(self):
        super(ChainOptions, self).__init__()
        self.id = None
        self.num = None
        self.sequence = None
        self.create_pdb = True
        self.overwrite_pdb = True
        self.instructions = None
        self.sequence_str = None

    def _create_options(self):
        super(ChainOptions, self)._create_options()
        self._OPTIONS['id'] = self._parse_id
        self._OPTIONS['num'] = self._parse_num
        self._OPTIONS['sequence'] = self._parse_sequence
        self._OPTIONS['createPDB'] = self._parse_create_pdb
        self._OPTIONS['overwritePDB'] = self._parse_overwrite_pdb
        self._OPTIONS['instructions'] = self._parse_instructions

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.num is None:
            self._incomplete_error('num')
        if self.sequence is None:
            self._incomplete_error('structure')

    # =========================================================================

    def _parse_id(self, *args):
        self.id = args[0]

    def _parse_num(self, *args):
        self.num = literal_eval(args[0])

    def _parse_sequence(self, *args):
        structure_parser = _StructureParser(args[0])
        self.sequence = structure_parser.get_chain_sequence()
        self.sequence_str = str(structure_parser)

    def _parse_create_pdb(self, *args):
        self.create_pdb = literal_eval(args[0])

    def _parse_overwrite_pdb(self, *args):
        self.overwrite_pdb = literal_eval(args[0])

    def _parse_instructions(self, *args):
        self.instructions = [instruction.strip() for instruction in args[0].split('/')]

    # =========================================================================

    def add_chain_to_topology(self, topology):

        # Map chain id to sequence
        id_to_sequence = {}

        # Add specified number of chains
        for _ in range(self.num):

            # Initialize topology and positions array for creating chain pdb
            topology_pdb = Topology()
            positions = []

            # Determine chain pdb
            if self.id is not None:
                chain_id = "{}-{}".format(topology.getNumChains() + 1, self.id)
            else:
                chain_id = "{}-{}".format(topology.getNumChains() + 1, self.sequence_str)

            # Create chain
            chain = topology.addChain(id=chain_id)
            chain_pdb = topology_pdb.addChain(id=chain_id)

            # Add chain id and sequence to dictionary
            id_to_sequence[chain.id] = self.sequence_str

            # Initialize atom from previous residue
            prev_res_atom = None
            prev_res_atom_pdb = None

            # Initialize pos of atom on previous residue
            prev_res_atom_pos = np.array([0.0, 0.0, 0.0])

            # Iterate through sequence to add residues to chain
            for j in range(len(self.sequence)):

                # String representation of monomer
                monomer = self.sequence[j]

                # Determine whether or not the residue is a terminal one
                left_ter = False
                right_ter = False
                if j == 0:
                    left_ter = True
                if j == len(self.sequence) - 1:
                    right_ter = True

                # Add residue to chain
                prev_res_atom, prev_res_atom_pdb, prev_res_atom_pos = self._add_residue_to_chain(topology, topology_pdb,
                                                                                                 chain, chain_pdb,
                                                                                                 prev_res_atom,
                                                                                                 prev_res_atom_pdb,
                                                                                                 positions,
                                                                                                 prev_res_atom_pos,
                                                                                                 monomer,
                                                                                                 left_ter, right_ter)

            # Create pdb file for chain
            if self.create_pdb:
                self._create_chain_pdb(topology_pdb, positions)

        return id_to_sequence

    def _add_residue_to_chain(self, topology, topology_pdb, chain, chain_pdb, prev_res_atom, prev_res_atom_pdb,
                              positions, prev_res_atom_pos, monomer, left_ter=False, right_ter=False):

        def _deg_to_rad(ang_deg):
            return ang_deg * 2 * np.pi / 360.0

        # Sines and cosines
        sin19 = np.sin(_deg_to_rad(19.5))
        cos19 = np.cos(_deg_to_rad(19.5))
        sin30 = np.sin(_deg_to_rad(30.0))
        cos30 = np.cos(_deg_to_rad(30.0))
        sin54 = np.sin(_deg_to_rad(54.75))
        cos54 = np.cos(_deg_to_rad(54.75))

        # Determine monomer type
        if monomer.startswith('mA'):
            monomer_type = 'mA'
            is_methyl = True
        else:
            monomer_type = 'A'
            is_methyl = False

        # Determine chain length
        end_chain_length = literal_eval(monomer.replace(monomer_type, ''))

        # Determine residue id
        if left_ter:
            residue_id = "LEFTTER"
        elif right_ter:
            residue_id = "RIGHTTER"
        else:
            residue_id = None

        # Add residue to topology
        residue = topology.addResidue(monomer, chain, id=residue_id)
        residue_pdb = topology_pdb.addResidue(monomer, chain_pdb, id=residue_id)

        # Add first two carbons
        if left_ter and right_ter:
            carbon_element = self.NITROGEN
            carbon_1_element = self.NITROGEN
        elif left_ter:
            carbon_element = self.BORON
            carbon_1_element = self.CARBON
        elif right_ter:
            carbon_element = self.NITROGEN
            carbon_1_element = self.CARBON
        else:
            carbon_element = self.CARBON
            carbon_1_element = self.CARBON
        carbon_pos = prev_res_atom_pos + 1.54*np.array([cos19*cos30, -cos19*sin30, sin19])
        carbon, carbon_pdb = self._add_atom_to_topology('C', carbon_element, residue, residue_pdb,
                                                        topology, topology_pdb,
                                                        positions, carbon_pos)
        carbon_1_pos = carbon_pos + 1.54*np.array([cos19*cos30, cos19*sin30, -sin19])
        carbon_1, carbon_1_pdb = self._add_atom_to_topology('C1', carbon_1_element, residue, residue_pdb,
                                                            topology, topology_pdb,
                                                            positions, carbon_1_pos)

        # Add bond to previous residue if applicable
        if prev_res_atom is not None:
            self._add_bond_to_topology(carbon, carbon_pdb, prev_res_atom, prev_res_atom_pdb, topology, topology_pdb)

        # Add bond between first two carbons
        self._add_bond_to_topology(carbon, carbon_pdb, carbon_1, carbon_1_pdb, topology, topology_pdb)

        # Add methyl group if methacrylate monomer
        if is_methyl:
            carbon_methyl_pos = carbon_1_pos + 1.54*np.array([0.0, 1.0, 0.0])
            carbon_methyl, carbon_methyl_pdb = self._add_atom_to_topology('Cm', self.CARBON, residue, residue_pdb,
                                                                          topology, topology_pdb,
                                                                          positions, carbon_methyl_pos)
            self._add_bond_to_topology(carbon_1, carbon_1_pdb, carbon_methyl, carbon_methyl_pdb,
                                       topology, topology_pdb)

        # Add ester group atoms
        carbon_2_pos = carbon_1_pos + 1.52*np.array([0.0, 0.0, -1.0])
        carbon_2, carbon_2_pdb = self._add_atom_to_topology('C2', self.CARBON, residue, residue_pdb,
                                                            topology, topology_pdb,
                                                            positions, carbon_2_pos)
        self._add_bond_to_topology(carbon_1, carbon_1_pdb, carbon_2, carbon_2_pdb, topology, topology_pdb)
        oxygen_carbonyl_pos = carbon_2_pos + 1.20*np.array([0.0, cos30, -sin30])
        oxygen_carbonyl, oxygen_carbonyl_pdb = self._add_atom_to_topology('O', self.OXYGEN, residue, residue_pdb,
                                                                          topology, topology_pdb,
                                                                          positions, oxygen_carbonyl_pos)
        self._add_bond_to_topology(carbon_2, carbon_2_pdb, oxygen_carbonyl, oxygen_carbonyl_pdb, topology, topology_pdb)
        oxygen_ether_pos = carbon_2_pos + 1.344*np.array([0.0, -cos30, -sin30])
        oxygen_ether, oxygen_ether_pdb = self._add_atom_to_topology('O', self.OXYGEN, residue, residue_pdb,
                                                                    topology, topology_pdb,
                                                                    positions, oxygen_ether_pos)
        self._add_bond_to_topology(carbon_2, carbon_2_pdb, oxygen_ether, oxygen_ether_pdb, topology, topology_pdb)

        # Add carbon chain
        prev_atom_pos = oxygen_ether_pos
        prev_atom = oxygen_ether
        prev_atom_pdb = oxygen_ether_pdb
        for i in range(end_chain_length):
            if i == 0:
                curr_atom_pos = prev_atom_pos + 1.41*np.array([0.0, (-1)**(i+1)*cos54, -sin54])
            else:
                curr_atom_pos = prev_atom_pos + 1.54*np.array([0.0, (-1)**(i+1)*cos54, -sin54])
            curr_atom, curr_atom_pdb = self._add_atom_to_topology('C{}'.format(i + 3), self.CARBON,
                                                                  residue, residue_pdb,
                                                                  topology, topology_pdb,
                                                                  positions, curr_atom_pos)
            self._add_bond_to_topology(prev_atom, prev_atom_pdb, curr_atom, curr_atom_pdb, topology, topology_pdb)
            prev_atom_pos = curr_atom_pos
            prev_atom = curr_atom
            prev_atom_pdb = curr_atom_pdb

        return carbon_1, carbon_1_pdb, carbon_1_pos

    @staticmethod
    def _add_atom_to_topology(name, element, residue, residue_pdb, topology, topology_pdb,
                              positions, atom_pos):

        # Add atom to topology
        atom = topology.addAtom(name, element, residue)
        atom_pdb = topology_pdb.addAtom(name, element, residue_pdb)

        # Add position to array
        positions.append(atom_pos)

        # Return atoms and position
        return atom, atom_pdb

    @staticmethod
    def _add_bond_to_topology(atom1, atom1_pdb, atom2, atom2_pdb, topology, topology_pdb):
        topology.addBond(atom1, atom2)
        topology_pdb.addBond(atom1_pdb, atom2_pdb)

    def _create_chain_pdb(self, topology_pdb, positions):
        dirname = os.path.dirname(__file__)
        file_path = os.path.join(dirname, "data/{}.pdb".format(self.sequence_str))
        if not os.path.isfile(file_path) or self.overwrite_pdb:
            self.overwrite_pdb = False
            PDBFile.writeFile(topology_pdb, positions, open(file_path, 'w'))
