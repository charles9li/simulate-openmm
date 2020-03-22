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
        self.overwrite_pdb = False
        self._sequence_str = None

    def _create_options(self):
        super(ChainOptions, self)._create_options()
        self._OPTIONS['id'] = self._parse_id
        self._OPTIONS['num'] = self._parse_num
        self._OPTIONS['sequence'] = self._parse_sequence
        self._OPTIONS['createPDB'] = self._parse_create_pdb
        self._OPTIONS['overwritePDB'] = self._parse_overwrite_pdb

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
        self._sequence_str = str(structure_parser)

    def _parse_create_pdb(self, *args):
        self.create_pdb = literal_eval(args[0])

    def _parse_overwrite_pdb(self, *args):
        self.overwrite_pdb = literal_eval(args[0])

    # =========================================================================

    def add_chain_to_topology(self, topology):
        id_to_sequence = {}
        for _ in range(self.num):
            if self.id is not None:
                chain = topology.addChain("{}-{}".format(topology.getNumChains() + 1, self.id))
            else:
                chain = topology.addChain("{}-{}".format(topology.getNumChains() + 1, self._sequence_str))
            id_to_sequence[chain.id] = self._sequence_str
            prev_residue_atom = None
            for j in range(len(self.sequence)):
                monomer = self.sequence[j]
                left_ter = False
                right_ter = False
                if j == 0:
                    left_ter = True
                if j == len(self.sequence) - 1:
                    right_ter = True
                prev_residue_atom = self._add_residue_to_chain(topology, chain, prev_residue_atom, monomer,
                                                               left_ter, right_ter)
            if self.create_pdb:
                self._create_chain_pdb(chain)
        return id_to_sequence

    def _add_residue_to_chain(self, topology, chain, prev_residue_atom, monomer, left_ter=False, right_ter=False):

        # Determine monomer type and end chain length
        is_methyl = False
        if monomer.startswith('mA'):
            monomer_type = 'mA'
            is_methyl = True
        else:
            monomer_type = 'A'
        end_chain_length = literal_eval(monomer.replace(monomer_type, ''))

        if left_ter:
            residue_id = "TER0"
        elif right_ter:
            residue_id = "TER1"
        else:
            residue_id = None
        residue = topology.addResidue(monomer, chain, id=residue_id)

        if left_ter:
            carbon = topology.addAtom('C', self.NITROGEN, residue)
        else:
            carbon = topology.addAtom('C', self.CARBON, residue)
        if right_ter:
            carbon_1 = topology.addAtom('C1', self.NITROGEN, residue)
        else:
            carbon_1 = topology.addAtom('C1', self.CARBON, residue)
        if prev_residue_atom is not None:
            topology.addBond(carbon, prev_residue_atom)
        topology.addBond(carbon, carbon_1)

        carbon_2 = topology.addAtom('C2', self.CARBON, residue)
        topology.addBond(carbon_1, carbon_2)
        oxygen_carbonyl = topology.addAtom('O', self.OXYGEN, residue)
        topology.addBond(carbon_2, oxygen_carbonyl)
        oxygen_ether = topology.addAtom('O1', self.OXYGEN, residue)
        topology.addBond(carbon_2, oxygen_ether)

        prev = oxygen_ether
        for i in range(end_chain_length):
            curr = topology.addAtom('C{}'.format(i + 3), self.CARBON, residue)
            topology.addBond(prev, curr)
            prev = curr

        # If methacrylate monomer
        if is_methyl:
            carbon_methyl = topology.addAtom('Cm', self.CARBON, residue)
            topology.addBond(carbon_1, carbon_methyl)
        return carbon_1

    def _create_chain_pdb(self, chain):
        dirname = os.path.dirname(__file__)
        file_path = os.path.join(dirname, "data/{}.pdb".format(self._sequence_str))
        if not os.path.isfile(file_path) or self.overwrite_pdb:
            self.overwrite_pdb = False
            chain_positions = None
            initial_position = np.array([0, 0, 0])
            for residue in chain.residues():
                residue_positions = self._create_residue_positions(residue.name, initial_position)
                initial_position = residue_positions[1]
                if chain_positions is None:
                    chain_positions = residue_positions
                else:
                    chain_positions = np.concatenate((chain_positions, residue_positions), axis=0)
            topology = Topology()
            topology._chains.append(chain)  # TODO: make this more elegant
            PDBFile.writeFile(topology, chain_positions, open(file_path, 'w'))

    @staticmethod
    def _create_residue_positions(monomer, initial_position):

        # Tetrahedral angle in radians
        tetra_angle = 109.5 / 2 * 2 * np.pi / 360

        # Determine monomer type and end chain length
        is_methyl = False
        if monomer.startswith('mA'):
            monomer_type = 'mA'
            is_methyl = True
        else:
            monomer_type = 'A'
        end_chain_length = literal_eval(monomer.replace(monomer_type, ''))

        # Create positions
        positions = []
        c0_pos = initial_position + 1.54*np.array([np.cos(np.pi/6), np.sin(np.pi/6), 0])
        positions.append(c0_pos)
        c1_pos = c0_pos + 1.54*np.array([np.cos(np.pi/6), -np.sin(np.pi/6), 0])
        positions.append(c1_pos)
        c2_pos = c1_pos + 1.52*np.array([0, -1, 0])
        positions.append(c2_pos)
        o_pos = c2_pos + 1.20*np.array([-np.cos(np.pi/6), -np.sin(np.pi/6), 0])
        positions.append(o_pos)
        o1_pos = c2_pos + 1.344*np.array([np.cos(tetra_angle), -np.sin(tetra_angle), 0])
        positions.append(o1_pos)

        prev_pos = o1_pos
        for i in range(end_chain_length):
            if i == 0:
                bond_length = 1.41
            else:
                bond_length = 1.54
            curr_pos = prev_pos + bond_length*np.array([(-1) ** (i+1) * np.cos(tetra_angle), -np.sin(tetra_angle), 0])
            positions.append(curr_pos)
            prev_pos = curr_pos

        # If methacrylate
        if is_methyl:
            cm_pos = c1_pos + 1.54*np.array([0, 1, 0])
            positions.append(cm_pos)

        return np.array(positions)

    # def _create_chain_pdb(self, chain):
    #     dirname = os.path.dirname(__file__)
    #     file_path = os.path.join(dirname, "data/{}.pdb".format(self._sequence_str))
    #     if not os.path.isfile(file_path):
    #         chain_positions = None
    #         prev_residue_positions = None
    #         for residue in chain.residues():
    #             residue_file_path = os.path.join(dirname, "data/{}.pdb".format(residue.name))
    #             residue_positions = md.load(residue_file_path).xyz[0]*10
    #             if prev_residue_positions is not None:
    #                 residue_positions = self._translated_residue_positions(residue_positions, prev_residue_positions)
    #             if chain_positions is None:
    #                 chain_positions = residue_positions
    #             else:
    #                 chain_positions = np.concatenate((chain_positions, residue_positions), axis=0)
    #             prev_residue_positions = residue_positions
    #         topology = Topology()
    #         topology._chains.append(chain)
    #         PDBFile.writeFile(topology, chain_positions, open(file_path, 'w'))
    #
    # def _translated_residue_positions(self, residue_positions, prev_residue_positions):
    #     a0 = prev_residue_positions[0]
    #     a1 = prev_residue_positions[1]
    #     a2 = prev_residue_positions[2]
    #     b0 = residue_positions[0]
    #     b1 = residue_positions[1]
    #     b2 = residue_positions[2]
    #
    #     # def func(x):
    #     #     # Rotate and translate matrices
    #     #     c0 = self._rotate(self._translate(b0, *x[3:]), *x[0:3])
    #     #     c1 = self._rotate(self._translate(b1, *x[3:]), *x[0:3])
    #     #     c2 = self._rotate(self._translate(b2, *x[3:]), *x[0:3])
    #     #     # Objectives
    #     #     obj_1 = 1.54 - np.linalg.norm(a1-c0)
    #     #     obj_2 = np.cos(2*np.pi*109.5/360) - np.dot(a0-a1, c0-a1)/(np.linalg.norm(a0-a1)*np.linalg.norm(c0-a1))
    #     #     obj_3 = np.cos(2*np.pi*109.5/360) - np.dot(a1-c0, c1-c0)/(np.linalg.norm(a1-c0)*np.linalg.norm(c1-c0))
    #     #     obj_4 = np.dot(a1-a0, c1-c0) - np.linalg.norm(a1-a0)*np.linalg.norm(c1-c0)
    #     #     obj_5 = np.dot(np.cross(a1-a0, c1-a0), c0-a0)
    #     #     obj_6 = np.dot(a2-a1, c2-c1) - np.linalg.norm(a2-a1)*np.linalg.norm(c2-c1)
    #     #     return obj_1, obj_2, obj_3, obj_4, obj_5, obj_6
    #     # sol = fsolve(func, np.array([0, 0, 0, 1, 1, 1]))
    #     # return self._rotate_positions(self._translate_positions(residue_positions, *sol[3:]), *sol[0:3])
    #
    #     def func(x):
    #         # Translate matrices
    #         c0 = self._translate(b0, *x)
    #         c1 = self._translate(b1, *x)
    #         # Objectives
    #         obj_1 = 1.54 - np.linalg.norm(a1-c0)
    #         obj_2 = np.cos(2*np.pi*109.5/360) - np.dot(a0-a1, c0-a1)/(np.linalg.norm(a0-a1)*np.linalg.norm(c0-a1))
    #         obj_3 = np.cos(2*np.pi*109.5/360) - np.dot(a1-c0, c1-c0)/(np.linalg.norm(a1-c0)*np.linalg.norm(c1-c0))
    #         return obj_1, obj_2, obj_3
    #     sol = fsolve(func, a0-b0)
    #     return self._translate_positions(residue_positions, *sol)
    #
    # def _rotate_positions(self, positions, alpha, beta, gamma):
    #     return np.array([self._rotate(position, alpha, beta, gamma) for position in positions])
    #
    # @staticmethod
    # def _rotate(position, alpha, beta, gamma):
    #     rotate_x = np.array([[1, 0, 0],
    #                          [0, np.cos(alpha), np.sin(alpha)],
    #                          [0, -np.sin(alpha), np.cos(alpha)]])
    #     rotate_y = np.array([[np.cos(beta), 0, -np.sin(beta)],
    #                          [0, 1, 0],
    #                          [np.sin(beta), 0, np.cos(beta)]])
    #     rotate_z = np.array([[np.cos(gamma), np.sin(gamma), 0],
    #                          [-np.sin(gamma), np.cos(gamma), 0],
    #                          [0, 0, 1]])
    #     return np.dot(rotate_z, np.dot(rotate_y, np.dot(rotate_x, position)))
    #
    # def _translate_positions(self, positions, del_x, del_y, del_z):
    #     return np.array([self._translate(position, del_x, del_y, del_z) for position in positions])
    #
    # @staticmethod
    # def _translate(position, del_x, del_y, del_z):
    #     return position + np.array([del_x, del_y, del_z])
    #
