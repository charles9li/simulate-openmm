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

import mdtraj as md
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
                level += 1
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
                elif isinstance(token, str) or isinstance(token, list):
                    for _ in range(multiplier):
                        self._parse_tokenized_string(token)
                    multiplier = 1
                else:
                    raise ValueError("Invalid sequence.")


class ChainOptions(_Options):

    _SECTION_NAME = "Chain"

    # =========================================================================

    HYDROGEN = Element.getBySymbol('H')
    BORON = Element.getBySymbol('B')
    CARBON = Element.getBySymbol('C')
    NITROGEN = Element.getBySymbol('N')
    OXYGEN = Element.getBySymbol('O')

    # =========================================================================

    def __init__(self, topology_options):
        super(ChainOptions, self).__init__()
        self.topology_options = topology_options
        self.forceField_str = topology_options.forceField_str
        self.id = None
        self.num = 1
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
        if self.sequence is None:
            self._incomplete_error('structure')

    # =====================================================================[====

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
                if self.forceField_str == 'TraPPE-UA':
                    add_residue_function = self._add_residue_to_chain_trappeua
                elif self.forceField_str == 'OPLS-AA':
                    add_residue_function = self._add_residue_to_chain_oplsaa
                else:
                    raise ValueError("Invalid force field.")
                prev_res_atom, prev_res_atom_pdb, prev_res_atom_pos = add_residue_function(topology, topology_pdb,
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

    def _add_residue_to_chain_oplsaa(self, topology, topology_pdb, chain, chain_pdb, prev_res_atom, prev_res_atom_pdb,
                                     positions, prev_res_atom_pos, monomer, left_ter=False, right_ter=False):
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

        # Function to add atoms to topology
        def _add_atom_to_topology_oplsaa(name, element, counter):
            # Add atom to topology
            name = "{}{}".format(name, counter)
            atom = topology.addAtom(name, element, residue)
            atom_pdb = topology_pdb.addAtom(name, element, residue_pdb)

            # Return atoms and position
            return atom, atom_pdb, counter + 1

        # Import positions from pdb
        residue_positions = md.load(os.path.join(self.topology_options.data_directory, "mA12_aa.pdb")).xyz[0]*10. + prev_res_atom_pos
        residue_positions = list(residue_positions)
        if end_chain_length < 12:
            residue_positions = residue_positions[:-(12-end_chain_length)*3]
        if not is_methyl:
            residue_positions.pop(9)
            residue_positions.pop(8)
            residue_positions.pop(7)
        if not right_ter or (left_ter and right_ter):
            residue_positions.pop(5)
        if (left_ter and right_ter) or (not left_ter):
            residue_positions.pop(2)
        positions += residue_positions

        # Add first carbon and attached hydrogen atoms
        C_counter = 0
        H_counter = 0
        C0, C0_pdb, C_counter = _add_atom_to_topology_oplsaa("C", self.CARBON, C_counter)
        if left_ter and not right_ter:
            num_H = 3
        else:
            num_H = 2
        for _ in range(num_H):
            H, H_pdb, H_counter = _add_atom_to_topology_oplsaa("H", self.HYDROGEN, H_counter)
            self._add_bond_to_topology(H, H_pdb, C0, C0_pdb, topology, topology_pdb)

        # Add bond to previous residue if applicable
        if prev_res_atom is not None:
            self._add_bond_to_topology(C0, C0_pdb, prev_res_atom, prev_res_atom_pdb, topology, topology_pdb)

        # Add second carbon
        C1, C1_pdb, C_counter = _add_atom_to_topology_oplsaa("C", self.CARBON, C_counter)
        self._add_bond_to_topology(C1, C1_pdb, C0, C0_pdb, topology, topology_pdb)
        if right_ter and not left_ter:
            H, H_pdb, H_counter = _add_atom_to_topology_oplsaa("H", self.HYDROGEN, H_counter)
            self._add_bond_to_topology(H, H_pdb, C1, C1_pdb, topology, topology_pdb)

        # Add methyl group or hydrogen if not methyl
        if is_methyl:
            Cm, Cm_pdb, C_counter = _add_atom_to_topology_oplsaa("C", self.CARBON, C_counter)
            for _ in range(3):
                H, H_pdb, H_counter = _add_atom_to_topology_oplsaa("H", self.HYDROGEN, H_counter)
                self._add_bond_to_topology(H, H_pdb, Cm, Cm_pdb, topology, topology_pdb)
            self._add_bond_to_topology(Cm, Cm_pdb, C1, C1_pdb, topology, topology_pdb)
        else:
            H, H_pdb, H_counter = _add_atom_to_topology_oplsaa("H", self.HYDROGEN, H_counter)
            self._add_bond_to_topology(H, H_pdb, C1, C1_pdb, topology, topology_pdb)

        # Add carbonyl carbon
        Cc, Cc_pdb, C_counter = _add_atom_to_topology_oplsaa("C", self.CARBON, C_counter)
        self._add_bond_to_topology(Cc, Cc_pdb, C1, C1_pdb, topology, topology_pdb)

        # Add carbonyl oxygen
        Oc, Oc_pdb, _, = _add_atom_to_topology_oplsaa("O", self.OXYGEN, 0)
        self._add_bond_to_topology(Oc, Oc_pdb, Cc, Cc_pdb, topology, topology_pdb)

        # Add ether oxygen
        Oe, Oe_pdb, _ = _add_atom_to_topology_oplsaa("O", self.OXYGEN, 1)
        self._add_bond_to_topology(Oe, Oe_pdb, Cc, Cc_pdb, topology, topology_pdb)

        # Add alkyl chain
        prev_atom = Oe
        prev_atom_pdb = Oe_pdb
        for _ in range(end_chain_length):
            curr_atom, curr_atom_pdb, C_counter = _add_atom_to_topology_oplsaa("C", self.CARBON, C_counter)
            for _ in range(2):
                H, H_pdb, H_counter = _add_atom_to_topology_oplsaa("H", self.HYDROGEN, H_counter)
                self._add_bond_to_topology(H, H_pdb, curr_atom, curr_atom_pdb, topology, topology_pdb)
            self._add_bond_to_topology(curr_atom, curr_atom_pdb, prev_atom, prev_atom_pdb, topology, topology_pdb)
            prev_atom = curr_atom
            prev_atom_pdb = curr_atom_pdb

        # Add hydrogen to last carbon
        H, H_pdb, H_counter = _add_atom_to_topology_oplsaa("H", self.HYDROGEN, H_counter)
        self._add_bond_to_topology(H, H_pdb, prev_atom, prev_atom_pdb, topology, topology_pdb)

        return C1, C1_pdb, prev_res_atom_pos + np.array([2.527, 0, 0])

    def _add_residue_to_chain_trappeua(self, topology, topology_pdb, chain, chain_pdb, prev_res_atom, prev_res_atom_pdb,
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
        carbon, carbon_pdb = self._add_atom_to_topology_trappeua('C', carbon_element, residue, residue_pdb,
                                                                 topology, topology_pdb,
                                                                 positions, carbon_pos)
        carbon_1_pos = carbon_pos + 1.54*np.array([cos19*cos30, cos19*sin30, -sin19])
        carbon_1, carbon_1_pdb = self._add_atom_to_topology_trappeua('C1', carbon_1_element, residue, residue_pdb,
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
            carbon_methyl, carbon_methyl_pdb = self._add_atom_to_topology_trappeua('Cm', self.CARBON, residue, residue_pdb,
                                                                                   topology, topology_pdb,
                                                                                   positions, carbon_methyl_pos)
            self._add_bond_to_topology(carbon_1, carbon_1_pdb, carbon_methyl, carbon_methyl_pdb,
                                       topology, topology_pdb)

        # Add ester group atoms
        carbon_2_pos = carbon_1_pos + 1.52*np.array([0.0, 0.0, -1.0])
        carbon_2, carbon_2_pdb = self._add_atom_to_topology_trappeua('C2', self.CARBON, residue, residue_pdb,
                                                                     topology, topology_pdb,
                                                                     positions, carbon_2_pos)
        self._add_bond_to_topology(carbon_1, carbon_1_pdb, carbon_2, carbon_2_pdb, topology, topology_pdb)
        oxygen_carbonyl_pos = carbon_2_pos + 1.20*np.array([0.0, cos30, -sin30])
        oxygen_carbonyl, oxygen_carbonyl_pdb = self._add_atom_to_topology_trappeua('O', self.OXYGEN, residue, residue_pdb,
                                                                                   topology, topology_pdb,
                                                                                   positions, oxygen_carbonyl_pos)
        self._add_bond_to_topology(carbon_2, carbon_2_pdb, oxygen_carbonyl, oxygen_carbonyl_pdb, topology, topology_pdb)
        oxygen_ether_pos = carbon_2_pos + 1.344*np.array([0.0, -cos30, -sin30])
        oxygen_ether, oxygen_ether_pdb = self._add_atom_to_topology_trappeua('O', self.OXYGEN, residue, residue_pdb,
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
            curr_atom, curr_atom_pdb = self._add_atom_to_topology_trappeua('C{}'.format(i + 3), self.CARBON,
                                                                           residue, residue_pdb,
                                                                           topology, topology_pdb,
                                                                           positions, curr_atom_pos)
            self._add_bond_to_topology(prev_atom, prev_atom_pdb, curr_atom, curr_atom_pdb, topology, topology_pdb)
            prev_atom_pos = curr_atom_pos
            prev_atom = curr_atom
            prev_atom_pdb = curr_atom_pdb

        return carbon_1, carbon_1_pdb, carbon_1_pos

    @staticmethod
    def _add_atom_to_topology_trappeua(name, element, residue, residue_pdb, topology, topology_pdb,
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
        filename = "{}.pdb".format(self.sequence_str)
        if self.forceField_str == 'OPLS-AA':
            filename = "{}_aa.pdb".format(self.sequence_str)
        file_path = os.path.join(dirname, "data/{}".format(filename))
        if not os.path.isfile(file_path) or self.overwrite_pdb:
            self.overwrite_pdb = False
            PDBFile.writeFile(topology_pdb, positions, open(file_path, 'w'))


_OPERATOR_PRECEDENCE = {'+': 2,
                        '*': 3}
_OPERATORS = list(_OPERATOR_PRECEDENCE.keys())
_DELIMITERS = ['(', ')']


def tokenize_expr(expr):
    token_list = expr.split()
    for char in _OPERATORS + _DELIMITERS:
        temp_list = []
        for token in token_list:
            token = [t.strip() for t in token.split(char)]
            token_split = [char] * (len(token) * 2 - 1)
            token_split[0::2] = token
            temp_list += token_split
        token_list = temp_list
    return [t for t in token_list if t]


def create_stack(token_list):
    output = []
    stack = []
    for token in token_list:
        if token in _OPERATORS:
            while len(stack) > 0 and stack[-1] != '(' and _OPERATOR_PRECEDENCE[stack[-1]] >= _OPERATOR_PRECEDENCE[token]:
                output.append(stack.pop())
            stack.append(token)
        elif token == '(':
            stack.append(token)
        elif token == ')':
            while stack[-1] != '(':
                output.append(stack.pop())
            stack.pop()
        else:
            output.append(token)
    while len(stack) > 0:
        output.append(stack.pop())
    return output


def evaluate_output(output):
    stack = []
    for item in output:
        if item in _OPERATORS:
            op2 = stack.pop()
            op1 = stack.pop()
            if item == '+':
                stack.append(op1 + op2)
            elif item == '*':
                stack.append(op1 * op2)
        else:
            try:
                stack.append(int(item))
            except ValueError:
                stack.append([item])
    if len(stack) != 1:
        raise ValueError("Malformed expression.")
    return stack.pop()


class BranchedChainOptions(ChainOptions):

    _SECTION_NAME = "BranchedChain"

    # =========================================================================

    PHOSPHORUS = Element.getBySymbol('P')

    # =========================================================================

    def __init__(self, topology_options):
        super(BranchedChainOptions, self).__init__(topology_options)
        self.backbone = None
        self.branches = []
        self.pdb = None

    def _create_options(self):
        super(BranchedChainOptions, self)._create_options()
        self._OPTIONS['backbone'] = self._parse_backbone

    def _create_sections(self):
        super(BranchedChainOptions, self)._create_sections()
        self._SECTIONS['Branch'] = self._parse_branch

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.backbone is None:
            self._incomplete_error('backbone')

    # =========================================================================

    def _parse_backbone(self, *args):
        self.backbone = evaluate_output(create_stack(tokenize_expr(args[0])))

    def _parse_branch(self, *args):
        line_deque = args[1]
        branch_options = BranchOptions()
        branch_options.parse(line_deque.popleft())
        self.branches.append((branch_options.sequence, branch_options.indices))

    # =========================================================================

    def add_chain_to_topology(self, topology):
        for _ in range(self.num):
            self._add_chain(topology)

    @staticmethod
    def _find_atom_by_name(name, residue):
        for atom in residue.atoms():
            if atom.name == name:
                return atom
        raise ValueError("No atom with name {} found.".format(name))

    def _add_chain(self, topology):

        # create chain
        chain = topology.addChain()

        # determine residue adding function
        if self.forceField_str == 'TraPPE-UA':
            add_residue_function = self._add_residue_trappeua
        elif self.forceField_str == 'OPLS-AA':
            add_residue_function = self._add_residue_to_chain_oplsaa
        else:
            raise ValueError("Invalid force field.")

        # initialize previous residue atom
        prev_res_atom = None

        # add backbone
        for i, monomer in enumerate(self.backbone):

            # determine whether or not the residue is a terminal one
            left_ter = False
            right_ter = False
            if i == 0:
                left_ter = True
            if i == len(self.backbone) - 1:
                right_ter = True

            # add residue
            prev_res_atom = add_residue_function(topology, chain, prev_res_atom, monomer, left_ter, right_ter)

        # add branches
        for sequence, indices in self.branches:
            residues = list(chain.residues())
            for index in indices:
                atom = self._find_atom_by_name('C1', residues[index])
                atom.element = self.PHOSPHORUS
                prev_res_atom = atom
                for i, monomer in enumerate(sequence):
                    right_ter = i == len(sequence) - 1
                    prev_res_atom = add_residue_function(topology, chain, prev_res_atom, monomer, right_ter=right_ter)

        # create pdb file
        if self.pdb is None:
            positions = self._create_positions_trappeua()
            topology_pdb = Topology()
            topology_pdb._chains.append(chain)
            bc_num = 0
            dirname = os.path.join(os.path.dirname(__file__), "data")
            while os.path.exists(os.path.join(dirname, "branched_chain_{}.pdb".format(bc_num))):
                bc_num += 1
            self.pdb = os.path.join(dirname, "branched_chain_{}.pdb".format(bc_num))
            PDBFile.writeFile(topology_pdb, positions, open(self.pdb, 'w'))

    def _create_positions_trappeua(self):

        positions = None

        # load positions for mA12
        mA12_pos = self._mA12_positions_trappeua()

        # initialize list mapping residue number to atom index
        res2atom = []

        # number of atoms in each monomer type
        Natoms_dict = {'A4': 9,
                       'A12': 17,
                       'mA12': 18}

        # add backbone positions
        for i, monomer in enumerate(self.backbone):

            # add atom indices to res2atom
            if res2atom:
                res2atom.append(np.arange(Natoms_dict[monomer]) + res2atom[-1][-1] + 1)
            else:
                res2atom.append(np.arange(Natoms_dict[monomer]))

            # create monomer positions
            pos = mA12_pos[:]
            if monomer.startswith('A'):
                pos = np.delete(pos, 2, axis=0)
            if monomer.endswith('4'):
                pos = pos[:-8]

            # add monomer positions to chain positions
            if positions is None:
                positions = pos
            else:
                cos19 = np.cos(np.deg2rad(19.5))
                cos30 = np.cos(np.deg2rad(30.0))
                pos += 1.54*np.array([cos19*cos30, 0., 0.]) + np.array([positions[res2atom[i-1][1]][0], 0., 0.]) - np.array([mA12_pos[0][0], 0., 0.])
                positions = np.vstack((positions, pos))

        # rotate methacrylate monomer positions
        Rz = np.array([[0., 1., 0.],
                       [-1., 0., 0.],
                       [0., 0., 1.]])
        mA12_pos_rot = mA12_pos.dot(Rz)

        # add side chain positions
        for sequence, indices in self.branches:
            for res_index in indices:
                prev_res_pos = positions[res2atom[res_index][1]]
                for monomer in sequence:

                    # create monomer positions
                    pos = mA12_pos_rot[:] + np.array([positions[res2atom[res_index][1]][0] - mA12_pos_rot[1][0], 0., 0.])
                    if monomer.startswith('A'):
                        pos = np.delete(pos, 2, axis=0)
                    if monomer.endswith('4'):
                        pos = pos[:-8]

                    # add monomer positions
                    cos19 = np.cos(np.deg2rad(19.5))
                    cos30 = np.cos(np.deg2rad(30.0))
                    pos += 1.54*np.array([0., cos19*cos30, 0.]) + np.array([0., prev_res_pos[1], 0.]) - np.array([0., pos[0][1], 0.])
                    positions = np.vstack((positions, pos))

                    prev_res_pos = pos[1]

        return positions

    @staticmethod
    def _mA12_positions_trappeua():

        # sines and cosines
        sin19 = np.sin(np.deg2rad(19.5))
        cos19 = np.cos(np.deg2rad(19.5))
        sin30 = np.sin(np.deg2rad(30.0))
        cos30 = np.cos(np.deg2rad(30.0))
        sin54 = np.sin(np.deg2rad(54.75))
        cos54 = np.cos(np.deg2rad(54.75))

        # methacrylate positions
        pos = 1.54*np.array([cos19*cos30, -cos19*sin30, sin19])
        pos = np.vstack((pos, pos + 1.54*np.array([cos19*cos30, cos19*sin30, -sin19])))
        pos = np.vstack((pos, pos[1] + 1.54*np.array([0., cos19, sin19])))
        pos = np.vstack((pos, pos[1] + 1.52*np.array([0., 0., -1.])))
        pos = np.vstack((pos, pos[3] + 1.20*np.array([0., cos30, -sin30])))
        pos = np.vstack((pos, pos[3] + 1.344*np.array([0., -cos30, -sin30])))

        # C12 tail positions
        for i in range(12):
            if i == 0:
                pos = np.vstack((pos, pos[-1] + 1.41*np.array([0.0, (-1)**(i+1)*cos54, -sin54])))
            else:
                pos = np.vstack((pos, pos[-1] + 1.54*np.array([0.0, (-1)**(i+1)*cos54, -sin54])))

        return pos

    def _add_residue_trappeua(self, topology, chain, prev_res_atom, monomer, left_ter=False, right_ter=False):

        # determine monomer type
        if monomer.startswith('mA'):
            monomer_type = 'mA'
            is_methyl = True
        else:
            monomer_type = 'A'
            is_methyl = False

        # determine chain length
        end_chain_length = literal_eval(monomer.replace(monomer_type, ''))

        # add residue to topology
        residue = topology.addResidue(monomer, chain)

        # add first two carbons
        if left_ter and right_ter:
            c_element = self.NITROGEN
            c1_element = self.NITROGEN
        elif left_ter:
            c_element = self.BORON
            c1_element = self.CARBON
        elif right_ter:
            c_element = self.NITROGEN
            c1_element = self.CARBON
        else:
            c_element = self.CARBON
            c1_element = self.CARBON
        c = topology.addAtom('C', c_element, residue)
        c1 = topology.addAtom('C1', c1_element, residue)

        # add bond to previous residue if applicable
        if prev_res_atom is not None:
            topology.addBond(prev_res_atom, c)

        # add bond between first two carbons
        topology.addBond(c, c1)

        # add methyl group if methacryalte
        if is_methyl:
            cm = topology.addAtom('Cm', self.CARBON, residue)
            topology.addBond(c1, cm)

        # add ester group atoms
        c2 = topology.addAtom('C2', self.CARBON, residue)
        topology.addBond(c1, c2)
        oc = topology.addAtom('O', self.OXYGEN, residue)
        topology.addBond(c2, oc)
        oe = topology.addAtom('O1', self.OXYGEN, residue)
        topology.addBond(c2, oe)

        # add carbon chain
        prev_atom = oe
        for i in range(end_chain_length):
            curr_atom = topology.addAtom('C{}'.format(i+3), self.CARBON, residue)
            topology.addBond(prev_atom, curr_atom)
            prev_atom = curr_atom

        return c1

    def _add_residue_oplsaa(self, topology, chain, prev_res_atom, monomer, left_ter=False, right_ter=False):
        raise NotImplementedError("OPLS-AA not implemented for branched chain.")


class BranchOptions(_Options):

    _SECTION_NAME = "Branch"

    def __init__(self):
        super(BranchOptions, self).__init__()
        self.indices = []
        self.sequence = None

    def _create_options(self):
        super(BranchOptions, self)._create_options()
        self._OPTIONS['sequence'] = self._parse_sequence
        self._OPTIONS['indices'] = self._parse_indices

    # =========================================================================
    
    def _check_for_incomplete_input(self):
        if self.sequence is None:
            self._incomplete_error('sequence')

    # =========================================================================

    def _parse_sequence(self, *args):
        self.sequence = evaluate_output(create_stack(tokenize_expr(args[0])))

    def _parse_indices(self, *args):
        self.indices = [literal_eval(i) for i in args[0].split()]
