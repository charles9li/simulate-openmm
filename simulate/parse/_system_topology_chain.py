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

from simtk.openmm.app import Element

from ._options import _Options


class _StructureParser(object):

    def __init__(self, structure_string):
        self._structure_string = "".join(structure_string.split())
        self._index = 0
        tokenized_string = self._parse_structure_string(0)
        self._chain_sequence = []
        self._parse_tokenized_string(tokenized_string)

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
        self.num = None
        self.sequence = None

    def _create_options(self):
        super(ChainOptions, self)._create_options()
        self._OPTIONS['num'] = self._parse_num
        self._OPTIONS['sequence'] = self._parse_sequence

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.num is None:
            self._incomplete_error('num')
        if self.sequence is None:
            self._incomplete_error('structure')

    # =========================================================================

    def _parse_num(self, *args):
        self.num = literal_eval(args[0])

    def _parse_sequence(self, *args):
        self.sequence = _StructureParser(args[0]).get_chain_sequence()

    # =========================================================================

    def add_chain_to_topology(self, topology):
        for _ in range(self.num):
            chain = topology.addChain()
            prev_residue_atom = None
            for i in range(self.N):
                left_ter = False
                right_ter = False
                if i == 0:
                    left_ter = True
                if i == self.N - 1:
                    right_ter = True
                prev_residue_atom = self._add_residue_to_chain(topology, chain, prev_residue_atom, left_ter, right_ter)

    def _add_residue_to_chain(self, topology, chain, prev_residue_atom, left_ter=False, right_ter=False):
        if left_ter:
            residue_id = "TER0"
        elif right_ter:
            residue_id = "TER1"
        else:
            residue_id = None
        residue = topology.addResidue(self.monomer, chain, id=residue_id)

        if left_ter:
            C = topology.addAtom('C', self.NITROGEN, residue)
        else:
            C = topology.addAtom('C', self.CARBON, residue)
        if right_ter:
            C1 = topology.addAtom('C1', self.NITROGEN, residue)
        else:
            C1 = topology.addAtom('C1', self.CARBON, residue)
        if prev_residue_atom is not None:
            topology.addBond(C, prev_residue_atom)
        topology.addBond(C, C1)

        # TODO: add option if methacrylate

        C2 = topology.addAtom('C2', self.CARBON, residue)
        topology.addBond(C1, C2)
        O = topology.addAtom('O', self.OXYGEN, residue)
        topology.addBond(C2, O)
        O1 = topology.addAtom('O1', self.OXYGEN, residue)
        topology.addBond(C2, O1)

        prev = O1
        for i in range(self.end_chain_length):
            curr = topology.addAtom('C{}'.format(i + 3), self.CARBON, residue)
            topology.addBond(prev, curr)
            prev = curr

        return C1
