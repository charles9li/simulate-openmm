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


class _ChainOptions(_Options):

    CARBON = Element.getBySymbol('C')
    NITROGEN = Element.getBySymbol('N')
    OXYGEN = Element.getBySymbol('O')

    # =========================================================================

    def __init__(self):
        super(_ChainOptions, self).__init__()

    def _add_chain_to_topology(self, topology):
        pass


class HomopolymerOptions(_ChainOptions):

    _SECTION_NAME = "Homopolymer"

    # =========================================================================

    def __init__(self):
        super(HomopolymerOptions, self).__init__()
        self.num = 0
        self.monomer = 'A1'
        self.N = 1

    def _create_options(self):
        super(HomopolymerOptions, self)._create_options()
        self._OPTIONS['num'] = self._parse_num
        self._OPTIONS['monomer'] = self._parse_monomer
        self._OPTIONS['N'] = self._parse_N

    # =========================================================================

    def _parse_num(self, *args):
        self.num = literal_eval(args[0])

    def _parse_monomer(self, *args):
        self.monomer = args[0]
        if self.monomer.startswith('mA'):
            monomer_type = 'mA'
            self.methyl = True
        else:
            monomer_type = 'A'
        self.end_chain_length = literal_eval(self.monomer.replace(monomer_type, ''))

    def _parse_N(self, *args):
        self.N = literal_eval(args[0])

    # =========================================================================

    def _add_chain_to_topology(self, topology):
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
