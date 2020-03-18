"""
_simulations_options.py: Parses position initialization options for simulations.

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

from simtk.openmm.app import PDBFile, Topology
from openmmtools.testsystems import subrandom_particle_positions
import mdtraj as md
import numpy as np
from scipy.optimize import fsolve

from ._options import _Options


class _PositionOptions(_Options):

    def __init__(self):
        super(_PositionOptions, self).__init__()

    def set_positions(self, simulation, *args):
        pass


class FileOptions(_PositionOptions):

    def __init__(self):
        super(FileOptions, self).__init__()
        self.file = None
        self.frame = 0

    def _create_options(self):
        super(FileOptions, self)._create_options()
        self._OPTIONS['file'] = self._parse_file
        self._OPTIONS['frame'] = self._parse_frame

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.file is None:
            self._incomplete_error('file')

    # =========================================================================

    def _parse_file(self, *args):
        self.file = args[0]

    def _parse_frame(self, *args):
        self.frame = literal_eval(args[0])

    # =========================================================================

    def set_positions(self, simulation, *args):
        t = md.load(self.file)
        simulation.context.setPositions(t.xyz[self.frame])


class SubrandomParticlePositions(_PositionOptions):

    _SECTION_NAME = "SubrandomParticlePositions"

    # =========================================================================

    def __init__(self):
        super(SubrandomParticlePositions, self).__init__()
        self.method = 'sobol'

    def _create_options(self):
        super(SubrandomParticlePositions, self)._create_options()
        self._OPTIONS['method'] = self._parse_method

    # =========================================================================

    def _parse_method(self, *args):
        self.method = args[0]

    # =========================================================================

    def set_positions(self, simulation, *args):
        topology = simulation.topology
        system = simulation.system
        num_residues = topology.getNumAtoms()
        box_vectors = system.getDefaultPeriodicBoxVectors()
        positions = subrandom_particle_positions(num_residues, box_vectors, method=self.method)
        simulation.context.setPositions(positions)


class DodecaneAcrylatePositionOptions(_PositionOptions):

    _SECTION_NAME = "DodecaneAcrylatePosition"

    # =========================================================================

    def __init__(self):
        super(DodecaneAcrylatePositionOptions, self).__init__()

    # =========================================================================

    _CHAIN_TO_PDB = {}
    _CHAIN_TO_NUM = {}

    def set_positions(self, simulation, *args):

        # Get topology and its user options
        topology = simulation.topology
        topology_options = args[0]
        id_to_sequence = topology_options.id_to_sequence

        # Create dictionary
        sequence_to_pdb = {}

        # Create PDB files
        sequences = []
        num = []
        index = -1
        prev_sequence = None
        for chain in topology:
            sequence = id_to_sequence[chain.id]
            if sequence != prev_sequence:
                sequences.append(sequence)
                index += 1
                num.append(1)
            else:
                num[index] += 1

    # TODO: call mdapackmol to produce initial positions
