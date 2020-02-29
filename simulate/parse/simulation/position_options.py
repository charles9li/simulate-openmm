from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from ast import literal_eval

from simtk.openmm.app import Simulation
from openmmtools.testsystems import subrandom_particle_positions
import mdtraj as md

from simulate.parse._options import _Options


class _PositionOptions(_Options):

    def __init__(self):
        super(_PositionOptions, self).__init__()

    def set_positions(self, simulation):
        pass


class FileOptions(_PositionOptions):

    def __init__(self):
        super(FileOptions, self).__init__()
        self.file = None
        self.frame = 0

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.file is None:
            self._incomplete_error('file')

    # =========================================================================

    def _parse_file(self, *args):
        self.file = args[0]

    def _parse_frame(self, *args):
        self.frame = literal_eval(args[0])

    OPTIONS = {'file': _parse_file,
               'frame': _parse_frame}

    # =========================================================================

    def set_positions(self, simulation):
        t = md.load(self.file)
        simulation.context.setPositions(t.xyz[self.frame])


class SubrandomParticlePositions(_PositionOptions):

    SECTION_NAME = "SubrandomParticlePositions"

    DATA_PATH = "../../data/"

    COMPOUNDS = {'A4': DATA_PATH + "A4.pdb"}

    # =========================================================================

    def __init__(self):
        super(SubrandomParticlePositions, self).__init__()
        self.method = 'sobol'

    # =========================================================================

    def _parse_method(self, *args):
        self.method = args[0]

    OPTIONS = {'method': _parse_method}

    # =========================================================================

    def set_positions(self, simulation):
        topology = simulation.topology
        system = simulation.system
        num_residues = topology.getNumAtoms()
        box_vectors = system.getDefaultPeriodicBoxVectors()
        positions = subrandom_particle_positions(num_residues, box_vectors, method=self.method)
        print(type(positions))
        simulation.context.setPositions(positions)
