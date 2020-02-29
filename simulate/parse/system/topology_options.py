from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from ast import literal_eval

import numpy as np
from simtk.openmm.app import ForceField, NoCutoff, Topology
from simtk.unit import nanometer
from parmed import gromacs

from simulate.parse._options import _Options
from simulate.parse.system.chain_options import *


class _TopologyOptions(_Options):

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

    SECTION_NAME = 'AmberTopologyOptions'

    def __init__(self):
        super(AmberTopologyOptions, self).__init__()
        raise NotImplementedError("'{}' is not supported yet.".format(self.SECTION_NAME))


class GromacsTopologyOptions(_TopologyOptions):

    SECTION_NAME = 'GromacsTopologyOptions'

    # =========================================================================

    def __init__(self):
        super(_Options, self).__init__()
        self.topFilename = None
        self.groFilename = None
        self._gromacs_topology = None

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

    OPTIONS = {'topFilename': _parse_top_filename,
               'groFilename': _parse_gro_filename}

    # =========================================================================

    def topology(self):
        self._create_gromacs_topology()
        return self._gromacs_topology.topology

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

    def _create_gromacs_topology(self):
        if self._gromacs_topology is None:
            gro = gromacs.GromacsGroFile.parse(self.groFilename)
            self._gromacs_topology = gromacs.GromacsTopologyFile(self.topFilename)
            self._gromacs_topology.box = gro.box
            
            
class DodecaneAcrylateTopologyOptions(_TopologyOptions):

    SECTION_NAME = "DodecaneAcrylateTopologyOptions"

    TRAPPEUA_FF_PATH = "/home/charlesli/lab/shell/simulate-openmm/simulate/data/trappeua-acrylates.xml"
    
    def __init__(self):
        super(DodecaneAcrylateTopologyOptions, self).__init__()
        self.force_field = ForceField(self.TRAPPEUA_FF_PATH)
        self.numDodecane = 0
        self.box_vectors = None
        self.chains = []

    # =========================================================================

    CHAIN_METHODS = {'Homopolymer': HomopolymerOptions}

    def _parse_num_dodecane(self, *args):
        self.numDodecane = literal_eval(args[0])

    def _parse_box(self, *args):
        a, b, c = args[0].split(' ')
        self.box_vectors = np.array([[literal_eval(a), 0.0, 0.0],
                                     [0.0, literal_eval(b), 0.0],
                                     [0.0, 0.0, literal_eval(c)]])*nanometer

    def _parse_chains(self, *args):
        line_deque = args[1].popleft()
        while len(line_deque) > 0:
            line = line_deque.popleft()
            chain_name = self._parse_option_name(line)
            chain_name = self._parse_option_value(line, chain_name)
            chain_options = self.CHAIN_METHODS[chain_name]()
            chain_options.parse(line_deque.popleft())
            self.chains.append(chain_options)

    OPTIONS = {'numDodecane': _parse_num_dodecane,
               'box': _parse_box}

    SECTIONS = {'chains': _parse_chains}

    # =========================================================================

    def topology(self):
        if self._topology is None:
            topology = Topology()
            if self.box_vectors is not None:
                topology.setPeriodicBoxVectors(self.box_vectors)
            for chain_option in self.chains:
                chain_option.add_chain_to_topology(topology)
            self._topology = topology
        return self._topology

    def create_system(self, nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*nanometer,
                      constraints=None, rigidWater=True, implicitSolvent=None,
                      soluteDielectric=1.0, solventDielectric=78.5,
                      ewaldErrorTolerance=0.0005, removeCMMotion=True,
                      hydrogenMass=None):
        return self.force_field.createSystem(self._topology, nonbondedMethod=nonbondedMethod,
                                             nonbondedCutoff=nonbondedCutoff,
                                             constraints=constraints, rigidWater=rigidWater,
                                             implicitSolvent=implicitSolvent, soluteDielectric=soluteDielectric,
                                             solventDielectric=solventDielectric, ewaldErrorTolerance=ewaldErrorTolerance,
                                             removeCMMotion=removeCMMotion, hydrogenMass=hydrogenMass)
