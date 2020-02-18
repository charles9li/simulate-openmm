from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from simtk.openmm.app import NoCutoff
from simtk.unit import nanometer
from parmed import gromacs

from simulate.parse._options import _Options


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
