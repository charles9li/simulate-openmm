from parmed import gromacs as gmx

from ..parse import OpenMMOptions


class SimulationInformation(object):

    def __init__(self, filename):
        self.openmm_options = OpenMMOptions(filename)
        self._create_topology_and_system()

    def _create_topology_and_system(self):
        topology_options = self.openmm_options.topology_options
        top = gmx.GromacsTopologyFile(topology_options.topFilename)
        gro = gmx.GromacsGroFile.parse(topology_options.groFilename)
        top.box = gro.box
        self.topology = top.topology

        system_options = self.openmm_options.system_options
        self.system = top.createSystem(nonbondedMethod=system_options.nonbondedMethod,
                                       nonbondedCutoff=system_options.nonbondedCutoff,
                                       ewaldErrorTolerance=system_options.ewaldErrorTolerance)
