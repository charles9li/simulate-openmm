from openmmtools.integrators import VelocityVerletIntegrator
from parmed import gromacs as gmx
from simtk.openmm.app import Simulation
from simtk.openmm import Context, NonbondedForce
import mdtraj as md

from .run_simulation import run_simulation
from ..parse import OpenMMOptions, SimulationOptions


class SimulationInformation(object):

    def __init__(self, filename):
        self.openmm_options = OpenMMOptions(filename)
        self._create_topologies()
        self._initialize_positions()
        self._initialize_velocities()

    def _create_topologies(self):
        topology_options = self.openmm_options.topology_options
        top = gmx.GromacsTopologyFile(topology_options.topFilename)
        gro = gmx.GromacsGroFile.parse(topology_options.groFilename)
        top.box = gro.box
        self.gmx_topology = top
        self.topology = top.topology

    def _create_system(self):
        system_options = self.openmm_options.system_options
        system = self.gmx_topology.createSystem(nonbondedMethod=system_options.nonbondedMethod,
                                                nonbondedCutoff=system_options.nonbondedCutoff,
                                                ewaldErrorTolerance=system_options.ewaldErrorTolerance)
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                force.setUseDispersionCorrection(system_options.useDispersionCorrection)
                force.setNonbondedMethod(system_options.nonbondedMethod)
        return system

    def _create_simulation(self, simulation_options):
        system = self._create_system()
        if simulation_options.barostat is not None:
            system.addForce(simulation_options.barostat)
        simulation = Simulation(self.topology, system, simulation_options.integrator)
        for reporter in simulation_options.reporters:
            simulation.reporters.append(reporter)
        return simulation

    def _initialize_positions(self):
        topology_options = self.openmm_options.topology_options
        self.positions = md.load(topology_options.coordFilename)

    def _initialize_velocities(self):
        system = self._create_system()
        integrator = VelocityVerletIntegrator()
        context = Context(system, integrator)
        initial_velocity_temperature = self.openmm_options.system_options.initialVelocityTemperature
        if initial_velocity_temperature is not None:
            context.setVelocitiesToTemperature(initial_velocity_temperature)
            state = context.getState(getVelocities=True)
            self.velocities = state.getVelocities(asNumpy=True)
        else:
            self.velocities = None

    def run(self):
        for simulation_options in self.openmm_options.simulation_options:
            simulation = self._create_simulation(simulation_options)
            simulation.context.setPositions(self.positions)
            if self.velocities is not None:
                simulation.context.setVelocities(self.velocities)
            if isinstance(simulation_options, SimulationOptions):
                run_simulation(simulation, simulation_options)
            state = simulation.context.getState(getPositions=True, getVelocities=True)
            self.positions = state.getVelocities(asNumpy=True)
            self.velocities = state.getVelocities(asNumpy=True)
