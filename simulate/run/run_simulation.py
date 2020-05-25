from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from simtk.openmm.app import PDBFile

from simulate.parse import InputOptions
from simulate.parse import RNEMDOptions
from simulate.run.run_rnemd import run_rnemd
from simulate.run.run_average import run_average
from simulate.run.run_temperature_ramp import run_temperature_ramp


class RunSimulation(object):

    def __init__(self, input_file):
        self.input_options = InputOptions(input_file)
        self.system_options = self.input_options.system_options
        self.topology_options = self.system_options.topology_options
        self.simulation_ensembles = self.input_options.simulation_ensembles
        self.ensembles = self.simulation_ensembles.ensembles
        self.positions = None
        self.velocities = None
        self.periodic_box_vectors = None

    def run(self):

        # initialize simulation
        simulation = None

        # iterate through and run ensembles
        for ensemble_options in self.ensembles:

            # # create topology and system
            # if simulation is None:
            #     topology = self.system_options.topology()
            #     system = self.system_options.create_system()
            # else:
            #     topology = simulation.topology
            #     system = simulation.system

            topology = self.system_options.topology()
            system = self.system_options.create_system()

            # modify periodic box vectors if necessary
            if self.periodic_box_vectors is not None:
                system.setDefaultPeriodicBoxVectors(*self.periodic_box_vectors)

            # create simulation
            simulation = ensemble_options.create_simulation(topology, system)

            # initialize positions
            if self.positions is None:
                self.simulation_ensembles.set_positions(simulation, self.topology_options)
            else:
                simulation.context.setPositions(self.positions)

            # initialize velocities
            if self.velocities is None:
                self.simulation_ensembles.set_velocities(simulation)
            else:
                simulation.context.setVelocities(self.velocities)

            # load save state if specified
            if ensemble_options.loadState is not None:
                simulation.loadState(ensemble_options.loadState)

            # minimize energy
            minimize_energy_options = ensemble_options.minimize_energy_options
            if minimize_energy_options is not None:
                tolerance = minimize_energy_options.tolerance
                max_iterations = minimize_energy_options.maxIterations
                simulation.minimizeEnergy(tolerance=tolerance, maxIterations=max_iterations)
                if minimize_energy_options.file is not None:
                    positions = simulation.context.getState(getPositions=True).getPositions()
                    PDBFile.writeFile(topology, positions, file=open(minimize_energy_options.file, 'w'))

            # create voids
            void_options = ensemble_options.void_options
            if void_options is not None:
                system_options = ensemble_options.simulations_options.input_options.system_options
                simulation = void_options.remove_molecules(simulation, system_options, ensemble_options)

            # apply constraints
            simulation.context.applyConstraints(simulation.integrator.getConstraintTolerance())

            # run simulation
            simulation = self._decide_run_type(topology, system, simulation, ensemble_options)

            # store positions and velocities for next ensemble
            state = simulation.context.getState(getPositions=True, getVelocities=True)
            self.positions = state.getPositions(asNumpy=True)
            self.velocities = state.getVelocities(asNumpy=True)
            self.periodic_box_vectors = state.getPeriodicBoxVectors(asNumpy=True)

            # save state if specified
            if ensemble_options.saveState is not None:
                simulation.saveState(ensemble_options.saveState)

    @staticmethod
    def _decide_run_type(topology, system, simulation, ensemble_options):
        if isinstance(ensemble_options, RNEMDOptions):
            run_rnemd(simulation, ensemble_options)
        else:
            if ensemble_options.average_options is not None:
                simulation = run_average(topology, system, simulation, ensemble_options)
            elif ensemble_options.temperature_ramp_options is not None:
                simulation = run_temperature_ramp(simulation, ensemble_options)
            else:
                simulation.step(ensemble_options.steps)
        return simulation
