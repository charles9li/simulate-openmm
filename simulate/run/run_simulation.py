from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from simulate.parse import InputOptions
from simulate.parse.simulation.ensemble_options import RNEMDOptions
from simulate.run.run_rnemd import run_rnemd
from simulate.run.run_average import run_average


class RunSimulation(object):

    def __init__(self, input_file):
        self.input_options = InputOptions(input_file)
        self.system_options = self.input_options.system_options
        self.simulation_ensembles = self.input_options.simulation_ensembles
        self.ensembles = self.simulation_ensembles.ensembles
        self.positions = None
        self.velocities = None
        self.periodic_box_vectors = None

    def run(self):
        for ensemble_options in self.ensembles:

            # create system and topology
            topology = self.system_options.topology()
            system = self.system_options.create_system()

            # modify periodic box vectors if necessary
            if self.periodic_box_vectors is not None:
                system.setDefaultPeriodicBoxVectors(*self.periodic_box_vectors)

            # create simulation
            simulation = ensemble_options.simulation(topology, system)

            # initialize positions
            if self.positions is None:
                self.simulation_ensembles.set_positions(simulation)
            else:
                simulation.context.setPositions(self.positions)

            # initialize velocities
            if self.velocities is None:
                self.simulation_ensembles.set_velocities(simulation)
            else:
                simulation.context.setVelocities(self.velocities)

            # minimize energy
            minimize_energy_options = ensemble_options.minimize_energy_options
            if minimize_energy_options is not None:
                tolerance = minimize_energy_options.tolerance
                max_iterations = minimize_energy_options.maxIterations
                simulation.minimizeEnergy(tolerance=tolerance, maxIterations=max_iterations)

            # run simulation
            simulation = self._decide_run_type(simulation, ensemble_options)

            # store positions and velocities for next ensemble
            state = simulation.context.getState(getPositions=True, getVelocities=True)
            self.positions = state.getPositions(asNumpy=True)
            self.velocities = state.getVelocities(asNumpy=True)
            self.periodic_box_vectors = state.getPeriodicBoxVectors(asNumpy=True)

    @staticmethod
    def _decide_run_type(simulation, ensemble_options):
        if isinstance(ensemble_options, RNEMDOptions):
            run_rnemd(simulation, ensemble_options)
        else:
            if ensemble_options.average_options is None:
                simulation.step(ensemble_options.steps)
            else:
                simulation = run_average(simulation, ensemble_options)
        return simulation
