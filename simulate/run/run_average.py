from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

import numpy as np
from simtk.openmm.app import Simulation
from simtk.openmm import MonteCarloBarostat, MonteCarloAnisotropicBarostat, MonteCarloMembraneBarostat
from simtk.unit import kilojoule_per_mole, nanometer

from simulate.utils.statistics import StatisticalInformation


def run_average(topology, system, simulation, ensemble_options):
    average_options = ensemble_options.average_options

    # initialize data
    average_data = {}
    if average_options.volume:
        average_data['volume'] = []
    if average_options.energy:
        average_data['energy'] = []

    # frequency at which data is stored
    store_frequency = 100

    # run simulation
    remaining_steps = ensemble_options.steps
    while remaining_steps > 0:
        if remaining_steps < store_frequency:
            simulation.step(remaining_steps)
            break
        simulation.step(store_frequency)
        remaining_steps -= store_frequency
        if average_options.volume or average_options.energy:
            state = simulation.context.getState(getEnergy=average_options.energy)
            if average_options.volume:
                average_data['volume'].append(state.getPeriodicBoxVolume().value_in_unit(nanometer**3))
            if average_options.energy:
                total_energy = state.getPotentialEnergy() + state.getKineticEnergy()
                average_data['energy'].append(total_energy.value_in_unit(kilojoule_per_mole))

    # average volume and energy
    tol = 1e-3              # TODO: add way to change this
    check_frequency = 25    # TODO: add way to change this
    if average_options.volume:
        volume_stats = StatisticalInformation(average_data['volume'])
        volume_mean = volume_stats.mean
        while True:
            simulation.step(check_frequency)
            volume_curr = simulation.context.getState().getPeriodicBoxVolume().value_in_unit(nanometer**3)
            if np.abs((volume_curr - volume_mean)/volume_mean) < tol:
                break
        simulation = _remove_barostat_from_simulation(topology, system, simulation, ensemble_options)
    if average_options.energy:
        energy_stats = StatisticalInformation(average_data['energy'])
        energy_mean = energy_stats.mean
        while True:
            simulation.step(check_frequency)
            state = simulation.context.getState(getEnergy=True)
            energy_curr = state.getPotentialEnergy() + state.getKineticEnergy()
            energy_curr = energy_curr.value_in_unit(kilojoule_per_mole)
            if np.abs((energy_curr - energy_mean)/energy_mean) < tol:
                break

    return simulation


def _remove_barostat_from_simulation(topology, system, simulation, ensemble_options):
    _remove_barostat_from_system(system)
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    positions = state.getPositions(asNumpy=True)
    velocities = state.getVelocities(asNumpy=True)
    periodic_box_vectors = state.getPeriodicBoxVectors(asNumpy=True)
    integrator = ensemble_options.create_integrator()
    system.setDefaultPeriodicBoxVectors(*periodic_box_vectors)
    simulation = Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    simulation.context.setVelocities(velocities)
    return simulation


def _remove_barostat_from_system(system):
    forces = system.getForces()
    for force in forces:
        if _is_barostat(force):
            barostat_index = forces.index(force)
            return system.removeForce(barostat_index)


def _is_barostat(force):
    barostats = [MonteCarloBarostat, MonteCarloAnisotropicBarostat, MonteCarloMembraneBarostat]
    for barostat in barostats:
        if isinstance(force, barostat):
            return True
    return False
