from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from simtk.openmm.app import Simulation
from simtk.openmm import MonteCarloBarostat, MonteCarloAnisotropicBarostat, MonteCarloMembraneBarostat
from simtk.unit import kelvin


def run_temperature_ramp(simulation, ensemble_options):
    system_options = ensemble_options.simulations_options.input_options.system_options
    temperature_ramp_options = ensemble_options.temperature_ramp_options
    start_temperature = temperature_ramp_options.startTemperature
    end_temperature = temperature_ramp_options.endTemperature
    step_temperature = temperature_ramp_options.stepTemperature
    steps_per_interval = temperature_ramp_options.stepsPerInterval
    current_temperature = start_temperature
    while not is_done(current_temperature, end_temperature, step_temperature):
        simulation = _create_new_simulation(simulation, current_temperature, system_options, ensemble_options)
        simulation.step(steps_per_interval)
        current_temperature += step_temperature
    return simulation


def is_done(current_temperature, end_temperature, step_temperature):
    if step_temperature < 0.0*kelvin:
        return current_temperature < end_temperature
    else:
        return current_temperature > end_temperature


def _create_new_simulation(simulation_old, current_temperature, system_options, ensemble_options):

    # retrieve information from old simulation
    topology = simulation_old.topology
    system_old = simulation_old.system
    state = simulation_old.context.getState(getPositions=True, getVelocities=True, getForces=True)
    positions = state.getPositions(asNumpy=True)
    velocities = state.getVelocities(asNumpy=True)
    periodic_box_vectors = state.getPeriodicBoxVectors(asNumpy=True)

    # create new system
    system_new = system_options.create_system_with_new_topology(topology)
    if _has_barostat(system_old):
        barostat_old = _get_barostat(system_old)
        default_pressure = barostat_old.getDefaultPressure()
        frequency = barostat_old.getFrequency()
        barostat_new = MonteCarloBarostat(default_pressure, current_temperature, frequency)
        system_new.addForce(barostat_new)
        barostat_new.setForceGroup(system_new.getNumForces() - 1)

    # create new integrator and set integrator temperature
    integrator_new = ensemble_options.create_integrator()
    integrator_new.setTemperature(current_temperature)

    # create new simulation
    simulation_new = Simulation(topology, system_new, integrator_new)

    # set positions, velocities, and periodic box vectors of new simulation
    simulation_new.context.setPositions(positions)
    simulation_new.context.setVelocities(velocities)
    simulation_new.context.setPeriodicBoxVectors(*periodic_box_vectors)
    simulation_new.currentStep = simulation_old.currentStep
    simulation_new.context.setTime(state.getTime())

    # add reporters from old simulation to new simulation
    for reporter in simulation_old.reporters:
        simulation_new.reporters.append(reporter)

    return simulation_new


def _get_barostat(system):
    forces = system.getForces()
    for force in forces:
        if _is_barostat(force):
            barostat_index = forces.index(force)
            return system.getForce(barostat_index)


def _has_barostat(system):
    forces = system.getForces()
    for force in forces:
        if _is_barostat(force):
            return True
    return False


def _is_barostat(force):
    barostats = [MonteCarloBarostat, MonteCarloAnisotropicBarostat, MonteCarloMembraneBarostat]
    for barostat in barostats:
        if isinstance(force, barostat):
            return True
    return False
