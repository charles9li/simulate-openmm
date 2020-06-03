from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from simtk.openmm import MonteCarloBarostat, MonteCarloAnisotropicBarostat, MonteCarloMembraneBarostat
from simtk.unit import kelvin


def run_temperature_ramp(simulation, ensemble_options):
    temperature_ramp_options = ensemble_options.temperature_ramp_options
    start_temperature = temperature_ramp_options.startTemperature
    end_temperature = temperature_ramp_options.endTemperature
    step_temperature = temperature_ramp_options.stepTemperature
    steps_per_interval = temperature_ramp_options.stepsPerInterval
    current_temperature = start_temperature
    while not is_done(current_temperature, end_temperature, step_temperature):
        simulation.integrator.setTemperature(current_temperature)
        _set_barostat_temperature(simulation, current_temperature)
        simulation.context.reinitialize(preserveState=True)
        simulation.step(steps_per_interval)
        current_temperature += step_temperature
    return simulation


def is_done(current_temperature, end_temperature, step_temperature):
    if step_temperature < 0.0*kelvin:
        return current_temperature < end_temperature
    else:
        return current_temperature > end_temperature


def _set_barostat_temperature(simulation, current_temperature):
    system = simulation.system
    barostat = _get_barostat(system)
    barostat.setDefaultTemperature(current_temperature)


def _get_barostat(system):
    forces = system.getForces()
    for force in forces:
        if _is_barostat(force):
            barostat_index = forces.index(force)
            return forces[barostat_index]


def _is_barostat(force):
    barostats = [MonteCarloBarostat, MonteCarloAnisotropicBarostat, MonteCarloMembraneBarostat]
    for barostat in barostats:
        if isinstance(force, barostat):
            return True
    return False
