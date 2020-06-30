from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from simtk.unit import amu, nanometer, picosecond
import numpy as np


# Units
MASS_UNIT = amu
VELOCITY_UNIT = nanometer / picosecond
MOMENTUM_UNIT = MASS_UNIT * VELOCITY_UNIT


def run_rnemd(simulation, ensemble_options):
    remaining_steps = ensemble_options.steps
    M = ensemble_options.numSlabs
    W = ensemble_options.swapFrequency
    while remaining_steps > 0:
        if remaining_steps < W:
            simulation.step(remaining_steps)
            break
        simulation.step(W)
        del_p = _exchange_momentum(simulation, M)
        simulation.context.totalMomentumExchanged += del_p
        simulation.context.xMomentumExchanged += del_p[0]
        simulation.context.yMomentumExchanged += del_p[1]
        simulation.context.zMomentumExchanged += del_p[2]
        remaining_steps -= W


def _exchange_momentum(simulation, M):

    # Get simulation information
    system = simulation.system
    context = simulation.context
    state = context.getState(getPositions=True, getVelocities=True)
    positions = state.getPositions(asNumpy=True)
    velocities = state.getVelocities(asNumpy=True)
    molecules = context.getMolecules()

    # Find which molecules are in upper and lower
    lower_particles, upper_particles = _separate_particles(molecules, system, positions, M)

    # Find molecule with most negative px in lower slab
    lower_molecule_indices = _find_most_negative_px(lower_particles, system, velocities)

    # Find molecule with most positive px in upper slab
    upper_molecule_indices = _find_most_positive_px(upper_particles, system, velocities)

    # Initialize change in momentum in each direction
    del_p = np.array([0.0] * 3) * MOMENTUM_UNIT

    # Swap momentum
    if lower_molecule_indices is not None and upper_molecule_indices is not None:

        # Compute molecule masses
        lower_molecule_mass = _compute_molecule_mass(lower_molecule_indices, system)
        upper_molecule_mass = _compute_molecule_mass(upper_molecule_indices, system)

        # Compute molecule momentum
        lower_molecule_momentum = _compute_molecule_p(lower_molecule_indices, system, velocities)
        upper_molecule_momentum = _compute_molecule_p(upper_molecule_indices, system, velocities)

        # Compute center of mass velocities
        lower_molecule_vel = lower_molecule_momentum / lower_molecule_mass
        upper_molecule_vel = upper_molecule_momentum / upper_molecule_mass

        # Swap velocities
        for atom in lower_molecule_indices:
            velocities[atom] = velocities[atom] - lower_molecule_vel + upper_molecule_vel
        for atom in upper_molecule_indices:
            velocities[atom] = velocities[atom] - upper_molecule_vel + lower_molecule_vel

        # Compute change in momentum
        del_p = upper_molecule_momentum - lower_molecule_momentum

        # Set new velocities
        simulation.context.setVelocities(velocities)
    return del_p


def _separate_particles(molecules, system, positions, M):
    lz = system.getDefaultPeriodicBoxVectors()[2][2]
    lower_particles = []
    upper_particles = []
    for atom_indices in molecules:
        cm_z_coord = _compute_cm_z_coord(atom_indices, system, positions)
        if 0.0*nanometer <= cm_z_coord < lz/M:
            lower_particles.append(atom_indices)
        elif lz/2 <= cm_z_coord < (0.5 + 1.0/M)*lz:
            upper_particles.append(atom_indices)
    return lower_particles, upper_particles


# def _set_momentum(px_target, atom_indices, velocities, system):
#     mass = _compute_molecule_mass(atom_indices, system)
#     vx = px_target/mass
#     for atom in atom_indices:
#         velocities[atom] = vx

def _swap_velocities(lower_molecule_indices, upper_molecule_indices, velocities):
    for i in range(len(lower_molecule_indices)):
        lower_atom_index = lower_molecule_indices[i]
        upper_atom_index = upper_molecule_indices[i]
        vel_temp = velocities[lower_atom_index][0]
        velocities[lower_atom_index][0] = velocities[upper_atom_index][0]
        velocities[upper_atom_index][0] = vel_temp


def _compute_cm_z_coord(atom_indices, system, positions):
    numerator = 0.0*amu*nanometer
    denominator = 0.0*amu
    for atom in atom_indices:
        mass = system.getParticleMass(atom)
        z = positions[atom][2]
        numerator += mass*z
        denominator += mass
    return numerator/denominator


def _find_most_positive_px(molecules, system, velocities):
    curr_px = 0.0*amu*nanometer/picosecond
    return_indices = None
    for atom_indices in molecules:
        px = _compute_molecule_px(atom_indices, system, velocities)
        if px > curr_px:
            curr_px = px
            return_indices = atom_indices
    return return_indices


def _find_most_negative_px(molecules, system, velocities):
    curr_px = 0.0*amu*nanometer/picosecond
    return_indices = None
    for atom_indices in molecules:
        px = _compute_molecule_px(atom_indices, system, velocities)
        if px < curr_px:
            curr_px = px
            return_indices = atom_indices
    return return_indices


def _compute_molecule_p(atom_indices, system, velocities):
    p = np.array([0.0] * 3) * MOMENTUM_UNIT
    for atom in atom_indices:
        mass = system.getParticleMass(atom)
        v = velocities[atom]
        p += mass * v
    return p


def _compute_molecule_px(atom_indices, system, velocities):
    px = 0.0*amu*nanometer/picosecond
    for atom in atom_indices:
        mass = system.getParticleMass(atom)
        vx = velocities[atom][0]
        px += mass*vx
    return px


def _compute_molecule_vx(atom_indices, system, velocities):
    numerator = 0.0*amu*nanometer/picosecond
    denominator = 0.0*amu
    for atom in atom_indices:
        mass = system.getParticleMass(atom)
        vx = velocities[atom][0]
        numerator += mass*vx
        denominator += mass
    return numerator/denominator


def _compute_molecule_mass(atom_indices, system):
    mass = 0.0*amu
    for atom in atom_indices:
        mass += system.getParticleMass(atom)
    return mass
