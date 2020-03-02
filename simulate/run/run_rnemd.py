from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from simtk.unit import amu, nanometer, picosecond


def run_rnemd(simulation, ensemble_options):
    remaining_steps = ensemble_options.steps
    M = ensemble_options.numSlabs
    W = ensemble_options.swapFrequency
    while remaining_steps > 0:
        if remaining_steps < W:
            simulation.step(remaining_steps)
            break
        simulation.step(W)
        del_px = _exchange_momentum(simulation, M)
        simulation.context.totalMomentumExchanged += del_px
        remaining_steps -= W


def _exchange_momentum(simulation, M):
    system = simulation.system
    context = simulation.context
    state = context.getState(getPositions=True, getVelocities=True)
    positions = state.getPositions(asNumpy=True)
    velocities = state.getVelocities(asNumpy=True)
    molecules = context.getMolecules()
    lower_particles, upper_particles = _separate_particles(molecules, system, positions, M)
    lower_molecule_indices, px_lower = _find_most_negative_px(lower_particles, system, velocities)
    upper_molecule_indices, px_upper = _find_most_positive_px(upper_particles, system, velocities)
    del_px = 0.0*amu*nanometer/picosecond
    if lower_molecule_indices is not None and upper_molecule_indices is not None:
        _set_momentum(px_upper, lower_molecule_indices, velocities, system)
        _set_momentum(px_lower, upper_molecule_indices, velocities, system)
        del_px = px_upper - px_lower
        simulation.context.setVelocities(velocities)
    return del_px


def _set_momentum(px_target, atom_indices, velocities, system):
    mass = _compute_molecule_mass(atom_indices, system)
    vx = px_target/mass
    for atom in atom_indices:
        velocities[atom] = vx


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
    return return_indices, curr_px


def _find_most_negative_px(molecules, system, velocities):
    curr_px = 0.0*amu*nanometer/picosecond
    return_indices = None
    for atom_indices in molecules:
        px = _compute_molecule_px(atom_indices, system, velocities)
        if px < curr_px:
            curr_px = px
            return_indices = atom_indices
    return return_indices, curr_px


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
