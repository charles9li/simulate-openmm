from __future__ import absolute_import, division, print_function
__author__ = "Charles Li"
__version__ = "1.0"

import numpy as np
from simtk import unit


class EndToEndDistanceReporter(object):

    def __init__(self, file, reporterInterval, step=False):
        self._reportInterval = reporterInterval
        self._openedFile = isinstance(file, str)
        if self._openedFile:
            self._out = open(file, 'w')
        else:
            self._out = file
        self._step = step
        self._separator = ','
        self._hasInitialized = False
        self._reportedMolecules = []
        self._reportedMoleculesMass = []

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False)

    def report(self, simulation, state):
        if not self._hasInitialized:
            headers = self._constructHeaders(simulation)
            print('#"{}"'.format(('"' + self._separator + '"').join(headers)), file=self._out)
            try:
                self._out.flush()
            except AttributeError:
                pass
            self._initialSteps = simulation.currentStep
            self._hasInitialized = True

        # Query for the values
        values = self._constructReportValues(simulation, state)

        # Write the values
        print(self._separator.join(str(v) for v in values), file=self._out)
        try:
            self._out.flush()
        except AttributeError:
            pass

    def _constructReportValues(self, simulation, state):
        values = []
        if self._step:
            values.append(simulation.currentStep)

        # Get values from simulation
        positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        values.append(state.getTime().value_in_unit(unit.picosecond))
        for molecule, mass in zip(self._reportedMolecules, self._reportedMoleculesMass):
            first_residue_atoms = molecule[0]
            last_residue_atoms = molecule[1]
            first_residue_masses = mass[0]
            last_residue_masses = mass[1]
            first_residue_cm = self._computeCenterOfMass(positions[first_residue_atoms], first_residue_masses)
            last_residue_cm = self._computeCenterOfMass(positions[last_residue_atoms], last_residue_masses)
            values.append(np.linalg.norm(first_residue_cm - last_residue_cm))
        return values

    @staticmethod
    def _computeCenterOfMass(positions, mass):
        return np.sum(positions*np.vstack((mass, mass, mass)).T, axis=0)/np.sum(mass)

    def _constructHeaders(self, simulation):
        headers = []
        if self._step:
            headers.append('Step')
        headers.append('Time (ps)')
        system = simulation.system
        topology = simulation.topology
        for molecule in simulation.context.getMolecules():
            if self._hasMultipleResidues(topology, molecule):
                headers.append('{} (nm)'.format(list(topology.atoms())[molecule[0]].residue.chain.id))
                first_residue_atoms = self._findFirstResidue(topology, molecule)
                last_residue_atoms = self._findLastResidue(topology, molecule)
                first_residue_masses = self._getResidueMasses(system, first_residue_atoms)
                last_residue_masses = self._getResidueMasses(system, last_residue_atoms)
                self._reportedMolecules.append((first_residue_atoms, last_residue_atoms))
                self._reportedMoleculesMass.append((first_residue_masses, last_residue_masses))
        return headers

    @staticmethod
    def _getResidueMasses(system, residue_atoms):
        masses = []
        for atom_index in residue_atoms:
            atom_index = int(atom_index)
            masses.append(system.getParticleMass(atom_index).value_in_unit(unit.amu))
        return np.array(masses)

    @staticmethod
    def _findFirstResidue(topology, molecule):
        first_residue = None
        first_residue_atoms = []
        atoms = list(topology.atoms())
        for atom_index in molecule:
            if first_residue is None:
                first_residue = atoms[atom_index].residue
            current_residue = atoms[atom_index].residue
            if current_residue == first_residue:
                first_residue_atoms.append(atom_index)
            else:
                return np.array(first_residue_atoms)

    @staticmethod
    def _findLastResidue(topology, molecule):
        current_residue = None
        current_residue_atoms = []
        atoms = list(topology.atoms())
        for atom_index in molecule:
            if current_residue is None:
                current_residue = atoms[atom_index].residue
            next_residue = atoms[atom_index].residue
            if next_residue == current_residue:
                current_residue_atoms.append(atom_index)
            else:
                current_residue = next_residue
                current_residue_atoms = []
        return np.array(current_residue_atoms)

    @staticmethod
    def _hasMultipleResidues(topology, molecule):
        previous_residue = None
        atoms = list(topology.atoms())
        for atom_index in molecule:
            if previous_residue is None:
                previous_residue = atoms[atom_index].residue
            current_residue = atoms[atom_index].residue
            if current_residue != previous_residue:
                return True
        return False

    def __del__(self):
        if self._openedFile:
            self._out.close()
