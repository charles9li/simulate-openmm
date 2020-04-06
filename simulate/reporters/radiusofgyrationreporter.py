from __future__ import absolute_import, division, print_function
__author__ = "Charles Li"
__version__ = "1.0"

import numpy as np
from simtk import unit


class RadiusOfGyrationReporter(object):

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
            values.append(self._computeRadiusOfGyration(positions[molecule], mass))

        return values

    def _computeRadiusOfGyration(self, positions, mass):
        total_mass = np.sum(mass)
        center_of_mass = np.sum(positions*np.vstack((mass, mass, mass)).T, axis=0)/total_mass
        relative_positions = positions - center_of_mass
        return np.sqrt(np.sum(mass*np.sum(relative_positions*relative_positions, axis=1))/total_mass)

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
                self._reportedMolecules.append(np.array(molecule))
                self._reportedMoleculesMass.append(self._getMoleculeMass(system, molecule))
        return headers

    @staticmethod
    def _getMoleculeMass(system, molecule):
        masses = []
        for atom_index in molecule:
            masses.append(system.getParticleMass(atom_index).value_in_unit(unit.amu))
        return np.array(masses)

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
