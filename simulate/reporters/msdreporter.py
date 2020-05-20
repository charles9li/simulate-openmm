from __future__ import print_function
__author__ = "Charles Li"
__version__ = "1.0"

from simtk import unit
import numpy as np


class MSDReporter(object):

    def __init__(self, file, reportInterval, step=False, enforcePeriodicBox=None):
        self._reportInterval = reportInterval
        self._openedFile = isinstance(file, str)
        if self._openedFile:
            self._out = open(file, 'w')
        else:
            self._out = file
        self._step = step
        self._enforcePeriodicBox = enforcePeriodicBox
        self._separator = ','
        self._hasInitialized = False
        self._initialPositions = None
        self._reportedMolecules = []
        self._reportedMoleculesMass = []

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False, self._enforcePeriodicBox)

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
            self._initialPositions = state.getPositions(asNumpy=True)

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
        values.append(state.getTime().value_in_unit(unit.picosecond))

        # Get values from simulation
        positions = state.getPositions(asNumpy=True)
        MSD_list = []
        for molecule, mass in zip(self._reportedMolecules, self._reportedMoleculesMass):
            MSD_list.append(self._computeMoleculeMSD(self._initialPositions[molecule], positions[molecule], mass))
        values.append(np.mean(MSD_list))

        return values

    def _computeMoleculeMSD(self, initial_positions, positions, mass):
        total_mass = np.sum(mass)
        mass_array = np.vstack((mass, mass, mass)).T
        initial_center_of_mass = np.sum(initial_positions*mass_array, axis=0) / total_mass
        center_of_mass = np.sum(positions*mass_array, axis=0) / total_mass
        relative_positions = center_of_mass - initial_center_of_mass
        return np.dot(relative_positions, relative_positions)

    def _constructHeaders(self, simulation):
        headers = []
        if self._step:
            headers.append('Step')
        headers.append('Time (ps)')
        headers.append("MSD (nm^2)")
        system = simulation.system
        for molecule in simulation.context.getMolecules():
            self._reportedMolecules.append(np.array(molecule))
            self._reportedMoleculesMass.append(self._getMoleculeMass(system, molecule))
        return headers

    @staticmethod
    def _getMoleculeMass(system, molecule):
        masses = []
        for atom_index in molecule:
            masses.append(system.getParticleMass(atom_index).value_in_unit(unit.amu))
        return np.array(masses)

    def __del__(self):
        if self._openedFile:
            self._out.close()
