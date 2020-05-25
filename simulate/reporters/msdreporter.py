"""
msdreporter.py: Outputs mean-squared displacement of molecules in simulation.

Copyright (c) 2020 Charles Li // UCSB, Department of Chemical Engineering

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
from __future__ import print_function
__author__ = "Charles Li"
__version__ = "1.0"

from collections import OrderedDict
import os

from simtk.openmm.app import DCDFile, PDBFile
from simtk import unit
import mdtraj as md
import numpy as np
import pandas as pd


class MSDReporter(object):
    """MSD"Reporter outputs the mean-squared displacement from a Simulation.

    To use it, create a MSDReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, step=False):
        """Create a MSDReporter.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        """

        # Set report interval
        self._reportInterval = reportInterval

        # Set file name
        self._fileName = file

        # Create temporary dcd file
        self._tempDCDFileName = file.split('.')[0] + "_tempDCD.dcd"
        self._tempOut = open(self._tempDCDFileName, 'wb')
        self._tempDCD = None

        # Create temporary pdb filepath
        self._tempPDBFileName = file.split('.')[0] + "_tempPDB.pdb"

        self._step = step
        self._separator = ','
        self._reportedMolecules = []
        self._reportedMoleculesMass = []

        self._times = []
        if self._step:
            self._steps = []

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A six element tuple. The first element is the number of steps
            until the next report. The next four elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.  The final element specifies whether
            positions should be wrapped to lie in a single periodic box.
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False, False)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if self._tempDCD is None:

            # Create DCDFile object
            self._tempDCD = DCDFile(
                self._tempOut, simulation.topology, simulation.integrator.getStepSize(),
                simulation.currentStep, self._reportInterval, False
            )

            # Save molecules and masses
            for molecule in simulation.context.getMolecules():
                self._reportedMolecules.append(np.array(molecule))
                self._reportedMoleculesMass.append(self._getMoleculeMass(simulation.system, molecule))

            # Create temporary pdb file
            PDBFile.writeFile(simulation.topology, state.getPositions(), open(self._tempPDBFileName, 'w'))

        # Get time and step from simulation
        self._times.append(state.getTime().value_in_unit(unit.picosecond))
        if self._step:
            self._steps.append(simulation.currentStep)

        # Write to dcd file
        self._tempDCD.writeModel(state.getPositions(), periodicBoxVectors=state.getPeriodicBoxVectors())

    @staticmethod
    def _getMoleculeMass(system, molecule):
        masses = []
        for atom_index in molecule:
            masses.append(system.getParticleMass(atom_index).value_in_unit(unit.amu))
        return np.array(masses)

    def __del__(self):

        # Get positions from pdb file
        traj = md.load_dcd(self._tempDCDFileName, top=self._tempPDBFileName)

        # Compute MSD for each frame
        MSD_list = []
        initial_positions = traj.xyz[0]
        for positions in traj.xyz:
            molecule_MSD_list = []
            for molecule, mass in zip(self._reportedMolecules, self._reportedMoleculesMass):
                molecule_MSD = self._computeMoleculeMSD(initial_positions[molecule],
                                                        positions[molecule],
                                                        mass)
                molecule_MSD_list.append(molecule_MSD)
            MSD_list.append(np.mean(molecule_MSD_list))

        # Save data to csv file
        odict = OrderedDict()
        odict['Time (ps)'] = self._times
        if self._step:
            odict['Step'] = self._steps
        odict['MSD (nm^2)'] = MSD_list
        df = pd.DataFrame(data=odict)
        df.to_csv(path_or_buf=self._fileName, index=False)

        # Close temporary dcd file
        self._tempOut.close()

        # Delete temporary dcd file
        try:
            os.remove(self._tempDCDFileName)
        except OSError:
            print("Temporary DCD file not found.")

        # Delete temporary pdb file
        try:
            os.remove(self._tempPDBFileName)
        except OSError:
            print("Temporary PDB file not found.")

    def _computeMoleculeMSD(self, initial_positions, positions, mass):
        total_mass = np.sum(mass)
        mass_array = np.vstack((mass, mass, mass)).T
        initial_center_of_mass = np.sum(initial_positions*mass_array, axis=0) / total_mass
        center_of_mass = np.sum(positions*mass_array, axis=0) / total_mass
        relative_positions = center_of_mass - initial_center_of_mass
        return np.dot(relative_positions, relative_positions)
