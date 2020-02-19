"""
rnemdreporter.py: Outputs x-component velocity profile in z-direction

Used for reporting velocity profiles when using the reverse
nonequilibrium molecular dynamics (RNEMD) technique to compute viscosity.

The Author claims no ownership over this software and the user is free to use
and distribute.
Authors: Charles Li

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import print_function
__author__ = "Charles Li"
__version__ = "1.0"

import simtk.openmm as mm
import simtk.unit as unit

import numpy as np


class RNEMDReporter(object):
    """RNEMDReprter outputs the x-component velocity profile, averaged in each slab, in the z direction to a csv file.

    To use it, create a RNEMDReporter, then add it to the Simulation's list of reporters. The number of slabs, M, can be
    specified in the constructor.
    """

    def __init__(self, file, reporterInterval, num_slabs, step=False):
        """Create a RNEMDReporter.

        Parameters
        ----------
        file : string or file
            The file to write to, specified as a file name or file object
        reporterInterval : int
            The interval (in time steps) at which to write frames
        num_slabs : int
            The number of slabs to divide the box into. Must be an even number
        """
        self._reportInterval = reporterInterval
        self._openedFile = isinstance(file, str)
        if self._openedFile:
            self._out = open(file, 'w')
        else:
            self._out = file
        if not isinstance(num_slabs, int):
            raise TypeError("The number of slabs must be an integer.")
        if num_slabs % 2 != 0:
            raise ValueError("The number of slabs must be an even integer.")
        self._num_slabs = num_slabs
        self._step = step
        self._separator = ','
        self._hasInitialized = False

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A five element tuple. The first element is the number of steps
            until the next report. The remaining elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, True, False, False)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if not self._hasInitialized:
            headers = self._constructHeaders()
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
        """Query the simulation for the current state of our observables of interest.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation

        Returns
        -------
        A list of values summarizing the current state of
        the simulation, to be printed or saved. Each element in the list
        corresponds to one of the columns in the resulting CSV file.
        """
        values = []
        if self._step:
            values.append(simulation.currentStep)

        # Get values from State
        Lz = state.getPeriodicBoxVectors()[2][2].value_in_unit(unit.nanometer)
        positions = state.getPositions(asNumpy=True)
        velocities = state.getVelocities(asNumpy=True)

        # Get z coordinates and x component velocities
        z_coordinates = [pos[2].value_in_unit(unit.nanometer) for pos in positions]
        x_velocities = [vel[0].value_in_unit(unit.nanometer/unit.picosecond) for vel in velocities]

        # Compute average vx for each slab
        x_velocities_by_slab = [[] for _ in range(self._num_slabs)]
        for i in range(len(z_coordinates)):
            z = z_coordinates[i]
            slab_num = int(z*self._num_slabs/Lz) % self._num_slabs
            x_velocities_by_slab[slab_num].append(x_velocities[i])
        for vx_list in x_velocities_by_slab:
            values.append(np.mean(vx_list))

        return values

    def _constructHeaders(self):
        """Construct the headers for the CSV output

        Returns: a list of strings giving the title of each observable being reported on.
        """
        headers = []
        if self._step:
            headers.append('Step')
        for i in range(self._num_slabs):
            headers.append('vx{} (nm/ps)'.format(i + 1))
        return headers

    def __del__(self):
        if self._openedFile:
            self._out.close()
