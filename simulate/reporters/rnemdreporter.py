"""
rnemdreporter.py: Outputs total exchanged momentum as a function of time

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

from simtk import unit


class RNEMDReporter(object):
    """RNEMDReporter outputs the total exchanged momentum and momentum flux
    used to compute viscosities in the RNMED method.

    To use it, create a RNMEDReporter, then add it to the Simulation's list of
    reporters. The number of slabs, M, can be specified in the constructor.
    """

    def __init__(self, file, reporterInterval, step=False):
        """Create a RNEMDVelocityReporter.

        Parameters
        ----------
        file : string or file
            The file to write to, specified as a file name or file object
        reporterInterval : int
            The interval (in time steps) at which to write frames
        step : bool
            True if also outputting current step
        """
        self._reportInterval = reporterInterval
        self._openedFile = isinstance(file, str)
        if self._openedFile:
            self._out = open(file, 'w')
        else:
            self._out = file
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
        return (steps, False, False, False, False)

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

        # Get values from simulation
        elapsed_ps = state.getTime().value_in_unit(unit.picosecond)
        Lx = state.getPeriodicBoxVectors()[0][0].value_in_unit(unit.nanometer)
        Ly = state.getPeriodicBoxVectors()[1][1].value_in_unit(unit.nanometer)
        total_momentum_exchanged = simulation.context.totalMomentumExchanged.value_in_unit(unit.amu*unit.nanometer/unit.picosecond)

        # Compute flux
        momentum_flux = total_momentum_exchanged/(2*Lx*Ly*elapsed_ps)

        # Add values
        values.append(elapsed_ps)
        values.append(total_momentum_exchanged)
        values.append(momentum_flux)

        return values

    def _constructHeaders(self):
        """Construct the headers for the CSV output.

        Returns: a list of strings giving the title of each observable being reported on.
        """
        headers = []
        if self._step:
            headers.append('Step')
        headers.append('Time (ps)')
        headers.append('Total Exchanged Momentum (amu*nm/ps)')
        headers.append('Momentum Flux (amu/nm*ps^2)')
        return headers

    def __del__(self):
        if self._openedFile:
            self._out.close()
