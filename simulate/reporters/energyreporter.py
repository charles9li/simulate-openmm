from __future__ import print_function
__author__ = "Charles Li"
__version__ = "1.0"


from simtk import unit


class EnergyReporter(object):

    def __init__(self, file, reporterInterval):
        self._reportInterval = reporterInterval
        self._openedFile = isinstance(file, str)
        if self._openedFile:
            self._out = open(file, 'w')
        else:
            self._out = file
        self._separator = ','
        self._hasInitialized = False

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, True)

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
        values.append(simulation.currentStep)
        values.append(state.getTime().value_in_unit(unit.picosecond))
        for i in range(simulation.system.getNumForces()):
            state = simulation.context.getState(getEnergy=True, groups={i})
            values.append(self._getEnergy(state).value_in_unit(unit.kilojoule_per_mole))
        return values

    @staticmethod
    def _getEnergy(state):
        return state.getPotentialEnergy() + state.getKineticEnergy()

    def _constructHeaders(self, simulation):
        headers = []
        headers.append('Step')
        headers.append('Time (ps)')
        system = simulation.system
        for i in range(system.getNumForces()):
            headers.append('{} (kJ/mol)'.format(system.getForce(i).__class__.__name__))
        return headers

    def __del__(self):
        if self._openedFile:
            self._out.close()


class PotentialEnergyReporter(EnergyReporter):

    @staticmethod
    def _getEnergy(state):
        return state.getPotentialEnergy()


class KineticEnergyReporter(EnergyReporter):

    @staticmethod
    def _getEnergy(state):
        return state.getKineticEnergy()
