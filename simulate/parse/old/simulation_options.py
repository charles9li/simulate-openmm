from ._options import _Options
from .integrator_options import *
from .barostat_options import *
from .reporter_options import *


class SimulationOptions(_Options):

    INTEGRATORS = {'LangevinIntegrator': LangevinIntegratorOptions,
                   'VelocityVerletIntegrator': VelocityVerletIntegratorOptions,
                   'VerletIntegrator': VerletIntegratorOptions}

    BAROSTATS = {'MonteCarloBarostat': MonteCarloBarostatOptions}

    REPORTERS = {'DCDReporter': DCDReporterOptions,
                 'StateDataReporter': StateDataReporterOptions}

    # =========================================================================

    def __init__(self):
        super(SimulationOptions, self).__init__()
        self.integrator = None
        self.barostat = None
        self.minimizeEnergy = False
        self.steps = 0
        self.reporters = []

    # =========================================================================

    def _parse_integrator(self, *args):
        option_value = args[0]
        line_deque = args[1]
        integrator_options = self.INTEGRATORS[option_value]()
        integrator_options.parse(line_deque.popleft())
        self.integrator = integrator_options.integrator()

    def _parse_barostat(self, *args):
        option_value = args[0]
        line_deque = args[1]
        barostat_options = self.BAROSTATS[option_value]()
        barostat_options.parse(line_deque.popleft())
        self.barostat = barostat_options.barostat()

    def _parse_minimize_energy(self, *args):
        self.minimizeEnergy = literal_eval(args[0])

    def _parse_steps(self, *args):
        self.steps = literal_eval(args[0])

    def _parse_reporter(self, *args):
        option_value = args[0]
        line_deque = args[1]
        reporter_options = self.REPORTERS[option_value]()
        reporter_options.parse(line_deque.popleft())
        reporter = reporter_options.reporter()
        self.reporters.append(reporter)

    OPTIONS = {'integrator': _parse_integrator,
               'barostat': _parse_barostat,
               'minimizeEnergy': _parse_minimize_energy,
               'steps': _parse_steps,
               'reporter': _parse_reporter}
