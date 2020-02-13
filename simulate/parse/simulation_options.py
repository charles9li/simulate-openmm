from ast import literal_eval

from ._options import _Options
from .integrator_options import *
from .reporter_options import *


class SimulationOptions(_Options):

    INTEGRATORS = {'LangevinIntegrator': LangevinIntegratorOptions,
                   'VelocityVerletIntegrator': VelocityVerletIntegratorOptions,
                   'VerletIntegrator': VerletIntegratorOptions}

    REPORTERS = {'DCDReporter': DCDReporterOptions}

    # =========================================================================

    def __init__(self):
        super(SimulationOptions, self).__init__()
        self.integrator = None
        self.minimizeEnergy = False
        self.steps = 0
        self.reporters = []

    # =========================================================================

    def _parse_integrator(self, *args):
        option_value = args[0]
        line_deque = args[1]
        integrator_options = self.INTEGRATORS[option_value]()
        integrator_options.parse(line_deque)
        self.integrator = integrator_options.integrator()

    def _parse_minimize_energy(self, *args):
        self.minimizeEnergy = literal_eval(args[0])

    def _parse_steps(self, *args):
        self.steps = literal_eval(args[0])

    def _parse_reporter(self, *args):
        pass

    OPTIONS = {'integrator': _parse_integrator,
               'minimizeEnergy': _parse_minimize_energy,
               'steps': _parse_steps,
               'reporter': _parse_reporter}
