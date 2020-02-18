from ast import literal_eval

from simtk.unit import femtosecond, picosecond

from ._options import _Options


class _IntegratorOptions(_Options):

    # =========================================================================

    def __init__(self):
        super(_IntegratorOptions, self).__init__()
        self.stepSize = None

    # =========================================================================

    def _parse_step_size(self, *args):
        self.stepSize = literal_eval(args[0])*femtosecond

    # =========================================================================

    def integrator(self):
        pass

    # =========================================================================

    def _no_option_specified_exception(self, option_name):
        raise ValueError("No {} specified for {}".format(option_name, self.SECTION_NAME))


class LangevinIntegratorOptions(_IntegratorOptions):

    SECTION_NAME = "LangevinIntegrator"

    def __init__(self):
        super(LangevinIntegratorOptions, self).__init__()
        self.temperature = None
        self.frictionCoeff = None
        self.stepSize = None

    # =========================================================================

    def _parse_temperature(self, *args):
        self.temperature = literal_eval(args[0])

    def _parse_friction_coeff(self, *args):
        self.frictionCoeff = literal_eval(args[0])/picosecond

    def _parse_step_size(self, *args):
        super(LangevinIntegratorOptions, self)._parse_step_size(*args)

    OPTIONS = {'temperature': _parse_temperature,
               'frictionCoeff': _parse_friction_coeff,
               'stepSize': _parse_step_size}

    # =========================================================================

    def integrator(self):
        if self.temperature is None:
            self._no_option_specified_exception("temperature")
        if self.frictionCoeff is None:
            self._no_option_specified_exception("friction coefficient")
        if self.stepSize is None:
            self._no_option_specified_exception("step size")
        from simtk.openmm import LangevinIntegrator
        return LangevinIntegrator(self.temperature, self.frictionCoeff, self.stepSize)


class VerletIntegratorOptions(_IntegratorOptions):

    SECTION_NAME = "VerletIntegrator"

    # =========================================================================

    def __init__(self):
        super(VerletIntegratorOptions, self).__init__()
        self.stepSize = None

    # =========================================================================

    def _parse_step_size(self, *args):
        super(VerletIntegratorOptions, self)._parse_step_size(*args)

    OPTIONS = {'stepSize': _parse_step_size}

    # =========================================================================

    def integrator(self):
        if self.stepSize is None:
            self._no_option_specified_exception("step size")
        from simtk.openmm import VerletIntegrator
        return VerletIntegrator(self.stepSize)


class VelocityVerletIntegratorOptions(VerletIntegratorOptions):

    SECTION_NAME = "VelocityVerletIntegrator"

    def __init__(self):
        super(VelocityVerletIntegratorOptions, self).__init__()
        self.stepSize = 1.0*femtosecond

    # =========================================================================

    def _parse_step_size(self, *args):
        super(VelocityVerletIntegratorOptions, self)._parse_step_size(*args)

    OPTIONS = {'stepSize', _parse_step_size}

    # =========================================================================

    def integrator(self):
        from openmmtools.integrators import VelocityVerletIntegrator
        return VelocityVerletIntegrator(self.stepSize)
