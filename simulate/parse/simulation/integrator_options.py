from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from ast import literal_eval

from simtk.unit import femtosecond, kelvin, picosecond

from simulate.parse._options import _Options


class _IntegratorOptions(_Options):

    # =========================================================================

    def __init__(self):
        super(_IntegratorOptions, self).__init__()
        self.stepSize = None

    # =========================================================================

    def _parse_step_size(self, *args):
        self.stepSize = literal_eval(args[0])*femtosecond

    OPTIONS = {'stepSize': _parse_step_size}

    # =========================================================================

    def integrator(self):
        pass


class VerletIntegratorOptions(_IntegratorOptions):

    SECTION_NAME = "VerletIntegrator"

    # =========================================================================

    def __init__(self):
        super(VerletIntegratorOptions, self).__init__()
        self.stepSize = None

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.stepSize is None:
            self._incomplete_error('stepSize')

    # =========================================================================

    def integrator(self):
        from simtk.openmm import VerletIntegrator
        return VerletIntegrator(self.stepSize)


class VelocityVerletIntegratorOptions(VerletIntegratorOptions):

    SECTION_NAME = "VelocityVerletIntegrator"

    def __init__(self):
        super(VelocityVerletIntegratorOptions, self).__init__()
        self.stepSize = 1.0*femtosecond

    # =========================================================================

    def integrator(self):
        from openmmtools.integrators import VelocityVerletIntegrator
        return VelocityVerletIntegrator(self.stepSize)


class LangevinIntegratorOptions(_IntegratorOptions):

    SECTION_NAME = "LangevinIntegrator"

    def __init__(self):
        super(LangevinIntegratorOptions, self).__init__()
        self.temperature = None
        self.frictionCoeff = None
        self.stepSize = None

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.temperature is None:
            self._incomplete_error('temperature')
        if self.frictionCoeff is None:
            self._incomplete_error('frictionCoeff')
        if self.stepSize is None:
            self._incomplete_error('stepSize')

    # =========================================================================

    def _parse_temperature(self, *args):
        self.temperature = literal_eval(args[0])*kelvin

    def _parse_friction_coeff(self, *args):
        self.frictionCoeff = literal_eval(args[0])/picosecond

    def _parse_step_size(self, *args):
        super(LangevinIntegratorOptions, self)._parse_step_size(*args)

    OPTIONS = {'temperature': _parse_temperature,
               'frictionCoeff': _parse_friction_coeff,
               'stepSize': _parse_step_size}

    # =========================================================================

    def integrator(self):
        from simtk.openmm import LangevinIntegrator
        return LangevinIntegrator(self.temperature, self.frictionCoeff, self.stepSize)
