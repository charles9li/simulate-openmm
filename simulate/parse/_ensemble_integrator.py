"""
_ensemble_integrator.py: Parses integrator options for an ensemble.

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
from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from ast import literal_eval

from simtk.unit import femtosecond, kelvin, picosecond

from ._options import _Options

__all__ = ['VerletIntegratorOptions', 'VelocityVerletIntegratorOptions',
           'LangevinIntegratorOptions', 'NoseHooverChainVelocityVerletIntegratorOptions']


class _IntegratorOptions(_Options):

    # =========================================================================

    def __init__(self):
        super(_IntegratorOptions, self).__init__()
        self.stepSize = None
        self.constraintTolerance = 1.0e-5

    def _create_options(self):
        super(_IntegratorOptions, self)._create_options()
        self._OPTIONS['stepSize'] = self._parse_step_size
        self._OPTIONS['constraintTolerance'] = self._parse_constraint_tolerance

    # =========================================================================

    def _parse_step_size(self, *args):
        self.stepSize = literal_eval(args[0])*femtosecond

    def _parse_constraint_tolerance(self, *args):
        self.constraintTolerance = literal_eval(args[0])

    # =========================================================================

    def _create_integrator(self):
        from simtk.openmm import VerletIntegrator
        return VerletIntegrator(self.stepSize)

    def integrator(self):
        integrator = self._create_integrator()
        integrator.setConstraintTolerance(self.constraintTolerance)
        return integrator


class VerletIntegratorOptions(_IntegratorOptions):

    _SECTION_NAME = "VerletIntegrator"

    # =========================================================================

    def __init__(self):
        super(VerletIntegratorOptions, self).__init__()

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.stepSize is None:
            self._incomplete_error('stepSize')


class VelocityVerletIntegratorOptions(VerletIntegratorOptions):

    _SECTION_NAME = "VelocityVerletIntegrator"

    def __init__(self):
        super(VelocityVerletIntegratorOptions, self).__init__()
        self.stepSize = 1.0*femtosecond

    # =========================================================================

    def _create_integrator(self):
        from openmmtools.integrators import VelocityVerletIntegrator
        return VelocityVerletIntegrator(self.stepSize)


class LangevinIntegratorOptions(_IntegratorOptions):

    _SECTION_NAME = "LangevinIntegrator"

    def __init__(self):
        super(LangevinIntegratorOptions, self).__init__()
        self.temperature = None
        self.frictionCoeff = None

    def _create_options(self):
        super(LangevinIntegratorOptions, self)._create_options()
        self._OPTIONS['temperature'] = self._parse_temperature
        self._OPTIONS['frictionCoeff'] = self._parse_friction_coeff

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

    # =========================================================================

    def _create_integrator(self):
        from simtk.openmm import LangevinIntegrator
        return LangevinIntegrator(self.temperature, self.frictionCoeff, self.stepSize)


class NoseHooverChainVelocityVerletIntegratorOptions(_IntegratorOptions):

    def __init__(self):
        super(NoseHooverChainVelocityVerletIntegratorOptions, self).__init__()
        self.temperature = 298*kelvin
        self.collision_frequency = 50/picosecond
        self.chain_length = 5
        self.num_mts = 5
        self.num_yoshidasuzuki = 5

    def _create_options(self):
        super(NoseHooverChainVelocityVerletIntegratorOptions, self)._create_options()
        self._OPTIONS['temperature'] = self._parse_temperature
        self._OPTIONS['collision_frequency'] = self._parse_collision_frequency
        self._OPTIONS['chain_length'] = self._parse_chain_length
        self._OPTIONS['num_mts'] = self._parse_num_mts
        self._OPTIONS['num_yoshidasuzuki'] = self._parse_num_yoshidasuzuki

    # =========================================================================

    def _parse_temperature(self, *args):
        self.temperature = literal_eval(args[0])*kelvin

    def _parse_collision_frequency(self, *args):
        self.collision_frequency = literal_eval(args[0])/picosecond

    def _parse_chain_length(self, *args):
        self.chain_length = literal_eval(args[0])

    def _parse_num_mts(self, *args):
        self.num_mts = literal_eval(args[0])

    def _parse_num_yoshidasuzuki(self, *args):
        self.num_yoshidasuzuki = literal_eval(args[0])

    # =========================================================================

    def _create_integrator(self):
        from openmmtools.integrators import NoseHooverChainVelocityVerletIntegrator
        return NoseHooverChainVelocityVerletIntegrator(self.temperature, self.collision_frequency, self.stepSize,
                                                       self.chain_length, self.num_mts, self.num_yoshidasuzuki)
