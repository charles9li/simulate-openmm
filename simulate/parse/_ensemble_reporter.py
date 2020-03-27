"""
_ensemble_reporter.py: Parses reporter options for an ensemble.

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

from ._options import _Options


class _ReporterOptions(_Options):

    _SECTION_NAME = "_ReporterOptions"

    def __init__(self):
        super(_ReporterOptions, self).__init__()
        self.file = None
        self.reportInterval = None

    def _create_options(self):
        super(_ReporterOptions, self)._create_options()
        self._OPTIONS['file'] = self._parse_file
        self._OPTIONS['reportInterval'] = self._parse_report_interval

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.file is None:
            self._incomplete_error('file')
        if self.reportInterval is None:
            self._incomplete_error('reportInterval')

    # =========================================================================

    def _parse_file(self, *args):
        self.file = args[0]

    def _parse_report_interval(self, *args):
        self.reportInterval = literal_eval(args[0])

    # =========================================================================

    def reporter(self):
        pass


class PDBReporterOptions(_ReporterOptions):

    _SECTION_NAME = "PDBReporter"

    # =========================================================================

    def __init__(self):
        super(PDBReporterOptions, self).__init__()
        self.enforcePeriodicBox = None

    def _create_options(self):
        super(PDBReporterOptions, self)._create_options()
        self._OPTIONS['enforcePeriodicBox'] = self._parse_enforce_periodic_box

    # =========================================================================

    def _parse_enforce_periodic_box(self, *args):
        self.enforcePeriodicBox = literal_eval(args[0])

    # =========================================================================

    def reporter(self):
        from simtk.openmm.app import PDBReporter
        return PDBReporter(self.file, self.reportInterval, enforcePeriodicBox=self.enforcePeriodicBox)


class DCDReporterOptions(PDBReporterOptions):

    _SECTION_NAME = "DCDReporter"

    # =========================================================================

    def __init__(self):
        super(DCDReporterOptions, self).__init__()
        self.append = False

    def _create_options(self):
        super(DCDReporterOptions, self)._create_options()
        self._OPTIONS['append'] = self._parse_append

    # =========================================================================

    def _parse_append(self, *args):
        self.append = literal_eval(args[0])

    # =========================================================================

    def reporter(self):
        from simtk.openmm.app import DCDReporter
        return DCDReporter(self.file, self.reportInterval,
                           append=self.append, enforcePeriodicBox=self.enforcePeriodicBox)


class StateDataReporterOptions(_ReporterOptions):

    _SECTION_NAME = "StateDataReporter"

    # =========================================================================

    def __init__(self):
        super(StateDataReporterOptions, self).__init__()

    # =========================================================================

    def reporter(self):
        from simtk.openmm.app import StateDataReporter
        return StateDataReporter(self.file, self.reportInterval, step=True, time=True,
                                 potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                 temperature=True, volume=True, density=True,
                                 speed=True, elapsedTime=True)


class EnergyReporterOptions(_ReporterOptions):

    _SECTION_NAME = "EnergyReporter"

    # =========================================================================

    def __init__(self):
        super(EnergyReporterOptions, self).__init__()

    # =========================================================================

    def reporter(self):
        from simulate.reporters import EnergyReporter
        return EnergyReporter(self.file, self.reportInterval)


class PotentialEnergyReporterOptions(_ReporterOptions):

    _SECTION_NAME = "PotentialEnergyReporter"

    # =========================================================================

    def __init__(self):
        super(PotentialEnergyReporterOptions, self).__init__()

    # =========================================================================

    def reporter(self):
        from simulate.reporters import PotentialEnergyReporter
        return PotentialEnergyReporter(self.file, self.reportInterval)


class KineticEnergyReporterOptions(_ReporterOptions):

    _SECTION_NAME = "KineticEnergyReporter"

    # =========================================================================

    def __init__(self):
        super(KineticEnergyReporterOptions, self).__init__()

    # =========================================================================

    def reporter(self):
        from simulate.reporters import KineticEnergyReporter
        return KineticEnergyReporter(self.file, self.reportInterval)


class RNEMDReporterOptions(_ReporterOptions):

    _SECTION_NAME = "RNEMDReporter"

    # =========================================================================

    def __init__(self):
        super(RNEMDReporterOptions, self).__init__()

    # =========================================================================

    def reporter(self):
        from simulate.reporters import RNEMDReporter
        return RNEMDReporter(self.file, self.reportInterval, step=True)


class RNEMDVelocityReporterOptions(_ReporterOptions):

    _SECTION_NAME = "RNEMDVelocityReporter"

    # =========================================================================

    def __init__(self):
        super(RNEMDVelocityReporterOptions, self).__init__()
        self.numSlabs = None

    def _create_options(self):
        super(RNEMDVelocityReporterOptions, self)._create_options()
        self._OPTIONS['numSlabs'] = self._parse_num_slabs

    # =========================================================================
    
    def _check_for_incomplete_input(self):
        super(RNEMDVelocityReporterOptions, self)._check_for_incomplete_input()
        if self.numSlabs is None:
            self._incomplete_error('numSlabs')

    # =========================================================================

    def _parse_num_slabs(self, *args):
        self.numSlabs = literal_eval(args[0])

    # =========================================================================

    def reporter(self):
        from simulate.reporters import RNEMDVelocityReporter
        return RNEMDVelocityReporter(self.file, self.reportInterval, self.numSlabs, step=True)


class CheckpointReporterOptions(_ReporterOptions):

    _SECTION_NAME = "CheckPointReporter"

    # =========================================================================

    def __init__(self):
        super(CheckpointReporterOptions, self).__init__()

    # =========================================================================

    def reporter(self):
        from simtk.openmm.app import CheckpointReporter
        return CheckpointReporter(self.file, self.reportInterval)
