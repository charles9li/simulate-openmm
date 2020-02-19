from ast import literal_eval
from copy import copy

from simulate.parse._options import _Options


class _ReporterOptions(_Options):

    SECTION_NAME = "_ReporterOptions"

    def __init__(self):
        super(_ReporterOptions, self).__init__()
        self.file = None
        self.reportInterval = None

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

    OPTIONS = {'file': _parse_file,
               'reportInterval': _parse_report_interval}

    # =========================================================================

    def reporter(self):
        pass


class DCDReporterOptions(_ReporterOptions):

    SECTION_NAME = "DCDReporter"

    def __init__(self):
        super(DCDReporterOptions, self).__init__()
        self.file = None
        self.reportInterval = None
        self.append = False
        self.enforcePeriodicBox = None

    # =========================================================================

    def _parse_file(self, *args):
        super(DCDReporterOptions, self)._parse_file(*args)

    def _parse_report_interval(self, *args):
        super(DCDReporterOptions, self)._parse_report_interval(*args)

    def _parse_append(self, *args):
        self.append = literal_eval(args[0])

    def _parse_enforce_periodic_box(self, *args):
        self.enforcePeriodicBox = args[0]

    OPTIONS = {'file': _parse_file,
               'reportInterval': _parse_report_interval,
               'append': _parse_append,
               'enforcePeriodicBox': _parse_enforce_periodic_box}

    # =========================================================================

    def reporter(self):
        from simtk.openmm.app import DCDReporter
        return DCDReporter(self.file, self.reportInterval,
                           append=self.append, enforcePeriodicBox=self.enforcePeriodicBox)


class StateDataReporterOptions(_ReporterOptions):

    SECTION_NAME = "StateDataReporter"

    # =========================================================================

    def __init__(self):
        super(StateDataReporterOptions, self).__init__()
        self.file = None
        self.reportInterval = None

    # =========================================================================

    def _parse_file(self, *args):
        super(StateDataReporterOptions, self)._parse_file(*args)

    def _parse_report_interval(self, *args):
        super(StateDataReporterOptions, self)._parse_report_interval(*args)

    OPTIONS = {'file': _parse_file,
               'reportInterval': _parse_report_interval}

    # =========================================================================

    def reporter(self):
        from simtk.openmm.app import StateDataReporter
        return StateDataReporter(self.file, self.reportInterval, step=True, time=True,
                                 potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                 temperature=True, volume=True, density=True,
                                 speed=True, elapsedTime=True)


class RNEMDReporterOptions(_ReporterOptions):

    SECTION_NAME = "RNEMDReporter"

    # =========================================================================

    def __init__(self):
        super(RNEMDReporterOptions, self).__init__()
        self.file = None
        self.reportInterval = None

    # =========================================================================

    def _parse_file(self, *args):
        super(RNEMDReporterOptions, self)._parse_file(*args)

    def _parse_report_interval(self, *args):
        super(RNEMDReporterOptions, self)._parse_report_interval(*args)

    OPTIONS = {'file': _parse_file,
               'reportInterval': _parse_report_interval}

    # =========================================================================

    def reporter(self):
        from simulate.reporters import RNEMDReporter
        return RNEMDReporter(self.file, self.reportInterval, step=True)


class RNEMDVelocityReporterOptions(_ReporterOptions):

    SECTION_NAME = "RNEMDVelocityReporter"

    # =========================================================================

    def __init__(self):
        super(RNEMDVelocityReporterOptions, self).__init__()
        self.file = None
        self.reportInterval = None
        self.numSlabs = None

    # =========================================================================
    
    def _check_for_incomplete_input(self):
        super(RNEMDVelocityReporterOptions, self)._check_for_incomplete_input()
        if self.numSlabs is None:
            self._incomplete_error('numSlabs')

    # =========================================================================

    def _parse_num_slabs(self, *args):
        self.numSlabs = literal_eval(args[0])

    OPTIONS = copy(_ReporterOptions.OPTIONS)
    OPTIONS['numSlabs'] = _parse_num_slabs

    # =========================================================================

    def reporter(self):
        from simulate.reporters import RNEMDVelocityReporter
        return RNEMDVelocityReporter(self.file, self.reportInterval, self.numSlabs, step=True)


class CheckpointReporterOptions(_ReporterOptions):

    SECTION_NAME = "CheckPointReporter"

    # =========================================================================

    def __init__(self):
        super(CheckpointReporterOptions, self).__init__()

    # =========================================================================

    def reporter(self):
        from simtk.openmm.app import CheckpointReporter
        return CheckpointReporter(self.file, self.reportInterval)
