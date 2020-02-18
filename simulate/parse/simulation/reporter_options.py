from ast import literal_eval

from simulate.parse._options import _Options


class _ReporterOptions(_Options):

    SECTION_NAME = "_ReporterOptions"

    def __init__(self):
        super(_ReporterOptions, self).__init__()
        self.file = None
        self.reportInterval = None

    # =========================================================================

    def _parse_file(self, *args):
        self.file = args[0]

    def _parse_report_interval(self, *args):
        self.reportInterval = literal_eval(args[0])

    # =========================================================================

    def reporter(self):
        pass

    # =========================================================================

    def _no_option_specified_exception(self, option_name):
        raise ValueError("No {} specified for {}".format(option_name, self.SECTION_NAME))


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
        if self.file is None:
            self._no_option_specified_exception('file')
        if self.reportInterval is None:
            self._no_option_specified_exception('report interval')
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
        if self.file is None:
            self._no_option_specified_exception('file')
        if self.reportInterval is None:
            self._no_option_specified_exception('report interval')
        from simtk.openmm.app import StateDataReporter
        return StateDataReporter(self.file, self.reportInterval, step=True, time=True,
                                 potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
                                 temperature=True, volume=True, density=True,
                                 speed=True, elapsedTime=True)
