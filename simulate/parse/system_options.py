from ast import literal_eval

from simtk.openmm.app import NoCutoff, CutoffPeriodic, CutoffNonPeriodic, Ewald, PME, LJPME
from simtk.unit import kelvin, nanometer

from ._options import _Options


class SystemOptions(_Options):

    SECTION_NAME = 'system'

    NONBONDED_METHODS = {'NoCutoff': NoCutoff,
                         'CutoffPeriodic': CutoffPeriodic,
                         'CutoffNonPeriodic': CutoffNonPeriodic,
                         'Ewald': Ewald,
                         'PME': PME,
                         'LJPME': LJPME}

    # =========================================================================

    def __init__(self):
        super(SystemOptions, self).__init__()
        self.nonbondedMethod = NoCutoff
        self.nonbondedCutoff = 0.9*nanometer
        self.ewaldErrorTolerance = 0.0005
        self.dispersionCorrection = True

    # =========================================================================

    def _parse_nonbonded_method(self, *args):
        option_value = args[0]
        try:
            self.nonbondedMethod = self.NONBONDED_METHODS[option_value]
        except KeyError:
            raise ValueError("{} is not a valid option for nonbondedMethod.".format(option_value))

    def _parse_nonbonded_cutoff(self, *args):
        self.nonbondedCutoff = literal_eval(args[0])

    def _parse_ewald_error_tolerance(self, *args):
        self.ewaldErrorTolerance = literal_eval(args[0])

    def _parse_dispersion_correction(self, *args):
        self.dispersionCorrection = literal_eval(args[0])

    # =========================================================================

    OPTIONS = {'nonbondedMethod': _parse_nonbonded_method,
               'nonbondedCutoff': _parse_nonbonded_cutoff,
               'ewaldErrorTolerance': _parse_ewald_error_tolerance,
               'dispersionCorrection': _parse_dispersion_correction}
