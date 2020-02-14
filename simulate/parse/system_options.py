from ast import literal_eval

from simtk.openmm import app
from simtk.unit import kelvin, nanometer

from ._options import _Options


class SystemOptions(_Options):

    SECTION_NAME = 'system'

    NONBONDED_METHODS = {'NoCutoff': app.NoCutoff,
                         'CutoffPeriodic': app.CutoffPeriodic,
                         'CutoffNonPeriodic': app.CutoffNonPeriodic,
                         'Ewald': app.Ewald,
                         'PME': app.PME,
                         'LJPME': app.LJPME}

    # =========================================================================

    def __init__(self):
        super(SystemOptions, self).__init__()
        self.nonbondedMethod = app.NoCutoff
        self.nonbondedCutoff = 0.9*nanometer
        self.ewaldErrorTolerance = 0.0005
        self.useDispersionCorrection = True
        self.initialVelocityTemperature = None

    # =========================================================================

    def _parse_nonbonded_method(self, *args):
        option_value = args[0]
        try:
            self.nonbondedMethod = self.NONBONDED_METHODS[option_value]
        except KeyError:
            raise ValueError("{} is not a valid option for nonbondedMethod.".format(option_value))

    def _parse_nonbonded_cutoff(self, *args):
        self.nonbondedCutoff = literal_eval(args[0])*nanometer

    def _parse_ewald_error_tolerance(self, *args):
        self.ewaldErrorTolerance = literal_eval(args[0])

    def _parse_use_dispersion_correction(self, *args):
        self.useDispersionCorrection = literal_eval(args[0])

    def _parse_initial_velocity_temperature(self, *args):
        self.initialVelocityTemperature = literal_eval(args[0])*kelvin

    # =========================================================================

    OPTIONS = {'nonbondedMethod': _parse_nonbonded_method,
               'nonbondedCutoff': _parse_nonbonded_cutoff,
               'ewaldErrorTolerance': _parse_ewald_error_tolerance,
               'useDispersionCorrection': _parse_use_dispersion_correction,
               'initialVelocityTemperature': _parse_initial_velocity_temperature}
