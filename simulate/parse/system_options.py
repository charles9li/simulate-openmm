from ast import literal_eval

from simtk.openmm.app import NoCutoff, CutoffPeriodic, CutoffNonPeriodic, Ewald, PME, LJPME
from simtk.unit import kelvin, nanometer

from ._options import _Options


class SystemOptions(_Options):

    NONBONDED_METHODS = {'NoCutoff': NoCutoff,
                         'CutoffPeriodic': CutoffPeriodic,
                         'CutoffNonPeriodic': CutoffNonPeriodic,
                         'Ewald': Ewald,
                         'PME': PME,
                         'LJPME': LJPME}

    def __init__(self):
        super(SystemOptions, self).__init__()
        self.nonbondedMethod = NoCutoff
        self.nonbondedCutoff = 0.9*nanometer
        self.ewaldErrorTolerance = 0.0005
        self.dispersionCorrection = True
        self.temperature = 298.0*kelvin

    def parse(self, line_deque):
        self._check_is_line_deque(line_deque)
        while len(line_deque) > 0:
            line = line_deque.popleft()
            option_name = self._parse_option_name(line)
            if option_name == 'nonbondedMethod':
                option_value = self._parse_option_value(line, option_name)
                self._parse_nonbonded_method(option_value)
            elif option_name == 'nonbondedCutoff':
                self.nonbondedCutoff = literal_eval(self._parse_option_value(line, option_name))*nanometer
            elif option_name == 'ewaldErrorTolerance':
                self.ewaldErrorTolerance = literal_eval(self._parse_option_value(line, option_name))
            elif option_name == 'dispersionCorrection':
                self.dispersionCorrection = literal_eval(self._parse_option_value(line, option_name))
            elif option_name == 'temperature':
                self.temperature = literal_eval(self._parse_option_value(line, option_name))
            else:
                raise ValueError("{} is not a valid option for the system section.".format(option_name))

    def _parse_nonbonded_method(self, option_value):
        try:
            self.nonbondedMethod = self.NONBONDED_METHODS[option_value]
        except KeyError:
            raise ValueError("{} is not a valid option for nonbondedMethod.".format(option_value))
