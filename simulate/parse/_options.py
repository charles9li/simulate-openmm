from simulate.utils import LineDeque


class _Options(object):

    SECTION_NAME = '_options'

    # =========================================================================

    def __init__(self):
        pass

    # =========================================================================

    OPTIONS = {}

    def parse(self, line_deque):
        self._check_is_line_deque(line_deque)
        while len(line_deque) > 0:
            line = line_deque.popleft()
            option_name = self._parse_option_name(line)
            if option_name in self.OPTIONS.keys():
                option_value = self._parse_option_value(line, option_name)
                self.OPTIONS[option_name](option_value, line_deque)
            else:
                raise ValueError("{} is not a valid option for the {} section.".format(option_name, self.SECTION_NAME))

    @staticmethod
    def _parse_option_name(line):
        return line.split('=')[0].strip()

    @staticmethod
    def _parse_option_value(line, option_name):
        try:
            option_value = line.split('=')[1].strip()
        except IndexError:
            option_value = ''
        if not option_value:
            raise ValueError("No value specified for {} option.".format(option_name))
        return option_value

    @staticmethod
    def _check_is_line_deque(line_deque):
        if not isinstance(line_deque, LineDeque):
            raise ValueError("Curly braces must be provided, even if no options are specified.")
