from simulate.utils import LineDeque


class _Options(object):

    def __init__(self):
        pass

    def parse(self, line_deque):
        pass

    @staticmethod
    def _parse_option_name(line):
        return line.split('=')[0].strip()

    @staticmethod
    def _parse_option_value(line, option_name):
        try:
            option_value = line[1].strip()
        except IndexError:
            option_value = ''
        if not option_value:
            raise ValueError("No value specified for {} option.".format(option_name))
        return option_value

    @staticmethod
    def _check_is_line_deque(line_deque):
        if not isinstance(line_deque, LineDeque):
            raise ValueError("Curly braces must be provided, even if no options are specified.")
