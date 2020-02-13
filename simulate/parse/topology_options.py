from ._options import _Options


class TopologyOptions(_Options):

    def __init__(self):
        super(_Options, self).__init__()
        self.topFilename = None
        self.groFilename = None
        self.coordFilename = None

    def parse(self, line_deque):
        self._check_is_line_deque(line_deque)
        while len(line_deque) > 0:
            line = line_deque.popleft()
            option_name = self._parse_option_name(line)
            if option_name == 'topFilename':
                self.topFilename = self._parse_option_value(line, option_name)
            elif option_name == 'groFilename':
                self.groFilename = self._parse_option_value(line, option_name)
            elif option_name == 'coordFilename':
                self.coordFilename = self._parse_option_value(line, option_name)
            else:
                raise ValueError("{} is not a valid option for the topology section.".format(option_name))
