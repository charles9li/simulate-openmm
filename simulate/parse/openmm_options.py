from .system_options import SystemOptions
from .topology_options import TopologyOptions
from .simulation_options import SimulationOptions
from simulate.utils import LineDeque, NestedParser


class OpenMMOptions(object):
    """ Class providing a parser for an input file containing OpenMM options.

    Parameters
    ----------
    input_filename : str
        The name of the file to be read

    """

    # =========================================================================

    COMMENT_CHAR = '#'
    BEGIN_SECTION_CHAR = '{'
    END_SECTION_CHAR = '}'

    # =========================================================================

    def __init__(self, input_filename=None):
        self.topology_options = TopologyOptions()
        self.system_options = SystemOptions()
        self.simulation_options = []
        self.read(input_filename)

    # =========================================================================

    def read(self, input_filename):
        """ Reads the input file and parses each section.

        Parameters
        ----------
        input_filename : str
            The name of the file to be read
        """
        if input_filename is None:
            raise ValueError("No input file was provided.")

        # remove comments
        line_deque = LineDeque()
        with open(input_filename) as f:
            for line in f:
                line = self._remove_comment(line)
                line_deque.append(line)

        # parse nested curly braces
        nested_parser = NestedParser(line_deque, self.BEGIN_SECTION_CHAR, self.END_SECTION_CHAR)
        parsed_lines = nested_parser.parsed_lines

        # parse options for each section
        current_section = None
        while len(parsed_lines) > 0:
            if current_section is None:
                current_section = parsed_lines.popleft().strip()
                continue
            elif current_section == 'topology':
                self.topology_options.parse(parsed_lines.popleft())
            elif current_section == 'system':
                self.system_options.parse(parsed_lines.popleft())
            elif current_section == 'simulation':
                simulation_options = SimulationOptions()
                simulation_options.parse(parsed_lines.popleft())
                self.simulation_options.append(simulation_options)
            else:
                raise ValueError("{} is not a valid section name.".format(current_section))
            current_section = None

    # =========================================================================

    # Private helper methods for parsing file

    def _remove_comment(self, line):
        """ Removes comment from line. """
        return line.split(self.COMMENT_CHAR)[0].strip()
