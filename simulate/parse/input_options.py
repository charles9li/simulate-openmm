from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from simulate.parse.simulation import SimulationEnsembles
from simulate.parse.system import SystemOptions
from simulate.utils import LineDeque, NestedParser


class InputOptions(object):
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
        """ Create an InputOptions instance. """
        self.system_options = SystemOptions()
        self.simulation_ensembles = SimulationEnsembles()
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
            elif current_section == 'system':
                self.system_options.parse(parsed_lines.popleft())
            elif current_section == 'simulations':
                self.simulation_ensembles.parse(parsed_lines.popleft())
                pass
            else:
                raise ValueError("{} is not a valid section name.".format(current_section))
            current_section = None

    # =========================================================================

    # Private helper methods for parsing file

    def _remove_comment(self, line):
        """ Removes comment from line. """
        return line.split(self.COMMENT_CHAR)[0].strip()