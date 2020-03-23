"""
input_options.py: Takes an input file and parses its contents to determine information about the simulation.

Copyright (c) 2020 Charles Li // UCSB, Department of Chemical Engineering

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
from __future__ import absolute_import

__author__ = "Charles Li"
__version__ = "1.0"

from ._simulations import SimulationsOptions
from ._system import SystemOptions
from simulate.utils import LineDeque, NestedParser


class InputOptions(object):
    """Class providing a parser for an input file containing OpenMM options."""

    # =========================================================================

    COMMENT_CHAR = '#'
    BEGIN_SECTION_CHAR = '{'
    END_SECTION_CHAR = '}'

    # =========================================================================

    def __init__(self, input_filename=None):
        """ Create an InputOptions instance.

        Parameters
        ----------
        input_filename : str
            The name of the file to be read
        """
        self.system_options = SystemOptions()
        self.simulation_ensembles = SimulationsOptions()
        self._read(input_filename)

    # =========================================================================

    # Read input file

    def _read(self, input_filename):
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
            else:
                raise ValueError("{} is not a valid section name.".format(current_section))
            current_section = None

    # =========================================================================

    # Private helper methods for parsing file

    def _remove_comment(self, line):
        """ Removes comment from line. """
        return line.split(self.COMMENT_CHAR)[0].strip()
