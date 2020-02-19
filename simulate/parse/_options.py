""" _options.py: A base class for other classes that store user options for
each section in the input script.

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

from simulate.utils import LineDeque


class _Options(object):
    """ Base class for parsing options for an OpenMM simulation.

    This class is inherited by classes that contain options for each section.
    The __init__, _parse_options, and _parse_sections methods are overridden by
    classes that inherit this class."""

    # Name of the section. Used for error messages. Overridden by child classes
    SECTION_NAME = '_options'

    # =========================================================================

    def __init__(self):
        """Creates an _Options instance.

        Options are parsed by calling the parse method().
        """
        pass

    # =========================================================================

    # Contains option/section name, private _parse method pairs. Used in the
    # parse() method to determine values for parameters.
    OPTIONS = {}
    SECTIONS = {}

    def parse(self, line_deque):
        """ Parses the line_deque into options and subsections for this
        section.

        """
        self._check_is_line_deque(line_deque)
        while len(line_deque) > 0:
            line = line_deque.popleft()
            name = self._parse_option_name(line)
            if self._is_option(line):
                if name in self.OPTIONS.keys():
                    option_value = self._parse_option_value(line, name)
                    self.OPTIONS[name](self, option_value, line_deque)
                else:
                    self._invalid_option_error(name)
            else:
                if name in self.SECTIONS.keys():
                    self.SECTIONS[name](self, name, line_deque)
                else:
                    self._invalid_section_error(name)
        self._check_for_incomplete_input()

    def _check_for_incomplete_input(self):
        """ Checks to see if any crucial information is missing after user
        input is parsed.

        """
        pass

    # =========================================================================

    # Private helper methods for parsing options

    @staticmethod
    def _check_is_line_deque(line_deque):
        """ Throws an error if argument is not a LineDeque instance. """
        if not isinstance(line_deque, LineDeque):
            raise ValueError("Curly braces must be provided, even if no options are specified.")

    @staticmethod
    def _parse_option_name(line):
        """ Parses an option for this section. """
        return line.split('=')[0].strip()

    @staticmethod
    def _parse_option_value(line, option_name):
        """ Parses a subsection for this section. """
        try:
            option_value = line.split('=')[1].strip()
        except IndexError:
            option_value = ''
        if not option_value:
            raise ValueError("No value specified for {} option.".format(option_name))
        return option_value

    @staticmethod
    def _is_option(line):
        """ Returns True if line contains a option, False if it is section
        name.

        """
        return '=' in line

    # =========================================================================

    # Private helper methods that throw errors

    def _invalid_option_error(self, option_name):
        """ Throws ValueError if the option is not valid for this section. """
        msg = "'{}' is not a valid option for the '{}' section.".format(option_name, self.SECTION_NAME)
        raise ValueError(msg)

    def _invalid_section_error(self, section_name):
        """ Throws ValueError if the section is not a valid subsection for
        this section.

        """
        msg = "'{}' is not a subsection for the '{}' section.".format(section_name, self.SECTION_NAME)
        raise ValueError(msg)

    def _incomplete_error(self, option_name):
        """ Throws ValueError if input is missing. """
        msg = "'{}' must be specified for the '{}' section.".format(option_name, self.SECTION_NAME)
        raise ValueError(msg)
