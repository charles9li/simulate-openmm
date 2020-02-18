from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from ast import literal_eval

import mdtraj as md

from simulate.parse._options import _Options


class _PositionOptions(_Options):

    def __init__(self):
        super(_PositionOptions, self).__init__()

    def set_positions(self, simulation):
        pass


class FileOptions(_PositionOptions):

    def __init__(self):
        super(FileOptions, self).__init__()
        self.file = None
        self.frame = 0

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.file is None:
            self._incomplete_error('file')

    # =========================================================================

    def _parse_file(self, *args):
        self.file = args[0]

    def _parse_frame(self, *args):
        self.frame = literal_eval(args[0])

    OPTIONS = {'file': _parse_file,
               'frame': _parse_frame}

    def set_positions(self, simulation):
        t = md.load(self.file)
        simulation.context.setPositions(t.xyz[self.frame])

