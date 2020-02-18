from simulate.parse._options import _Options


class TopologyOptions(_Options):

    SECTION_NAME = 'topology'

    # =========================================================================

    def __init__(self):
        super(_Options, self).__init__()
        self.topFilename = None
        self.groFilename = None
        self.coordFilename = None

    # =========================================================================

    def _parse_top_filename(self, *args):
        self.topFilename = args[0]

    def _parse_gro_filename(self, *args):
        self.groFilename = args[0]

    def _parse_coord_filename(self, *args):
        self.coordFilename = args[0]

    # =========================================================================

    OPTIONS = {'topFilename': _parse_top_filename,
               'groFilename': _parse_gro_filename,
               'coordFilename': _parse_coord_filename}
