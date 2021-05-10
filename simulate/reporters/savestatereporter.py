from __future__ import absolute_import, print_function

import os


class SaveStateReporter(object):

    def __init__(self, file, reportInterval):
        self._reportInterval = reportInterval
        self._file = file

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False)

    def report(self, simulation, state):
        tempFilename1 = self._file + ".backup1"
        tempFilename2 = self._file + ".backup2"
        simulation.saveState(tempFilename1)
        exists = os.path.exists(self._file)
        if exists:
            os.rename(self._file, tempFilename2)
        os.rename(tempFilename1, self._file)
        if exists:
            os.rename(tempFilename2)
