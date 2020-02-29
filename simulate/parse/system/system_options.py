from __future__ import absolute_import
__author__ = "Charles Li"
__version__ = "1.0"

from ast import literal_eval

from simtk.openmm import app, NonbondedForce
from simtk.unit import nanometer

from simulate.parse._options import _Options
from simulate.parse.system.topology_options import GromacsTopologyOptions, DodecaneAcrylateTopologyOptions


class SystemOptions(_Options):

    SECTION_NAME = 'system'

    # =========================================================================

    def __init__(self):
        super(SystemOptions, self).__init__()
        self.topology_options = None
        self.nonbondedMethod = app.NoCutoff
        self.nonbondedCutoff = 0.9*nanometer
        self.constraints = None
        self.rigidWater = True
        self.implicitSolvent = None
        self.soluteDielectric = 1.0
        self.solventDielectric = 78.5
        self.ewaldErrorTolerance = 0.0005
        self.removeCMMotion = True
        self.hydrogenMass = None
        self.useDispersionCorrection = True

    # =========================================================================

    def _check_for_incomplete_input(self):
        if self.topology_options is None:
            self._incomplete_error('topology')

    # =========================================================================

    NONBONDED_METHODS = {'NoCutoff': app.NoCutoff,
                         'CutoffPeriodic': app.CutoffPeriodic,
                         'CutoffNonPeriodic': app.CutoffNonPeriodic,
                         'Ewald': app.Ewald,
                         'PME': app.PME,
                         'LJPME': app.LJPME}
    TOPOLOGY_METHODS = {'GromacsTopology': GromacsTopologyOptions,
                        'DodecaneAcrylateTopology': DodecaneAcrylateTopologyOptions}

    def _parse_topology(self, *args):
        option_value = args[0]
        line_deque = args[1]
        topology_options = self.TOPOLOGY_METHODS[option_value]()
        topology_options.parse(line_deque.popleft())
        self.topology_options = topology_options

    def _parse_nonbonded_method(self, *args):
        option_value = args[0]
        try:
            self.nonbondedMethod = self.NONBONDED_METHODS[option_value]
        except KeyError:
            raise ValueError("{} is not a valid option for nonbondedMethod.".format(option_value))

    def _parse_nonbonded_cutoff(self, *args):
        self.nonbondedCutoff = literal_eval(args[0])*nanometer

    def _parse_ewald_error_tolerance(self, *args):
        self.ewaldErrorTolerance = literal_eval(args[0])

    def _parse_use_dispersion_correction(self, *args):
        self.useDispersionCorrection = literal_eval(args[0])

    OPTIONS = {'topology': _parse_topology,
               'nonbondedMethod': _parse_nonbonded_method,
               'nonbondedCutoff': _parse_nonbonded_cutoff,
               'ewaldErrorTolerance': _parse_ewald_error_tolerance,
               'useDispersionCorrection': _parse_use_dispersion_correction}

    # =========================================================================

    def topology(self):
        return self.topology_options.topology()

    def create_system(self):
        system = self.topology_options.create_system(nonbondedMethod=self.nonbondedMethod,
                                                     nonbondedCutoff=self.nonbondedCutoff,
                                                     constraints=self.constraints,
                                                     rigidWater=self.rigidWater,
                                                     implicitSolvent=self.implicitSolvent,
                                                     soluteDielectric=self.soluteDielectric,
                                                     solventDielectric=self.solventDielectric,
                                                     ewaldErrorTolerance=self.ewaldErrorTolerance,
                                                     removeCMMotion=self.removeCMMotion,
                                                     hydrogenMass=self.hydrogenMass)
        for force in system.getForces():
            if isinstance(force, NonbondedForce):
                force.setUseDispersionCorrection(self.useDispersionCorrection)
        return system
