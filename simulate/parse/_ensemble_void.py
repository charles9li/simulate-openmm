"""
_ensemble_void.py: Parses options for creating voids.

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

from ast import literal_eval
import os
import random

from simtk.openmm.app import PDBFile, Simulation, Topology
from simtk.openmm import MonteCarloBarostat, MonteCarloAnisotropicBarostat, MonteCarloMembraneBarostat
import numpy as np

from simulate.parse._options import _Options

__all__ = ['VoidOptions']


class VoidOptions(_Options):

    def __init__(self, ensemble_options):
        super(VoidOptions, self).__init__()
        self.ensemble_options = ensemble_options
        self.numVoid = 0
        self.file = None

    def _create_options(self):
        super(VoidOptions, self)._create_options()
        self._OPTIONS['numVoid'] = self._parse_num_void
        self._OPTIONS['file'] = self._parse_file

    # =========================================================================

    def _parse_num_void(self, *args):
        self.numVoid = literal_eval(args[0])

    def _parse_file(self, *args):
        self.file = self._create_filepath(args[0])

    # =========================================================================

    def _create_filepath(self, filepath):
        directory = self.ensemble_options.simulations_options.input_options.directory
        return os.path.join(directory, filepath)

    # =========================================================================

    def remove_molecules(self, simulation_old, system_options, ensemble_options):
        topology_old = simulation_old.topology
        system_old = simulation_old.system
        state = simulation_old.context.getState(getPositions=True, getVelocities=True)
        positions_old = state.getPositions(asNumpy=True)
        velocities_old = state.getVelocities(asNumpy=True)
        periodic_box_vectors = state.getPeriodicBoxVectors(asNumpy=True)

        # randomly determine which molecules are removed
        removed_molecules = random.sample(range(topology_old.getNumChains()), self.numVoid)

        # create new topology and dictionary mapping old atom indices to new atom indices
        removed_atoms = []
        topology_new = Topology()
        old_to_new = {}
        atom_index_new = 0
        for chain_index, chain_old in enumerate(topology_old.chains()):
            if chain_index not in removed_molecules:
                chain_id = chain_old.id.split('-', 1)[-1]
                chain_new = topology_new.addChain(id="{}-{}".format(topology_new.getNumChains()+1, chain_id))
                for residue_old in chain_old.residues():
                    residue_new = topology_new.addResidue(residue_old.name, chain_new)
                    for atom_old in residue_old.atoms():
                        topology_new.addAtom(atom_old.name, atom_old.element, residue_new)
                        old_to_new[atom_old.index] = atom_index_new
                        atom_index_new += 1
            else:
                for atom_old in chain_old.atoms():
                    removed_atoms.append(atom_old.index)

        # add bonds to new topology
        atoms_new = list(topology_new.atoms())
        for bond_old in topology_old.bonds():
            atom1_index_old = bond_old[0].index
            atom2_index_old = bond_old[1].index
            try:
                atom1_new = atoms_new[old_to_new[atom1_index_old]]
                atom2_new = atoms_new[old_to_new[atom2_index_old]]
                topology_new.addBond(atom1_new, atom2_new)
            except KeyError:
                pass

        # set box vectors for topology
        topology_new.setPeriodicBoxVectors(periodic_box_vectors)

        # create system
        system_new = system_options.create_system_with_new_topology(topology_new)

        # check if old system had barostat and add to new system if applicable
        if self._has_barostat(system_old):
            barostat_old = self._get_barostat(system_old)
            defaultPressure = barostat_old.getDefaultPressure()
            defaultTemperature = barostat_old.getDefaultTemperature()
            frequency = barostat_old.getFrequency()
            if isinstance(barostat_old, MonteCarloAnisotropicBarostat):
                scaleX = barostat_old.getScaleX()
                scaleY = barostat_old.getScaleY()
                scaleZ = barostat_old.getScaleZ()
                barostat_new = MonteCarloAnisotropicBarostat(defaultPressure, defaultTemperature,
                                                             scaleX, scaleY, scaleZ, frequency)
            else:
                barostat_new = MonteCarloBarostat(defaultPressure, defaultTemperature, frequency)
            system_new.addForce(barostat_new)
            barostat_new.setForceGroup(system_new.getNumForces() - 1)

        # create integrator
        integrator = ensemble_options.create_integrator()

        # create new positions and velocities arrays
        positions_new = np.delete(positions_old, removed_atoms, axis=0)
        velocities_new = np.delete(velocities_old, removed_atoms, axis=0)

        # create new simulation
        simulation_new = Simulation(topology_new, system_new, integrator)
        simulation_new.context.setPositions(positions_new)
        simulation_new.context.setVelocities(velocities_new)
        simulation_new.context.setPeriodicBoxVectors(*periodic_box_vectors)

        # move reporters from old simulation to new simulation
        while simulation_old.reporters:
            simulation_new.reporters.append(simulation_old.reporters.pop(0))

        # create pdb file if specified
        if self.file is not None:
            PDBFile.writeFile(topology_new, positions_new, open(self.file, 'w'))

        return simulation_new

    def _get_barostat(self, system):
        forces = system.getForces()
        for force in forces:
            if self._is_barostat(force):
                barostat_index = forces.index(force)
                return system.getForce(barostat_index)

    def _has_barostat(self, system):
        forces = system.getForces()
        for force in forces:
            if self._is_barostat(force):
                return True
        return False

    @staticmethod
    def _is_barostat(force):
        barostats = [MonteCarloBarostat, MonteCarloAnisotropicBarostat, MonteCarloMembraneBarostat]
        for barostat in barostats:
            if isinstance(force, barostat):
                return True
        return False
