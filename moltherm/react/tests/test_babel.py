# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Created on Apr 28, 2012
"""


__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Apr 28, 2012"

import unittest
import os
import copy
import warnings
from pymatgen.core.structure import Molecule
from pymatgen.io.xyz import XYZ
# from pymatgen.io.babel import BabelMolAdaptor
from ..babel import BabelMolAdaptor

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files", "molecules")

try:
    import openbabel as ob
    import pybel as pb
except ImportError:
    pb = None
    ob = None


@unittest.skipIf(not (pb and ob), "OpenBabel not present. Skipping...")
class BabelMolAdaptorTest(unittest.TestCase):

    def setUp(self):
        coords = [[0.000000, 0.000000, 0.000000],
                  [0.000000, 0.000000, 1.089000],
                  [1.026719, 0.000000, -0.363000],
                  [-0.513360, -0.889165, -0.363000],
                  [-0.513360, 0.889165, -0.363000]]
        self.mol = Molecule(["C", "H", "H", "H", "H"], coords)
        warnings.simplefilter("ignore")

    def tearDown(self):
        warnings.resetwarnings()

    def test_init(self):
        adaptor = BabelMolAdaptor(self.mol)
        obmol = adaptor.openbabel_mol
        self.assertEqual(obmol.NumAtoms(), 5)

        adaptor = BabelMolAdaptor(adaptor.openbabel_mol)
        self.assertEqual(adaptor.pymatgen_mol.formula, "H4 C1")

    # def test_from_file(self):
    #     adaptor = BabelMolAdaptor.from_file(
    #         os.path.join(test_dir, "Ethane_e.pdb"), "pdb")
    #     mol = adaptor.pymatgen_mol
    #     self.assertEqual(mol.formula, "H6 C2")

    # def test_from_molecule_graph(self):
    #     pass

    def test_from_string(self):
        xyz = XYZ(self.mol)
        adaptor = BabelMolAdaptor.from_string(str(xyz), "xyz")
        mol = adaptor.pymatgen_mol
        self.assertEqual(mol.formula, "H4 C1")

    def test_localopt(self):
        mol = copy.deepcopy(self.mol)
        mol[1] = "H", [0, 0, 1.05]
        adaptor = BabelMolAdaptor(mol)
        adaptor.localopt()
        optmol = adaptor.pymatgen_mol
        for site in optmol[1:]:
            self.assertAlmostEqual(site.distance(optmol[0]), 1.09216, 2)

    def test_make3D(self):
        mol_0d = pb.readstring("smi", "CCCC").OBMol
        adaptor = BabelMolAdaptor(mol_0d)
        adaptor.make3D()
        self.assertEqual(mol_0d.GetDimension(), 3)
        # for atom in ob.OBMolAtomIter(mol_0d):
        #     print(atom.GetAtomicNum())
        #     print([atom.GetX(), atom.GetY(), atom.GetZ()])

    def add_hydrogen(self):
        mol_0d = pb.readstring("smi", "CCCC").OBMol
        self.assertEqual(len(pb.Molecule(mol_0d).atoms), 2)
        adaptor = BabelMolAdaptor(mol_0d)
        adaptor.add_hydrogen()
        self.assertEqual(len(adaptor.pymatgen_mol.sites), 14)

    def test_rotor_search_wrs(self):
        mol = copy.deepcopy(self.mol)
        mol[1] = "H", [0, 0, 1.05]
        adaptor = BabelMolAdaptor(mol)
        rotor_args = (250, 50)
        adaptor.rotor_conformer(*rotor_args, algo="WeightedRotorSearch")
        self.assertNotEquals(adaptor.pymatgen_mol[1].coords[2], 1.05)

    def test_rotor_search_srs(self):
        mol = copy.deepcopy(self.mol)
        mol[1] = "H", [0, 0, 1.05]
        adaptor = BabelMolAdaptor(mol)
        adaptor.rotor_conformer(200, algo="SystematicRotorSearch")
        self.assertNotEquals(adaptor.pymatgen_mol[1].coords[2], 1.05)

    def test_rotor_search_frs(self):
        mol = copy.deepcopy(self.mol)
        mol[1] = "H", [0, 0, 1.05]
        adaptor = BabelMolAdaptor(mol)
        adaptor.rotor_conformer(algo="FastRotorSearch")
        self.assertNotEquals(adaptor.pymatgen_mol[1].coords[2], 1.05)

    def test_rotor_search_rrs(self):
        mol = copy.deepcopy(self.mol)
        mol[1] = "H", [0, 0, 1.05]
        adaptor = BabelMolAdaptor(mol)
        adaptor.rotor_conformer(250, 50, algo="RandomRotorSearch")
        self.assertNotEquals(adaptor.pymatgen_mol[1].coords[2], 1.05)

    def test_gen3d_conformer(self):
        mol = copy.deepcopy(self.mol)
        mol[1] = "H", [0, 0, 1.05]
        adaptor = BabelMolAdaptor(mol)
        # adaptor.gen3d_conformer()
        # temporally auto-pass this test as sthe demand for memory is large
        pass

    def test_confab_conformers(self):
        mol = copy.deepcopy(self.mol)
        mol[1] = "H", [0, 0, 1.05]
        adaptor = BabelMolAdaptor(mol)
        conformers = adaptor.confab_conformers()
        self.assertEquals(BabelMolAdaptor(mol).openbabel_mol.NumRotors(), 0)
        self.assertAlmostEquals(BabelMolAdaptor(conformers[0]).pymatgen_mol[1].
                                coords[2], 1.05)

if __name__ == "__main__":
    unittest.main()
