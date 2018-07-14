# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals

"""
Created on July 13, 2018
"""


__author__ = "Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__date__ = "July 13, 2018"

import unittest
import os
import copy
import warnings
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.io.babel import BabelMolAdaptor
from moltherm.react.func_groups import FunctionalGroupExtractor

test_dir = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                        "test_files")

try:
    import openbabel as ob
    import pybel as pb
    import networkx as nx
except ImportError:
    pb = None
    ob = None
    nx = None


@unittest.skipIf(not (pb and ob and nx), "OpenBabel or NetworkX not present. Skipping...")
class FunctionalGroupExtractorTest(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter("ignore")

        self.file = os.path.join(test_dir, "func_group_test.mol")
        self.mol = Molecule.from_file(self.file)
        self.strat = OpenBabelNN()
        self.mg = MoleculeGraph.with_local_env_strategy(self.mol, self.strat,
                                                        reorder=False,
                                                        extend_structure=False)
        self.extractor = FunctionalGroupExtractor(self.mg)

    def tearDown(self):
        warnings.resetwarnings()
        del self.extractor
        del self.mg
        del self.strat
        del self.mol
        del self.file

    def test_init(self):
        # Ensure that instantiation is equivalent for all valid input types
        extractor_str = FunctionalGroupExtractor(self.file)
        extractor_mol = FunctionalGroupExtractor(self.mol)
        extractor_mg = self.extractor

        self.assertEqual(extractor_str.molgraph, extractor_mol.molgraph)
        self.assertEqual(extractor_str.molgraph, extractor_mg.molgraph)
        self.assertEqual(extractor_str.species, extractor_mol.species)
        self.assertEqual(extractor_str.species, extractor_mg.species)

        # Test optimization
        file_no_h = os.path.join(test_dir, "func_group_test_no_h.mol")
        mol_no_h = Molecule.from_file(file_no_h)
        extractor_no_h = FunctionalGroupExtractor(mol_no_h, optimize=True)

        self.assertEqual(len(extractor_no_h.molecule), len(extractor_mol.molecule))
        # Not sure this one will work
        self.assertEqual(extractor_no_h.species, extractor_mol.species)

    def test_get_heteroatoms(self):
        heteroatoms = self.extractor.get_heteroatoms()
        hetero_species = [self.extractor.species[x] for x in heteroatoms]

        self.assertEqual(len(heteroatoms), 3)
        self.assertEqual(sorted(hetero_species), ["N", "O", "O"])

        # Test with limitation
        hetero_no_o = self.extractor.get_heteroatoms(elements=["N"])
        self.assertEqual(len(hetero_no_o), 1)

    def test_get_special_carbon(self):
        special_cs = self.extractor.get_special_carbon()

        self.assertEqual(len(special_cs), 4)

        # Test with limitation
        special_cs_no_o = self.extractor.get_special_carbon(elements=["N"])
        self.assertEqual(len(special_cs_no_o), 2)

    def test_link_marked_atoms(self):
        heteroatoms = self.extractor.get_heteroatoms()
        special_cs = self.extractor.get_special_carbon()

        link = self.extractor.link_marked_atoms(heteroatoms.union(special_cs))

        self.assertEqual(len(link), 1)
        self.assertEqual(len(link[0]), 9)

        # Exclude Oxygen-related functional groups
        heteroatoms_no_o = self.extractor.get_heteroatoms(elements=["N"])
        special_cs_no_o = self.extractor.get_special_carbon(elements=["N"])
        all_marked = heteroatoms_no_o.union(special_cs_no_o)

        link_no_o = self.extractor.link_marked_atoms(all_marked)

        self.assertEqual(len(link_no_o), 2)

    def test_get_basic_functional_groups(self):
        pass

    def test_get_all_functional_groups(self):
        pass

    def test_categorize_functional_groups(self):
        pass

if __name__ == "__main__":
    unittest.main()
