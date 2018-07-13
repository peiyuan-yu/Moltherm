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

    def tearDown(self):
        warnings.resetwarnings()

    def test_init(self):
        pass

    def test_get_heteroatoms(self):
        pass

    def test_get_special_carbon(self):
        pass

    def test_link_marked_atoms(self):
        pass

    def test_get_basic_functional_groups(self):
        pass

    def test_get_all_functional_groups(self):
        pass

    def test_categorize_functional_groups(self):
        pass

if __name__ == "__main__":
    unittest.main()
