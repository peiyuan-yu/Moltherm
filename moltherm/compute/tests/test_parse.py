import os
from os.path import join, dirname, abspath
import shutil
import unittest

import bs4

from pymatgen.io.babel import BabelMolAdaptor

from moltherm.compute.parse import ReaxysParser, EPISuiteParser


__author__ = "Evan Spotte-Smith"
__version__ = "0.2"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Beta"
__date__ = "June 2018"

module_dir = join(dirname(abspath(__file__)))
files_dir = join(module_dir, "..", "..", "..", "test_files")

class TestReaxysParser(unittest.TestCase):
    def setUp(self):
        self.filename = "reaxys.xml"
        self.parser = ReaxysParser(files_dir)

    def tearDown(self):
        del self.parser
        del self.filename

    def test_parse_reaxys_xml(self):

        parsed = self.parser.parse_reaxys_xml(self.filename)

        # Ensure that all proper reactions are present
        self.assertEqual(len(parsed), 3)

        indices = [x["meta"]["index"] for x in parsed]
        for index in [1, 9, 10]:
            self.assertTrue(index in indices)

        rxn_ids = [x["meta"]["rxn_id"] for x in parsed]
        for rxn_id in [28100547, 3553459, 8633298]:
            self.assertTrue(rxn_id in rxn_ids)

        # Ensure that all reactions are actually appropriate
        for rxn in parsed:
            pro_species = list()
            rct_species = list()

            for pro in rxn["pros"]:
                molecule = BabelMolAdaptor.from_string(pro, "mol")
                molecule.add_hydrogen()
                mol = molecule.pymatgen_mol
                pro_species += [str(site.specie) for site in mol]
            for rct in rxn["rcts"]:
                molecule = BabelMolAdaptor.from_string(rct, "mol")
                molecule.add_hydrogen()
                mol = molecule.pymatgen_mol
                rct_species += [str(site.specie) for site in mol]

            self.assertSequenceEqual(sorted(pro_species), sorted(rct_species))

        # Test one reaction to make sure metadata is recorded correctly
        rxn = parsed[-1]["meta"]
        self.assertEqual(rxn["solvents"], {'dichloromethane'})
        self.assertEqual(len(rxn["pro_meta"]), 1)
        self.assertEqual(len(rxn["rct_meta"]), 2)
        self.assertSequenceEqual(rxn["pro_meta"], [(5424566, '(4,5-dimethylcyclohexa-1,4-dienyl)trimethylsilane')])
        self.assertSequenceEqual(rxn["rct_meta"], [(605285, '2,3-dimethyl-buta-1,3-diene'),
                                                   (906752, 'trimethylsilylacetylene')])

    def test_get_unique_reactions(self):
        pass

    def test_store_reaxys_reactions_files(self):
        pass

    def test_store_reaxys_reactions_db(self):
        #TODO: How do I test this?
        pass

class TestEPISuiteParser(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_parse_epi_suite_summary(self):
        pass

    def test_parse_epi_suite_complete(self):
        pass

    def test_store_epi_suite_db(self):
        #TODO: How do I test this?
        pass