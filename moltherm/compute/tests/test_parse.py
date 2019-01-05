import os
from os.path import join, dirname, abspath
import shutil
import unittest

from bs4 import BeautifulSoup

from pymatgen.io.babel import BabelMolAdaptor

from monty.serialization import loadfn, #dumpfn

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

        if "tmp" in os.listdir(files_dir):
            shutil.rmtree(join(files_dir, "tmp"))

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

        files = [join("reaxys_xml", "diels_alder_1.xml"),
                 join("reaxys_xml", "diels_alder_2.xml")]

        res_1 = self.parser.parse_reaxys_xml(files[0])
        res_2 = self.parser.parse_reaxys_xml(files[1])

        self.assertEqual(len(res_1), 207)
        self.assertEqual(len(res_2), 231)

        unique = self.parser.get_unique_reactions(files)
        self.assertEqual(len(unique[0]), len(unique[1]))
        self.assertEqual(len(unique[0]), 328)

        for i in range(len(unique[0])):
            self.assertEqual(unique[0][i], unique[1][i]["meta"]["rxn_id"])

    def test_store_reaxys_reactions_files(self):

        parsed = self.parser.parse_reaxys_xml(self.filename)

        rxns = [28100547, 3553459, 8633298]
        mols = [15815515, 471171, 15815519, 636190, 969158,
                605285, 906752, 4992354, 5424566]

        # Test with base_path already created
        os.mkdir(join(files_dir, "tmp"))

        self.parser.store_reaxys_reactions_files(parsed,
                                                 base_path=join(files_dir,
                                                                "tmp"))

        files_rxn = os.listdir(join(files_dir, "tmp", "reactions"))
        for rxn_id in rxns:
            self.assertTrue("{}.xml".format(rxn_id) in files_rxn)

        files_mol = os.listdir(join(files_dir, "tmp", "molecules"))
        for mol_id in mols:
            self.assertTrue("{}".format(mol_id) in files_mol)
            files_this_mol = os.listdir(join(files_dir, "tmp", "molecules",
                                             str(mol_id)))
            self.assertTrue("{}.mol".format(mol_id) in files_this_mol)

        shutil.rmtree(join(files_dir, "tmp"))

        # Test with base_path not created

        self.parser.store_reaxys_reactions_files(parsed,
                                                 base_path=join(files_dir,
                                                                "tmp"))

        files_rxn = os.listdir(join(files_dir, "tmp", "reactions"))
        for rxn_id in rxns:
            self.assertTrue("{}.xml".format(rxn_id) in files_rxn)

        files_mol = os.listdir(join(files_dir, "tmp", "molecules"))
        for mol_id in mols:
            self.assertTrue("{}".format(mol_id) in files_mol)
            files_this_mol = os.listdir(join(files_dir, "tmp", "molecules",
                                             str(mol_id)))
            self.assertTrue("{}.mol".format(mol_id) in files_this_mol)

        for rxn in parsed:
            rxn_path = join(files_dir, "tmp", "reactions",
                            "{}.xml".format(rxn["meta"]["rxn_id"]))
            with open(rxn_path, 'r') as fileobj:
                xml = fileobj.read()
                xml_parsed = BeautifulSoup(xml, "lxml-xml")

                self.assertEqual(str(rxn["meta"]["index"]),
                                 xml_parsed.find("index").text)
                self.assertEqual(str(rxn["meta"]["rxn_id"]),
                                 xml_parsed.find("reaxysid").text)
                self.assertSequenceEqual(rxn["meta"]["solvents"],
                                         set(xml_parsed.find("solvents").text.split(" || ")))

                pro_ids = [x[0] for x in rxn["meta"]["pro_meta"]]
                pro_names = [x[1] for x in rxn["meta"]["pro_meta"]]
                rct_ids = [x[0] for x in rxn["meta"]["rct_meta"]]
                rct_names = [x[1] for x in rxn["meta"]["rct_meta"]]

                self.assertSequenceEqual([str(e) for e in pro_ids],
                                         [e.text for e in xml_parsed.find_all("proid")])
                self.assertSequenceEqual([e for e in pro_names],
                                         [e.text for e in xml_parsed.find_all("proname")])
                self.assertSequenceEqual([str(e) for e in rct_ids],
                                         [e.text for e in xml_parsed.find_all("rctid")])
                self.assertSequenceEqual([e for e in rct_names],
                                         [e.text for e in xml_parsed.find_all("rctname")])

                for i, e in enumerate(pro_ids):
                    with open(join(files_dir, "tmp", "molecules", str(e), "{}.mol".format(e)), 'r') as pro_file:
                        pro_entry = pro_file.read()
                        self.assertEqual(pro_entry, rxn["pros"][i])

                for i, e in enumerate(rct_ids):
                    with open(join(files_dir, "tmp", "molecules", str(e), "{}.mol".format(e)), 'r') as rct_file:
                        rct_entry = rct_file.read()
                        self.assertEqual(rct_entry, rxn["rcts"][i])

        shutil.rmtree(join(files_dir, "tmp"))

    def test_store_reaxys_reactions_db(self):
        #TODO: How do I test this?
        pass


class TestEPISuiteParser(unittest.TestCase):
    def setUp(self):

        self.parser = EPISuiteParser()

        self.summary_file = join(files_dir, "episuite_summary.out")
        self.complete_file = join(files_dir, "episuite_full.out")

    def tearDown(self):

        del self.complete_file
        del self.summary_file
        del self.parser

    def test_parse_epi_suite_summary(self):

        parsed = self.parser.parse_epi_suite_summary(self.summary_file)
        #dumpfn(parsed, join(files_dir, "episuite_summary.json"))
        reference = loadfn(join(files_dir, "episuite_summary.json"))
        self.assertSequenceEqual(parsed, reference)

    def test_parse_epi_suite_complete(self):

        parsed = self.parser.parse_epi_suite_complete(self.complete_file)
        #dumpfn(parsed, join(files_dir, "episuite_complete.json"))
        reference = loadfn(join(files_dir, "episuite_complete.json"))
        self.assertSequenceEqual(parsed, reference)

    def test_store_epi_suite_db(self):
        #TODO: How do I test this?
        pass