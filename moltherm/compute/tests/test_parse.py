import os
from os.path import join, dirname, abspath
import shutil
import unittest

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
        pass

    def tearDown(self):
        pass

    def test_parse_reaxys_xml(self):
        pass

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