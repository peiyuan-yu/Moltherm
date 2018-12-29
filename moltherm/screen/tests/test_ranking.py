from os.path import join, dirname, abspath
import unittest
import operator

import numpy as np

from moltherm.screen.ranking import ReactionRanker

module_dir = join(dirname(abspath(__file__)))
files_dir = join(module_dir, "..", "..", "..", "test_files")

__author__ = "Evan Spotte-Smith"
__version__ = "0.2"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Beta"
__date__ = "December 2018"


