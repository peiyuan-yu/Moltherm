from os.path import join, dirname, abspath
import unittest
import operator
import json

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

module_dir = join(dirname(abspath(__file__)))
files_dir = join(module_dir, "..", "..", "..", "test_files")


class TestReactionRanker(unittest.TestCase):

    def setUp(self):
        with open(join(files_dir, "rxn_set.json"), 'rt') as file:
            self.reactions = json.load(file)
        with open(join(files_dir, "mol_set.json"), 'rt') as file:
            self.molecules = json.load(file)

    def tearDown(self):
        del self.molecules
        del self.reactions

    def test_defaults(self):
        ranker = ReactionRanker(self.reactions, self.molecules)
        self.assertEqual(ranker.parameters_min, ["mp", "vp", "solubility"])
        self.assertEqual(ranker.parameters_max, ["bp", "log_kow"])
        self.assertEqual(ranker.parameters, ranker.parameters_min + ranker.parameters_max)
        self.assertTrue(all([x in ranker.goals.keys() for x in ranker.parameters]))
        for param in ranker.parameters:
            self.assertEqual(ranker.goals[param], 0)
        self.assertEqual(ranker.molecules_raw, self.molecules)
        self.assertEqual(ranker.reactions, self.reactions)
        self.assertEqual(len(ranker.molecules.index), 10)
        self.assertEqual(len(ranker.all_reactants.index), 5)

        # Make sure that minimizing and maximizing parameters are treated differently
        mol_set = [e for e in self.molecules if e["mol_id"] == "9932269"][0]
        mol_df = ranker.molecules.loc["9932269"]
        self.assertEqual(mol_set["bp"], mol_df["bp"])
        self.assertEqual(mol_set["log_kow"], mol_df["log_kow"])
        self.assertEqual(mol_set["vp"], -1 * mol_df["vp"])
        self.assertEqual(mol_set["mp"], -1 * mol_df["mp"])
        self.assertEqual(mol_set["solubility"], -1 * mol_df["solubility"])

        rcts = ["108550", "9932269"]
        rxn_all_reactants = ranker.all_reactants.loc["9765072"]
        rxn_molecules = ranker.molecules.loc[ranker.molecules.index.isin(rcts)].sum()
        for param in ranker.parameters:
            self.assertEqual(rxn_all_reactants[param], rxn_molecules[param])

    def test_exclusion(self):
        # Test to make sure that no reactions are added to the test set if their
        # reactants are not present

        molecules = [m for m in self.molecules if m["mol_id"] != "9932269"]
        ranker = ReactionRanker(self.reactions, molecules)
        self.assertTrue("9765072" not in [r["rxn_id"] for r in ranker.reactions])

    def test_pareto_ranking(self):

        goals = {"bp": 257, "mp": 12, "log_kow": 5}
        ranker = ReactionRanker(self.reactions, self.molecules, goals=goals)

        # For this simple case, all rankings are Pareto optimal, so there is
        # only one class
        ranking = ranker.pareto_ranking(group_by_class=True)
        self.assertEqual(len(ranking[0]), 5)
        # At most, one reaction could be ranked differently in this test case
        ranking = ranker.pareto_ranking(by_worst=True)
        self.assertTrue(len(ranking[0]) >= 4)

        # Check that multiple groups form when only some parameters are used
        ranking = ranker.pareto_ranking(parameters=["vp", "log_kow"],
                                        group_by_class=True)
        self.assertEqual(len(ranking.keys()), 4)
        self.assertEqual(ranking[0], ['3942060'])
        self.assertEqual(ranking[1], ['1918886'])
        self.assertEqual(ranking[2], ['28714022', '2298180'])
        self.assertEqual(ranking[3], ['9765072'])

        # Check that list properly forms
        ranking = ranker.pareto_ranking(parameters=["log_kow"])
        self.assertEqual(ranking, ['3942060', '1918886', '2298180', '28714022',
                                   '9765072'])
        self.assertSequenceEqual(ranking,
                                 ranker.all_reactants.sort_values(by="log_kow").index.tolist())

        # Check that ranking works properly when reaction set is limited
        ranking = ranker.pareto_ranking(reactions=["9765072", "1918886",
                                                   "2298180"])
        self.assertTrue("3942060" not in ranking)

    def test_heuristic_ranking(self):
        pass

    def test_tiered_ranking(self):
        pass
