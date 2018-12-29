from os.path import join, dirname, abspath
import unittest
import operator

import numpy as np

from pymatgen.core.structure import Molecule

from moltherm.compute.utils import (get_molecule, find_common_solvents,
                                    get_reactions_common_solvent, extract_id,
                                    get_smiles)

__author__ = "Evan Spotte-Smith"
__version__ = "0.2"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Beta"
__date__ = "December 2018"

module_dir = join(dirname(abspath(__file__)))
files_dir = join(module_dir, "..", "..", "..", "test_files")


class TestUtils(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_get_molecule(self):
        benzene_file = join(files_dir, "benzene.mol")

        benzene_pmg = Molecule.from_file(benzene_file)
        benzene_moltherm = get_molecule(benzene_file)

        species_no_h = [e for e in benzene_moltherm.species
                        if str(e).upper() != "H"]
        self.assertEqual(species_no_h, benzene_pmg.species)

        # Test that implicit hydrogens are added appropriately
        species = sorted([str(e) for e in benzene_moltherm.species])
        self.assertEqual(species, ["C", "C", "C", "C",
                                                            "C", "C", "H", "H",
                                                            "H", "H", "H", "H"])

        # Ensure that get_molecule is deterministic, always produces the same
        # molecule
        self.assertEqual(benzene_moltherm, get_molecule(benzene_file))

        coords = np.array([[-7.59858151e-01,  1.16908119e+00, -1.61105859e-03],
                           [-1.39065495e+00, -7.49583582e-02, -9.63095317e-04],
                           [-6.28683825e-01, -1.24326334e+00,  6.60526465e-04],
                           [ 7.64084196e-01, -1.16752892e+00,  1.63614012e-03],
                           [ 1.39488100e+00,  7.65106237e-02,  9.87135623e-04],
                           [ 6.32909871e-01,  1.24481561e+00, -6.36441101e-04],
                           [-1.35352141e+00,  2.07932532e+00, -2.87509442e-03],
                           [-2.47578162e+00, -1.33964201e-01, -1.72330284e-03],
                           [-1.12014718e+00, -2.21251339e+00,  1.16530208e-03],
                           [ 1.35774746e+00, -2.07777305e+00,  2.90204589e-03],
                           [ 2.48000766e+00,  1.35516465e-01,  1.74638272e-03],
                           [ 1.12437322e+00,  2.21406566e+00, -1.14215271e-03]])

        self.assertTrue(np.allclose(benzene_moltherm.cart_coords, coords))

    def test_find_common_solvents(self):
        common_solvents = find_common_solvents(join(files_dir, "reaction_xml"))

        result = [('N,N-dimethyl-formamide', 1), ('', 1), ('chloroform', 1),
                  ('tetrahydrofuran', 1), ('1,2-dichloro-ethane', 1),
                  ('acetonitrile', 1), ('chlorobenzene', 1), ('benzene', 1),
                  ('toluene', 2)]

        # Ensure that result is expected
        self.assertSequenceEqual(common_solvents, result)

        # Ensure that sorting is correct
        self.assertEqual(common_solvents, sorted(common_solvents,
                                                 key=operator.itemgetter(1)))

    def test_get_reactions_common_solvent(self):
        toluene = get_reactions_common_solvent(join(files_dir, "reaction_xml"),
                                               ["toluene"])

        rxn_ids = ["30460605", "46313497"]
        self.assertEqual(sorted([e["rxn_id"] for e in toluene]), rxn_ids)

        # Attempt to catch multiple solvents
        tol_benzene = get_reactions_common_solvent(join(files_dir,
                                                        "reaction_xml"),
                                                   ["toluene", "benzene"])

        rxn_ids.append("4833625")
        self.assertEqual(sorted([e["rxn_id"] for e in tol_benzene]), rxn_ids)

        rct_ids = ['108550', '1747205', '1919880', '21267099', '605285',
                   '7876230']
        pro_ids = ['21267114', '31776431', '7886204']

        res_rct_ids = []
        res_pro_ids = []
        for entry in tol_benzene:
            for pro in entry["pro_ids"]:
                res_pro_ids.append(pro)
            for rct in entry["rct_ids"]:
                res_rct_ids.append(rct)

        # Check that parsing is accurate for other ids as well as rxn_ids
        # Could also check names, but these are far less important
        self.assertEqual(sorted(res_rct_ids), rct_ids)
        self.assertEqual(sorted(res_pro_ids), pro_ids)

    def test_extract_id(self):
        no_filename = "/home/moltherm/"
        no_path = "1234.mol"
        no_mol = "44511"
        filename_path = "/data/moltherm/112358.mol"
        underscores = "/global/homes/m/moltherm/1_90210.mol"

        self.assertEqual(extract_id(no_filename), "")
        self.assertEqual(extract_id(no_path), "1234")
        self.assertEqual(extract_id(no_mol), "44511")
        self.assertEqual(extract_id(filename_path), "112358")
        self.assertEqual(extract_id(underscores), "90210")

    def test_get_smiles(self):
        single_molecule = get_smiles(join(files_dir, "molecules"), ["1453094"])
        self.assertEqual(single_molecule[0], 'c12c(cc(c(=O)n1C)C=O)cccc2')

        all_mols = ["1453094", "1738108", "1873402", "2045554", "21925165",
                    "22125071", "28599994", "31695576", "5078635", "6657763"]

        for mol in all_mols:
            smiles = get_smiles(join(files_dir, "molecules"), [mol])[0]
            file = join(files_dir, "molecules", mol, "{}.mol".format(mol))

            mol_smiles = Molecule.from_str(smiles, "smi")
            mol_file = Molecule.from_file(file)

            smiles_species = sorted([str(e) for e in mol_smiles.species
                                     if str(e) != "H"])
            file_species = sorted([str(e) for e in mol_file.species
                                   if str(e) != "H"])

            self.assertSequenceEqual(smiles_species, file_species)