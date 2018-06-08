import os
import numpy as np
import pandas as pd
from pymatgen.core.structure import Molecule
from pymatgen.analysis.molecule_matcher import MoleculeMatcher


# TODO: add bonding information to pymatgen molecule object? -qw
class ProductSplit:
    """
    Split product molecule to reactants by breaking bonds.
    Currently supports breaking one bond, breaking two bonds for circular
    molecules to be added.
    Args:

    """
    def __init__(self, mol, mol_name=None, write_output=True,
                 output_path="", output_type='sdf'):
        self.mol = mol
        self.mol_name = mol_name if mol_name else mol.formula
        self.write_output = write_output
        self.output_path = output_path
        self.output_type = output_type

    @classmethod
    def from_file(cls, input_path, input_file, mol_name=None, write_output=True,
                  output_path="", output_type='sdf'):
        f = os.path.join(input_path, input_file)
        print(f)
        mol = Molecule.from_file(f)
        mol_name = mol_name if mol_name else os.path.splitext(input_file)[0]
        return cls(mol, mol_name=mol_name, write_output=write_output,
                   output_path=output_path, output_type=output_type)

    def split(self, break_bond_sites):
        list_like = (list, tuple, np.ndarray, pd.Series)
        if isinstance(break_bond_sites, list_like) \
                and not isinstance(break_bond_sites[0], list_like):
            break_bond_sites = [list(break_bond_sites)]
        mol = self.mol
        reac1 = None
        reac2 = None
        # TODO: current break_bond method in pymatgen use distances to separate reactants, due to lost bonding message in Molecule object -qw
        # TODO: breaking circular mols can return same reac1 & reac2 but with different site sequence -qw
        # TODO: if reac1 and reac2 are same reactants and setting wrong break_bond_sites can have problems -qw
        for sites in break_bond_sites:
            reac1, reac2 = split_chain_mol(mol, sites)
            if MoleculeMatcher(reac1, reac2):
                mol = reac1
            else:
                break

        if self.write_output:
            reac1.to(filename=os.path.join(
                self.output_path, "{}_1_{}.{}".format(self.mol_name,
                                                      reac1.formula,
                                                      self.output_type)))
            reac2.to(filename=os.path.join(
                self.output_path, "{}_2_{}.{}".format(self.mol_name,
                                                      reac2.formula,
                                                      self.output_type)))


# TODO: add split_circular_mol -qw
def split_chain_mol(molecule, break_bond_sites):
    reac1, reac2 = molecule.break_bond(break_bond_sites[0],
                                       break_bond_sites[1])
    print("Split reactant 1: {}".format(reac1.formula))
    print("Split reactant 2: {}".format(reac2.formula))
    return reac1, reac2


# class ReactantJoin:





# def mol_from_file(path, file):
#     """
#     Get pymatgen molecule object from input file.
#     Args:
#         path:
#         file:
#
#     Returns:
#         pymatgen molecule object
#     """
#     input_file = os.path.join(path, file)
#     try:
#         molecule = IMolecule.from_file(input_file)
#         return molecule
#     except AttributeError:
#         print("This file {} cannot be parsed by pymatgen!".format(file))
#
#
# def replace_number_of_bonds():
#     with open(filename, 'r') as input_file, open('a.sdf', 'w') as output_file:
#         for line in input_file:
#             if line.strip() == '26 26  0  0  0  0            999 V2000':
#                 output_file.write(' 26 24  0  0  0  0            999 V2000\n')
#             else:
#                 output_file.write(line)
#
# def replace_1st_bond():
#     with open('a.sdf', 'r') as input_file, open('b.sdf', 'w') as output_file:
#         for line in input_file:
#             if line.strip() == '2  3  2  0  0  0  0':
#                 output_file.write('  2  3  1  0  0  0  0\n')
#             else:
#                 output_file.write(line)
#
# def replace_2nd_bond():
#     with open('b.sdf', 'r') as input_file, open('c.sdf', 'w') as output_file:
#         for line in input_file:
#             if line.strip() == '3  4  1  0  0  0  0':
#                 output_file.write('  3  4  2  0  0  0  0\n')
#             else:
#                 output_file.write(line)
#
# def replace_3rd_bond():
#     with open('c.sdf', 'r') as input_file, open('d.sdf', 'w') as output_file:
#         for line in input_file:
#             if line.strip() == '5  6  1  0  0  0  0':
#                 output_file.write('  5  6  2  0  0  0  0\n')
#             else:
#                 output_file.write(line)
#
# def replace_4th_bond():
#     with open('d.sdf', 'r') as input_file, open('e.sdf', 'w') as output_file:
#         for line in input_file:
#             if line.strip() == '1  2  1  0  0  0  0':
#                 output_file.write('  1  2  2  0  0  0  0\n')
#             else:
#                 output_file.write(line)
#
# def delete_5th_bond():
#     with open('e.sdf', 'r') as input_file, open('f.sdf', 'w') as output_file:
#         for line in input_file:
#             if line.strip() == '4  5  1  0  0  0  0':
#                 output_file.write('')
#             else:
#                 output_file.write(line)
#
# def delete_6th_bond():
#     with open('f.sdf', 'r') as input_file, open(filename, 'w') as output_file:
#         for line in input_file:
#             if line.strip() == '1  6  1  0  0  0  0':
#                 output_file.write('')
#             else:
#                 output_file.write(line)
#
# replace_number_of_bonds()
#
# replace_1st_bond()
#
# replace_2nd_bond()
#
# replace_3rd_bond()
#
# replace_4th_bond()
#
# delete_5th_bond()
#
# delete_6th_bond()
#
# os.remove('a.sdf')
# os.remove('b.sdf')
# os.remove('c.sdf')
# os.remove('d.sdf')
# os.remove('e.sdf')
# os.remove('f.sdf')
