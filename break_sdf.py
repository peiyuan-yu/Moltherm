from sys import argv
import os
from pymatgen.core.structure import Molecule

# script, filename = argv


class ProductSplit:

    def __init__(self, input_path, input_file, output_path,
                 breaking_bond_sites):
        self.input_path = input_path
        self.input_file = input_file
        self.output_path = output_path
        self.breaking_bond_sites = breaking_bond_sites

    def split(self):
        f = os.path.join(self.input_path, self.input_file)
        print(f)
        # try:
        product = Molecule.from_file(f)
        # except AttributeError:
        #     print("This file {} cannot be parsed by pymatgen!".
        #           format(self.input_file))
        (reac1, reac2) = product.break_bond(self.breaking_bond_sites[0],
                                            self.breaking_bond_sites[1])
        print("Splitted reactant 1: {}".format(reac1.formula))
        print("Splitted reactant 2: {}".format(reac2.formula))

        (header, suffix) = os.path.splitext(self.input_file)
        reac1.to(filename=os.path.join(self.output_path,
                                       header + "_1_" + reac1.formula + suffix))
        reac2.to(filename=os.path.join(self.output_path,
                                       header + "_2_" + reac2.formula + suffix))


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
