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
        # TODO: current break_bond method use distances to separate reactants, due to lost bonding message in Molecule object. -qw
        # TODO: only supports break_chain_mol now. -qw
        # TODO: output filename of react1 and react2 needs improvement as one molecule can be split many times. -qw
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

