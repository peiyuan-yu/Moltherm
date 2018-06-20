import re
import openbabel as ob
import pybel as pb
from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph

__author__ = "Qi Wang, Peiyuan Yu"
__version__ = "0.1"
__maintainer__ = "Qi Wang"
__email__ = "wqthu11@gmail.com"
__status__ = "Beta"
__date__ = "June 2018"


class OBMolecule:
    __hash__ = None

    def __init__(self, mol):
        """
        Creates a OBMolecule.

        """
        if isinstance(mol, ob.OBMol):
            self.obmol = mol
        elif isinstance(mol, Molecule):
            self.obmol = BabelMolAdaptor(mol).openbabel_mol

    @classmethod
    def from_file(cls, filename, fmt=None):
        fmt = filename.strip().split('.')[-1] if fmt is None else fmt
        fmt_std = re.search(r"(pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv)",
                            fmt.lower()).group(1)
        obmol = list(pb.readfile(str(fmt_std), str(filename)))[0].OBMol
        return cls(obmol)

    @classmethod
    def from_mol(cls, mol):
        if isinstance(mol, Molecule):
            obmol = BabelMolAdaptor(mol).openbabel_mol
        elif isinstance(mol, MoleculeGraph):
            obmol = BabelMolAdaptor(mol.molecule).openbabel_mol
        else:
            raise TypeError("This input mol format {} is not supported."
                            "Please input a pymatgen Molecule object or"
                            "MoleculeGraph".format(type(mol)))
        return cls(obmol)

    def to(self, filename=None, fmt=None):
        """
        Outputs the OBMolecule to a file.

        Args:
            fmt (str): Format to output to.
            filename (str): If provided, output will be written to a file. If
                fmt is not specified, the format is determined from the
                filename.

        Returns:
            (str) if filename is None. None otherwise.
        """
        fmt = filename.strip().split('.')[-1] if fmt is None else fmt
        m = re.search(r"\.(pdb|mol|mdl|sdf|sd|ml2|sy2|mol2|cml|mrv)",
                      fmt.lower())
        if m:
            pb.Molecule(self.obmol).write(fmt, filename)
        else:
            BabelMolAdaptor(self.obmol).pymatgen_mol.to(fmt, filename)
