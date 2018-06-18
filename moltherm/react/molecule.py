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


# Todo: consider whether to move confm_search to BabelMolAdaptor. -qw
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

    def confm_search(self, forcefield="mmff94", freeze_atoms=None,
                     rmsd_cutoff=0.5, energy_cutoff=50.0,
                     conf_cutoff=100000, verbose=False, make_3d=False,
                     add_hydrogens=False):
        """
        Perform conformer search.
        Args:
            forcefield:
            freeze_atoms:
            rmsd_cutoff:
            energy_cutoff:
            conf_cutoff:
            verbose:
            make_3d:
            add_hydrogens:
        Returns:

        """

        # Had to remove the copying, as apparently self.obmol cannot be copied
        obmol = self.obmol

        # TODO: Does this belong here? Is this functionality general enough?
        if make_3d:
            builder = ob.OBBuilder()
            builder.Build(obmol)

        if add_hydrogens:
            obmol.AddHydrogens()


        ff = ob.OBForceField_FindType(forcefield)
        if ff == 0:
            print("Could not find forcefield {} in openbabel, the forcefield "
                  "will be reset as default 'mmff94'".format(forcefield))
            ff = ob.OBForceField_FindType("mmff94")
        # Make sure setup works
        assert (ff.Setup(obmol))

        if freeze_atoms:
            print('{} atoms will be freezed'.format(len(freeze_atoms)))
            constraints = ob.OBFFConstraints()
            for atom in ob.OBMolAtomIter(obmol):
                atom_id = atom.GetIndex() + 1
                if id in freeze_atoms:
                    constraints.AddAtomConstraint(atom_id)
            ff.SetConstraints(constraints)

        # To improve 3D coordinates
        # TODO: If we integrate this into BabelMolAdapter, we can use localopt
        # to optimize 3D structure
        ff.WeightedRotorSearch(500, 25)

        # Run Confab conformer generation
        ff.DiverseConfGen(rmsd_cutoff, conf_cutoff, energy_cutoff,
                          verbose)

        ff.GetConformers(obmol)
        print("Generated {} conformers total".format(obmol.NumConformers()))
        # return the best conformer.
        return OBMolecule(obmol)

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
