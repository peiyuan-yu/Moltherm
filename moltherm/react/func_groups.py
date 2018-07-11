import openbabel as ob
import pybel as pb

from pymatgen.core.strucuture import Molecule, FunctionalGroups
from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

from moltherm.compute.utils import get_molecule

from

try:
    import networkx as nx
except ImportError:
    raise ImportError("moltherm.react.func_groups requires the NetworkX "
                      "graph library to be installed.")


__author__ = "Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "ewcspottesmith@lbl.gov"
__status__ = "Alpha"
__date__ = "July 2018"


class FunctionalGroupExtractor:
    """
    This class is used to algorithmically parse a molecule (represented by an
    instance of pymatgen.analysis.graphs.MoleculeGraph) and determine arbitrary
    functional groups.
    """

    def __init__(self, molecule, optimize=False):
        """
        Instantiation method for FunctionalGroupExtractor.

        :param molecule: Either a filename, a pymatgen.core.structure.Molecule
            object, or a pymatgen.analysis.graphs.MoleculeGraph object.
        :param optimize: Default False. If True, then the input molecule will be
            modified, adding Hydrogens, performing a simple conformer search,
            etc.
        """

        self.molgraph = None

        if isinstance(molecule, str):
            try:
                if optimize:
                    self.molecule = get_molecule(molecule)
                else:
                    self.molecule = Molecule.from_file(molecule)
            except OSError:
                raise ValueError("Input must be a valid molecule file, a "
                                 "Molecule object, or a MoleculeGraph object.")

        elif isinstance(molecule, Molecule):
            if optimize:
                obmol = BabelMolAdaptor(molecule)
                obmol.add_hydrogen()
                obmol.make3d()
                obmol.localopt()

                self.molecule = obmol.pymatgen_mol
            else:
                self.molecule = molecule

        elif isinstance(molecule, MoleculeGraph):
            if optimize:
                obmol = BabelMolAdaptor(molecule.molecule)
                obmol.add_hydrogen()
                obmol.make3d()
                obmol.localopt()

                self.molecule = obmol.pymatgen_mol
                molecule.molecule = obmol.pymatgen_mol

                self.molgraph = molecule
            else:
                self.molecule = molecule.molecule
                self.molgraph = molecule

        if self.molgraph is None:
            self.molgraph = MoleculeGraph.with_local_env_strategy(self.molecule,
                                                                  OpenBabelNN())