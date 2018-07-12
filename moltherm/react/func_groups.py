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

            else:
                self.molecule = molecule.molecule
                self.molgraph = molecule

        if self.molgraph is None:
            self.molgraph = MoleculeGraph.with_local_env_strategy(self.molecule,
                                                                  OpenBabelNN())

        # Assign a specie and coordinates to each node in the graph,
        # corresponding to the Site in the Molecule object
        self.molgraph.set_node_attributes()

    def get_heteroatoms(self, elements=None):
        """
        Identify non-H, non-C atoms in the MoleculeGraph, returning a list of
        their node indices.

        :param elements: List of elements to identify (if only certain
            functional groups are of interest).
        :return: list of ints representing node indices
        """

        species = nx.get_node_attributes(self.molgraph.graph, "specie")

        heteroatoms = []

        for node in self.molgraph.graph.nodes():
            if elements is not None:
                if str(species[node]) in elements:
                    heteroatoms.append(node)
            else:
                if str(species[node]) not in ["C", "H"]:
                    heteroatoms.append(node)

        return heteroatoms

    def get_special_carbon(self, elements=None):
        """
        Identify Carbon atoms in the MoleculeGraph that fit the characteristics
        defined Ertl (2017), returning a list of their node indices.

        :param elements: List of elements that will qualify a carbon as special
            (if only certain functional groups are of interest).
            Default None.
        :return: list of ints representing node indices.
        """
        pass

    def link_marked_atoms(self, atoms):
        """
        Take a list of marked "interesting" atoms (heteroatoms, special carbons)
        and attempt to connect them, returning a list of disjoint groups of
        special atoms (and their connected hydrogens).

        :param atoms: List of marked "interesting" atoms, presumably identified
            using other functions in this class.
        :return: list of lists of ints, representing groups of connected atoms.
        """
        pass

    def get_basic_functional_groups(self, func_groups=None):
        """
        Identify functional groups that cannot be identified by the Ertl method
        of get_special_carbon and get_heteroatoms, such as benzene rings, methyl
        groups, and ethyl groups.

        :param func_groups: List of strs representing the functional groups of
            interest. Default to None, meaning that all of the functional groups
            defined in this function will be sought.
        :return:
        """
        pass

    def get_all_functional_groups(self, elements=None, func_groups=None):
        """
        Identify all functional groups (or all within a certain subset) in the
        molecule, combining the methods described above.

        :param elements: List of elements that will qualify a carbon as special
            (if only certain functional groups are of interest).
            Default None.
        :param func_groups: List of strs representing the functional groups of
            interest. Default to None, meaning that all of the functional groups
            defined in this function will be sought.
        :return: List of lists of ints, representing groups of connected atoms.
        """
        pass

    def categorize_functional_groups(self, groups):
        """
        Determine classes of functional groups present in a set.

        :param groups: Set of functional groups.
        :return: dict containing representations of the groups, the indices of
            where the group occurs in the MoleculeGraph, and how many of each
            type of group there is.
        """
        pass

    def make_smarts(self, group):
        """
        Transforms a functional group into a SMARTS string representation.

        :param group: List of indices representing a functional group.
        :return: str, a SMARTS string.
        """
        pass