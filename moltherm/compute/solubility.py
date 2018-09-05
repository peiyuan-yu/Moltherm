from unifac.facade import Facade

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.io.babel import BabelMolAdaptor


class MixtureCalculator:
    """
    This class uses the UNIFAC group contribution method of activity coefficient
    calculation to estimate the solubility of molecules in solutions which may
    involve a number of components.
    """

    def __init__(self, components, temperature=273.15):
        """
        :param components: dict ({Molecule: float}, {MoleculeGraph: float}
            or {str: float}, where the str represents a molecule SMILES
            string and float represents the quantity of the component in the
            mixture.
        :param temperature: float representing the temperature of the system in
            Kelvin (273.15 K by default)
        """

        self.calculator = Facade()
        self.components = {}

        for component, quantity in components.items():
            if isinstance(component, str):
                self.add_smiles(component, quantity)
            elif isinstance(component, Molecule):
                self.add_molecule(component, quantity)
            elif isinstance(component, MoleculeGraph):
                self.add_molgraph(component, quantity)

        self.calculator.set_temperature(temperature)

    def set_temperature(self, temperature):
        """
        Wrapper function for Facade.set_temperature

        :param temperature: float representing the temperature of the system in
            Kelvin (273.15 K by default)
        :return:
        """

        self.calculator.set_temperature(temperature)

    def add_smiles(self, smiles, quantity):
        """
        Wrapper function for Facade.add_molecule_smiles.

        :param smiles: str representing molecule SMILES string
        :param quantity: float representing the quantity of the component in the
            mixture.
        :return:
        """

        self.calculator.add_molecule_smiles(smiles, quantity)
        self.components[smiles] = quantity

    def add_molecule(self, mol, quantity):
        """
        Wrapper function for Facade.add_molecule_smiles using pymatgen Molecule
        class.

        :param mol: pymatgen.core.structure.Molecule object
        :param quantity: float representing the quantity of the component in the
            mixture.
        :return:
        """

        adaptor = BabelMolAdaptor(mol)
        smiles = adaptor.pybel_mol.write("smi")

        self.calculator.add_molecule_smiles(smiles, quantity)
        self.components[smiles] = quantity

    def add_molgraph(self, mg, quantity):
        """
        Wrapper function for Facade.add_molecule_smiles using pymatgen
        MoleculeGraph class.

        :param mg: pymatgen.analysis.graphs.MoleculeGraph object
        :param quantity: float representing the quantity of the component in the
            mixture.
        :return:
        """

        adaptor = BabelMolAdaptor(mg.molecule)
        smiles = adaptor.pybel_mol.write("smi")

        self.calculator.add_molecule_smiles(smiles, quantity)
        self.components[smiles] = quantity

    def edit_component(self, smiles, quantity):
        """
        Change the quantity of a component in the mixture.

        :param smiles: str representing molecule SMILES string to be adjusted
        :param quantity: float representing the quantity of the component in the
            mixture.
        :return:
        """

        self.components[smiles] = quantity

        self.calculator.reset_solution()

        for c, q in self.components.items():
            self.add_smiles(c, q)

    def remove_component(self, smiles):
        """
        Remove a component completely from the mixture.

        :param smiles: str representing molecule SMILES string to be removed
        :return:
        """

        del self.components[smiles]

        self.calculator.reset_solution()

        for c, q in self.components.items():
            self.add_smiles(c, q)

    def get_coeff(self, smiles):
        """
        Wrapper for Facade.getcoeff.

        :param smiles: str representing molecule SMILES string of molecule to be
            considered
        :return: Molecule activity coefficient in the current mixture
        """

        return self.calculator.get_coeff(smiles)