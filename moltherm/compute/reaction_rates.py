import logging
import operator
from functools import reduce

import numpy as np
from scipy.constants import h, k, R, N_A, pi

from pymatgen.core.structure import Molecule
from pymatgen.core.units import amu_to_kg
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.entries.mol_entry import MoleculeEntry
from pymatgen.analysis.reaction_calculator import Reaction, ReactionError

from atomate.qchem.database import QChemCalcDb


__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__email__ = "espottesmith@gmail.com"

logger = logging.getLogger(__name__)


def product(iterable):
    return reduce(operator.mul, iterable, 1)


class ReactionRateCalculator:

    """
    An object which represents a chemical reaction (in terms of reactants, transition state,
    and products) and which can, from the energetics of those individual molecules, predict the
    rate constant, rate law, and thus the chemical kinetics of the reaction.

    NOTE: It is assumed that only one transition state is present.

    Args:
        reactants (list): list of MoleculeEntry objects
        products (list): list of MoleculeEntry objects
        transition_state (MoleculeEntry): MoleculeEntry representing the transition state between
            the reactants and the products
        reaction (dict, or None): optional. If None (default), the "reactants" and
        "products" lists will serve as the basis for a Reaction object which represents the
        balanced stoichiometric reaction. Otherwise, this dict will show the number of molecules
        present in the reaction for each reactant and each product in the reaction.

    Returns:
        None
    """

    def __init__(self, reactants, products, transition_state, reaction=None):
        """

        """

        self.reactants = reactants
        self.products = products
        self.transition_state = transition_state

        if reaction is None:
            rct_mols = [r.mol_graph.molecule for r in self.reactants]
            pro_mols = [p.mol_graph.molecule for p in self.reactions]
            try:
                self.reaction = self.generate_reaction(rct_mols, pro_mols)
            except ReactionError:
                # Reaction cannot be balanced
                self.reaction = None
        else:
            self.reaction = reaction

    @classmethod
    def from_atomate_tasks(cls, reactants, products, transition_state):
        """
        Constructor using task documents (MSON-dicts) from atomate workflows.

        Note: This constructor is NOT FLEXIBLE, and requires that task docs at least have access to
        an "output" field (as well as subfields "optimized_geometry", "final_energy", "enthalpy",
        and "entropy").

        Args:
            reactants (list): list of dicts representing task docs for reactant molecules
            products (list): list of dicts representing task docs for product molecules
            transition_state (dict): dict representing task doc for transition state molecule

        Returns:
            ReactionRateCalculator
        """

        rcts = list()
        pros = list()

        # Construct main dicts with molecular properties
        for rct_doc in reactants:
            try:
                rcts.append(MoleculeEntry(Molecule.from_dict(rct_doc["output"]["optimized_molecule"]),
                                          rct_doc["output"]["final_energy"],
                                          enthalpy=rct_doc["output"]["enthalpy"],
                                          entropy=rct_doc["output"]["entropy"],
                                          entry_id=rct_doc["task_id"]))
            except KeyError:
                raise ValueError("Reactant task doc does not follow schema! Docs must contain"
                                 " an output field with subfields optimized_geometry, final_energy,"
                                 " enthalpy, and entropy.")

        for pro_doc in products:
            try:
                pros.append(MoleculeEntry(Molecule.from_dict(pro_doc["output"]["optimized_molecule"]),
                                          pro_doc["output"]["final_energy"],
                                          enthalpy=pro_doc["output"]["enthalpy"],
                                          entropy=pro_doc["output"]["entropy"],
                                          entry_id=pro_doc["task_id"]))

            except KeyError:
                raise ValueError("Product task doc does not follow schema! Docs must contain"
                                 " an output field with subfields optimized_geometry, final_energy,"
                                 " enthalpy, and entropy.")

        try:
            transition = MoleculeEntry(Molecule.from_dict(transition_state["output"]["optimized_molecule"]),
                                       transition_state["output"]["final_energy"],
                                       enthalpy=transition_state["output"]["enthalpy"],
                                       entropy=transition_state["output"]["entropy"],
                                       entry_id=transition_state["task_id"])
        except KeyError:
            raise ValueError("Transition state task doc does not follow schema! Docs must contain"
                             " an output field with subfields optimized_geometry, final_energy,"
                             " enthalpy, and entropy.")

        # Calculate stoichiometry
        rct_mols = [r.mol_graph.molecule for r in rcts]
        pro_mols = [p.mol_graph.molecule for p in pros]
        reaction = cls.generate_reaction(rct_mols, pro_mols)

        return cls(rcts, pros, transition, reaction=reaction)

    @classmethod
    def generate_reaction(cls, rct_mols, pro_mols):
        """
        Generate a pymatgen Reaction object, given lists of reactant and product Molecule objects.

        Args:
            rct_mols (list): list of Molecule objects representing reactants
            pro_mols (list): list of Molecule objects representing products

        Note: if reaction cannot be balanced, this method will raise a ReactionError

        Returns:
            reaction (Reaction): pymatgen Reaction object representing the reaction between the
                reactants and products.
        """

        rct_comps = [r.composition for r in rct_mols]
        pro_comps = [p.composition for p in pro_mols]

        return Reaction(rct_comps, pro_comps)

    @property
    def net_energy(self):
        """
        Calculate net reaction energy (in Hartree).
        """
        rct_energies = [r.energy for r in self.reactants]
        pro_energies = [p.energy for p in self.products]

        return sum(pro_energies) - sum(rct_energies)

    @property
    def net_enthalpy(self):
        """
        Calculate net reaction enthalpy (in kcal/mol).
        """
        rct_enthalpies = [r.enthalpy for r in self.reactants]
        pro_enthalpies = [p.enthalpy for p in self.products]

        return sum(pro_enthalpies) - sum(rct_enthalpies)

    @property
    def net_entropy(self):
        """
        Calculate net reaction entropy (in cal/mol-K).
        """
        rct_entropies = [r.entropy for r in self.reactants]
        pro_entropies = [p.entropy for p in self.products]

        return sum(pro_entropies) - sum(rct_entropies)

    def calculate_net_gibbs(self, temperature):
        """
        Calculate net reaction Gibbs free energy at a given temperature.

        ΔG = ΔH - T ΔS

        Args:
            temperature (float): absolute temperature in Kelvin

        Returns:
            float: net Gibbs free energy (in kcal/mol)
        """

        return self.net_enthalpy - temperature * self.net_entropy / 1000

    def calculate_net_thermo(self, temperature=300.0):
        """
        Calculate net energy, enthalpy, and entropy for the reaction.

        Args:
            temperature (float): absolute temperature in Kelvin (default 300.0K)

        Returns:
            net_thermo: dict with relevant net thermodynamic variables
        """

        net_thermo = {"energy": self.net_energy,
                      "enthalpy": self.net_enthalpy,
                      "entropy": self.net_entropy,
                      "gibbs": self.calculate_net_gibbs(temperature)}

        return net_thermo

    def calculate_act_energy(self, reverse=False):
        """
        Calculate energy of activation.

        Args:
            reverse (bool): if True (default False), consider the reverse reaction; otherwise,
                consider the forwards reaction

        Returns:
            float: energy of activation (in Hartree)

        """

        trans_energy = self.transition_state.energy

        if reverse:
            pro_energies = [p.energy for p in self.products]
            return trans_energy - sum(pro_energies)
        else:
            rct_energies = [r.energy for r in self.reactants]
            return trans_energy - sum(rct_energies)

    def calculate_act_enthalpy(self, reverse=False):
        """
        Calculate enthalpy of activation.

        Args:
            reverse (bool): if True (default False), consider the reverse reaction; otherwise,
                consider the forwards reaction

        Returns:
            float: enthalpy of activation (in kcal/mol)

        """

        trans_enthalpy = self.transition_state.enthalpy

        if reverse:
            pro_enthalpies = [p.enthalpy for p in self.products]
            return trans_enthalpy - sum(pro_enthalpies)
        else:
            rct_enthalpies = [r.enthalpy for r in self.reactants]
            return trans_enthalpy - sum(rct_enthalpies)

    def calculate_act_entropy(self, reverse=False):
        """
        Calculate entropy of activation.

        Args:
            reverse (bool): if True (default False), consider the reverse reaction; otherwise,
                consider the forwards reaction

        Returns:
            float: entropy of activation (in cal/mol-K)

        """

        trans_entropy = self.transition_state.entropy

        if reverse:
            pro_entropies = [p.entropy for p in self.products]
            return trans_entropy - sum(pro_entropies)
        else:
            rct_entropies = [r.entropy for r in self.reactants]
            return trans_entropy - sum(rct_entropies)

    def calculate_act_gibbs(self, temperature, reverse=False):
        """
        Calculate Gibbs free energy of activation at a given temperature.

        ΔG = ΔH - T ΔS

        Args:
            temperature (float): absolute temperature in Kelvin
            reverse (bool): if True (default False), consider the reverse reaction; otherwise,
                consider the forwards reaction

        Returns:
            float: Gibbs free energy of activation (in kcal/mol)
        """

        return self.calculate_act_enthalpy(reverse=reverse) - temperature * self.calculate_act_entropy(reverse=reverse) / 1000

    def calculate_activation_thermo(self, temperature=300.0, reverse=False):
        """
        Calculate thermodynamics of activation for the reaction.

        Args:
            temperature (float): absolute temperature in Kelvin (default 300.0K)
            reverse (bool): if True (default False), consider the reverse reaction; otherwise,
                consider the forwards reaction

        Returns:
            act_thermo: dict with relevant activation thermodynamic variables
        """

        act_thermo = {"energy": self.calculate_act_energy(reverse=reverse),
                      "enthalpy": self.calculate_act_enthalpy(reverse=reverse),
                      "entropy": self.calculate_act_entropy(reverse=reverse),
                      "gibbs": self.calculate_act_gibbs(temperature, reverse=reverse)}

        return act_thermo

    def calculate_rate_constant(self, temperature=300.0, reverse=False, kappa=1.0):
        """
        Calculate the rate constant k by the Eyring-Polanyi equation of transition state theory.

        Args:
            temperature (float): absolute temperature in Kelvin
            reverse (bool): if True (default False), consider the reverse reaction; otherwise,
                consider the forwards reaction
            kappa (float): transmission coefficient (by default, we assume the assumptions of
                transition-state theory are valid, so kappa = 1.0

        Returns:
            k_rate (float): temperature-dependent rate constant
        """

        gibbs = self.calculate_act_gibbs(temperature=temperature, reverse=reverse)

        k_rate = kappa * k * temperature / h * np.exp(-gibbs * 4184 / (R * temperature))
        return k_rate

    def calculate_rate(self, concentrations, temperature=300.0, reverse=False, kappa=1.0):
        """
        Calculate the based on the reaction stoichiometry.

        NOTE: Here, we assume that the reaction is an elementary step.

        Args:
            concentrations (list): concentrations of reactant molecules. Order of the reactants
                DOES matter.
            temperature (float): absolute temperature in Kelvin
            reverse (bool): if True (default False), consider the reverse reaction; otherwise,
                consider the forwards reaction
            kappa (float): transmission coefficient (by default, we assume the assumptions of
                transition-state theory are valid, so kappa = 1

        Returns:
            rate (float): reaction rate, based on the stoichiometric rate law and the rate constant
        """

        rate_constant = self.calculate_rate_constant(temperature=temperature, reverse=reverse,
                                                     kappa=kappa)

        if reverse:
            pro_comps = [p.mol_graph.molecule.composition for p in self.products]
            exponents = np.array([self.reaction.get_coeff(comp) for comp in pro_comps])
            rate = rate_constant * np.array(concentrations) ** exponents
        else:
            rct_comps = [r.mol_graph.molecule.composition for r in self.reactants]
            exponents = np.array([self.reaction.get_coeff(comp) for comp in rct_comps])
            rate = rate_constant * np.array(concentrations) ** exponents

        return rate


class BellEvansPolanyiRateCalculator(ReactionRateCalculator):
    """
    A modified reaction rate calculator that uses the Bell-Evans-Polanyi principle to predict the
    activation energies (and, thus, the rate constants and reaction rates) of chemical reactions).

    The Bell-Evans-Polanyi principle states that, for reactions within a similar class or family,
    the difference in activation energy between the reactions is proportional to the difference in
    reaction enthalpy. That is,

    E_a = E_a,0 + alpha * ΔH, where

    E_a = the activation energy of the reaction of interest
    E_a,0 = the activation energy of some reference reaction in the same reaction family
    alpha = the location of the transition state along the reaction coordinate (0 <= alpha <= 1)
    ΔH = the enthalpy of reaction

    Whereas ReactionRateCalculator uses the Eyring equation, here we are forced to use collision
    theory to estimate reaction rates.

    Args:
        reactants (list): list of MoleculeEntry objects
        products (list): list of MoleculeEntry objects
        ea_reference (float): activation energy reference point (in kcal/mol)
        delta_h_reference (float): reaction enthalpy reference point (in kcal/mol)
        reaction (dict, or None): optional. If None (default), the "reactants" and
        "products" lists will serve as the basis for a Reaction object which represents the
        balanced stoichiometric reaction. Otherwise, this dict will show the number of molecules
        present in the reaction for each reactant and each product in the reaction.
        alpha (float): the reaction coordinate (must between 0 and 1)
    """

    def __init__(self, reactants, products, ea_reference, delta_h_reference, reaction=None,
                 alpha=0.5):
        """

        """

        self.ea_reference = ea_reference
        self.delta_h_reference = delta_h_reference
        self.alpha = alpha

        super().__init__(reactants, products, None, reaction=reaction)

    def calculate_act_energy(self, reverse=False):
        """
        Use the Bell-Evans-Polanyi principle to calculate the activation energy of the reaction.

        Args:
            reverse (bool): if True (default False), consider the reverse reaction; otherwise,
                consider the forwards reaction

        Returns:
            ea (float): the predicted energy of activation in kcal/mol
        """

        if reverse:
            enthalpy = -self.net_enthalpy
        else:
            enthalpy = self.net_enthalpy

        ea = self.ea_reference + self.alpha * (enthalpy - self.delta_h_reference)
        return ea

    def calculate_act_enthalpy(self, reverse=False):
        raise NotImplementedError("Method calculate_act_enthalpy is not valid for "
                                  "BellEvansPolanyiRateCalculator,")

    def calculate_act_entropy(self, reverse=False):
        raise NotImplementedError("Method calculate_act_entropy is not valid for "
                                  "BellEvansPolanyiRateCalculator,")

    def calculate_act_gibbs(self, temperature, reverse=False):
        raise NotImplementedError("Method calculate_act_gibbs is not valid for "
                                  "BellEvansPolanyiRateCalculator,")

    def calculate_activation_thermo(self, temperature=300.0, reverse=False):
        raise NotImplementedError("Method calculate_activation_thermo is not valid for "
                                  "BellEvansPolanyiRateCalculator,")

    def calculate_rate_constant(self, temperature=300.0, reverse=False, kappa=None):
        """
        Calculate the rate constant predicted by collision theory.

        Args:
            temperature (float): absolute temperature in Kelvin
            reverse (bool): if True (default False), consider the reverse reaction; otherwise,
                consider the forwards reaction
            kappa (None): not used for BellEvansPolanyiRateCalculator

        Returns:
            k_rate (float): temperature-dependent rate constant
        """

        ea = self.calculate_act_energy(reverse=reverse)

        k_rate = np.exp(-ea * 4184 / (R * temperature))

        return k_rate

    def calculate_rate(self, concentrations, temperature=300.0, reverse=False, kappa=1.0):
        """
        Calculate the rate using collision theory.

        Args:
            concentrations (list): concentrations of reactant molecules. Order of the reactants
                DOES matter.
            temperature (float): absolute temperature in Kelvin
            reverse (bool): if True (default False), consider the reverse reaction; otherwise,
                consider the forwards reaction
            kappa (float): here, kappa represents the steric factor (default 1.0, meaning that all
                collisions lead to appropriate conditions for a reaction

        Returns:
            rate (float): reaction rate, based on the stoichiometric rate law and the rate constant
        """

        k_rate = self.calculate_rate_constant(temperature=temperature, reverse=reverse)

        if reverse:
            mols = [p.mol_graph.molecule for p in self.products]
        else:
            mols = [r.mol_graph.molecule for r in self.reactants]

        masses = [m.composition.weight for m in mols]
        exponents = np.array([self.reaction.get_coeff(mol.composition) for mol in mols])

        radius_factor = pi * sum([(max(mol.distance_matrix) / 2) ** 2 for mol in mols])

        total_exponent = sum(exponents)
        number_prefactor = (1000 * N_A) ** total_exponent
        concentration_factor = np.array(concentrations) ** exponents
        mass_factor = product(masses)/sum(masses) * amu_to_kg
        root_factor = np.sqrt(8 * k * temperature / (pi * mass_factor))

        z = number_prefactor * concentration_factor * radius_factor * root_factor

        rate = z * kappa * k_rate
        return rate

