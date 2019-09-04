import logging

import numpy as np

from pymatgen.core.structure import Molecule
from pymatgen.analysis.reaction_calculator import Reaction, ReactionError

from atomate.qchem.database import QChemCalcDb


__author__ = "Evan Spotte-Smith"
__copyright__ = "Copyright 2019, The Materials Project"
__version__ = "0.1"
__email__ = "espottesmith@gmail.com"

logger = logging.getLogger(__name__)


class ReactionRateCalculator:

    """
    An object which represents a chemical reaction (in terms of reactants, transition state,
    and products) and which can, from the energetics of those individual molecules, predict the
    rate constant, rate law, and thus the chemical kinetics of the reaction.
    """

    def __init__(self, reactants, products, transition_state, reaction=None):
        """
        If constructing this class manually, the use of one of the constructor methods is strongly
        preferred.

        NOTE: It is assumed that only one transition state is present.

        Args:
            reactants (list): list of dicts with reactant thermodynamics (energy, enthalpy, entropy)
                and other parameters (including Molecule object representing structure)
            products (list): list of dicts with product thermodynamics (energy, enthalpy, entropy)
                and other parameters (including Molecule object representing structure)
            transition_state (dict): dict with transition state thermodynamics (energy, enthalpy,
                entropy) and other parameters (including Molecule object representing structure)
            reaction (dict, or None): optional. If None (default), the "reactants" and
            "products" lists will serve as the basis for a Reaction object which represents the
            balanced stoichiometric reaction. Otherwise, this dict will show the number of molecules
            present in the reaction for each reactant and each product in the reaction.

        Returns:
            None
        """

        self.reactants = reactants
        self.products = products
        self.transition_state = transition_state

        if reaction is None:
            rct_mols = [r["molecule"] for r in self.reactants]
            pro_mols = [p["molecule"] for p in self.reactions]
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
            None
        """

        rcts = list()
        pros = list()

        # Construct main dicts with molecular properties
        for rct_doc in reactants:
            reactant = dict()
            try:
                reactant["task_id"] = rct_doc["task_id"]
                reactant["label"] = rct_doc["task_label"]
                reactant["smiles"] = rct_doc["smiles"]
                # convert all energies from Hartree to kcal/mol
                reactant["energy"] = rct_doc["output"]["final_energy"] * 627.509474
                reactant["enthalpy"] = rct_doc["output"]["enthalpy"]
                reactant["entropy"] = rct_doc["output"]["entropy"]
                reactant["molecule"] = Molecule.from_dict(rct_doc["output"]["optimized_molecule"])

                rcts.append(reactant)

            except KeyError:
                raise ValueError("Reactant task doc does not follow schema! Docs must contain"
                                 " an output field with subfields optimized_geometry, final_energy,"
                                 " enthalpy, and entropy.")

        for pro_doc in products:
            product = dict()

            try:
                product["task_id"] = pro_doc["task_id"]
                product["label"] = pro_doc["task_label"]
                product["smiles"] = pro_doc["smiles"]
                # convert all energies from Hartree to kcal/mol
                product["energy"] = pro_doc["output"]["final_energy"] * 627.509474
                product["enthalpy"] = pro_doc["output"]["enthalpy"]
                product["entropy"] = pro_doc["output"]["entropy"]
                product["molecule"] = Molecule.from_dict(pro_doc["output"]["optimized_molecule"])

                pros.append(product)

            except KeyError:
                raise ValueError("Product task doc does not follow schema! Docs must contain"
                                 " an output field with subfields optimized_geometry, final_energy,"
                                 " enthalpy, and entropy.")

        transition = dict()
        try:
            transition["task_id"] = transition_state["task_id"]
            transition["label"] = transition_state["task_label"]
            transition["smiles"] = transition_state["smiles"]
            # convert all energies from Hartree to kcal/mol
            transition["energy"] = transition_state["output"]["final_energy"] * 627.509474
            transition["enthalpy"] = transition_state["output"]["enthalpy"]
            transition["entropy"] = transition_state["output"]["entropy"]
            transition["molecule"] = Molecule.from_dict(transition_state["output"]["optimized_molecule"])

        except KeyError:
            raise ValueError("Transition state task doc does not follow schema! Docs must contain"
                             " an output field with subfields optimized_geometry, final_energy,"
                             " enthalpy, and entropy.")

        # Calculate stoichiometry
        rct_mols = [r["molecule"] for r in rcts]
        pro_mols = [p["molecule"] for p in pros]
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