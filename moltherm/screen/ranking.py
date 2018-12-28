import numpy as np
import pandas as pd

# TODO: Do we/will we need these?
# from pymatgen.core.structure import Molecule
# from pymatgen.io.babel import BabelMolAdaptor

from atomate.qchem.database import QChemCalcDb


__author__ = "Evan Spotte-Smith"
__version__ = "0.2"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Beta"
__date__ = "December 2018"


class ReactionRanker:
    """
    This class includes several ranking algorithms which can distinguish between
    different reactions. For each algorithm, flexibility in parameters is
    allowed, and in general, the same set of parameters can be used for each
    algorithm.

    At present, the allowed set of parameters includes the following:
    bp - normal boiling point
    mp - normal melting point
    vp - vapor pressure
    solubility - aqueous solubility
    low_kow - octanol-water partition coefficient
    """

    def __init__(self, reactions, db_file="db.json", collection="episuite",
                 parameters_min=None, parameters_max=None, goals=None):
        """
        :param reactions: list of dicts representing reactions to be ranked
        :param db_file: str representing a path to a database config file
        :param collection: Collection in database which contains molecule
            information.
        :param parameters_min: list of strings, each representing a parameter
            which can be used in ranking. High-ranking reactions will minimize
             these paramaters. Default is None, which means that the full
            allowed set of parameters will be used.
        :param parameters_max: list of strings, each representing a parameter
            which can be used in ranking. High-ranking reactions will maximize
             these paramaters. Default is None, which means that the full
            allowed set of parameters will be used.
        :param goals: dict with {parameter: value} pairs. These key-value pairs
            represent ideal values which the molecule parameters will be
            compared against. If goals is None (default), all parameters will be
            compared to 0. If within goals a particular paramter is not present,
            or is set to None, then this parameter will be set to 0.
        """

        self.db_file = db_file
        try:
            #TODO: Make a custom database class, instead of relying on the
            # atomate classes
            self.db = QChemCalcDb.from_db_file(self.db_file)
        except:
            self.db = None

        if self.db is not None:
            self.collection = self.db.db[collection]

        if parameters_min is None:
            self.parameters_min = ["mp", "vp", "solubility"]
        else:
            self.parameters_min = parameters_min

        if parameters_max is None:
            self.parameters_max = ["bp", "log_kow"]
        else:
            self.parameters_max = parameters_max

        self.parameters = self.parameters_min + self.parameters_max

        if goals is None:
            self.goals = {key: 0 for key in self.parameters}
        else:
            self.goals = {}
            for parameter in self.parameters:
                if parameter in goals:
                    if goals[parameter] is None:
                        self.goals[parameter] = 0
                    else:
                        self.goals[parameter] = goals[parameter]
                else:
                    self.goals[parameter] = 0

        # Filter reactions to ensure that all molecules are in DB
        self.reactions = []
        for rxn in reactions:
            rcts = rxn["rct_ids"]
            entries = [self.collection.find_one({"mol_id": r}) for r in rcts]
            if all([e is not None for e in entries]):
                self.reactions.append(rxn)

        # Information about individual molecules
        self.molecules = self._construct_molecules()
        # Information about all reactants in a reaction
        self.all_reactants = self._construct_all_reactants()

    def _construct_molecules(self):
        """
        Uses data from database to construct pandas DataFrame with molecule
            data, relative to goal values.

        Note: ALL molecule entries in the database (and the collection
            self.collection, specifically) must have values for ALL keys in
            self.parameters.

        :return: molecules (pd.DataFrame)
        """

        all_molecules = [e for e in self.collection.find()]

        mol_ids = np.array([e["mol_id"] for e in all_molecules])

        data = {key: [] for key in self.parameters}
        for mol in all_molecules:
            for key in self.parameters:
                val = mol[key]
                if key in self.parameters_min:
                    val = self.goals[key] - mol[key]
                elif key in self.parameters_max:
                    val = mol[key] - self.goals[key]
                data[key].append(val)

        return pd.DataFrame(data=data, index=mol_ids)

    def _construct_all_reactants(self):
        """
        Uses data from database to construct pandas DataFrame with data compiled
            from all reactants relevant to each reaction.

        :return: all_reactants (pd.DataFrame)
        """

        rxn_ids = np.array([e["rxn_id"] for e in self.reactions])

        data = {key: [] for key in self.parameters}
        for rxn in self.reactions:
            rcts = rxn["rct_ids"]
            entries = [self.collection.find_one({"mol_id": r}) for r in rcts]

            for key in self.parameters:
                total = 0
                for entry in entries:
                    if key in self.parameters_min:
                        total += self.goals[key] - entry[key]
                    elif key in self.parameters_max:
                        total += entry[key] - self.goals[key]

                data[key].append(total)

        return pd.DataFrame(data=data, index=rxn_ids)

    def pareto_ranking(self, reactions=None, parameters=None, by_worst=False):
        pass

    def heuristic_ranking(self, reactions=None, parameters=None,
                          by_worst=False):
        pass

    def tiered_ranking(self, reactions=None, parameters=None, by_worst=False):
        pass