from random import choice

import numpy as np
import pandas as pd


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

    def __init__(self, reactions, molecules, parameters_min=None,
                 parameters_max=None, goals=None):
        """
        :param reactions: list of dicts representing reactions to be ranked
        :param molecules: list of dicts representing the molecules relevant
            to the reactions list
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

        # Filter reactions to ensure that all molecules are in molecule set
        mol_ids = [e["mol_id"] for e in molecules]
        self.reactions = []
        for rxn in reactions:
            rcts = rxn["rct_ids"]
            entries = [r if r in mol_ids else None for r in rcts]
            if all([e is not None for e in entries]):
                self.reactions.append(rxn)

        self.molecules_raw = molecules

        # Information about all individual molecules
        self.molecules = self._construct_molecules()
        # Information about all reactants in a reaction
        self.all_reactants = self._construct_all_reactants()

    def _construct_molecules(self):
        """
        Uses data from database to construct pandas DataFrame with molecule
            data, relative to goal values.

        Note: ALL molecule entries must have values for ALL keys in
            self.parameters.

        :return: molecules (pd.DataFrame)
        """

        mol_ids = [e["mol_id"] for e in self.molecules_raw]

        data = {key: [] for key in self.parameters}
        for mol in self.molecules_raw:
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

        data = {key: list() for key in self.parameters}
        for rxn in self.reactions:
            rcts = rxn["rct_ids"]
            entries = [e for e in self.molecules_raw if e["mol_id"] in rcts]

            for key in self.parameters:
                total = 0
                for entry in entries:
                    if key in self.parameters_min:
                        total += self.goals[key] - entry[key]
                    elif key in self.parameters_max:
                        total += entry[key] - self.goals[key]

                data[key].append(total)

        return pd.DataFrame(data=data, index=rxn_ids)

    def pareto_ranking(self, reactions=None, parameters=None, by_worst=False,
                       group_by_class=False):
        """
        Rank using a Pareto optimization scheme - that is, find all reactions
            that cannot be beat in all categories by another reaction. This
            algorithm effectively creates "classes", rather than a true ranking
            (meaning that while the first and second rank will be
            indistinguishable, the first and the fiftieth may be.

        NOTE: When using by_worst with this algorithm, the result is not
            deterministic. If the multiple choices of reactants are equivalent
            (in terms of their Pareto ranking), then one will randomly be
            chosen.

        :param reactions: list of strings representing molecule IDs. By default
            this is None, meaning that all reactions will be considered.
        :param parameters: list of strings representing parameters to be
            considered in the ranking. By default, this is None, meaning that
            all parameters are considered
        :param by_worst: If True (default False), then reactions should be
            ranked by a single reactant which is ranked lowest, rather than by
            the combination of all reactants.
        :param group_by_class: Boolean value. If True (default False), return
            a dict {class: [items]}, rather than a sorted list.
        :return: ranking, either a sorted list of reaction IDs, or a dict
            containing lists separated by Pareto class
        """

        if parameters is None:
            parameters = ["bp", "log_kow", "mp", "vp", "solubility"]

        if by_worst:
            if reactions is None:
                of_interest = self.reactions
            else:
                of_interest = [r for r in self.reactions if
                               r["rxn_id"] in reactions]

            sortings = list()
            for parameter in parameters:
                sortings.append(self.molecules.sort_values(by=parameter,
                                                           ascending=False))

            mol_rankings = dict()
            for index, rxn in self.molecules.iterrows():
                rankings = list()
                for sorting in sortings:
                    rankings.append(sorting.index.get_loc(index))
                mol_rankings[index] = np.array(rankings)

            rxn_rankings = dict()
            for rxn in of_interest:
                rxn_id = rxn["rxn_id"]
                rct_ids = rxn["rct_ids"]

                choices = list()
                for rct_1 in rct_ids:
                    if all(rct_1 != rct_2
                           and all(np.greater(mol_rankings[rct_1],
                                              mol_rankings[rct_2]))
                           for rct_2 in rct_ids):
                        choices.append(rct_1)
                if len(choices) == 1:
                    rxn_rankings[rxn_id] = mol_rankings[choices[0]]
                elif len(choices) == 0:
                    rxn_rankings[rxn_id] = mol_rankings[choice(rct_ids)]
                else:
                    rxn_rankings[rxn_id] = mol_rankings[choice(choices)]

        else:
            if reactions is None:
                of_interest = self.all_reactants
            else:
                of_interest = self.all_reactants.loc[self.all_reactants.index.isin(reactions)]

            sortings = list()
            for parameter in parameters:
                sortings.append(of_interest.sort_values(by=parameter,
                                                        ascending=False))

            rxn_rankings = dict()
            for index, rxn in of_interest.iterrows():
                rankings = list()
                for sorting in sortings:
                    rankings.append(sorting.index.get_loc(index))
                rxn_rankings[index] = np.array(rankings)

        classes = dict()
        class_num = 0
        while len(rxn_rankings) != 0:
            classes[class_num] = list()
            for ii, ri in rxn_rankings.items():
                superior = list()
                for jj, rj in rxn_rankings.items():
                    if jj != ii and all(np.greater(rj, ri)):
                        superior.append(jj)

                if len(superior) == 0:
                    classes[class_num].append(ii)

            for rxn in classes[class_num]:
                del rxn_rankings[rxn]

            class_num += 1

        if group_by_class:
            ranking = classes
        else:
            ranking = list()
            classes_sorted = sorted(classes.keys())

            for class_index in classes_sorted:
                ranking += classes[class_index]

        return ranking

    def heuristic_ranking(self, reactions=None, parameters=None,
                          by_worst=False):
        """
        Rank using a custom (and somewhat arbitrary) weighting and scaling
            scheme.

        :param reactions: list of strings representing molecule IDs. By default
            this is None, meaning that all reactions will be considered.
        :param parameters: list of strings representing parameters to be
            considered in the ranking. By default, this is None, meaning that
            all parameters are considered
        :param by_worst: If True (default False), then reactions should be
            ranked by a single reactant which is ranked lowest, rather than by
            the combination of all reactants.
        :return: ranking, a sorted list of reaction IDs
        """

        if parameters is None:
            parameters = ["bp", "log_kow", "mp", "vp", "solubility"]

        if reactions is None:
            of_interest = self.reactions
        else:
            of_interest = [r for r in self.reactions if
                           r["rxn_id"] in reactions]

        scores = dict()

        if by_worst:
            for rxn in of_interest:
                rcts = [e for e in self.molecules_raw if e in rxn["rct_ids"]]

                rct_scores = list()
                for rct in rcts:
                    score = 0
                    if "bp" in parameters:
                        if rct["bp"] >= self.goals["bp"]:
                            score += 50
                        else:
                            score += rct["bp"] - self.goals["bp"]
                    if "mp" in parameters:
                        if rct["mp"] <= self.goals["mp"]:
                            score += 25
                        else:
                            score += self.goals["mp"] - rct["mp"]
                    if "vp" in parameters:
                        score -= np.log10(rct["vp"])
                    if "log_kow" in parameters:
                        score += rct["log_kow"]
                    if "solubility" in parameters:
                        score -= np.log10(rct["solubility"] / 1000)

                    rct_scores.append(score)

                scores[rxn] = min(rct_scores)

        else:
            for rxn in of_interest:
                rcts = [e for e in self.molecules_raw if e in rxn["rct_ids"]]

                score = 0
                for rct in rcts:
                    if "bp" in parameters:
                        if rct["bp"] >= self.goals["bp"]:
                            score += 50
                        else:
                            score += rct["bp"] - self.goals["bp"]
                    if "mp" in parameters:
                        if rct["mp"] <= self.goals["mp"]:
                            score += 25
                        else:
                            score += self.goals["mp"] - rct["mp"]
                    if "vp" in parameters:
                        score -= np.log10(rct["vp"])
                    if "log_kow" in parameters:
                        score += rct["log_kow"]
                    if "solubility" in parameters:
                        score -= np.log10(rct["solubility"] / 1000)

                scores[rxn] = score

        ranking = sorted(scores.keys(), key=lambda x: scores[x])
        return ranking


    def tiered_ranking(self, reactions=None, parameters=None, by_worst=False):
        """
        Rank by several parameters in order (for instance, give preference to
            boiling point, and in case of ties sort by melting point).

        :param reactions: list of strings representing molecule IDs. By default
            this is None, meaning that all reactions will be considered.
        :param parameters: list of strings representing parameters to be
            considered in the ranking. These should be in order of preference.
            For instance, ["bp", "mp"] would mean that boiling point would be
            the main factor of consideration, melting point would be the second
            factor, and no other factor will be considered. By default, this is
            None, meaning that all parameters are considered with a default
            ranking
        :param by_worst: If True (default False), then reactions should be
            ranked by a single reactant which is ranked lowest, rather than by
            the combination of all reactants.
        :return: ranking, a sorted list of reaction IDs
        """

        if parameters is None:
            parameters = ["bp", "log_kow", "mp", "vp", "solubility"]

        if by_worst:
            sorted_mols = self.molecules.sort_values(by=parameters,
                                                     ascending=False)
            sorted_ids = sorted_mols.index.tolist()

            if reactions is None:
                of_interest = self.reactions
            else:
                of_interest = [r for r in self.reactions if r["rxn_id"] in reactions]

            highest_index = {}
            for rxn in of_interest:
                max_index = 0
                for rct in rxn["rct_ids"]:
                    index = sorted_ids.index(rct)
                    if index > max_index:
                        max_index = index
                highest_index[rxn] = max_index

            rxn_ids = [r["rxn_id"] for r in of_interest]
            return sorted(rxn_ids, key=lambda r: highest_index[r])

        else:
            if reactions is None:
                of_interest = self.all_reactants
            else:
                of_interest = self.all_reactants.loc[self.all_reactants.index.isin(reactions)]

            sorted = of_interest.sort_values(by=parameters, ascending=False)

            return sorted.index.tolist()
