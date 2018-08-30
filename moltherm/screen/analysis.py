import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import seaborn
from seaborn import stripplot, regplot

__author__ = "Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "July 2018"


class MolThermAnalyzer:
    """
    This class performs analysis based on the data obtained from
    MolThermWorkflow and extracted via MolThermDataProcessor.
    """

    def __init__(self, dataset, setup=True, in_features=None, dep_features=None,
                 func_groups=None, species=None):
        """
        :param dataset: A list of dicts representing all data necessary to
            represent a reaction.
        :param setup: If True (default), then clean the data (put it in a format
            that is appropriate for analysis).
        :param in_features: list of feature/descriptor names for independent
            variables (number_atoms, molecular_weight, etc.).
        :param dep_features: list of feature/descriptor names for dependent
            variables (enthalpy, entropy)
        :param func_groups: list of str representations of functional groups
        :param species: list of str representations of species present in the
            dataset
        """

        if in_features is None:
           self.in_features = ["number_atoms", "molecular_weight", "tpsa",
                               "double_bonds", "triple_bonds", "species",
                               "functional_groups"]
        else:
            self.in_features = in_features

        if dep_features is None:
            self.dep_features = ["enthalpy", "entropy", "t_star"]
        else:
            self.dep_features = dep_features

        if species is None:
            if setup:
                self.species = self._setup_species(dataset)
            else:
                self.species = ["C", "H", "O", "N", "S"]

        if func_groups is None:
            if setup:
                self.func_groups = self._setup_func_groups(dataset)
            else:
                self.func_groups = np.arange(len(dataset["reactants"]["functional_groups"][0]))
        else:
            self.func_groups = func_groups

        if setup:
            self.dataset = self._setup_dataset(dataset)
        else:
            self.dataset = dataset

    @staticmethod
    def _setup_func_groups(dataset):
        """
        Construct a numpy array with labels corresponding to each functional
        group that appears in the dataset.

        :param dataset: list of dicts representing all data necessary to
            represent a reaction.
        :return: np.ndarray
        """

        func_groups = {}

        for datapoint in dataset:
            molecules = [datapoint["product"]] + datapoint["reactants"]
            for molecule in molecules:
                func_groups.update(molecule["functional_groups"])

        return np.array(list(func_groups.keys()))

    @staticmethod
    def _setup_species(dataset):
        """
        Construct a numpy array with str representations of each species that
        appears in the dataset.

        :param dataset: list of dicts representing all data necessary to
            represent a reaction.
        :return: np.ndarray
        """

        species = {}

        for datapoint in dataset:
            molecules = [datapoint["product"]] + datapoint["reactants"]
            for molecule in molecules:
                species.update(molecule["species"])

        return np.array(list(species.keys()))

    def _setup_dataset(self, dataset):
        """
        Alter dataset to make it appropriate for analysis.

        :param dataset: list of dicts, with each dict representing a reaction
        :return: dict containing individual molecule information as well as
            overall reaction information.
        """

        new_dset = {"molecules": {}, "reactions": {}}

        all_molecules = []

        num_reactions = len(dataset)

        for datapoint in dataset:
            if datapoint["product"] not in all_molecules:
                all_molecules.append(datapoint["product"])

            for rct in datapoint["reactants"]:
                if rct not in all_molecules:
                    all_molecules.append(rct)

        new_dset["molecules"]["ids"] = np.array([m["mol_id"] for m in all_molecules])
        new_dset["reactions"]["ids"] = np.array([p["mol_ids"] for p in dataset])
        new_dset["reactions"]["dirs"] = np.array([p["dir_name"] for p in dataset])

        num_molecules = len(all_molecules)

        # Vectorize molecule and reaction features, including thermodynamic
        # properties, surface area, etc.
        for marker in (self.in_features + self.dep_features):
            # Fill dataset for molecules
            if marker == "species":
                new_dset["molecules"][marker] = np.zeros((num_molecules,
                                                         len(self.species)))
            elif marker == "functional_groups":
                new_dset["molecules"][marker] = np.zeros((num_molecules,
                                                          len(self.func_groups)))
            else:
                new_dset["molecules"][marker] = np.zeros(num_molecules)

            for i, mol in enumerate(all_molecules):
                if marker == "enthalpy":
                    new_dset["molecules"][marker][i] = mol["enthalpy"] + mol["energy"]
                elif marker == "t_star":
                    # Turning temperature is not defined for individual molecules
                    continue
                elif marker == "species":
                    for j, spe in enumerate(self.species):
                        if spe in mol["species"].keys():
                            new_dset["molecules"]["species"][i, j] = mol["species"][spe]
                elif marker == "functional_groups":
                    for j, grp in enumerate(self.func_groups):
                        if grp in mol["functional_groups"].keys():
                            new_dset["molecules"]["functional_groups"][i, j] = \
                            mol["functional_groups"][grp]["count"]
                else:
                    new_dset["molecules"][marker][i] = mol[marker]

            # Now fill dataset for reactions

            if marker == "species":
                new_dset["reactions"][marker] = np.zeros((num_reactions,
                                                          len(self.species)))
            elif marker == "functional_groups":
                new_dset["reactions"][marker] = np.zeros((num_reactions,
                                                          len(self.func_groups)))
            else:
                new_dset["reactions"][marker] = np.zeros(num_reactions)

            for i, react in enumerate(dataset):
                thermo = react["thermo"]
                pro = react["product"]
                rcts = react["reactants"]
                if marker in thermo.keys():
                    new_dset["reactions"][marker][i] = thermo[marker]
                elif marker == "species":
                    pro_species = np.zeros(len(self.species))
                    rct_species = np.zeros(len(self.species))

                    for j, spe in enumerate(self.func_groups):
                        if spe in react["product"]["species"].keys():
                            pro_species[j] = react["product"]["species"][spe]

                        for mol in react["reactants"]:
                            if spe in mol["species"].keys():
                                rct_species[j] += mol["species"][spe]

                    new_dset["reactions"]["species"][i] = pro_species - rct_species
                elif marker == "functional_groups":
                    pro_grps = np.zeros(len(self.func_groups))
                    rct_grps = np.zeros(len(self.func_groups))

                    for j, grp in enumerate(self.func_groups):
                        if grp in react["product"]["functional_groups"].keys():
                            pro_grps[j] = react["product"]["functional_groups"][grp]["count"]

                        for mol in react["reactants"]:
                            if grp in mol["functional_groups"].keys():
                                rct_grps[j] += mol["functional_groups"][grp]["count"]

                    new_dset["reactions"]["functional_groups"][i] = pro_grps - rct_grps
                else:
                    pro_data = pro[marker]
                    rct_data = sum(rct[marker] for rct in rcts)
                    new_dset["reactions"][marker][i] = pro_data - rct_data

        return new_dset

    def analyze_features(self, in_features, dep_feature, molecules=False):
        """
        Perform a regression analysis to determine the effect of various
        parameters (molecular weight, for instance) on a particular dependent
        feature (for instance, enthalpy)

        :param in_features: list of strs representing independent variables to
            be analyzed
        :param dep_feature: str representing a dependent variable to be
            analyzed
        :param molecules: If True, perform analysis on an individual molecule
            basis, rather than on a reaction basis
        :return: dict of statistical values
        """

        if molecules:
            in_dataset = {feat: self.dataset["molecules"][feat] for feat in
                          in_features if not (feat == "species" or
                                              feat == "functional_groups")}
            if "species" in in_features:
                species = {s: self.dataset["molecules"]["species"][:, i]
                           for i, s in enumerate(self.species)}
                in_dataset.update(species)
            if "functional_groups" in in_features:
                func_grps = {f: self.dataset["molecules"]["functional_groups"][:, i]
                           for i, f in enumerate(self.func_groups)}
                in_dataset.update(func_grps)

            in_frame = pd.DataFrame(data=in_dataset)
            dep_frame = pd.DataFrame(self.dataset["molecules"][dep_feature],
                                     columns=[dep_feature])

        else:
            in_dataset = {feat: self.dataset["reactions"][feat] for feat in
                          in_features if not (feat == "species" or
                                              feat == "functional_groups")}

            if "species" in in_features:
                species = {s: self.dataset["reactions"]["species"][:, i]
                           for i, s in enumerate(self.species)}
                in_dataset.update(species)
            if "functional_groups" in in_features:
                func_grps = {f: self.dataset["reactions"]["functional_groups"][:, i]
                           for i, f in enumerate(self.func_groups)}
                in_dataset.update(func_grps)

            in_frame = pd.DataFrame(data=in_dataset)
            dep_frame = pd.DataFrame(self.dataset["reactions"][dep_feature],
                                     columns=[dep_feature])

        lm = LinearRegression()
        lm.fit(in_frame, dep_frame)

        score = lm.score(in_frame, dep_frame)
        coefficients = lm.coef_
        intercept = lm.intercept_

        coefficients = {e: coefficients[0][i] for i, e in
                        enumerate(in_frame.columns.values.tolist())}

        return {"r_squared": score,
                "coefficients": coefficients,
                "intercept": intercept}

    def plot_relation(self, in_feature, dep_feature, categorical=False,
                      molecules=False):
        """

        :param in_feature: Independent feature to be evaluated. Must be a member
            of self.in_features
        :param dep_feature: Dependent feature to be evaluated. Must be a member
            of self.dep_features
        :param categorical: If True (default False), use a strip plot rather
            than a regular plot to display data.
        :param molecules: If true, plot on an individual molecule basis, rather
            than on a reaction basis
        :return:
        """

        seaborn.set(style="ticks", color_codes=True)

        if in_feature in self.species:
            col = np.where(self.species == in_feature)[0][0]
            if molecules:
                in_data = self.dataset["molecules"]["species"][:, col]
                dep_data = self.dataset["molecules"][dep_feature]
            else:
                in_data = self.dataset["reactions"]["species"][:, col]
                dep_data = self.dataset["reactions"][dep_feature]
        elif in_feature in self.func_groups:
            col = np.where(self.func_groups == in_feature)[0][0]
            if molecules:
                in_data = self.dataset["molecules"]["functional_groups"][:, col]
                dep_data = self.dataset["molecules"][dep_feature]
            else:
                in_data = self.dataset["reactions"]["functional_groups"][:, col]
                dep_data = self.dataset["reactions"][dep_feature]
        else:
            if molecules:
                in_data = self.dataset["molecules"][in_feature]
                dep_data = self.dataset["molecules"][dep_feature]
            else:
                in_data = self.dataset["reactions"][in_feature]
                dep_data = self.dataset["reactions"][dep_feature]

        dframe = pd.DataFrame(data={in_feature: in_data, dep_feature: dep_data})

        if categorical:
            stripplot(x=in_feature, y=dep_feature, data=dframe)
        else:
            regplot(x=in_feature, y=dep_feature, data=dframe)
        plt.show()

    def analyze_heat_capacity(self, reaction_index=None, reaction_dir=None,
                              max_conc=1.0, solvent_capacity=0.0,
                              plot_dest=None):
        """
        Determine the effective heat capacity of a thermochemical reaction
        system as a function of temperature.

        :param reaction_index: Index in the dataset for the reaction in
            question. Default None
        :param reaction_dir: Name of reaction directory. Default None
        :param max_conc: Maximum concentration of the reactants in the solution.
            Default is 1.0 mol/L
        :param solvent_capacity: Heat capacity of the solvent, which will be
            added to the effective heat capacity of the reaction. Default is
            0.0, which means that there is no solvent
        :param plot_dest: File path pointing to destination for plot. Default is
            None, meaning no plot will be made
        :return: tuple (temp_at_max, maximum) representing peak of heat capacity
            curve
        """

        plt.rc('font', size=14)
        R = 8.31446

        if reaction_index is not None:
            index = reaction_index
        elif reaction_dir is not None:
            index = np.where(self.dataset["reactions"]["dirs"] == reaction_dir)[0][0]
        else:
            raise ValueError("User must supply either a reaction index or "
                             "directory name.")

        delta_h = self.dataset["reactions"]["enthalpy"][index]
        delta_s = self.dataset["reactions"]["entropy"][index]

        fig, ax = plt.subplots()
        ax.set(xlabel='T (°C)', ylabel='C (J/g·K)', title='Specific heat')

        conc_i = 0.99 * max_conc
        conc_f = 0.01 * max_conc

        K_eq_i = conc_i / (max_conc - conc_i) ** 2
        K_eq_f = conc_f / (max_conc - conc_f) ** 2

        temp_i = delta_s / (delta_s - R * np.log(K_eq_i))
        temp_f = delta_s / (delta_s - R * np.log(K_eq_f))
        T = np.linspace(temp_i, temp_f, 1000)
        T_in_deg_C = T - 273.15

        K_eq = np.exp(delta_s / R - delta_h / (R * T))
        Cp = delta_h ** 2 / (R * (T ** 2) * K_eq) * (((2 * max_conc + 1 / K_eq) / (2 * np.sqrt((2 * max_conc + 1 / K_eq) ** 2 - 4 * max_conc ** 2))) - 0.5) + solvent_capacity

        max_cp = np.argmax(Cp)
        maximum = max(Cp) / 1000
        temp_at_max = T_in_deg_C[max_cp]
        print('The highest specific heat is ' + str(maximum) + ' J/g·K at ' + str(temp_at_max) + ' C°.')

        if plot_dest is not None:
            ax.plot(T_in_deg_C, Cp / 1000)
            fig.savefig(plot_dest, dpi=300)

        return temp_at_max, maximum

    def analyze_energy_density(self, reaction_index=None, reaction_dir=None,
                               max_conc=1.0, solvent_capacity=0.0,
                               plot_dest=None):
        """
        Determine the effective energy density of a thermochemical reaction
        system as a function of temperature.

        :param reaction_index: Index in the dataset for the reaction in
            question. Default None
        :param reaction_dir: Name of reaction directory. Default None
        :param max_conc: Maximum concentration of the reactants in the solution.
            Default is 1.0 mol/L
        :param solvent_capacity: Heat capacity of the solvent, which will be
            added to the effective heat capacity of the reaction. Default is
            0.0, which means that there is no solvent.
        :param plot_dest: File path pointing to destination for plot. Default is
            None, meaning no plot will be made
        :return: Maximum energy density over the range of interest
        """

        plt.rc('font', size=14)
        R = 8.31446

        if reaction_index is not None:
            index = reaction_index
        elif reaction_dir is not None:
            index = np.where(self.dataset["reactions"]["dirs"] == reaction_dir)[0][0]
        else:
            raise ValueError("User must supply either a reaction index or "
                             "directory name.")

        delta_h = self.dataset["reactions"]["enthalpy"][index]
        delta_s = self.dataset["reactions"]["entropy"][index]

        fig, ax = plt.subplots()
        ax.set(xlabel='T (°C)', ylabel='Energy density (MJ/kg)',
                 title='Energy density')

        conc_i = 0.99 * max_conc
        conc_f = 0.01 * max_conc

        K_eq_i = conc_i / (max_conc - conc_i) ** 2
        K_eq_f = conc_f / (max_conc - conc_f) ** 2

        temp_i = delta_s / (delta_s - R * np.log(K_eq_i))
        temp_f = delta_s / (delta_s - R * np.log(K_eq_f))
        T = np.linspace(temp_i, temp_f, 1000)
        T_in_deg_C = T - 273.15

        K_eq = np.exp(delta_s / R - delta_h / (R * T))
        H = delta_h * (((2 * max_conc + 1 / K_eq) - np.sqrt((2 * max_conc + 1 / K_eq) ** 2 - 4 * max_conc ** 2)) / 2 - conc_i) + solvent_capacity * (T - temp_i)

        if plot_dest is not None:
            ax.plot(T_in_deg_C, H / 1000000)
            fig.savefig(plot_dest, dpi=300)

        return max(H)
