    from os import listdir
from os.path import join, isfile, isdir, abspath
from collections import Counter
from itertools import chain

import networkx as nx

from pymatgen.core.structure import Molecule
from pymatgen.analysis.functional_groups import FunctionalGroupExtractor
from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.io.qchem.outputs import QCOutput

from atomate.qchem.database import QChemCalcDb
from atomate.qchem.drones import QChemDrone

from moltherm.compute.utils import (associate_qchem_to_mol,
                                    calculate_solubility,
                                    extract_id, get_molecule)

__author__ = "Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "August 2018"


class MolThermDataProcessor:
    """
    This class can be used to extract data from MolThermWorkflow workflows,
    including extracting thermo data from calculations and generating predicted
    boiling and melting points.
    """

    def __init__(self, base_dir, mol_dir="molecules", rxn_dir="reactions",
                 db_file="db.json", mol_coll="molecules", rxn_coll="reactions",
                 thermo_coll="thermo", sol_coll="solubility",
                 epi_coll="episuite"):
        """
        :param base_dir: Directory where all data should be stored.
        :param mol_dir: Subdirectory where molecule files and calculations are
            stored. Default is "molecules".
        :param rxn_dir: Subdirectory where reaction metadata files (and
            potentially calculations) are stored. Default is "reactions".
        :param db_file: Path to database config file.
        :param mol_coll: Database collection in which to store molecule data.
            Default is "molecules".
        :param rxn_coll: Database collection in which reaction metadata is
            stored. Default is "reaxys".
        :param thermo_coll: Database collection in which to store thermo data.
            Default is "thermo".
        :param sol_coll: Database collection in which to store solubility data.
            Default is "solubility".
        :param epi_coll: Database collection from which to extract data
            regarding physical parameters (for instance, vapor pressure and
            boiling point). Default is "episuite".
        """

        self.base_dir = base_dir
        self.mol_dir = join(self.base_dir, mol_dir)
        self.rxn_dir = join(self.base_dir, rxn_dir)
        self.db_file = db_file

        try:
            self.db = QChemCalcDb.from_db_file(self.db_file)
        except:
            self.db = None

        if self.db is None:
            self.mol_coll = mol_coll
            self.rxn_coll = rxn_coll
            self.thermo_coll = thermo_coll
            self.sol_coll = sol_coll
            self.epi_coll = epi_coll
        else:
            self.mol_coll = self.db.db[mol_coll]
            self.rxn_coll = self.db.db[rxn_coll]
            self.thermo_coll = self.db.db[thermo_coll]
            self.sol_coll = self.db.db[sol_coll]
            self.epi_coll = self.db.db[epi_coll]

    def extract_reaction_thermo_files(self, rxn_id, files_mol_prefix=False,
                                      runs_pattern=None):
        """
        Naively scrape thermo data from QChem output files.

        :param rxn_id: str representing a reaction ID
        :param files_mol_prefix: QChemDrone.assimilate() requires a template for
            QChem input and output file names. If True, the prefix for these
            files will be the mol_id. If False (default), the prefix will be
            "mol".
        :param runs_pattern: QChemDrone.assimilate() requires a template for
            what calculations have completed. This should be represented by a
            link. Default None

        :return: dict {prop: value}, where properties are enthalpy, entropy,
            etc.
        """

        base_path = self.mol_dir

        rxn_entry = self.rxn_coll.find_one({"rxn_id": rxn_id})

        if rxn_entry is None:
            raise ValueError("Invalid rxn_id. Reaction not found in "
                             "collection {}".format(self.rxn_coll))

        rct_ids = rxn_entry["rct_ids"]
        pro_ids = rxn_entry["pro_ids"]

        rct_thermo = {"enthalpy": 0, "entropy": 0, "energy": 0}
        pro_thermo = {"enthalpy": 0, "entropy": 0, "energy": 0}

        drone = QChemDrone(runs=runs_pattern)

        for mol_id in rct_ids:
            if files_mol_prefix:
                calc_doc = drone.assimilate(path=join(base_path, mol_id),
                                            input_file="{}.qin".format(mol_id),
                                            output_file="{}.qout".format(mol_id),
                                            multirun=False)
            else:
                calc_doc = drone.assimilate(path=join(base_path, mol_id),
                                            input_file="mol.qin",
                                            output_file="mol.qout",
                                            multirun=False)

            rct_thermo["entropy"] += calc_doc["output"]["entropy"]
            rct_thermo["enthalpy"] += calc_doc["output"]["enthalpy"]
            rct_thermo["energy"] += calc_doc["output"]["final_energy"]

        for mol_id in pro_ids:
            if files_mol_prefix:
                calc_doc = drone.assimilate(path=join(base_path, mol_id),
                                            input_file="{}.qin".format(mol_id),
                                            output_file="{}.qout".format(mol_id),
                                            multirun=False)
            else:
                calc_doc = drone.assimilate(path=join(base_path, mol_id),
                                            input_file="mol.qin",
                                            output_file="mol.qout",
                                            multirun=False)

            pro_thermo["entropy"] += calc_doc["output"]["entropy"]
            pro_thermo["enthalpy"] += calc_doc["output"]["enthalpy"]
            pro_thermo["energy"] += calc_doc["output"]["final_energy"]

        thermo_data = {}

        # Generate totals as ∆H = H_pro - H_rct, ∆S = S_pro - S_rct
        # Also ensures that units are appropriate (Joules/mol,
        # rather than cal/mol or kcal/mol, or hartree for energy)
        energy = (pro_thermo["energy"] - rct_thermo["energy"]) * 627.509
        enthalpy = (pro_thermo["enthalpy"] - rct_thermo["enthalpy"])
        thermo_data["enthalpy"] = (energy + enthalpy) * 1000 * 4.184
        thermo_data["entropy"] = (pro_thermo["entropy"] - rct_thermo["entropy"]) * 4.184
        try:
            thermo_data["t_critical"] = thermo_data["enthalpy"] / thermo_data["entropy"]
        except ZeroDivisionError:
            thermo_data["t_critical"] = None

        result = {"thermo": thermo_data,
                  "rxn_id": rxn_id,
                  "reactant_ids": rct_ids,
                  "product_ids": pro_ids}

        return result

    def extract_reaction_thermo_db(self, rxn_id, opt=None, freq=None, sp=None):
        """
        Gathers all relevant reaction parameters, including references to
        each job performed.

        :param rxn_id: str representing a reaction ID
        :param opt: dict containing information about the optimization jobs. By
            default, this is None, and that information will be obtained by
            querying the self.db.tasks collection.
        :param freq: dict containing information about the frequency jobs. By
            default, this is None, and that information will be obtained by
            querying the self.db.tasks collection.
        :param sp: dict containing information about the single-point jobs. By
            default, this is None, and that information will be obtained by
            querying the self.db.tasks collection.

        :return: dict
        """

        if self.db is None:
            raise RuntimeError("Could not connect to database. Check db_file"
                               "and try again later.")

        # To extract enthalpy and entropy from calculation results
        def get_thermo(job, sp=False):
            enthalpy = job["output"]["enthalpy"]
            entropy = job["output"]["entropy"]
            if sp:
                try:
                    energy = job["output"]["final_energy_sp"]
                except KeyError:
                    energy = job["output"]["final_energy"]
            else:
                energy = job["output"]["final_energy"]

            if enthalpy is None:
                enthalpy = 0.0
            if entropy is None:
                entropy = 0.0
            if energy is None:
                energy = 0.0

            return {"enthalpy": enthalpy,
                    "entropy": entropy,
                    "energy": energy}

        rxn_entry = self.rxn_coll.find_one({"rxn_id": rxn_id})

        # Sort files for if they are reactants or products
        reactants = [self.mol_coll.find_one({"mol_id": m}) for m in rxn_entry["rct_ids"]]
        products = [self.mol_coll.find_one({"mol_id": m}) for m in rxn_entry["pro_ids"]]

        records = reactants + products

        for i, record in enumerate(records):
            if opt is None:
                for calc in record["calcs_reversed"]:
                    if calc["input"]["rem"]["job_type"] == "opt" or \
                            calc["input"]["rem"]["job_type"] == "optimization":
                        method = calc["input"]["rem"]["method"]
                        basis = calc["input"]["rem"]["basis"]
                        solvent_method = calc["input"]["rem"].get(
                            "solvent_method", None)
                        if solvent_method == "smd":
                            if calc["input"]["smx"] is None:
                                solvent = None
                            else:
                                solvent = calc["input"]["smx"]["solvent"]
                        elif solvent_method == "pcm":
                            solvent = calc["input"]["solvent"]
                        else:
                            solvent = None

                        opt = {"method": method,
                               "basis": basis,
                               "solvent_method": solvent_method,
                               "solvent": solvent}
                        break
            if freq is None:
                for calc in record["calcs_reversed"]:
                    if calc["input"]["rem"]["job_type"] == "freq" or \
                            calc["input"]["rem"]["job_type"] == "frequency":
                        method = calc["input"]["rem"]["method"]
                        basis = calc["input"]["rem"]["basis"]
                        solvent_method = calc["input"]["rem"].get(
                            "solvent_method", None)
                        if solvent_method == "smd":
                            if calc["input"]["smx"] is None:
                                solvent = None
                            else:
                                solvent = calc["input"]["smx"]["solvent"]
                        elif solvent_method == "pcm":
                            solvent = calc["input"]["solvent"]
                        else:
                            solvent = None

                        freq = {"method": method,
                                "basis": basis,
                                "solvent_method": solvent_method,
                                "solvent": solvent}
                        break
            if sp is None:
                for calc in record["calcs_reversed"]:
                    if calc["input"]["rem"]["job_type"] == "sp":
                        method = calc["input"]["rem"]["method"]
                        basis = calc["input"]["rem"]["basis"]
                        solvent_method = calc["input"]["rem"].get(
                            "solvent_method", None)
                        if solvent_method == "smd":
                            if calc["input"]["smx"] is None:
                                solvent = None
                            else:
                                solvent = calc["input"]["smx"]["solvent"]
                        elif solvent_method == "pcm":
                            solvent = calc["input"]["solvent"]
                        else:
                            solvent = None

                        sp = {"method": method,
                              "basis": basis,
                              "solvent_method": solvent_method,
                              "solvent": solvent}
                        break

        # Get thermo data
        if sp is None:
            rct_thermo = [get_thermo(r) for r in reactants]
            pro_thermo = [get_thermo(p) for p in products]
        else:
            rct_thermo = [get_thermo(r, sp=True) for r in reactants]
            pro_thermo = [get_thermo(p, sp=True) for p in products]

        # Compile reaction thermo from reactant and product thermos
        delta_e = sum(p["energy"] for p in pro_thermo) - sum(r["energy"] for r in rct_thermo)
        delta_e *= 627.509
        delta_h = sum(p["enthalpy"] for p in pro_thermo) - sum(r["enthalpy"] for r in rct_thermo) + delta_e
        delta_h *= 1000 * 4.184
        delta_s = sum(p["entropy"] for p in pro_thermo) - sum(r["entropy"] for r in rct_thermo)
        delta_s *= 4.184
        thermo = {
            "enthalpy": delta_h,
            "entropy": delta_s
        }

        try:
            thermo["t_star"] = delta_h / delta_s
        except ZeroDivisionError:
            thermo["t_star"] = 0

        result = {"rxn_id": rxn_id,
                  "opt": opt,
                  "freq": freq,
                  "sp": sp,
                  "rct_ids": rxn_entry["rct_ids"],
                  "pro_ids": rxn_entry["pro_ids"],
                  "thermo": thermo}

        return result

    def record_molecule_data_db(self, mol_id, input_file="mol.qin",
                                output_file="mol.qout", runs_pattern=None):
        """
        Compile calculation information for a single molecule and record it in
        the molecules collection.

        :param mol_id: Unique identifier for molecule (str)
        :param input_file: Basic format for input files. The Drone which
            compiles the molecule information will use this to pattern-match
            files. Default is "mol.qin"
        :param output_file: Basic format for output files. The Drone which
            compiles the molecule information will use this to pattern-match
            files. Default is "mol.qout"
        :param runs_pattern: QChem drone assimilation requires a template for
            what calculations have completed. This should be represented by a
            link. Default None
        :return:
        """

        calc_dir = join(self.mol_dir, mol_id)

        drone = QChemDrone(runs=runs_pattern)

        task_doc = drone.assimilate(
            path=calc_dir,
            input_file=input_file,
            output_file=output_file,
            multirun=False)

        task_doc["mol_id"] = mol_id

        if "enthalpy" not in task_doc["output"] or "enthalpy" not in task_doc or "frequencies" not in task_doc:
            if task_doc["calcs_reversed"][1]["input"]["rem"]["job_type"] in ["freq", "frequency"]:
                freq = task_doc["calcs_reversed"][1]
                task_doc["output"]["enthalpy"] = freq["total_enthalpy"]
                task_doc["output"]["entropy"] = freq["total_entropy"]
                task_doc["output"]["frequencies"] = freq["frequencies"]


        if self.db is None:
            raise RuntimeError("Cannot record data to db without valid database"
                               " connection!")

        self.mol_coll.insert_one(task_doc)

    def record_reaction_data_db(self, rxn_id, use_db=False, opt=None,
                                freq=None, sp=None):
        """
        Record thermo data in thermo collection.

        :param rxn_id: str representing a reaction ID
        :param use_db: If set to True (default False), use
            extract_reaction_data to gather data
        :param opt: dict containing information about the optimization jobs. By
            default, this is None, and that information will be obtained by
            querying the self.db.tasks collection.
        :param freq: dict containing information about the frequency jobs. By
            default, this is None, and that information will be obtained by
            querying the self.db.tasks collection.
        :param sp: dict containing information about the single-point jobs. By
            default, this is None, and that information will be obtained by
            querying the self.db.tasks collection.

        :return:
        """

        if self.db is None:
            raise RuntimeError("Could not connect to database. Check db_file"
                               "and try again later.")

        if use_db:
            self.thermo_coll.insert_one(self.extract_reaction_thermo_db(rxn_id,
                                                                        opt=opt,
                                                                        freq=freq,
                                                                        sp=sp))
        else:
            self.thermo_coll.insert_one(self.extract_reaction_thermo_files(rxn_id))

    def populate_collections(self, thermo=False, solubility=False,
                             overwrite=False, files_mol_prefix=False,
                             runs_pattern=None, solvent=None,
                             vacuum_directory=None):
        """
        Mass insert into db collections using above data extraction and
        recording methods.

        :param thermo: If True (default False), store thermodynamics information
        :param solubility: If True (default False), store solubility information
        :param overwrite: If True (default False), overwrite molecules that are
            already included in the db.
         :param files_mol_prefix: QChemDrone.assimilate() requires a template
            for QChem input and output file names. If True, the prefix for these
            files will be the mol_id. If False (default), the prefix will be
            "mol".
        :param runs_pattern: QChemDrone.assimilate() requires a template for
            what calculations have completed. This should be represented by a
            link. Default None
        :param solvent: Solvent of interest. Only needed if solubility is True.
        :param vacuum_directory: FOR SOLVATION ONLY: Path to search for energy
            data in vacuum, to be compared against the solvated data in the
            molecule collection.

        :return:
        """

        if self.db is None:
            raise RuntimeError("Cannot connect to database. Check configuration"
                               " file and try again.")

        completed_mols = self.get_completed_molecules(extra=True,
                                                      files_mol_prefix=files_mol_prefix,
                                                      runs_pattern=runs_pattern)
        mols_in_db = [mol for mol in self.mol_coll.find({}, {"_id": 0})]

        for mol_id, path in completed_mols:
            if mol_id not in [m["mol_id"] for m in mols_in_db]:
                if files_mol_prefix:
                    self.record_molecule_data_db(mol_id,
                                                 "{}.qin".format(mol_id),
                                                 "{}.qout".format(mol_id),
                                                 runs_pattern=runs_pattern)
                else:
                    self.record_molecule_data_db(mol_id, "mol.qin", "mol.qout",
                                                 runs_pattern=runs_pattern)
            elif overwrite:
                drone = QChemDrone(runs=runs_pattern)

                base_path = self.mol_dir

                if files_mol_prefix:
                    task_doc = drone.assimilate(path=join(base_path, mol_id),
                                                input_file="{}.qin".format(
                                                    mol_id),
                                                output_file="{}.qout".format(
                                                    mol_id),
                                                multirun=False)
                else:
                    task_doc = drone.assimilate(path=join(base_path, mol_id),
                                                input_file="mol.qin",
                                                output_file="mol.qout",
                                                multirun=False)

                task_doc["mol_id"] = mol_id

                self.mol_coll.update_one({"mol_id": mol_id},
                                    {"$set": task_doc})

        if thermo:
            completed_rxns = self.get_completed_reactions()
            rxns_in_db = [rxn for rxn in self.thermo_coll.find()]

            for rxn in completed_rxns:
                if rxn["rxn_id"] not in [r["rxn_id"] for r in rxns_in_db]:
                    try:
                        self.record_reaction_data_db(rxn["rxn_id"], use_db=True)
                    except:
                        print("Failed to add reaction: {}".format(rxn.get("rxn_id", rxn)))

                elif overwrite:
                    self.thermo_coll.update_one({"rxn_id": rxn["rxn_id"]},
                                                {"$set": self.extract_reaction_thermo_db(rxn["rxn_id"])})
        if solubility:
            # Update list, presuming that there were some additions
            mols_in_db = [mol for mol in self.mol_coll.find({}, {"_id": 0})]
            mols_in_sol_coll = [mol for mol in self.sol_coll.find({}, {"_id": 0})]

            for mol in mols_in_db:
                mol_id = mol["mol_id"]

                try:
                    # Assume that vacuum data uses mol.qin format
                    # This should probably be changed in the future
                    data = self.get_solubility_data(mol_id, solvent,
                                                    vacuum_directory=vacuum_directory)
                    entry = {"mol_id": data["mol_id"],
                             "solubilities":
                                 {data["solvent"]: data["solubility"]}}

                    if mol_id not in [e["mol_id"] for e in mols_in_sol_coll]:
                        self.sol_coll.insert_one(entry)
                    elif overwrite:
                        new = self.sol_coll.find_one({"mol_id": mol_id})
                        new["solubilities"][data["solvent"]] = data["solubility"]
                        self.sol_coll.update_one({"mol_id": mol_id},
                                                 {"$set": new})

                except RuntimeError:
                    print("Could not process {}".format(mol_id))
                    continue

    def update_molecules(self, files_mol_prefix=False):
        """
        Update molecules collection with data from the subdirectories in
        self.mol_dir.

        :param files_mol_prefix: QChemDrone.assimilate() requires a template
            for QChem input and output file names. If True, the prefix for these
            files will be the mol_id. If False (default), the prefix will be
            "mol".

        :return:
        """

        if self.db is None:
            raise RuntimeError("Cannot access database. Check configuration"
                               " settings and try again.")

        dirs = [d for d in listdir(self.mol_dir) if
                isdir(join(self.mol_dir, d)) and not d.startswith("block")]

        drone = QChemDrone()

        for d in dirs:
            calc_dir = join(self.mol_dir, d)

            if files_mol_prefix:
                new_calc = drone.assimilate(path=calc_dir,
                                            input_file="{}.qin".format(d),
                                            output_file="{}.qout".format(d),
                                            multirun=False)
            else:
                new_calc = drone.assimilate(path=calc_dir,
                                            input_file="mol.qin",
                                            output_file="mol.qout",
                                            multirun=False)

            old_calc = self.mol_coll.find_one({"mol_id": d})

            if old_calc is not None:
                calcs_reversed = new_calc["calcs_reversed"] + old_calc["calcs_reversed"]

                output = new_calc["output"]

                walltime = new_calc["walltime"] + old_calc["walltime"]
                cputime = new_calc["cputime"] + old_calc["cputime"]

                self.mol_coll.update_one({"mol_id": d},
                                         {"$set": {"calcs_reversed": calcs_reversed,
                                                   "output": output,
                                                   "walltime": walltime,
                                                   "cputime": cputime}})

    def update_thermo(self):
        """
        Update thermo collection with data from the subdirectories in
        self.base_dir.

        :return:
        """

        if self.db is None:
            raise RuntimeError("Cannot access database. Check configuration"
                               " settings and try again.")

        to_update = []

        all_rxns = [r for r in self.thermo_coll.find({}, {"_id": 0})]

        for rxn in all_rxns:
            rxn_entry = self.rxn_coll.find_one({"rxn_id": rxn["rxn_id"]})

            pro_thermo = [self.get_molecule_data(m) for m in rxn_entry["pro_ids"]]
            rct_thermo = [self.get_molecule_data(m) for m in rxn_entry["rct_ids"]]

            # Compile reaction thermo from reactant and product thermos
            delta_e = sum(p["energy"] for p in pro_thermo) - sum(
                r["energy"] for r in rct_thermo)
            delta_h = sum(p["enthalpy"] for p in pro_thermo) - sum(
                r["enthalpy"] for r in rct_thermo) + delta_e
            delta_s = sum(p["entropy"] for p in pro_thermo) - sum(
                r["entropy"] for r in rct_thermo)
            thermo = {
                "enthalpy": delta_h,
                "entropy": delta_s
            }

            try:
                thermo["t_star"] = delta_h / delta_s
            except ZeroDivisionError:
                thermo["t_star"] = 0

            to_update.append((rxn["rxn_id"], thermo))

        for update in to_update:
            self.thermo_coll.update_one({"rxn_id": update[0]},
                                        {"$set": {"thermo": update[1]}})

    def update_solubility(self, solvent="water", vacuum_directory=None,
                          files_mol_prefix=False):
        """
        Update solubility collection with data from the subdirectories in
        self.base_dir.

        :param solvent: Solvent of interest. Default is "water"
        :param vacuum_directory: Path to search for energy data in vacuum, to be
            compared against the solvated data in the molecule collection.
        :param files_mol_prefix: QChemDrone.assimilate() requires a template
            for QChem input and output file names. If True, the prefix for these
            files will be the mol_id. If False (default), the prefix will be
            "mol".

        :return:
        """

        mols_in_db = [mol for mol in self.mol_coll.find({}, {"_id": 0})]

        for mol in mols_in_db:
            mol_id = mol["mol_id"]

            try:
                data = self.get_solubility_data(mol_id, solvent,
                                                vacuum_directory=vacuum_directory,
                                                files_mol_prefix=files_mol_prefix)

                new = self.sol_coll.find_one({"mol_id": mol_id})
                new["solubilities"][data["solvent"]] = data["solubility"]
                self.sol_coll.update_one({"mol_id": mol_id},
                                         {"$set": new})

            except RuntimeError:
                continue

    def find_reactions_common_reactant(self, mol_id):
        """
        Queries reactions to identify all those that share a common reactant or
        product.

        :param mol_id: ID (str) of reactant molecule to be searched.
        :return: list of dicts representing reaction metadata.
        """

        # Should return all reactions where mol_id is in rct_ids or pro_ids
        res_rct = self.rxn_coll.find({"rct_ids": mol_id})
        res_pro = self.rxn_coll.find({"pro_ids": mol_id})

        results = set()

        if res_rct is not None:
            for res in res_rct:
                results.add(res["rxn_id"])
        if res_pro is not None:
            for res in res_pro:
                results.add(res["rxn_id"])

        return results

    def map_reactants_to_reactions(self):
        """
        Construct a dict showing which directories share each reactant.

        This is useful for analysis of common reactants, and to identify the
        "source" of a given reactant (in which directory the calculation
        actually took place).

        :return: dict mapping molecule IDs to sets of reaction IDs.
        """

        molecules = [d for d in listdir(self.mol_dir)
                     if isdir(join(self.mol_dir, d))]

        mapping = {}

        for mol_id in molecules:
            res = self.find_reactions_common_reactant(mol_id)
            if len(res) != 0:
                mapping[mol_id] = res

        return mapping

    def get_completed_molecules(self, mols=None, extra=False,
                                files_mol_prefix=False, runs_pattern=None):
        """
        Returns a list of molecules with completed opt, freq, and sp output
        files.

        :param mols: List of molecule ids (equivalently, molecule directories)
            to search for completion.
        :params extra: If True, include directory of completed reaction and name
            of molfile along with mol_id
        :param files_mol_prefix: QChemDrone.assimilate() requires a template
            for QChem input and output file names. If True, the prefix for these
            files will be the mol_id. If False (default), the prefix will be
            "mol".
        :param runs_pattern: QChem drone assimilation requires a template for
            what calculations have completed. This should be represented by a
            link. Default None
        :return: set of completed molecules
        """

        completed = set()

        all_mols = [m for m in listdir(self.mol_dir)
                    if isdir(join(self.mol_dir, m))]

        if mols is not None:
            all_mols = [m for m in all_mols if m in mols]

        drone = QChemDrone(runs=runs_pattern)

        for m in all_mols:
            path = join(self.mol_dir, m)

            try:
                if files_mol_prefix:
                    result = drone.assimilate(path=path,
                                              input_file="{}.qin".format(m),
                                              output_file="{}.qout".format(m),
                                              multirun=False)
                else:
                    result = drone.assimilate(path=path,
                                              input_file="mol.qin",
                                              output_file="mol.qout",
                                              multirun=False)

                sp = False
                if runs_pattern is not None:
                    if "sp" in runs_pattern:
                        sp = True

                calcs = result["calcs_reversed"]
                parts = {"opt": False, "freq": False, "sp": False}

                freqs = None

                for calc in calcs:
                    if freqs is None:
                        freqs = calc.get("frequencies", None)

                    task_type = calc["task"]["type"]
                    part = bool(calc["completion"])
                    if "sp" in task_type and part:
                        parts["sp"] = True
                    elif "freq" in task_type and part:
                        parts["freq"] = True
                    elif "opt" in task_type and part:
                        parts["opt"] = True

                if sp:
                    if parts["opt"] and parts["freq"] and parts["sp"]:
                        if any([x < 0 for x in freqs]):
                            complete = False
                        else:
                            complete = True
                    else:
                        complete = False
                else:
                    if parts["opt"] and parts["freq"]:
                        if any([x < 0 for x in freqs]):
                            complete = False
                        else:
                            complete = True
                    else:
                        complete = False
            except (ValueError, IndexError):
                complete = False

            if complete:
                if extra:
                    completed.add((m, path))
                else:
                    completed.add(m)

        return completed

    def get_completed_reactions(self):
        """
        Returns a list of directories (reactions) where all molecules are
        completed.

        :return: list of directories with complete information.
        """

        if self.db is None:
            raise RuntimeError("Could not connect to database. Check db_file"
                               "and try again later.")

        completed_molecules = [x["mol_id"] for x in self.mol_coll.find()]
        completed_reactions = list()

        all_reactions = [r for r in self.rxn_coll.find({}, {"_id": 0})]

        for r in all_reactions:
            all_ids = r["pro_ids"] + r["rct_ids"]

            are_completed = [True if m in completed_molecules else False
                             for m in all_ids]

            if all(are_completed):
                completed_reactions.append(r)

        return completed_reactions

    def get_molecule_data(self, mol_id):
        """
        Compile all useful molecular data for analysis, including molecule size
        (number of atoms), molecular weight, enthalpy, entropy, and functional
        groups.

        NOTE: This function automatically converts energy, enthalpy, and entropy
        into SI units (J/mol and J/mol*K)

        :param mol_id: Unique ID associated with the molecule.
        :return: dict of relevant molecule data.
        """

        mol_data = {"mol_id": mol_id}

        if self.db is None:
            raise RuntimeError("Cannot query database; connection is invalid."
                               " Try to connect again.")

        mol_entry = self.mol_coll.find_one({"mol_id": mol_id})

        mol_data["enthalpy"] = mol_entry["output"]["enthalpy"] * 4.184 * 1000
        mol_data["entropy"] = mol_entry["output"]["entropy"] * 4.184

        try:
            final_energy = mol_entry["output"]["final_energy_sp"]
        except KeyError:
            final_energy = mol_entry["output"]["final_energy"]

        mol_data["energy"] = final_energy * 627.509 * 4.184 * 1000

        mol_dict = mol_entry["output"]["optimized_molecule"]
        mol_data["molecule"] = Molecule.from_dict(mol_dict)

        adaptor = BabelMolAdaptor(mol_data["molecule"])
        pbmol = adaptor.pybel_mol

        mol_data["number_atoms"] = len(mol_data["molecule"])
        mol_data["molecular_weight"] = pbmol.molwt
        mol_data["tpsa"] = pbmol.calcdesc()["TPSA"]

        extractor = FunctionalGroupExtractor(mol_data["molecule"])
        molgraph = extractor.molgraph
        func_grps = extractor.get_all_functional_groups()

        mol_data["functional_groups"] = extractor.categorize_functional_groups(func_grps)

        weights = nx.get_edge_attributes(molgraph.graph, "weight")
        bonds_checked = set()
        double_bonds = 0
        triple_bonds = 0
        for bond, weight in weights.items():
            # Remove index from multidigraph
            bond = (bond[0], bond[1])
            if int(weight) == 2 and bond not in bonds_checked:
                double_bonds += 1
            elif int(weight) == 3 and bond not in bonds_checked:
                triple_bonds += 1
            bonds_checked.add(bond)

        mol_data["double_bonds"] = double_bonds
        mol_data["triple_bonds"] = triple_bonds

        species = [str(s.specie) for s in mol_data["molecule"].sites]
        mol_data["species"] = dict(Counter(species))

        return mol_data

    def get_reaction_data(self, rxn_id):
        """
        Compile all useful data for a set of molecules associated with a
        particular reaction. This data will be compiled on a reaction basis
        (difference between reactants and products) as well as an individual
        molecule basis.

        :param rxn_id: str representing a reaction ID
        :return: dict of relevant reaction data.
        """

        reaction_data = {}

        rxn_entry = self.rxn_coll.find_one({"rxn_id": rxn_id})

        reaction_data["rxn_id"] = rxn_id
        reaction_data["mol_ids"] = rxn_entry["pro_ids"] + rxn_entry["rct_ids"]
        reaction_data["product"] = self.get_molecule_data(rxn_entry["pro_ids"][0])
        reaction_data["reactants"] = [self.get_molecule_data(m)
                                      for m in rxn_entry["rct_ids"]]

        reaction_data["thermo"] = {}
        pro_h = reaction_data["product"]["enthalpy"] + reaction_data["product"]["energy"]
        rct_h = sum(r["enthalpy"] + r["energy"]
                    for r in reaction_data["reactants"])
        reaction_data["thermo"]["enthalpy"] = pro_h - rct_h

        pro_s = reaction_data["product"]["entropy"]
        rct_s = sum(r["entropy"] for r in reaction_data["reactants"])
        reaction_data["thermo"]["entropy"] = pro_s - rct_s

        try:
            reaction_data["thermo"]["t_star"] = reaction_data["thermo"]["enthalpy"] / reaction_data["thermo"]["entropy"]
        except ZeroDivisionError:
            reaction_data["thermo"]["t_star"] = 0

        return reaction_data

    def get_solubility_data(self, mol_id, solvent, molecule_collection=None,
                            vacuum_directory=None, files_mol_prefix=False):
        """
        Extract solubility data from molecule data.

        :param mol_id: Unique ID associated with the molecule.
        :param solvent: Name of the solvent
        :param molecule_collection: Collection to search for solubility data.
            By default, this will be None, and self.mol_coll will be used.
        :param vacuum_directory: Path to search for energy data in vacuum, to be
            compared against the solvated data in the molecule collection.
        :param files_mol_prefix: QChemDrone.assimilate() requires a template
            for QChem input and output file names. If True, the prefix for these
            files will be the mol_id. If False (default), the prefix will be
            "mol".

        :return: sol_data
        """

        if molecule_collection is None:
            coll = self.mol_coll
        else:
            coll = self.db.db[molecule_collection]

        entry = coll.find_one({"mol_id": mol_id})

        try:
            vp = self.epi_coll.find_one({"mol_id": mol_id})["vp"]
        except TypeError:
            raise RuntimeError("{} is not in epi_coll!".format(mol_id))

        sol_data = dict()
        sol_data["mol_id"] = mol_id
        sol_data["solvent"] = solvent

        g_liq = None
        g_vac = None

        for calc in entry["calcs_reversed"]:
            job_type = calc["input"]["rem"]["job_type"]

            if job_type == "sp" and g_liq is None:
                try:
                    g_liq = calc["solvent_data"]["smd6"]
                except KeyError:
                    continue
            elif job_type in ["opt", "optimization"] and g_vac is None\
                    and vacuum_directory is None:
                try:
                    g_vac = calc["final_energy"]
                except KeyError:
                    continue

        if g_vac is None and vacuum_directory is not None:
            drone = QChemDrone(runs=[])
            calc_dir = join(vacuum_directory, mol_id)

            if files_mol_prefix:
                vac_calc = drone.assimilate(path=calc_dir,
                                            input_file="{}.qin".format(mol_id),
                                            output_file="{}.qout".format(mol_id),
                                            multirun=False)
            else:
                vac_calc = drone.assimilate(path=calc_dir,
                                            input_file="mol.qin",
                                            output_file="mol.qout",
                                            multirun=False)

            g_vac = vac_calc["output"].get("final_energy_sp",
                                           vac_calc["output"]["final_energy"])

        if g_vac is None or g_liq is None:
            raise RuntimeError("Could not find energy values from"
                               " molecule calculations!")
        else:
            delta_g_solv = (g_liq - g_vac) * 627.509471 * 4184
            sol_data["solubility"] = calculate_solubility(vp, delta_g_solv)

            return sol_data


class MolThermDataProcessorOld:
    """
    This class can be used to extract data from MolThermWorkflow workflows,
    including extracting thermo data from calculations and generating predicted
    boiling and melting points.
    """

    def __init__(self, base_dir, mol_dir="molecules", rxn_dir="reactions",
                 db_file="db.json"):
        """
        :param base_dir: Directory where all data should be stored.
        :param mol_dir: Subdirectory where molecule files and calculations are
        stored. Default is "molecules".
        :param rxn_dir: Subdirectory where reaction metadata files (and
        potentially calculations) are stored. Default is "reactions".
        :param db_file: Path to database config file.
        """

        self.reactant_pre = "rct"
        self.product_pre = "pro"

        self.base_dir = base_dir
        self.mol_dir = join(self.base_dir, mol_dir)
        self.rxn_dir = join(self.base_dir, rxn_dir)
        self.db_file = db_file

        try:
            self.db = QChemCalcDb.from_db_file(self.db_file)
        except:
            self.db = None

    def check_appropriate_dirs(self, dirs):
        """
        Returns only those reactions which have appropriate products and
        reactants (products, reactants have same number of atoms).

        This is not a sophisticated checking mechanism, and could probably be
        easily improved upon.

        :return:
        """

        add_up = []

        for d in dirs:
            path = join(self.base_dir, d)
            files = [f for f in listdir(path) if isfile(join(path, f))]
            rcts = [f for f in files if f.startswith(self.reactant_pre) and f.endswith(".mol")]
            pros = [f for f in files if f.startswith(self.product_pre) and f.endswith(".mol")]

            rct_mols = [get_molecule(join(self.base_dir, d, r)) for r in rcts]
            pro_mols = [get_molecule(join(self.base_dir, d, p)) for p in pros]

            total_pro_length = sum([len(p) for p in pro_mols])
            total_rct_length = sum([len(r) for r in rct_mols])

            if total_pro_length == total_rct_length:
                add_up.append(d)

        return add_up

    def extract_reaction_thermo_files(self, path, runs_pattern=None):
        """
        Naively scrape thermo data from QChem output files.

        :param path: Path to a subdirectory.
        :param runs_pattern: QChem drone assimilation requires a template for
            what calculations have completed. This should be represented by a
            link. Default None

        :return: dict {prop: value}, where properties are enthalpy, entropy.
        """

        base_path = join(self.base_dir, path)

        rct_files = [f for f in listdir(base_path) if f.endswith(".mol") and
                     f.startswith(self.reactant_pre)]
        pro_files = [f for f in listdir(base_path) if f.endswith(".mol") and
                     f.startswith(self.product_pre)]

        rct_ids = [extract_id(f) for f in rct_files]
        pro_ids = [extract_id(f) for f in pro_files]

        rct_thermo = {"enthalpy": 0, "entropy": 0, "energy": 0}
        pro_thermo = {"enthalpy": 0, "entropy": 0, "energy": 0}

        drone = QChemDrone(runs=runs_pattern)

        for mol in rct_files:
            prefix = mol.replace(".mol", "")

            calc_doc = drone.assimilate(path=base_path,
                                        input_file=prefix + ".in",
                                        output_file=prefix + ".out",
                                        multirun=False)

            rct_thermo["entropy"] += calc_doc["output"]["entropy"]
            rct_thermo["enthalpy"] += calc_doc["output"]["enthalpy"]
            rct_thermo["energy"] += calc_doc["output"]["final_energy"]

        for mol in pro_files:
            prefix = mol.replace(".mol", "")

            calc_doc = drone.assimilate(path=base_path,
                                        input_file=prefix + ".in",
                                        output_file=prefix + ".out",
                                        multirun=False)

            pro_thermo["entropy"] += calc_doc["output"]["entropy"]
            pro_thermo["enthalpy"] += calc_doc["output"]["enthalpy"]
            pro_thermo["energy"] += calc_doc["output"]["final_energy"]

        thermo_data = {}

        # Generate totals as ∆H = H_pro - H_rct, ∆S = S_pro - S_rct
        # Also ensures that units are appropriate (Joules/mol,
        # rather than cal/mol or kcal/mol, or hartree for energy)
        energy = (pro_thermo["energy"] - rct_thermo["energy"]) * 627.509
        enthalpy = (pro_thermo["enthalpy"] - rct_thermo["enthalpy"])
        print(path, energy, enthalpy)
        thermo_data["enthalpy"] = (energy + enthalpy) * 1000 * 4.184
        thermo_data["entropy"] = (pro_thermo["entropy"] - rct_thermo["entropy"]) * 4.184
        try:
            thermo_data["t_critical"] = thermo_data["enthalpy"] / thermo_data["entropy"]
        except ZeroDivisionError:
            thermo_data["t_critical"] = None

        result = {"thermo": thermo_data,
                  "directory": path,
                  "reactant_ids": rct_ids,
                  "product_ids": pro_ids}

        return result

    def extract_reaction_thermo_db(self, directory, collection="molecules",
                                   opt=None, freq=None, sp=None):
        """
        Gathers all relevant reaction parameters, including references to
        each job performed.

        :param directory: Directory name where the reaction is stored. Right
            now, this is the easiest way to identify the reaction. In the
            future, more sophisticated searching should be used.
        :param collection: Database collection from which thermo data should be
            extracted. Default is "molecules".
        :param opt: dict containing information about the optimization jobs. By
            default, this is None, and that information will be obtained by
            querying the self.db.tasks collection.
        :param freq: dict containing information about the frequency jobs. By
            default, this is None, and that information will be obtained by
            querying the self.db.tasks collection.
        :param sp: dict containing information about the single-point jobs. By
            default, this is None, and that information will be obtained by
            querying the self.db.tasks collection.

        :return: dict
        """

        if self.db is None:
            raise RuntimeError("Could not connect to database. Check db_file"
                               "and try again later.")

        # To extract enthalpy and entropy from calculation results
        def get_thermo(job):
            enthalpy = job["output"]["enthalpy"]
            entropy = job["output"]["entropy"]

            try:
                energy = job["output"]["final_energy_sp"]
            except KeyError:
                energy = job["output"]["final_energy"]

            if enthalpy is None:
                enthalpy = 0.0
            if entropy is None:
                entropy = 0.0
            if energy is None:
                energy = 0.0

            return {"enthalpy": enthalpy,
                    "entropy": entropy,
                    "energy": energy}

        if abspath(directory) != directory:
            directory = join(self.base_dir, directory)

        mol_files = [f for f in listdir(directory) if f.endswith(".mol")]

        dir_ids = [extract_id(f) for f in mol_files]

        db_collection = self.db.db[collection]
        records = []

        for mol_id in dir_ids:
            record = db_collection.find_one({"mol_id": str(mol_id)})
            records.append(record)

        # Sort files for if they are reactants or products
        reactants = []
        products = []
        for i, record in enumerate(records):
            filename = mol_files[i]
            if opt is None:
                for calc in record["calcs_reversed"][::-1]:
                    if calc["task"]["type"] == "opt" or \
                            calc["task"]["type"] == "optimization":
                        method = calc["input"]["rem"]["method"]
                        basis = calc["input"]["rem"]["basis"]
                        solvent_method = calc["input"]["rem"].get(
                            "solvent_method", None)
                        if solvent_method == "smd":
                            if calc["input"]["smx"] is None:
                                solvent = None
                            else:
                                solvent = calc["input"]["smx"]["solvent"]
                        elif solvent_method == "pcm":
                            solvent = calc["input"]["solvent"]
                        else:
                            solvent = None

                        opt = {"method": method,
                               "basis": basis,
                               "solvent_method": solvent_method,
                               "solvent": solvent}
                        break
            if freq is None:
                for calc in record["calcs_reversed"][::-1]:
                    if calc["task"]["type"] == "freq" or \
                            calc["task"]["type"] == "frequency":
                        method = calc["input"]["rem"]["method"]
                        basis = calc["input"]["rem"]["basis"]
                        solvent_method = calc["input"]["rem"].get(
                            "solvent_method", None)
                        if solvent_method == "smd":
                            if calc["input"]["smx"] is None:
                                solvent = None
                            else:
                                solvent = calc["input"]["smx"]["solvent"]
                        elif solvent_method == "pcm":
                            solvent = calc["input"]["solvent"]
                        else:
                            solvent = None

                        freq = {"method": method,
                                "basis": basis,
                                "solvent_method": solvent_method,
                                "solvent": solvent}
                        break
            if sp is None:
                for calc in record["calcs_reversed"][::-1]:
                    if calc["task"]["type"] == "sp":
                        method = calc["input"]["rem"]["method"]
                        basis = calc["input"]["rem"]["basis"]
                        solvent_method = calc["input"]["rem"].get(
                            "solvent_method", None)
                        if solvent_method == "smd":
                            if calc["input"]["smx"] is None:
                                solvent = None
                            else:
                                solvent = calc["input"]["smx"]["solvent"]
                        elif solvent_method == "pcm":
                            solvent = calc["input"]["solvent"]
                        else:
                            solvent = None

                        sp = {"method": method,
                              "basis": basis,
                              "solvent_method": solvent_method,
                              "solvent": solvent}
                        break

            if filename.startswith(self.reactant_pre):
                reactants.append(record)
            elif filename.startswith(self.product_pre):
                products.append(record)
            else:
                print("Skipping {} because it cannot be determined if it is"
                      "reactant or product.".format(filename))
                continue

        # Get ids
        reactant_ids = [r["mol_id"] for r in reactants]
        product_ids = [p["mol_id"] for p in products]

        # Get thermo data
        rct_thermo = [get_thermo(r) for r in reactants]
        pro_thermo = [get_thermo(p) for p in products]

        # Compile reaction thermo from reactant and product thermos
        delta_e = sum(p["energy"] for p in pro_thermo) - sum(r["energy"] for r in rct_thermo)
        delta_e *= 627.509
        delta_h = sum(p["enthalpy"] for p in pro_thermo) - sum(r["enthalpy"] for r in rct_thermo) + delta_e
        delta_h *= 1000 * 4.184
        delta_s = sum(p["entropy"] for p in pro_thermo) - sum(r["entropy"] for r in rct_thermo)
        delta_s *= 4.184
        thermo = {
            "enthalpy": delta_h,
            "entropy": delta_s
        }

        try:
            thermo["t_star"] = delta_h / delta_s
        except ZeroDivisionError:
            thermo["t_star"] = 0

        result = {"dir_name": directory,
                  "opt": opt,
                  "freq": freq,
                  "sp": sp,
                  "reactant_ids": reactant_ids,
                  "product_ids": product_ids,
                  "thermo": thermo}

        return result

    def record_molecule_data_db(self, mol_id, calc_dir, input_file, output_file,
                                collection="molecules"):
        """
        Compile calculation information for a single molecule and record it in
        the molecules collection.

        :param mol_id: Unique identifier for molecule (str)
        :param calc_dir: Directory where molecule information is stored.
        :param input_file: Basic format for input files. The Drone which
            compiles the molecule information will use this to pattern-match
            files.
        :param output_file: Basic format for output files. The Drone which
            compiles the molecule information will use this to pattern-match
            files.
        :param collection: Collection in which to store molecule information.
            Default is "molecules".
        :return:
        """

        runs = list(chain.from_iterable([["opt_" + str(ii), "freq_" + str(ii)] for ii in range(3)]))
        runs.append("sp")

        drone = QChemDrone(runs=runs)

        task_doc = drone.assimilate(
            path=calc_dir,
            input_file=input_file,
            output_file=output_file,
            multirun=False)

        task_doc["mol_id"] = mol_id

        if self.db is None:
            raise RuntimeError("Cannot record data to db without valid database"
                               " connection!")

        db_collection = self.db.db[collection]

        db_collection.insert_one(task_doc)

    def record_reaction_data_db(self, directory, molecules="molecules",
                                thermo="thermo", use_files=True, use_db=False,
                                opt=None, freq=None, sp=None):
        """
        Record thermo data in thermo collection.

        :param directory: Directory name where the reaction is stored. Right
            now, this is the easiest way to identify the reaction. In the
            future, more sophisticated searching should be used.
        :param molecules: Database collection from which to extract molecule
            data. Default is "molecules".
        :param thermo: Database collection in which to store thermo data.
            Default is "thermo".
        :param use_files: If set to True (default True), use
            get_reaction_thermo_files to gather data
        :param use_db: If set to True (default False), use
            extract_reaction_data to gather data
        :param opt: dict containing information about the optimization jobs. By
            default, this is None, and that information will be obtained by
            querying the self.db.tasks collection.
        :param freq: dict containing information about the frequency jobs. By
            default, this is None, and that information will be obtained by
            querying the self.db.tasks collection.
        :param sp: dict containing information about the single-point jobs. By
            default, this is None, and that information will be obtained by
            querying the self.db.tasks collection.

        :return:
        """

        if self.db is None:
            raise RuntimeError("Could not connect to database. Check db_file"
                               "and try again later.")

        db_collection = self.db.db[thermo]

        if use_db:
            db_collection.insert_one(self.extract_reaction_thermo_db(directory,
                                                                     collection=molecules,
                                                                     opt=opt,
                                                                     freq=freq,
                                                                     sp=sp))
        elif use_files:
            db_collection.insert_one(self.extract_reaction_thermo_files(directory))
        else:
            raise RuntimeError("Either database or files must be used to "
                               "extract thermo data.")

    def populate_collections(self, molecules="molecules", thermo=None,
                             overwrite=False, sp_job=True):
        """
        Mass insert into db collections using above data extraction and
        recording methods.

        :param molecules: Database collection in which to store molecule data.
            Default is "molecules".
        :param thermo: Database collection in which to store thermo data.
            Default is None, meaning that no thermo data will be stored.
        :param overwrite: If True (default False), overwrite molecules that are
            already included in the db.
        :param sp_job: If True (default), then ensure that all completed
            molecules have a sp output file.
        :return:
        """

        if self.db is None:
            raise RuntimeError("Cannot connect to database. Check configuration"
                               " file and try again.")

        mol_coll = self.db.db[molecules]
        completed_mols = self.get_completed_molecules(extra=True, sp_job=sp_job)
        mols_in_db = [mol for mol in mol_coll.find()]

        for mol_id, d, molfile in completed_mols:
            molfile = molfile.replace(".mol", "")
            if mol_id not in [m["mol_id"] for m in mols_in_db]:
                self.record_molecule_data_db(mol_id, join(self.base_dir, d),
                                             molfile+".qin", molfile+".qout",
                                             collection=molecules)
            elif overwrite:
                drone = QChemDrone()

                task_doc = drone.assimilate(
                    path=join(self.base_dir, d),
                    input_file=molfile+".qin",
                    output_file=molfile+".qout",
                    multirun=False)

                task_doc["mol_id"] = mol_id

                mol_coll.update_one({"mol_id": mol_id},
                                    {"$set": task_doc})

        if thermo is not None:
            thermo_coll = self.db.db[thermo]
            completed_rxns = self.get_completed_reactions(molecules=molecules)
            rxns_in_db = [rxn for rxn in thermo_coll.find()]

            for rxn in completed_rxns:
                if rxn not in [r["dir_name"].split("/")[-1]
                               for r in rxns_in_db]:
                    try:
                        self.record_reaction_data_db(rxn, molecules=molecules,
                                                     thermo=thermo,
                                                     use_files=False,
                                                     use_db=True)
                    except:
                        print(rxn)

                elif overwrite:
                    thermo_coll.update_one({"dir_name": join(self.base_dir, rxn)},
                                           {"$set": self.extract_reaction_thermo_db(rxn)})

    def update_molecules(self, collection="molecules"):
        """
        Update molecules collection with data from the subdirectories in
        self.base_dir.

        :param collection: Database collection in which to store molecule data.
            Default is "molecules".

        :return:
        """

        if self.db is None:
            raise RuntimeError("Cannot access database. Check configuration"
                               " settings and try again.")

        dirs = [d for d in listdir(self.base_dir) if
                isdir(join(self.base_dir, d)) and not d.startswith("block")]

        drone = QChemDrone()
        mol_coll = self.db.db[collection]

        for d in dirs:
            calc_dir = join(self.base_dir, d)

            files = [f for f in listdir(calc_dir) if isfile(join(calc_dir, f))]

            mol_names = set()

            for file in files:
                mol_names.add(file.split(".")[0])

            for mol_name in mol_names:

                new_calc = drone.assimilate(path=calc_dir,
                                            input_file=mol_name+".in",
                                            output_file=mol_name+".out",
                                            multirun=False)

                old_calc = mol_coll.find_one({"mol_id": mol_name})

                calcs_reversed = new_calc["calcs_reversed"] + old_calc["calcs_reversed"]

                output = new_calc["output"]

                walltime = new_calc["walltime"] + old_calc["walltime"]
                cputime = new_calc["cputime"] + old_calc["cputime"]

                mol_coll.update_one({"mol_id": mol_name},
                                    {"$set": {"calcs_reversed": calcs_reversed,
                                              "output": output,
                                              "walltime": walltime,
                                              "cputime": cputime}})

    def update_thermo(self, collection="thermo"):
        """
        Update thermo collection with data from the subdirectories in
        self.base_dir.

        :param collection: Database collection in which to store thermo data.
            Default is "thermo".

        :return:
        """

        if self.db is None:
            raise RuntimeError("Cannot access database. Check configuration"
                               " settings and try again.")

        dirs = [d for d in listdir(self.base_dir) if
                isdir(join(self.base_dir, d)) and not d.startswith("block")]

        thermo_coll = self.db.db[collection]

        to_update = []

        all_rxns = [r for r in thermo_coll.find()]

        for d in dirs:
            for rxn in all_rxns:
                if rxn["dir_name"].split("/")[-1] == d:
                    old_data = rxn
                    break

            pro_thermo = [self.get_molecule_data(m) for m in old_data["product_ids"]]
            rct_thermo = [self.get_molecule_data(m) for m in old_data["reactant_ids"]]

            # Compile reaction thermo from reactant and product thermos
            delta_e = sum(p["energy"] for p in pro_thermo) - sum(
                r["energy"] for r in rct_thermo)
            delta_h = sum(p["enthalpy"] for p in pro_thermo) - sum(
                r["enthalpy"] for r in rct_thermo) + delta_e
            delta_s = sum(p["entropy"] for p in pro_thermo) - sum(
                r["entropy"] for r in rct_thermo)
            thermo = {
                "enthalpy": delta_h,
                "entropy": delta_s
            }

            try:
                thermo["t_star"] = delta_h / delta_s
            except ZeroDivisionError:
                thermo["t_star"] = 0

            to_update.append((rxn["dir_name"], thermo))

        for update in to_update:
            thermo_coll.update_one({"dir_name": update[0]},
                                   {"$set": {"thermo": update[1]}})

    def find_reactions_common_reactant(self, mol_id, collection="reaxys"):
        """
        Queries reactions to identify all those that share a common reactant or
        product.

        :param mol_id: ID (str) of reactant molecule to be searched.
        :param collection: Database collection to search for reactions. Default
        is "reaxys".
        :return: list of dicts representing reaction metadata.
        """

        # Should return all reactions where mol_id is in rct_ids or pro_ids
        res_rct = self.db.db[collection].find({"rct_ids": mol_id})
        res_pro = self.db.db[collection].find({"pro_ids": mol_id})

        results = []

        if res_rct is not None:
            for res in res_rct:
                results.append(res["rxn_id"])
        if res_pro is not None:
            for res in res_pro:
                results.append(res["rxn_id"])

        return results

    def map_reactants_to_reactions(self, collection="reaxys"):
        """
        Construct a dict showing which directories share each reactant.

        This is useful for analysis of common reactants, and to identify the
        "source" of a given reactant (in which directory the calculation
        actually took place).

        :param collection: Database collection to search for reactions. Default
        is "reaxys".
        :return: dict mapping molecule IDs to lists of reaction IDs.
        """

        molecules = [d for d in listdir(self.mol_dir)
                     if isdir(join(self.mol_dir, d))]

        mapping = {}

        for mol_id in molecules:
            res = self.find_reactions_common_reactant(mol_id,
                                                      collection=collection)

            if len(res) != 0:
                mapping[mol_id] = res

        return mapping

    def get_completed_molecules(self, dirs=None, extra=False, sp_job=True):
        """
        Returns a list of molecules with completed opt, freq, and sp output
        files.

        :param dirs: List of directories to search for completed molecules.
        :params extra: If True, include directory of completed reaction and name
            of molfile along with mol_id
        :return: set of completed molecules
        """

        completed = set()

        all_dirs = [d for d in listdir(self.base_dir)
                    if isdir(join(self.base_dir, d)) and not d.startswith("block")]

        if dirs is not None:
            all_dirs = [d for d in all_dirs if d in dirs]

        for d in all_dirs:
            path = join(self.base_dir, d)
            mapping = associate_qchem_to_mol(self.base_dir, d)

            for molfile, qcfiles in mapping.items():
                mol_id = extract_id(molfile)

                for outfile in qcfiles["out"]:
                    if sp_job:
                        if "sp" in outfile:
                            spfile = QCOutput(join(path, outfile))

                            completion = spfile.data.get("completion", False)

                            # Currently will catch iefpcm or smd
                            if completion:
                                if extra:
                                    completed.add((mol_id, d, molfile))
                                else:
                                    completed.add(mol_id)
                    else:
                        if "freq" in outfile:
                            freqfile = QCOutput(join(path, outfile))

                            completion = freqfile.data.get("completion", False)

                            if completion:
                                if extra:
                                    completed.add((mol_id, d, molfile))
                                else:
                                    completed.add(mol_id)

        return completed

    def get_completed_reactions(self, molecules="molecules"):
        """
        Returns a list of directories (reactions) where all molecules are
        completed.

        :param molecules: Database collection in which to store molecule data.
            Default is "molecules".

        :return: list of directories with complete information.
        """

        if self.db is None:
            raise RuntimeError("Could not connect to database. Check db_file"
                               "and try again later.")

        collection = self.db.db[molecules]

        completed_molecules = [x["mol_id"] for x in collection.find()]

        completed_reactions = set()

        dirs = [d for d in listdir(self.base_dir) if isdir(join(self.base_dir, d)) and not d.startswith("block")]

        for d in dirs:
            path = join(self.base_dir, d)

            mols = [extract_id(f) for f in listdir(path) if isfile(join(path, f)) and f.endswith(".mol")]

            are_completed = [True if m in completed_molecules else False for m in mols]

            if all(are_completed):
                completed_reactions.add(d)

        return completed_reactions

    def get_molecule_data(self, mol_id, collection="molecules"):
        """
        Compile all useful molecular data for analysis, including molecule size
        (number of atoms), molecular weight, enthalpy, entropy, and functional
        groups.

        NOTE: This function automatically converts energy, enthalpy, and entropy
        into SI units (J/mol and J/mol*K)

        :param mol_id: Unique ID associated with the molecule.
        :param collection: Collection from which to extract molecule
            information. Default is "molecules".
        :return: dict of relevant molecule data.
        """

        mol_data = {"mol_id": mol_id}

        if self.db is None:
            raise RuntimeError("Cannot query database; connection is invalid."
                               " Try to connect again.")

        db_collection = self.db.db[collection]

        mol_entry = db_collection.find_one({"mol_id": mol_id})

        if mol_entry is None:
            print(mol_id)

        mol_data["enthalpy"] = mol_entry["output"]["enthalpy"] * 4.184 * 1000
        mol_data["entropy"] = mol_entry["output"]["entropy"] * 4.184

        try:
            final_energy = mol_entry["output"]["final_energy_sp"]
        except KeyError:
            final_energy = mol_entry["output"]["final_energy"]

        mol_data["energy"] = final_energy * 627.509 * 4.184 * 1000

        mol_dict = mol_entry["output"]["optimized_molecule"]
        mol_data["molecule"] = Molecule.from_dict(mol_dict)

        adaptor = BabelMolAdaptor(mol_data["molecule"])
        pbmol = adaptor.pybel_mol

        mol_data["number_atoms"] = len(mol_data["molecule"])
        mol_data["molecular_weight"] = pbmol.molwt
        mol_data["tpsa"] = pbmol.calcdesc()["TPSA"]

        extractor = FunctionalGroupExtractor(mol_data["molecule"])
        molgraph = extractor.molgraph
        func_grps = extractor.get_all_functional_groups()

        mol_data["functional_groups"] = extractor.categorize_functional_groups(func_grps)

        weights = nx.get_edge_attributes(molgraph.graph, "weight")
        bonds_checked = set()
        double_bonds = 0
        triple_bonds = 0
        for bond, weight in weights.items():
            # Remove index from multidigraph
            bond = (bond[0], bond[1])
            if int(weight) == 2 and bond not in bonds_checked:
                double_bonds += 1
            elif int(weight) == 3 and bond not in bonds_checked:
                triple_bonds += 1
            bonds_checked.add(bond)

        mol_data["double_bonds"] = double_bonds
        mol_data["triple_bonds"] = triple_bonds

        species = [str(s.specie) for s in mol_data["molecule"].sites]
        mol_data["species"] = dict(Counter(species))

        return mol_data

    def get_reaction_data(self, directory=None, mol_ids=None):
        """
        Compile all useful data for a set of molecules associated with a
        particular reaction. This data will be compiled on a reaction basis
        (difference between reactants and products) as well as an individual
        molecule basis.

        :param directory: Subdirectory where molecule data is located.
        :param mol_ids: List of unique IDs for molecules associated with the
            reaction
        :return: dict of relevant reaction data.
        """

        reaction_data = {}
        reaction_data["thermo"] = None

        if directory is not None:
            mol_ids = [extract_id(f) for f in listdir(join(self.base_dir, directory))
                        if f.endswith(".mol")]

            component_data = [self.get_molecule_data(m) for m in mol_ids]

            reaction_data["thermo"] = self.extract_reaction_thermo_db(directory)["thermo"]

        elif mol_ids is not None:
            component_data = [self.get_molecule_data(m) for m in mol_ids]


        else:
            raise ValueError("get_reaction_data requires either a directory or "
                             "a set of molecule ids.")

        component_data = sorted(component_data, key=lambda x: len(x["molecule"]))

        reaction_data["dir_name"] = directory
        reaction_data["mol_ids"] = mol_ids
        reaction_data["product"] = component_data[-1]
        reaction_data["reactants"] = component_data[:-1]

        if reaction_data["thermo"] is None:
            reaction_data["thermo"] = {}
            pro_h = reaction_data["product"]["enthalpy"] + reaction_data["product"]["energy"]
            rct_h = sum(r["enthalpy"] + r["energy"]
                        for r in reaction_data["reactants"])
            reaction_data["thermo"]["enthalpy"] = pro_h - rct_h

            pro_s = reaction_data["product"]["entropy"]
            rct_s = sum(r["entropy"] for r in reaction_data["reactants"])
            reaction_data["thermo"]["entropy"] = pro_s - rct_s

            try:
                reaction_data["thermo"]["t_star"] = reaction_data["thermo"]["enthalpy"] / reaction_data["thermo"]["entropy"]
            except ZeroDivisionError:
                reaction_data["thermo"]["t_star"] = 0

        return reaction_data
