from os import listdir
from os.path import join, isfile, isdir, abspath
import shutil
from collections import Counter
import re

import networkx as nx

from pymatgen.core.structure import Molecule
from pymatgen.analysis.functional_groups import FunctionalGroupExtractor
from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.io.qchem.outputs import QCOutput

from atomate.qchem.database import QChemCalcDb
from atomate.qchem.drones import QChemDrone

from moltherm.compute.drones import MolThermDrone
from moltherm.compute.utils import get_molecule, extract_id, associate_qchem_to_mol

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

    def __init__(self, base_dir, reactant_pre="rct_", product_pre="pro_",
                 db_file="db.json"):
        """
        :param base_dir: Directory where input and output data should be stored.
        :param reactant_pre: Prefix for reactant files.
        :param product_pre: Prefix for product files.
        :param db_file: Path to database config file.
        """

        self.base_dir = base_dir
        self.reactant_pre = reactant_pre
        self.product_pre = product_pre
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
        # Note: After all sp jobs are finished, it should be unnecessary to use
        # energy_opt
        def get_thermo(job):
            enthalpy = None
            entropy = None
            energy = None

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

        drone = MolThermDrone()

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

    def record_reaction_data_db(self, directory, collection="thermo",
                                use_files=True, use_db=False, opt=None,
                                freq=None, sp=None):
        """
        Record thermo data in thermo collection.

        :param directory: Directory name where the reaction is stored. Right
            now, this is the easiest way to identify the reaction. In the
            future, more sophisticated searching should be used.
        :param collection: Database collection in which to store thermo data.
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

        db_collection = self.db.db[collection]

        if use_db:
            db_collection.insert_one(self.extract_reaction_thermo_db(directory,
                                                                  opt=opt,
                                                                  freq=freq,
                                                                  sp=sp))
        elif use_files:
            db_collection.insert_one(self.extract_reaction_thermo_files(directory))
        else:
            raise RuntimeError("Either database or files must be used to "
                               "extract thermo data.")

    def record_reaction_data_file(self, directory, filename="thermo.txt",
                                  use_files=True, use_db=False, opt=None,
                                  freq=None, sp=None):
        """
        Record thermo data in thermo.txt file.

        Note: This function does NOT store the reactant and product IDs

        :param directory: Directory name where the reaction is stored. Right
            now, this is the easiest way to identify the reaction. In the
            future, more sophisticated searching should be used.
        :param filename: File (within directory) where data should be stored.
            By default, it will be stored in thermo.txt.
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

        if abspath(directory) != directory:
            directory = join(self.base_dir, directory)

        with open(join(directory, filename), "w+") as file:
            if use_db:
                data = self.extract_reaction_data(directory, opt=opt, freq=freq,
                                                  sp=sp)
            elif use_files:
                data = self.get_reaction_thermo_files(directory)
            else:
                raise RuntimeError("Either database or files must be used to "
                                   "extract thermo data.")

            file.write("Directory: {}\n".format(data["dir_name"]))
            file.write("Optimization Input: {}\n".format(data.get("opt", "")))
            file.write("Frequency Input: {}\n".format(data.get("freq", "")))
            file.write("Single-Point Input: {}\n".format(data.get("sp", "")))
            file.write("Reaction Enthalpy: {}\n".format(data["thermo"]["enthalpy"]))
            file.write("Reaction Entropy: {}\n".format(data["thermo"]["entropy"]))
            file.write("Turning Temperature: {}\n".format(data["thermo"]["t_critical"]))

    def populate_collections(self, molecules="molecules", thermo=None, overwrite=False):
        """
        Mass insert into db collections using above data extraction and
        recording methods.

        :param molecules: Database collection in which to store molecule data.
            Default is "molecules".
        :param thermo: Database collection in which to store thermo data.
            Default is None, meaning that no thermo data will be stored.
        :param overwrite: If True (default False), overwrite molecules that are
        already included in the db.
        :return:
        """

        if self.db is None:
            raise RuntimeError("Cannot connect to database. Check configuration"
                               " file and try again.")

        mol_coll = self.db.db[molecules]
        completed_mols = self.get_completed_molecules(extra=True)
        mols_in_db = [mol for mol in mol_coll.find()]

        for mol_id, d, molfile in completed_mols:
            molfile = molfile.replace(".mol", "")
            if mol_id not in [m["mol_id"] for m in mols_in_db]:
                self.record_molecule_data_db(mol_id, join(self.base_dir, d),
                                             molfile+".in", molfile+".out",
                                             collection=molecules)
            elif overwrite:
                drone = MolThermDrone()

                task_doc = drone.assimilate(
                    path=join(self.base_dir, d),
                    input_file=molfile+".in",
                    output_file=molfile+".out",
                    multirun=False)

                task_doc["mol_id"] = mol_id

                mol_coll.update_one({"mol_id": mol_id},
                                    {"$set": task_doc})

        if thermo is not None:
            thermo_coll = self.db.db[thermo]
            completed_rxns = self.get_completed_reactions()
            rxns_in_db = [rxn for rxn in thermo_coll.find()]

            for rxn in completed_rxns:
                if rxn not in [r["dir_name"].split("/")[-1]
                               for r in rxns_in_db]:
                    self.record_reaction_data_db(rxn, collection=thermo,
                                                 use_files=False, use_db=True)

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

        drone = MolThermDrone()
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

    def copy_outputs_across_directories(self):
        """
        Copy output files between subdirectories to ensure that all reaction
        directories that need outputs of a given molecule will have them.

        Note: This function should not be used unless necessary. It was written
        because for each directory, only a single database entry was being made
        (because db entries were being overwritten by default.

        :return:
        """

        files_copied = 0

        dirs = [d for d in listdir(self.base_dir) if
                isdir(join(self.base_dir, d)) and not d.startswith("block")]
        print("Number of directories: {}".format(len(dirs)))

        for start_d in dirs:
            start_p = join(self.base_dir, start_d)
            mol_files = [f for f in listdir(start_p) if isfile(join(start_p, f)) and f.endswith(".mol")]
            out_files = [f for f in listdir(start_p) if isfile(join(start_p, f)) and ".out" in f]

            for mf in mol_files:
                is_covered = False
                mol_id = extract_id(mf)

                mol_obj = get_molecule(join(start_p, mf))

                for out in out_files:
                    qcout = QCOutput(join(start_p, out))
                    if sorted(qcout.data["initial_molecule"].species) == sorted(mol_obj.species):
                        # If there is already output, do not copy any files
                        is_covered = True

                if is_covered:
                    continue

                for other_d in dirs:
                    if other_d == start_d:
                        continue
                    if is_covered:
                        break

                    other_p = join(self.base_dir, other_d)
                    # Check if this id is present
                    other_mol_files = [f for f in listdir(other_p) if isfile(join(other_p, f)) and f.endswith(".mol") and mol_id in f]
                    other_out_files = [f for f in listdir(other_p) if isfile(join(other_p, f)) and ".out" in f]
                    to_copy = []
                    for other_mol in other_mol_files:
                        if other_mol.startswith(self.product_pre):
                            to_copy = [f for f in other_out_files if
                                       f.startswith(self.product_pre)]
                        elif other_mol.startswith(self.reactant_pre):
                            to_check = [f for f in other_out_files if f.startswith(self.reactant_pre)]
                            to_copy = []
                            for file in to_check:
                                qcout = QCOutput(join(other_p, file))
                                if qcout.data["initial_molecule"].species == mol_obj.species:
                                    to_copy.append(file)
                        else:
                            to_copy = []
                    for file in to_copy:
                        shutil.copyfile(join(other_p, file), join(start_p, file + "_copy"))
                        files_copied += 1

                    if files_copied > 0:
                        is_covered = True
        print("Number of files copied: {}".format(files_copied))

    def find_common_reactants(self, rct_id):
        """
        Searches all subdirectories for those that have reactant .mol files with
        unique id rct_id.

        :param rct_id: String representing unique identifier for Reaxys
            molecules.
        :return: List of reaction directories containing the given reactant.
        """
        results = []
        for d in listdir(self.base_dir):
            if isdir(join(self.base_dir, d)) and not d.startswith("block"):
                for f in listdir(join(self.base_dir, d)):
                    if rct_id in f:
                        results.append(d)
        return results

    def map_reactants_to_reactions(self):
        """
        Construct a dict showing which directories share each reactant.

        This is useful for analysis of common reactants, and to identify the
        "source" of a given reactant (in which directory the calculation
        actually took place).

        :return:
        """

        mapping = {}
        dirs = [d for d in listdir(self.base_dir)
                if isdir(join(self.base_dir, d)) and not d.startswith("block")]

        for d in dirs:
            if isdir(join(self.base_dir, d)) and not d.startswith("block"):
                molfiles = [f for f in listdir(join(self.base_dir, d))
                            if f.endswith(".mol")
                            and f.startswith(self.reactant_pre)]
                for file in molfiles:
                    f_id = extract_id(file)
                    if f_id in mapping:
                        mapping[f_id].append(d)
                    else:
                        mapping[f_id] = [d]

        return mapping

    def get_completed_molecules(self, dirs=None, extra=False):
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
                    if "sp" in outfile:
                        spfile = QCOutput(join(path, outfile))

                        completion = spfile.data.get("completion", False)

                        # Currently will catch iefpcm or smd
                        if completion:
                            if extra:
                                completed.add((mol_id, d, molfile))
                            else:
                                completed.add(mol_id)

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

        collection = self.db.db["molecules"]

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

        collection = self.db.db["molecules"]

        mol_entry = collection.find_one({"mol_id": mol_id})

        mol_data["enthalpy"] = mol_entry["output"]["enthalpy"] * 4.184 * 1000
        mol_data["entropy"]  = mol_entry["output"]["entropy"] * 4.184

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

    @staticmethod
    def parse_epi_suite_data(file):
        """
        Parse predicted data from the US EPA's EPI Suite batch mode.

        Currently, this function only extracts the predicted boiling point,
        melting point, and solubility.

        :param file: Path to EPI Suite output file.
        :return: list of dicts with predicted molecular data
        """

        parsed_results = []

        with open(file, 'r') as file:
            entries = file.read().split("\n\n========================\n\n")[1:-1]

            for entry in entries:
                smiles = re.search(r"SMILES\s+:\s+([A-Za-z0-9=\(\)#\[\]\+\-@]+)",
                                   entry)
                if smiles:
                    smiles = smiles.group(1)
                else:
                    smiles = None
                name = re.search(r"CHEM\s+:\s+([A-Z/_a-z0-9]+)\s+:\s+([A-Za-z0-9]+\n*\s*[A-Za-z0-9]+)", entry)
                if name:
                    dir_name = name.group(1)
                    mol_id = name.group(2).replace("\n         ", "").replace("\nMOL", "")
                else:
                    name = re.search(
                        r"CHEM\s+:\s+([A-Z/_a-z0-9]+\n\s+[A-Z/_a-z0-9]+)\s+:\s+([A-Za-z0-9]+)",
                        entry)
                    if name:
                        dir_name = name.group(1).replace("\n         ", "")
                        mol_id = name.group(2)
                    else:
                        print(entry)
                        dir_name = None
                        mol_id = None
                bp = re.search(r"\s+Boiling Pt \(deg C\):\s+([0-9]+\.[0-9]+)\s+\(Adapted Stein & Brown method\)", entry)
                if bp:
                    bp = float(bp.group(1))
                else:
                    bp = None
                mp = re.search(r"\s+Melting Pt \(deg C\):\s+([0-9]+\.[0-9]+)\s+\(Mean or Weighted MP\)", entry)
                if mp:
                    mp = float(mp.group(1))
                else:
                    mp = None
                solubility = re.search(r"\s+Water Solubility at 25 deg C \(mg/L\):\s+([e0-9\+\-\.]+)", entry)
                if solubility:
                    solubility = float(solubility.group(1))
                else:
                    solubility = None

                if dir_name is not None:
                    dir_name.replace("\n         ", "")

                parsed_results.append({"mol_id": mol_id,
                                       "dir_name": dir_name,
                                       "smiles": smiles,
                                       "bp": bp,
                                       "mp": mp,
                                       "solubility": solubility})
        return parsed_results