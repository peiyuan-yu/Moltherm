from os import listdir
from os.path import join, isfile, isdir, abspath
import operator
import shutil

from bs4 import BeautifulSoup

from pymatgen.io.qchem_io.sets import OptSet, FreqSet, SinglePointSet
from pymatgen.io.babel import BabelMolAdaptor

from fireworks import Workflow, LaunchPad

from atomate.qchem.database import QChemCalcDb

from moltherm.compute.fireworks import OptFreqSPFW
from moltherm.compute.outputs import QCOutput

__author__ = "Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "June 2018"


""" This file contains functions necessary to actually perform our workflow
What constitutes the workflow will naturally evolve over time. Right now, this
is nothing more than a script describing basic steps.
For now, we want to:
    - Find the ideal starting conformer for reactants & products
    - Run through QChem (for now using pymatgen, later using atomate/custodian)
    - Determine reaction enthalpy, entropy
    - Calculate heat capacity as function of temperature
    - Determine working temperature range using QSPR (or database, if available?)
    - Perform some analysis to rank candidate reactions

TODO list:
    - Learn how to use Drones and Queens to parallelize for large sets of data
    - Figure out how to query with pymatgen-db and pymongo
"""


def get_molecule(molfile):
    """
    Create pymatgen Molecule object from molecule data file.

    In addition to parsing the input, this function also performs a conformer
    search to get a reasonable starting structure.

    :param molfile: Absolute path to structure file (.mol, .sdf, etc.)
    :return: Molecule.
    """

    obmol = BabelMolAdaptor.from_file(molfile, file_format="mol")
    # OBMolecule does not contain pymatgen Molecule information
    # So, we need to wrap the obmol in a BabelMolAdapter and extract
    obmol.add_hydrogen()
    obmol.make3d()
    obmol.localopt()

    return obmol.pymatgen_mol


def generate_opt_input(molfile, qinfile, basis_set="6-311++G*",
                       pcm_dielectric=None, overwrite_inputs=None):
    """
    Generates a QChem input file from Molecule after conformer search.

    :param molfile: Absolute path to the input file (.mol, .sdf, etc.)
    :param qinfile: Absolute path to the output file (.in)
    :param basis_set: To overwrite default basis.
    :param pcm_dielectric: To use solvent
    :param overwrite_inputs: To overwrite any set defaults
    :return:

    """
    mol = get_molecule(molfile)

    qcinput = OptSet(mol, basis_set=basis_set, pcm_dielectric=pcm_dielectric,
                     overwrite_inputs=overwrite_inputs)

    qcinput.write_file(qinfile)


def generate_freq_input(qoutfile, qinfile, basis_set="6-311++G*",
                       pcm_dielectric=None, overwrite_inputs=None):
    """
    Parses a QChem output file for ideal structure and then returns a QChem
    input file for frequency calculations (to determine enthalpy and entropy).

    :param qoutfile: Absolute path to the QChem output file (.out)
    :param qinfile: Absolute path to the QChem input file (.in)
    :return:
    """

    output = QCOutput(qoutfile)

    if len(output.data.get("molecule_from_optimized_geometry", [])) > 0:
        mol = output.data["molecule_from_optimized_geometry"]
    else:
        try:
            mol = output.data["molecule_from_last_geometry"]
        except KeyError:
            raise RuntimeError("No molecule to use as input")

    qcinput = FreqSet(mol, basis_set=basis_set, pcm_dielectric=pcm_dielectric,
                      overwrite_inputs=overwrite_inputs)

    qcinput.write_file(qinfile)


def generate_single_point_input(qoutfile, qinfile, basis_set="6-311++G*",
                       pcm_dielectric=None, overwrite_inputs=None):
    """
    Parse QChem output file for ideal structure and then returns a QChem
    input file for single-point calculations.

    :param qoutfile:
    :param qinfile:
    :return:
    """

    output = QCOutput(qoutfile)

    if len(output.data.get("molecule_from_optimized_geometry", [])) > 0:
        mol = output.data["molecule_from_optimized_geometry"]
    else:
        try:
            mol = output.data["molecule_from_last_geometry"]
        except KeyError:
            raise RuntimeError("No molecule to use as input")

    qcinput = SinglePointSet(mol, basis_set=basis_set,
                             pcm_dielectric=pcm_dielectric,
                             overwrite_inputs=overwrite_inputs)

    qcinput.write_file(qinfile)


def find_common_solvents(base_dir):
    """
    Iteratively scrape through list of reaction subdirectories to create a
    {solvent: occurrence} mapping.

    :param base_dir: Directory to begin search in.
    :return: dict {solvent: occurrence}
    """

    rct_dirs = [d for d in listdir(base_dir) if isdir(join(base_dir, d))]

    solvent_occurrence = {}

    for rct_dir in rct_dirs:
        if "meta.xml" not in listdir(join(base_dir, rct_dir)):
            # If metadata has not been recorded, solvent cannot be determined
            continue

        with open(join(base_dir, rct_dir, "meta.xml"), "r") as file:
            parsed = BeautifulSoup(file.read(), "lxml-xml")

            solvents = parsed.find("solvents").text.split(",")

            for solvent in solvents:
                current_value = solvent_occurrence.get(solvent, 0)
                solvent_occurrence[solvent] = current_value + 1

    return sorted(solvent_occurrence.items(), key=operator.itemgetter(1))


def get_reactions_common_solvent(base_dir, outdir, solvent):
    """
    Identify all reactions from a set of reactions which share a solvent.

    :param base_dir: Directory to begin search in.
    :param outdir: Directory to put all reactions with the common solvent
    :param solvent: Solvent of interest
    :return:
    """

    rct_dirs = [d for d in listdir(base_dir) if isdir(join(base_dir, d))]

    common_solvent = []

    for rct_dir in rct_dirs:
        if "meta.xml" not in listdir(join(base_dir, rct_dir)):
            # If metadata has not been recorded, solvent cannot be determined
            continue

        with open(join(base_dir, rct_dir, "meta.xml"), "r") as file:
            parsed = BeautifulSoup(file.read(), "lxml-xml")

            solvents = parsed.find("solvents").text.split(",")
            solvents = [s.lower() for s in solvents]

            if solvent.lower() in solvents:
                common_solvent.append(rct_dir)


    num_copied = 0
    for rct_dir in common_solvent:
        shutil.copytree(join(base_dir, rct_dir), join(outdir, rct_dir))
        num_copied += 1

    print("{} reactions with solvent {}".format(str(num_copied), solvent))


class MolThermWorkflow:
    """
    This class contains all functionality needed to perform Atomate
    workflows on molecular data.
    """

    def __init__(self, base_dir, subdirs=False, reactant_pre="rct_",
                 product_pre="pro_", db_file="db.json"):
        """
        :param base_dir: Directory where input and output data should be stored.
        :param subdirs: Is data all stored in one directory (False), or is it
        separated among subdirectories (True)?
        :param reactant_pre: Prefix for reactant files.
        :param product_pre: Prefix for product files.
        :param db_file: Path to database config file.
        """

        self.base_dir = base_dir
        self.subdirs = subdirs
        self.reactant_pre = reactant_pre
        self.product_pre = product_pre
        self.db_file = db_file

        try:
            self.db = QChemCalcDb.from_db_file(self.db_file)
        except:
            self.db = None

    def quick_check(self, dirs):
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

    def get_single_reaction_workflow(self, name_pre="opt_freq_sp", path=None,
                                     filenames=None, max_cores=64,
                                     qchem_input_params=None,
                                     sp_params=None):
        """
        Generates a Fireworks Workflow to find the structures and energies of
        the reactants and products of a single reaction.

        :param name_pre: str indicating the prefix which should be used for all
        Firework names
        :param path: Specified (sub)path in which to run the reaction. By
        default, this is None, and the Fireworks will run in self.base_dir
        :param filenames: Specified files within the path (if self.base_dir or
        a subdirectory) that should be considered a part of this reaction. If
        None, assume all files in the directory are to be involved.
        :param max_cores: int specifying number of processes/threads that can
        be used for this workflow.
        :param qchem_input_params: dict
        :param sp_params: For OptFreqSPFW, single-point calculations can be
        treated differently from Opt and Freq. In this case, another dict
        for sp must be used.
        :return: Workflow
        """

        fws = []

        if self.subdirs:
            base_path = join(self.base_dir, path)
        else:
            base_path = self.base_dir

        if filenames:
            rcts = [f for f in filenames if f.startswith(self.reactant_pre) and
                    f.endswith(".mol")]
            pros = [f for f in filenames if f.startswith(self.product_pre) and
                    f.endswith(".mol")]
            print(rcts)
            print(pros)
        else:
            # Assume that every file in the directory is part of the reaction
            files = [f for f in listdir(base_path) if isfile(join(base_path, f))
                     and f.endswith(".mol")]
            rcts = [f for f in files if f.startswith(self.reactant_pre)]
            pros = [f for f in files if f.startswith(self.product_pre)]

        for i, rct in enumerate(rcts):
            mol = get_molecule(join(base_path, rct))

            infile = join(base_path, self.reactant_pre + str(i) + ".in")
            outfile = join(base_path, self.reactant_pre + str(i) + ".out")

            fw = OptFreqSPFW(molecule=mol,
                             name="{}: {}/{}".format(name_pre, path, rct),
                             qchem_cmd="qchem -slurm",
                             input_file=infile,
                             output_file=outfile,
                             qclog_file=join(base_path, self.reactant_pre + str(i) + ".qclog"),
                             max_cores=max_cores,
                             qchem_input_params=qchem_input_params,
                             sp_params=sp_params,
                             db_file=self.db_file)

            fws.append(fw)

        for i, pro in enumerate(pros):
            mol = get_molecule(join(base_path, pro))

            infile = join(base_path, self.product_pre + str(i) + ".in")
            outfile = join(base_path, self.product_pre + str(i) + ".out")

            fw = OptFreqSPFW(molecule=mol,
                             name="{}: {}/{}".format(name_pre, path, pro),
                             qchem_cmd="qchem -slurm",
                             input_file=infile,
                             output_file=outfile,
                             qclog_file=join(base_path, self.product_pre + str(i) + ".qclog"),
                             max_cores=max_cores,
                             qchem_input_params=qchem_input_params,
                             sp_params=sp_params,
                             db_file=self.db_file)

            fws.append(fw)

        return Workflow(fws)

    def get_reaction_set_workflow(self, name_pre="opt_freq_sp", max_cores=64,
                                  qchem_input_params=None,
                                  sp_params=None):
        """Generates a Fireworks Workflow to find the structures and energies of
        the reactants and products of a single reaction.

        Note: as written now, this function will only work if self.subdirs is
        True; that is, only if each reaction is in a separate subdirectory.
        Later additions could allow for some other means of specifying the
        separate reactions within a single directory.

        :param name_pre: str indicating the prefix which should be used for all
        Firework names
        :param max_cores: int specifying number of processes/threads that can
        be used for this workflow.
        :param qchem_input_params: dict
        :param sp_params: For OptFreqSPFW, single-point calculations can be
        treated differently from Opt and Freq. In this case, another dict
        for sp must be used.

        :return: Workflow
        """

        if not self.subdirs:
            raise RuntimeError("Cannot run get_reaction_set_workflow();"
                               "Need reactions components to be isolated in"
                               "different subdirectories.")

        def extract_id(string):
            return string.split("/")[-1].rstrip(".mol").split("_")[-1]

        fws = []

        dirs = [d for d in listdir(self.base_dir) if isdir(join(self.base_dir, d))]

        # Only set up a workflow if it is worthwhile (the reaction actually
        # proceeds as written, and all atoms add up)
        appropriate_dirs = self.quick_check(dirs)

        if self.db is not None:
            all_fws = self.db.collection.find()

            # Keep track of which molecules have already been run as jobs before
            molecules_registered = [extract_id(fw["task_label"])
                                    for fw in all_fws]
        else:
            molecules_registered = []

        for d in appropriate_dirs:
            path = join(self.base_dir, d)
            files = [f for f in listdir(path) if isfile(join(path, f)) and f.endswith(".mol")]
            rcts = [f for f in files if f.startswith(self.reactant_pre)]
            pros = [f for f in files if f.startswith(self.product_pre)]

            for i, rct in enumerate(rcts):
                mol_id = rct.rstrip(".mol").split("_")[-1]

                if mol_id in molecules_registered:
                    continue
                else:
                    molecules_registered.append(mol_id)

                mol = get_molecule(join(self.base_dir, d, rct))

                infile = join(path, self.reactant_pre + str(i) + ".in")
                outfile = join(path, self.reactant_pre + str(i) + ".out")

                fw = OptFreqSPFW(molecule=mol,
                                 name="{}: {}/{}".format(name_pre, d, rct),
                                 qchem_cmd="qchem -slurm",
                                 input_file=infile,
                                 output_file=outfile,
                                 qclog_file=join(path, self.reactant_pre + str(i) + ".qclog"),
                                 max_cores=max_cores,
                                 qchem_input_params=qchem_input_params,
                                 sp_params=sp_params,
                                 db_file=self.db_file)

                fws.append(fw)

            for i, pro in enumerate(pros):
                mol_id = pro.rstrip(".mol").split("_")[-1]

                if mol_id in molecules_registered:
                    continue
                else:
                    molecules_registered.append(mol_id)

                mol = get_molecule(join(self.base_dir, d, pro))

                infile = join(path, self.product_pre + str(i) + ".in")
                outfile = join(path, self.product_pre + str(i) + ".out")

                fw = OptFreqSPFW(molecule=mol,
                                 name="{}: {}/{}".format(name_pre, d, pro),
                                 qchem_cmd="qchem -slurm",
                                 input_file=infile,
                                 output_file=outfile,
                                 qclog_file=join(path, self.product_pre + str(i) + ".qclog"),
                                 max_cores=max_cores,
                                 qchem_input_params=qchem_input_params,
                                 sp_params=sp_params,
                                 db_file=self.db_file)

                fws.append(fw)

        return Workflow(fws)

    @staticmethod
    def add_workflow(workflow):
        """
        Use Fireworks to add a generated workflow.

        :param workflow: a Workflow object (should have been generated by one
        of the workflow-generating functions about).
        :return:
        """

        launchpad = LaunchPad.auto_load()
        launchpad.add_wf(workflow)


class MolThermAnalysis:
    """
    This class can be used to analyze data from MolThermWorkflow workflows,
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

    def get_reaction_thermo_files(self, path=None):
        """
        Naively scrape thermo data from QChem output files.

        :param path: Path to a subdirectory.

        :return: dict {prop: value}, where properties are enthalpy, entropy.
        """

        if path is not None:
            base_path = join(self.base_dir, path)
        else:
            base_path = self.base_dir

        def extract_id(string):
            return string.split("/")[-1].rstrip(".mol").split("_")[-1]

        rct_ids = [extract_id(f) for f in listdir(base_path) if
                   f.endswith(".mol") and f.startswith(self.reactant_pre)]

        pro_ids = [extract_id(f) for f in listdir(base_path) if
                   f.endswith(".mol") and f.startswith(self.product_pre)]

        rct_map = {mol: [f for f in listdir(base_path) if mol in f and ".out" in f] for mol in rct_ids}
        pro_map = {mol: [f for f in listdir(base_path) if mol in f and ".out" in f] for mol in pro_ids}

        rct_thermo = {"enthalpy": 0, "entropy": 0}
        pro_thermo = {"enthalpy": 0, "entropy": 0}

        for mol in rct_map.keys():
            enthalpy = 0
            entropy = 0
            energy_opt = 0
            energy_sp = 0

            for out in rct_map[mol]:
                qcout = QCOutput(join(base_path, mol))

                enthalpy += qcout.data.get("enthalpy", 0)
                entropy += qcout.data.get("entropy", 0)
                energy_opt += qcout.data.get("final_energy", 0)
                energy_sp += qcout.data.get("final_energy_sp", 0)

            rct_thermo["enthalpy"] += (enthalpy - energy_opt) + energy_sp
            rct_thermo["entropy"] += entropy

        for mol in pro_map.keys():
            enthalpy = 0
            entropy = 0
            energy_opt = 0
            energy_sp = 0

            for out in pro_map[mol]:
                qcout = QCOutput(join(base_path, mol))

                enthalpy += qcout.data.get("enthalpy", 0)
                entropy += qcout.data.get("entropy", 0)
                energy_opt += qcout.data.get("final_energy", 0)
                energy_sp += qcout.data.get("final_energy_sp", 0)

            pro_thermo["enthalpy"] += (enthalpy - energy_opt) + energy_sp
            pro_thermo["entropy"] += entropy

        thermo_data = {}

        # Generate totals as ∆H = H_pro - H_rct, ∆S = S_pro - S_rct
        thermo_data["enthalpy"] = pro_thermo["enthalpy"] - rct_thermo["enthalpy"]
        thermo_data["entropy"] = pro_thermo["entropy"] - pro_thermo["entropy"]
        thermo_data["t_critical"] = thermo_data["enthalpy"] / thermo_data["entropy"]

        return thermo_data

    def extract_reaction_data(self, directory, opt=None, freq=None, sp=None):
        """
        Gathers all relevant reaction parameters, including references to
        each job performed.

        :param directory: Directory name where the reaction is stored. Right
            now, this is the easiest way to identify the reaction. In the
            future, more sophisticated searching should be used.
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
            enthalpy = 0
            entropy = 0
            energy_opt = 0
            energy_sp = 0

            for calc in job["calcs_reversed"]:
                if calc["task"]["type"] == "opt" or calc["task"]["type"] == "optimization":
                    energy_opt = calc["final_energy"]
                if calc["task"]["type"] == "freq" or calc["task"]["type"] == "frequency":
                    enthalpy = calc["enthalpy"]
                    entropy = calc["entropy"]
                if calc["task"]["type"] == "sp":
                    energy_sp = calc["final_energy_sp"]

            return {"enthalpy": (enthalpy - energy_opt) + energy_sp,
                    "entropy": entropy}

        if abspath(directory) != directory:
            directory = join(self.base_dir, directory)

        # This is horribly inefficient. Should change the how data is stored
        # to allow for nicer queries
        all_records = self.db.collection.find()
        records = []

        dir_ids = [extract_id(f) for f in listdir(directory) if
                   f.endswith(".mol")]

        for record in all_records:
            molecule_id = extract_id(record["task_label"])
            if record["dir_name"] == directory or molecule_id in dir_ids:
                records.append(record)

        # Sort files for if they are reactants or products
        reactants = []
        products = []
        for record in records:
            filename = record["task_label"].split(" : ")[-1]

            if opt is None:
                for calc in record["calcs_reversed"]:
                    if calc["task"]["type"] == "opt" or \
                            calc["task"]["type"] == "optimization":
                        method = calc["input"]["rem"]["method"]
                        basis = calc["input"]["rem"]["basis"]
                        solvent_method  = calc["input"]["rem"].get(
                            "solvent_method", None)
                        if solvent_method == "smd" or solvent_method == "sm12":
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
                    if calc["task"]["type"] == "freq" or \
                            calc["task"]["type"] == "frequency":
                        method = calc["input"]["rem"]["method"]
                        basis = calc["input"]["rem"]["basis"]
                        solvent_method  = calc["input"]["rem"].get(
                            "solvent_method", None)
                        if solvent_method == "smd" or solvent_method == "sm12":
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
                    if calc["task"]["type"] == "sp":
                        method = calc["input"]["rem"]["method"]
                        basis = calc["input"]["rem"]["basis"]
                        solvent_method  = calc["input"]["rem"].get(
                            "solvent_method", None)
                        if solvent_method == "smd" or solvent_method == "sm12":
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
        reactant_ids = [r["_id"] for r in reactants]
        product_ids = [p["_id"] for p in products]

        # Get thermo data
        rct_thermo = [get_thermo(r) for r in reactants]
        pro_thermo = [get_thermo(p) for p in products]

        # Compile reaction thermo from reactant and product thermos
        delta_h = sum(p["enthalpy"] for p in pro_thermo) - sum(r["enthalpy"] for r in rct_thermo)
        delta_s = sum(p["entropy"] for p in pro_thermo) - sum(r["entropy"] for r in rct_thermo)
        thermo = {
            "enthalpy": delta_h,
            "entropy": delta_s,
            "t_critical": delta_h / delta_s
        }

        result = {"dir_name": directory,
                  "opt": opt,
                  "freq": freq,
                  "sp": sp,
                  "reactant_ids": reactant_ids,
                  "product_ids": product_ids,
                  "thermo": thermo}

        return result

    def record_data_db(self, directory, opt=None, freq=None, sp=None):
        """
        Record thermo data in thermo collection.

        :param directory: Directory name where the reaction is stored. Right
            now, this is the easiest way to identify the reaction. In the
            future, more sophisticated searching should be used.
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

        collection = self.db.db["thermo"]

        collection.insert_one(self.extract_reaction_data(directory, opt=opt,
                                                         freq=freq, sp=sp))

    def record_data_file(self, directory, filename="thermo.txt", opt=None, freq=None, sp=None):
        """
        Record thermo data in thermo.txt file.

        Note: This function does NOT store the reactant and product IDs

        :param directory: Directory name where the reaction is stored. Right
            now, this is the easiest way to identify the reaction. In the
            future, more sophisticated searching should be used.
        :param filename: File (within directory) where data should be stored.
            By default, it will be stored in thermo.txt.
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
            data = self.extract_reaction_data(directory, opt=opt, freq=freq,
                                              sp=sp)

            file.write("Directory: {}\n".format(data["dir_name"]))
            file.write("Optimization Input: {}\n".format(data["opt"]))
            file.write("Frequency Input: {}\n".format(data["freq"]))
            file.write("Single-Point Input: {}\n".format(data["sp"]))
            file.write("Reaction Enthalpy: {}\n".format(data["thermo"]["enthalpy"]))
            file.write("Reaction Entropy: {}\n".format(data["thermo"]["entropy"]))
            file.write("Critical/Switching Temperature: {}\n".format(data["thermo"]["t_critical"]))