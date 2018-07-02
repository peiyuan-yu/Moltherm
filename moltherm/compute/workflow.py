from os import listdir
from os.path import join, isfile, isdir
import operator
import shutil

from bs4 import BeautifulSoup

from pymatgen.io.qchem_io.sets import OptSet, FreqSet, SinglePointSet
from pymatgen.io.babel import BabelMolAdaptor

from fireworks import Workflow, LaunchPad

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


class MolTherm:
    """
    This workflow contains all functionality needed to perform Atomate
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

    def get_reaction_thermo_files(self, path=None):
        """
        Naively scrape thermo data from QChem output files.

        :param path: Path to a subdirectory. Will be ignored if self.subdirs
        is False.

        :return: dict {prop: value}, where properties are enthalpy, entropy.
        """
        thermo_data = {}

        rct_thermo = {"enthalpy": 0, "entropy": 0}
        pro_thermo = {"enthalpy": 0, "entropy": 0}

        if self.subdirs and path is not None:
            base_path = join(self.base_dir, path)
        else:
            base_path = self.base_dir

        files = [f for f in listdir(base_path) if isfile(join(base_path, f))
                 and f.endswith(".out")]
        rcts = [f for f in files if f.startswith(self.reactant_pre)]
        pros = [f for f in files if f.startswith(self.product_pre)]

        for rct in rcts:
            qcout = QCOutput(join(base_path, rct))
            rct_thermo["enthalpy"] += qcout.data["enthalpy"]
            rct_thermo["entropy"] += qcout.data["entropy"]

        for pro in pros:
            qcout = QCOutput(join(base_path, pro))
            pro_thermo["enthalpy"] += qcout.data["enthalpy"]
            pro_thermo["entropy"] += qcout.data["entropy"]

        # Generate totals as ∆H = H_pro - H_rct, ∆S = S_pro - S_rct
        thermo_data["enthalpy"] = pro_thermo["enthalpy"] - rct_thermo["enthalpy"]
        thermo_data["entropy"] = pro_thermo["entropy"] - pro_thermo["entropy"]

        return thermo_data

    def get_reaction_thermo_db(self, db_file):
        pass

    def record_data(self):
        pass

    def get_single_reaction_workflow(self, path=None, filenames=None,
                                     max_cores=64,
                                     qchem_input_params=None,
                                     sp_params=None):
        """
        Generates a Fireworks Workflow to find the structures and energies of
        the reactants and products of a single reaction.

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

        if path is not None:
            fw_pre = path
        else:
            fw_pre = "opt_freq_sp_"
            path = ""

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

            infile = join(self.base_dir, path, self.reactant_pre + str(i) + ".in")
            outfile = join(self.base_dir, path, self.reactant_pre + str(i) + ".out")

            fw = OptFreqSPFW(molecule=mol,
                             name=(fw_pre + " : " + rct),
                             qchem_cmd="qchem -slurm",
                             input_file=infile,
                             output_file=outfile,
                             qclog_file=join(self.base_dir, path, self.reactant_pre + str(i) + ".qclog"),
                             max_cores=max_cores,
                             qchem_input_params=qchem_input_params,
                             sp_params=sp_params,
                             db_file=self.db_file)

            fws.append(fw)

        for i, pro in enumerate(pros):
            mol = get_molecule(join(base_path, pro))

            infile = join(self.base_dir, path, self.product_pre + str(i) + ".in")
            outfile = join(self.base_dir, path, self.product_pre + str(i) + ".out")

            fw = OptFreqSPFW(molecule=mol,
                             name=(fw_pre + " : " + pro),
                             qchem_cmd="qchem -slurm",
                             input_file=infile,
                             output_file=outfile,
                             qclog_file=join(self.base_dir, path, self.reactant_pre + str(i) + ".qclog"),
                             max_cores=max_cores,
                             qchem_input_params=qchem_input_params,
                             sp_params=sp_params,
                             db_file=self.db_file)

            fws.append(fw)

        return Workflow(fws)

    def get_reaction_set_workflow(self, max_cores=64,
                                  qchem_input_params=None,
                                  sp_params=None):
        """Generates a Fireworks Workflow to find the structures and energies of
        the reactants and products of a single reaction.

        Note: as written now, this function will only work if self.subdirs is
        True; that is, only if each reaction is in a separate subdirectory.
        Later additions could allow for some other means of specifying the
        separate reactions within a single directory.

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

        fws = []

        dirs = [d for d in listdir(self.base_dir) if isdir(join(self.base_dir, d))]

        for d in dirs:
            path = join(self.base_dir, d)
            files = [f for f in listdir(path) if isfile(join(path, f))]
            rcts = [f for f in files if f.startswith(self.reactant_pre)]
            pros = [f for f in files if f.startswith(self.product_pre)]

            for i, rct in enumerate(rcts):
                mol = get_molecule(join(self.base_dir, d, rct))

                infile = join(self.base_dir, self.reactant_pre + str(i) + ".in")
                outfile = join(self.base_dir,
                               self.reactant_pre + str(i) + ".out")

                fw = OptFreqSPFW(molecule=mol,
                                 name=("opt_freq_sp_: " + rct),
                                 qchem_cmd="qchem -slurm",
                                 input_file=infile,
                                 output_file=outfile,
                                 max_cores=max_cores,
                                 qchem_input_params=qchem_input_params,
                                 sp_params=sp_params,
                                 db_file=self.db_file)

                fws.append(fw)

            for i, pro in enumerate(pros):
                mol = get_molecule(join(self.base_dir, d, pro))

                infile = join(self.base_dir, self.product_pre + str(i) + ".in")
                outfile = join(self.base_dir,
                               self.product_pre + str(i) + ".out")

                fw = OptFreqSPFW(molecule=mol,
                                 name=("opt_freq_sp_: " + rct),
                                 qchem_cmd="qchem -slurm",
                                 input_file=infile,
                                 output_file=outfile,
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
