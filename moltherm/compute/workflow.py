from os import listdir
from os.path import join, isfile, isdir

from moltherm.react.molecule import OBMolecule

from pymatgen.core.structure import Molecule, FunctionalGroups
from pymatgen.analysis.graphs import MoleculeGraph

from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.io.qchem_io.outputs import QCOutput
from pymatgen.io.qchem_io.sets import OptSet, FreqSet

from fireworks import Firework, Workflow, LaunchPad, FWorker

from custodian.qchem.new_handlers import QChemErrorHandler
from custodian.qchem.new_jobs import QCJob

from atomate.qchem.database import QChemCalcDb
from atomate.qchem.drones import QChemDrone
from atomate.qchem.firetasks.parse_outputs import *
from atomate.qchem.firetasks.run_calc import *
from atomate.qchem.firetasks.write_inputs import *
from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW

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
"""


def get_molecule(molfile):
    """
    Create pymatgen Molecule object from molecule data file.

    In addition to parsing the input, this function also performs a conformer
    search to get a reasonable starting structure.

    :param molfile: Absolute path to structure file (.mol, .sdf, etc.)
    :return: Molecule.
    """

    obmol = OBMolecule.from_file(molfile)
    # OBMolecule does not contain pymatgen Molecule information
    # So, we need to wrap the obmol in a BabelMolAdapter and extract
    obmol = obmol.confm_search(make_3d=True, add_hydrogens=True).obmol
    return BabelMolAdaptor(obmol).pymatgen_mol


def generate_opt_input(molfile, qinfile):
    """
    Generates a QChem input file from Molecule after conformer search.

    :param molfile: Absolute path to the input file (.mol, .sdf, etc.)
    :param qinfile: Absolute path to the output file (.in)
    :return:

    """
    mol = get_molecule(molfile)

    # Right now, just use all defaults.
    qcinput = OptSet(mol)

    qcinput.write_file(qinfile)


def generate_freq_input(qoutfile, qinfile):
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

    qcinput = FreqSet(mol)

    qcinput.write_file(qinfile)


class MolTherm:
    """
    This workflow contains all functionality needed to perform Atomate
    workflows on molecular data.
    """

    def __init__(self, base_dir, subdirs=False, indexing=False,
                 reactant_pre="rct_", product_pre="pro_", with_freq=True,
                 db_file="db.json"):
        """
        :param base_dir: Directory where input and output data should be stored.
        :param subdirs: Is data all stored in one directory (False), or is it
        separated among subdirectories (True)?
        :param indexing: If True, then assume that files (or subdirectories, if
        applicable) are indexed by number.
        :param reactant_pre: Prefix for reactant files.
        :param product_pre: Prefix for product files.
        :param with_freq: Should the workflow only optimize the structure of
        the molecules (False), or should it also perform frequency calculations
        (True)?
        :param db_file: Path to database config file.
        """

        self.base_dir = base_dir
        self.subdirs = subdirs
        self.indexing = indexing
        self.reactant_pre = reactant_pre
        self.product_pre = product_pre
        self.with_freq = with_freq
        self.db_file = db_file

    def get_reaction_enthalpy_files(self, path=None, index=None):
        pass

    def get_reaction_enthalpy_db(self, db_file):
        pass

    def record_data(self):
        pass

    def get_single_reaction_workflow(self, path=None, filenames=None):
        """
        Runs an atomate workflow on a single reaction.

        :param path: Specified (sub)path in which to run the reaction. By default,
        this is None, and the Fireworks will run in self.base_dir
        :param filenames: Specified files within the path (if self.base_dir or
        a subdirectory) that should be considered a part of this reaction. If
        None, assume all files in the directory are to be involved.
        :return: Workflow
        """

        fws = []

        if path is not None and self.subdirs:
            base_path = join(self.base_dir, path)
        else:
            base_path = self.base_dir

        if filenames:
            rcts = [f for f in filenames if f.startswith(self.reactant_pre)]
            pros = [f for f in filenames if f.startswith(self.product_pre)]
        else:
            files = [f for f in listdir(base_path) if isfile(join(base_path, f))]
            rcts = [f for f in files if f.startswith(self.reactant_pre)]
            pros = [f for f in files if f.startswith(self.product_pre)]

        for i, rct in enumerate(rcts):
            mol = get_molecule(join(base_path, rct))

            fw = FrequencyFlatteningOptimizeFW(molecule=mol,
                                               name=("opt+freq: " + rct),
                                               input_file=(self.reactant_pre + str(i) + ".in"),
                                               output_file=(self.reactant_pre + str(i) + ".out"),
                                               db_file=self.db_file)

            fws.append(fw)

        for pro in pros:
            mol = get_molecule(join(base_path, pro))

            fw = FrequencyFlatteningOptimizeFW(molecule=mol,
                                               name=("opt+freq: " + pro),
                                               input_file=(self.product_pre + str(i) + ".in"),
                                               output_file=(self.product_pre + str(i) + ".out"),
                                               db_file=self.db_file)

            fws.append(fw)

        return Workflow(fws)

    def get_n_reaction_workflow(self, n):
        pass
