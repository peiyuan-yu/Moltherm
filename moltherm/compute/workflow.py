from moltherm.compute.scrape_compounds import ReaxysScraper
from moltherm.react.molecule import OBMolecule

from pymatgen.core.structure import Molecule, FunctionalGroups
from pymatgen.analysis.graphs import MoleculeGraph

from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.io.qchem_io.outputs import QCOutput
from pymatgen.io.qchem_io.sets import OptSet, FreqSet

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


def generate_opt_input(molfile, qinfile):
    """
    Generates a QChem input file from Molecule after conformer search.

    :param molfile: Absolute path to the input file (.mol, .sdf, etc.)
    :param qinfile: Absolute path to the output file (.in)
    :return:

    """

    obmol = OBMolecule.from_file(molfile)
    # OBMolecule does not contain pymatgen Molecule information
    # So, we need to wrap the obmol in a BabelMolAdapter and extract
    obmol = obmol.confm_search(make_3d=True, add_hydrogens=True).obmol
    mol = BabelMolAdaptor(obmol).pymatgen_mol

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

    if len(output.data["molecule_from_optimized_geometry"]) > 0:
        mol = output.data["molecule_from_optimized_geometry"]
    else:
        mol = output.data["molecule_from_last_geometry"]

    qcinput = FreqSet(mol)

    qcinput.write_file(qinfile)