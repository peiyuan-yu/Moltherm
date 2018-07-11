from os import listdir, remove
from os.path import join, isfile, isdir
import operator
import shutil

from bs4 import BeautifulSoup

from pymatgen.io.qchem_io.sets import OptSet, FreqSet, SinglePointSet
from pymatgen.io.babel import BabelMolAdaptor

from moltherm.compute.outputs import QCOutput

"""
TODO list:
    - Add docstrings to util files
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


def extract_id(string):
    return string.split("/")[-1].rstrip(".mol").split("_")[-1]


def remove_copies(base_dir):
    for d in listdir(base_dir):
        if isdir(join(base_dir, d)) and not d.startswith("block"):
            for f in listdir(join(base_dir, d)):
                if f.endswith("_copy") and isfile(join(base_dir, d, f)):
                    remove(join(base_dir, d, f))
                elif f == "copies" and isdir(join(base_dir, d, f)):
                    shutil.rmtree(join(base_dir, d, f))


def mass_copy(base_dir, from_dir, files, directories):
    for file in files:
        for directory in directories:
            if file in listdir(join(base_dir, directory)):
                filename = file + "_copy"
            else:
                filename = file

            shutil.copy(join(base_dir, from_dir, file), join(base_dir, directory, filename))