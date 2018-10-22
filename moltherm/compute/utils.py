from os import listdir, remove, rename
from os.path import join, isfile, isdir
import operator
import shutil

from bs4 import BeautifulSoup

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.io.babel import BabelMolAdaptor


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


def get_reactions_common_solvent(base_dir, solvents, outdir=None):
    """
    Identify all reactions from a set of reactions which share a solvent.

    :param base_dir: Directory to begin search in.
    :param solvents: list of strings representing solvents
    :param outdir: Directory to put all reactions with the common solvent. If
        None, do not copy relevant files.
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

            this_solv = parsed.find("solvents").text.split(",")
            this_solv = [s.lower() for s in this_solv]

            for solvent in solvents:

                if solvent.lower() in this_solv and rct_dir not in [d['dir'] for d in common_solvent]:
                    rct_names = [d.text for d in parsed.find_all('rctname')]
                    rct_ids = [d.text for d in parsed.find_all('rctid')]
                    pro_name = [d.text for d in parsed.find_all('proname')]
                    pro_ids = [d.text for d in parsed.find_all('proid')]

                    common_solvent.append({'dir': rct_dir,
                                           'rct_names': rct_names,
                                           'rct_ids': rct_ids,
                                           'pro_name': pro_name,
                                           'pro_ids': pro_ids})

    if outdir is not None:
        num_copied = 0
        for rct_dir in common_solvent:
            shutil.copytree(join(base_dir, rct_dir), join(outdir, rct_dir))
            num_copied += 1

        print("{} reactions with solvents {}".format(str(num_copied), solvents))

    else:
        return common_solvent


def associate_qchem_to_mol(base_dir, directory):
    """
    Assign all .in and .out files in a directory to one of the .mol files in that
    directory, based on the non-H atoms in those molecules.

    :param directory:
    :return:
    """

    base_path = join(base_dir, directory)

    mol_files = [f for f in listdir(base_path) if isfile(join(base_path, f))
                 and f.endswith(".mol") and not f.startswith(".")]
    # Note: This will catch .in and .out files for incomplete computations
    # TODO: What's the best way to filter these out?
    in_files = [f for f in listdir(base_path) if isfile(join(base_path, f))
                and ".in" in f and not f.startswith("atomate")]
    out_files = [f for f in listdir(base_path) if isfile(join(base_path, f))
                 and ".out" in f and not f.startswith("atomate")]

    mapping = {mol: {"in": [], "out": []} for mol in mol_files}

    for file in in_files:
        qcin = QCInput.from_file(join(base_path, file))
        file_mol = qcin.molecule
        # Remove H because mol files may not begin with H included
        file_species = [str(s) for s in file_mol.species if str(s) != "H"]

        for mf in mol_files:
            mol_mol = Molecule.from_file(join(base_path, mf))
            mol_species = [str(s) for s in mol_mol.species if str(s) != "H"]
            # Preserve initial order because that gives a better guarantee
            # That the two are actually associated
            if mol_species == file_species:
                mapping[mf]["in"].append(file)
                break

    for file in out_files:
        qcout = QCOutput(join(base_path, file))
        file_mol = qcout.data["initial_molecule"]
        file_species = [str(s) for s in file_mol.species if str(s) != "H"]

        for mf in mol_files:
            mol_mol = Molecule.from_file(join(base_path, mf))
            mol_species = [str(s) for s in mol_mol.species if str(s) != "H"]

            # Preserve initial order because that gives a better guarantee
            # That the two are actually associated
            if mol_species == file_species:
                mapping[mf]["out"].append(file)
                break

    return mapping


def extract_id(string):
    """
    Extract unique molecule ID from a filepath.

    :param string: Path or filename.
    :return: str representing unique ID
    """

    return string.split("/")[-1].replace(".mol", "").split("_")[-1]
