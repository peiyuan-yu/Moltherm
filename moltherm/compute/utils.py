from os import listdir, remove, rename
from os.path import join, isfile, isdir
import operator
import shutil

from bs4 import BeautifulSoup

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.io.qchem.sets import OptSet, FreqSet, SinglePointSet
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
                 and f.endswith(".mol")]
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


def remove_copies(base_dir):
    """
    Move through subdirectories in base_dir to eliminate all *_copy files.

    :param base_dir: str representing a path to a directory.
    :return:
    """

    for d in listdir(base_dir):
        if isdir(join(base_dir, d)) and not d.startswith("block"):
            for f in listdir(join(base_dir, d)):
                if f.endswith("_copy") and isfile(join(base_dir, d, f)):
                    remove(join(base_dir, d, f))
                elif f == "copies" and isdir(join(base_dir, d, f)):
                    shutil.rmtree(join(base_dir, d, f))


def mass_copy(base_dir, from_dir, files, directories):
    """
    Copy files from one directory into many directories.

    :param base_dir: str representing a path to a the base directory (all other
        directories should be subdirectories of this base directory)
    :param from_dir: str representing a path to the original subdirectory
    :param files: list of filenames within from_dir
    :param directories: list of directories to copy into
    :return:
    """
    for file in files:
        for directory in directories:
            if file in listdir(join(base_dir, directory)):
                filename = file + "_copy"
            else:
                filename = file

            shutil.copy(join(base_dir, from_dir, file), join(base_dir, directory, filename))


def rename_products(base_dir, product_prefix, exclude=None):
    """
    Add additional information to files associated with products.

    This function exists because with the currently used id scheme for
    functional-group substituted molecules, there could be ambiguity between
    products with similar substitutions.

    :param base_dir: str representing a path to a the base directory (all other
        directories should be subdirectories of this base directory)
    :param product_prefix: str which all product files should begin with
    :param exclude: list of subdirectories which should not be changed. Default
        is None.
    :return:
    """

    dirs = [d for d in listdir(base_dir) if isdir(join(base_dir, d)) and not
            d.startswith("block")]
    if exclude is not None:
        dirs = [d for d in dirs if d not in exclude]

    for d in dirs:
        path = join(base_dir, d)

        product_files = [f for f in listdir(path) if isfile(join(path, f)) and
                         f.startswith(product_prefix)]

        if d.endswith("furan"):
            id_suffix = "x103221"
        elif d.endswith("maleimide"):
            id_suffix = "x106910"
        elif d.endswith("acrylonitrile"):
            id_suffix = "x605310"
        else:
            id_suffix = ""

        for file in product_files:
            before_suffix = file.split(".")[0]
            new_name = file.replace(before_suffix, before_suffix+id_suffix)
            rename(join(path, file), join(path, new_name))
