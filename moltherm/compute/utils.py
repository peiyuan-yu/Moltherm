from os import listdir
from os.path import join, isfile, basename
import operator
from difflib import SequenceMatcher
from statistics import mean
import itertools

from bs4 import BeautifulSoup

import numpy as np
import networkx.algorithms.isomorphism as iso

from pymatgen.core.structure import Molecule
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.io.babel import BabelMolAdaptor

__author__ = "Evan Spotte-Smith"
__version__ = "0.2"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Beta"
__date__ = "July 2018"


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

    meta = [f for f in listdir(base_dir) if isfile(join(base_dir, f)) and f.endswith(".xml")]

    solvent_occurrence = {}

    for f in meta:
        with open(join(base_dir, f), "r") as file:
            parsed = BeautifulSoup(file.read(), "lxml-xml")

            solvents = parsed.find("solvents").text.split(" || ")

            for solvent in solvents:
                current_value = solvent_occurrence.get(solvent, 0)
                solvent_occurrence[solvent] = current_value + 1

    return sorted(solvent_occurrence.items(), key=operator.itemgetter(1))


def get_reactions_common_solvent(base_dir, solvents):
    """
    Identify all reactions from a set of reactions which share a solvent.

    :param base_dir: Directory to begin search in.
    :param solvents: list of strings representing solvents
    :return: list of dicts with IDS and names for compounds in relevant
        reactions.
    """

    meta = [f for f in listdir(base_dir) if
            isfile(join(base_dir, f)) and f.endswith(".xml")]

    common_solvent = []
    common_ids = []

    for f in meta:
        with open(join(base_dir, f), "r") as file:
            rxn_id = f.replace(".xml", "")
            parsed = BeautifulSoup(file.read(), "lxml-xml")

            this_solv = parsed.find("solvents").text.split(" || ")
            this_solv = [s.lower() for s in this_solv]

            for solvent in solvents:

                if solvent.lower() in this_solv and rxn_id not in common_ids:
                    rct_names = [d.text for d in parsed.find_all('rctname')]
                    rct_ids = [d.text for d in parsed.find_all('rctid')]
                    pro_name = [d.text for d in parsed.find_all('proname')]
                    pro_ids = [d.text for d in parsed.find_all('proid')]

                    common_solvent.append({'rxn_id': rxn_id,
                                           'rct_names': rct_names,
                                           'rct_ids': rct_ids,
                                           'pro_name': pro_name,
                                           'pro_ids': pro_ids})

                    common_ids.append(rxn_id)

    return common_solvent


def extract_id(string):
    """
    Extract unique molecule ID from a filepath.

    :param string: Path or filename.
    :return: str representing unique ID
    """

    return basename(string).replace(".mol", "").split("_")[-1]


def get_smiles(base_path, molecules, extra=False):
    """
    Returns SMILES strings for a set of molecules

    :param base_path: Directory to begin search in.
    :param molecules: List of strings representing molecule IDs.
    :param extra: If True (default False), give molecule ID with SMILES
    :return: list of SMILES strings
    """

    smiles = []

    if base_path is not None:
        for mol in molecules:
            path = join(base_path, mol)
            molfile = join(path, "{}.mol".format(mol))

            obmol = BabelMolAdaptor.from_file(molfile, file_format="mol")

            smi = obmol.pybel_mol.write("smi").replace("\t\n", "")

            if extra:
                smiles.append((smi, mol))
            else:
                smiles.append(smi)

        return smiles
    else:
        raise ValueError("No path given.")


def associate_qchem_to_mol(base_dir, directory):
    """
    Assign all .in and .out files in a directory to one of the .mol files in that
    directory, based on the non-H atoms in those molecules.

    :param base_dir: Directory to begin search in.
    :param directory: Subdirectory of interes
    :return:
    """

    base_path = join(base_dir, directory)

    all_files = [e for e in listdir(base_path) if not e.startswith(".") and not e.endswith("swp")]

    mol_files = [f for f in all_files  if isfile(join(base_path, f))
                 and f.endswith(".mol") and not f.startswith(".")]
    # Note: This will catch .in and .out files for incomplete computations
    in_files = [f for f in all_files if isfile(join(base_path, f))
                and (".in" in f or ".qin" in f)
                and not f.startswith("atomate")]
    out_files = [f for f in all_files if isfile(join(base_path, f))
                 and (".out" in f or ".qout" in f)
                 and not f.startswith("atomate")]

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


def calculate_solubility(vp, delta_g_solv, temp_kelvin=298):
    """
    Get the solubility of an arbitrary material in an arbitrary solvent.

    Based on the method described in:

    Thompson, J.D., Cramer, C.J. and Truhlar, D.G., 2003. Predicting aqueous
    solubilities from aqueous free energies of solvation and experimental or
    calculated vapor pressures of pure substances.
    The Journal of chemical physics, 119(3), pp.1661-1670.

    :param vp: Vapor pressure of the solute molecules under standard temperature
    :param delta_g_solv: Free energy of solvation, given by (for instance) the
        SMD universal solvation method.
    :param temp_kelvin: Temperature in Kelvin

    :return: S, the equilibrium solubility in units of mols per unit volume
    """

    # Ideal gas pressure at 1 molar concentration and 298K
    p_std = 2477396.2
    R = 8.314

    S = vp / p_std * np.exp(-1 * delta_g_solv / (R * temp_kelvin))

    return S


def map_atoms_reaction(reactants, product):
    """
    Create a mapping of atoms between a set of reactant Molecules and a product
    Molecule.

    :param reactants: list of Molecule objects representing the reaction
        reactants
    :param product: Molecule object representing the reaction product

    NOTE: This currently only works with one product

    :return: dict with 'reactants' and 'product' keys
    """

    def get_ranked_atom_dists(mol):
        dist_matrix = mol.distance_matrix

        result = dict()
        for num, row in enumerate(dist_matrix):
            ranking = np.argsort(row)
            # The first member will always be the atom itself, which should be excluded
            result[num] = ranking[1:]
        return result

    rct_mgs = list()
    for reactant in reactants:
        rct_mgs.append(MoleculeGraph.with_local_env_strategy(reactant, OpenBabelNN(), reorder=False,
                                                             extend_structure=False))
    pro_mg = MoleculeGraph.with_local_env_strategy(product, OpenBabelNN(), reorder=False,
                                                   extend_structure=False)

    # Next, try to construct isomorphisms between reactant and product
    nm = iso.categorical_node_match("specie", "ERROR")
    # Prefer undirected graphs
    rct_graphs = [rct_mg.graph.to_undirected() for rct_mg in rct_mgs]
    pro_graph = pro_mg.graph.to_undirected()

    pro_dists = get_ranked_atom_dists(pro_mg.molecule)
    ranking_by_reactant = list()
    for e, rct_graph in enumerate(rct_graphs):

        meta_iso = {e: set() for e in range(len(pro_mg.molecule))}
        matcher = iso.GraphMatcher(pro_graph, rct_graph, node_match=nm)

        # Compile all isomorphisms
        isomorphisms = [i for i in matcher.subgraph_isomorphisms_iter()]
        for isomorphism in isomorphisms:
            for pro_node, rct_node in isomorphism.items():
                meta_iso[pro_node].add(rct_node)

        meta_iso = {k: v for (k, v) in meta_iso.items() if v != set()}

        # Determine which nodes need to be checked
        disputed_nodes = set()
        for pro_node, rct_nodes in meta_iso.items():
            if len(rct_nodes) > 1:
                disputed_nodes.add(pro_node)

        average_ratios = list()
        rct_dists = get_ranked_atom_dists(reactants[e])
        for isomorphism in isomorphisms:
            ratios = list()

            for node in disputed_nodes:
                if node in isomorphism:
                    rct_dist = rct_dists[isomorphism[node]]
                    pro_dist_old = pro_dists[node]

                    pro_dist = list()
                    for n in pro_dist_old:
                        if n in isomorphism:
                            pro_dist.append(isomorphism[n])

                    matcher = SequenceMatcher(None, rct_dist, pro_dist)
                    ratios.append(matcher.ratio())

            average_ratios.append(mean(ratios))

        # Rank isomorphisms by the average SequenceMatch ratio
        ranking_by_reactant.append([isom for _, isom in sorted(zip(average_ratios, isomorphisms),
                                                         key=lambda pair: pair[0])])

    combinations = list(itertools.product(*ranking_by_reactant))

    mapping = None
    for combination in combinations:
        index_total = 0
        this_mapping = dict()

        for part in combination:
            for key, value in part.items():
                this_mapping[key] = value + index_total
            index_total += len(part.keys())

        # If this is a valid mapping - all atoms accounted for
        if len(this_mapping) == len(product):
            mapping = this_mapping
            break

    return mapping