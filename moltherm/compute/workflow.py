from os import listdir
from os.path import join, isfile, isdir

from fireworks import Workflow, LaunchPad

from atomate.qchem.database import QChemCalcDb

from moltherm.compute.fireworks import OptFreqSPFW, SinglePointFW
from moltherm.compute.inputs import QCInput
from moltherm.compute.outputs import QCOutput
from moltherm.compute.utils import get_molecule, extract_id, associate_qchem_to_mol

__author__ = "Evan Spotte-Smith"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Alpha"
__date__ = "June 2018"


"""
TODO list:
    - Figure out how to query with pymatgen-db and pymongo
"""


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

    def get_unfinished_jobs(self, sp_params, name_pre="single_point", dirs=None, max_cores=24):
        """
        Look for jobs where optimization and frequency calculations have
        successfully completed, but single-point has not. Then, for these cases,
        construct a workflow which will only run the sp job.

        :param sp_params: dict containing input parameters for single-point job
        :param name_pre: str representing prefix for all jobs.
        :param dirs: list of subdirectories to check for unfinished jobs.
            Default None, meaning that all subdirectories will be checked.
        :param max_cores: max_cores (int): Maximum number of cores to
            parallelize over. Defaults to 24.
        :return:
        """

        if not self.subdirs:
            raise RuntimeError("Cannot run get_reaction_set_workflow();"
                               "Need reactions components to be isolated in"
                               "different subdirectories.")

        fws = []

        all_dirs = [d for d in listdir(self.base_dir)
                    if isdir(join(self.base_dir, d))]

        molecules_cleared = []

        appropriate_dirs = all_dirs

        if dirs is not None:
            appropriate_dirs = [d for d in appropriate_dirs if d in dirs]

        for d in appropriate_dirs:
            path = join(self.base_dir, d)
            file_map = associate_qchem_to_mol(self.base_dir, d)

            for key, values in file_map.items():
                mol_id = extract_id(key)

                if mol_id in molecules_cleared:
                    continue

                freq_complete = False
                sp_complete = False

                in_files = values["in"]
                out_files = values["out"]

                # Check if this molecule has finished freq, sp
                # If there is no sp output file, or if the sp output file did
                # not complete, then we may proceed
                for out_file in out_files:
                    if "freq" in out_file:
                        freq_out = QCOutput(join(path, out_file))

                        if freq_out.data.get("completion", []):
                           freq_complete = True
                    elif "sp" in out_file:
                        sp_out = QCOutput(join(path, out_file))

                        if sp_out.data.get("completion", []):
                            sp_complete = True

                if freq_complete and not sp_complete:
                    # Check if there is already an sp input file
                    freq_in_file = None
                    sp_in_file = None

                    for in_file in in_files:
                        if "freq" in in_file:
                            freq_in_file = in_file
                        if "sp" in in_file:
                            sp_in_file = in_file

                    if sp_in_file is None:
                        if freq_in_file is None:
                            # We could parse output files to get previous input
                            # information, but we should try to keep all input
                            # files in the same directory
                            continue
                        else:
                            infile = join(path, key.replace(".mol", "") + ".in")
                            outfile = join(path, key.replace(".mol", "") + ".out")
                            qclogfile = join(path, key.replace(".mol", "") + ".qclog")

                            freq_in_file = QCInput.from_file(join(path,
                                                                  freq_in_file))
                            mol = freq_in_file.molecule

                            fw = SinglePointFW(molecule=mol,
                                               name="{}: {}/{}".format(name_pre, d, mol_id),
                                               qchem_cmd="qchem -slurm",
                                               multimode="openmp",
                                               input_file=infile,
                                               output_file=outfile,
                                               qclog_file=qclogfile,
                                               max_cores=max_cores,
                                               sp_params=sp_params)

                            fws.append(fw)
                            molecules_cleared.append(mol_id)

                    else:
                        infile = sp_in_file
                        outfile = join(path, infile.replace(".in", ".out"))
                        qclogfile = join(path, infile.replace(".in", ".qclog"))

                        sp_in_file = QCInput.from_file(join(path, infile))

                        mol = sp_in_file.molecule

                        fw = SinglePointFW(molecule=mol,
                                           name="{}: {}/{}".format(name_pre, d, mol_id),
                                           qchem_cmd="qchem -slurm",
                                           multimode="openmp",
                                           input_file=join(path, infile),
                                           output_file=outfile,
                                           qclog_file=qclogfile,
                                           input_exists=True,
                                           max_cores=max_cores,
                                           sp_params=sp_params)

                        fws.append(fw)
                        molecules_cleared.append(mol_id)

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