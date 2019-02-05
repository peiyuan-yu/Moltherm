from __future__ import division, print_function, unicode_literals, absolute_import

import os

from pymatgen.io.qchem.inputs import QCInput

from fireworks import FiretaskBase, explicit_serialize
from fireworks import Firework

from custodian import Custodian
from custodian.qchem.handlers import QChemErrorHandler

from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.firetasks.run_calc import RunQChemCustodian as RunQChemCustodianAtomate
from atomate.qchem.firetasks.parse_outputs import QChemToDb
from atomate.utils.utils import env_chk, get_logger

from moltherm.compute.jobs import QCJob

logger = get_logger(__name__)

__author__ = "Evan Spotte-Smith, Samuel Blau"
__version__ = "0.2"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Beta"
__date__ = "Jan 2019"


class SinglePointFW(Firework):
    def __init__(self,
                 molecule=None,
                 name="single point",
                 qchem_cmd=">>qchem_cmd<<",
                 multimode=">>multimode<<",
                 max_cores=">>max_cores<<",
                 qchem_input_params=None,
                 directory=None,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """

        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Supports env_chk.
            multimode (str): Parallelization scheme, either openmp or mpi. Supports env_chk.
            max_cores (int): Maximum number of cores to parallelize over. Supports env_chk.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                       Basic uses would be to modify the default inputs of the set,
                                       such as dft_rung, basis_set, pcm_dielectric, scf_algorithm,
                                       or max_scf_cycles. See pymatgen/io/qchem/sets.py for default
                                       values of all input parameters. For instance, if a user wanted
                                       to use a more advanced DFT functional, include a pcm with a
                                       dielectric of 30, and use a larger basis, the user would set
                                       qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30,
                                       "basis_set": "6-311++g**"}. However, more advanced customization
                                       of the input is also possible through the overwrite_inputs key
                                       which allows the user to directly modify the rem, pcm, smd, and
                                       solvent dictionaries that QChemDictSet passes to inputs.py to
                                       print an actual input file. For instance, if a user wanted to
                                       set the sym_ignore flag in the rem section of the input file
                                       to true, then they would set qchem_input_params = {"overwrite_inputs":
                                       "rem": {"sym_ignore": "true"}}. Of course, overwrite_inputs
                                       could be used in conjuction with more typical modifications,
                                       as seen in the test_double_FF_opt workflow test.
            db_file (str): Path to file specifying db credentials to place output parsing.
            directory (str): path to directory where calculations should take place. Default is None,
                             meaning current working directory will be used.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        if directory is None:
            input_file = "mol.qin"
            output_file = "mol.qout"
            qclog_file = "mol.qclog"
        else:
            input_file = os.path.join(directory, "mol.qin")
            output_file = os.path.join(directory, "mol.qout")
            qclog_file = os.path.join(directory, "mol.qclog")
        t = []
        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="SinglePointSet",
                input_file=input_file,
                qchem_input_params=qchem_input_params))
        t.append(
            RunQChemCustodian(
                qchem_cmd=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                qclog_file=qclog_file,
                max_cores=max_cores,
                job_type="normal",
                gzipped_output=False))
        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                additional_fields={"task_label": name}))
        super(SinglePointFW, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)


class OptFreqSPFW(Firework):
    def __init__(self, molecule=None,
                 name="opt_freq_sp",
                 qchem_cmd="qchem",
                 multimode="openmp",
                 input_file="mol.qin",
                 output_file="mol.qout",
                 qclog_file="mol.qclog",
                 max_cores=24,
                 max_iterations=1,
                 qchem_input_params=None,
                 sp_params=None,
                 reversed_direction=False,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """
        Performs a QChem workflow with three steps: structure optimization,
        frequency calculation, and single-point calculation.

        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Defaults to qchem.
            multimode (str): Parallelization scheme, either openmp or mpi.
            input_file (str): Name of the QChem input file. Defaults to mol.qin.
            output_file (str): Name of the QChem output file. Defaults to mol.qout.
            qclog_file (str): Name of the QChem log file. Defaults to mol.qclog.
            max_cores (int): Maximum number of cores to parallelize over. Defaults to 24.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                       For example, if you want to change the DFT_rung, you should
                                       provide: {"DFT_rung": ...}. Defaults to None.
            sp_params (dict): Specify inputs for single-point calculation.
            max_iterations (int): Number of perturbation -> optimization -> frequency
                                  iterations to perform. Defaults to 10.
            max_molecule_perturb_scale (float): The maximum scaled perturbation that can be
                                                applied to the molecule. Defaults to 0.3.
            reversed_direction (bool): Whether to reverse the direction of the vibrational
                                       frequency vectors. Defaults to False.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        t = []
        t.append(
            WriteCustomInput(rem=qchem_input_params["rem"],
                             molecule=molecule,
                             opt=qchem_input_params.get("opt", None),
                             pcm=qchem_input_params.get("pcm", None),
                             solvent=qchem_input_params.get("solvent", None),
                             smx=qchem_input_params.get("smx", None),
                             input_file=input_file))
        t.append(
            RunQChemCustodian(
                qchem_cmd=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                qclog_file=qclog_file,
                max_cores=max_cores,
                max_iterations=max_iterations,
                sp_params=sp_params,
                job_type="opt_with_frequency_flattener",
                gzipped_output=False,
                handler_group="no_handler",
                reversed_direction=reversed_direction
            ))

        calc_dir, input_file = os.path.split(input_file)
        output_file = os.path.basename(output_file)

        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                calc_dir=calc_dir,
                additional_fields={
                    "task_label": name,
                    "special_run_type": "opt_freq_sp"
                }))
        super(OptFreqSPFW, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)


class FrequencyFlatteningOptimizeFW(Firework):
    def __init__(self,
                 molecule=None,
                 name="frequency flattening structure optimization",
                 qchem_cmd=">>qchem_cmd<<",
                 multimode=">>multimode<<",
                 max_cores=">>max_cores<<",
                 directory=None,
                 qchem_input_params=None,
                 max_iterations=10,
                 max_molecule_perturb_scale=0.3,
                 reversed_direction=False,
                 db_file=None,
                 parents=None,
                 **kwargs):
        """
        Iteratively optimize the given structure and flatten imaginary frequencies to ensure that
        the resulting structure is a true minima and not a saddle point.

        Args:
            molecule (Molecule): Input molecule.
            name (str): Name for the Firework.
            qchem_cmd (str): Command to run QChem. Supports env_chk.
            multimode (str): Parallelization scheme, either openmp or mpi. Supports env_chk.
            max_cores (int): Maximum number of cores to parallelize over. Supports env_chk.
            directory (str): Location where calculation should take place. Default is current working directory.
            qchem_input_params (dict): Specify kwargs for instantiating the input set parameters.
                                       Basic uses would be to modify the default inputs of the set,
                                       such as dft_rung, basis_set, pcm_dielectric, scf_algorithm,
                                       or max_scf_cycles. See pymatgen/io/qchem/sets.py for default
                                       values of all input parameters. For instance, if a user wanted
                                       to use a more advanced DFT functional, include a pcm with a
                                       dielectric of 30, and use a larger basis, the user would set
                                       qchem_input_params = {"dft_rung": 5, "pcm_dielectric": 30,
                                       "basis_set": "6-311++g**"}. However, more advanced customization
                                       of the input is also possible through the overwrite_inputs key
                                       which allows the user to directly modify the rem, pcm, smd, and
                                       solvent dictionaries that QChemDictSet passes to inputs.py to
                                       print an actual input file. For instance, if a user wanted to
                                       set the sym_ignore flag in the rem section of the input file
                                       to true, then they would set qchem_input_params = {"overwrite_inputs":
                                       "rem": {"sym_ignore": "true"}}. Of course, overwrite_inputs
                                       could be used in conjuction with more typical modifications,
                                       as seen in the test_double_FF_opt workflow test.
            max_iterations (int): Number of perturbation -> optimization -> frequency
                                  iterations to perform. Defaults to 10.
            max_molecule_perturb_scale (float): The maximum scaled perturbation that can be
                                                applied to the molecule. Defaults to 0.3.
            reversed_direction (bool): Whether to reverse the direction of the vibrational
                                       frequency vectors. Defaults to False.
            db_file (str): Path to file specifying db credentials to place output parsing.
            parents ([Firework]): Parents of this particular Firework.
            **kwargs: Other kwargs that are passed to Firework.__init__.
        """

        qchem_input_params = qchem_input_params or {}
        if directory is None:
            input_file = "mol.qin"
            output_file = "mol.qout"
            qclog_file = "mol.qclog"
        else:
            input_file = os.path.join(directory, "mol.qin")
            output_file = os.path.join(directory, "mol.qout")
            qclog_file = os.path.join(directory, "mol.qclog")
        t = []
        t.append(
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="OptSet",
                input_file=input_file,
                qchem_input_params=qchem_input_params))
        t.append(
            RunQChemCustodianAtomate(
                qchem_cmd=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                qclog_file=qclog_file,
                max_cores=max_cores,
                job_type="opt_with_frequency_flattener",
                max_iterations=max_iterations,
                max_molecule_perturb_scale=max_molecule_perturb_scale,
                reversed_direction=reversed_direction,
                gzipped_output=False))
        if directory is None:
            t.append(
                QChemToDb(
                    db_file=db_file,
                    input_file=input_file,
                    output_file=output_file,
                    additional_fields={
                        "task_label": name,
                        "special_run_type": "frequency_flattener"
                    }))
        else:
            t.append(
                QChemToDb(
                    db_file=db_file,
                    calc_dir=directory,
                    input_file="mol.qin",
                    output_file="mol.qout",
                    additional_fields={
                        "task_label": name,
                        "special_run_type": "frequency_flattener"
                    }))

        super(FrequencyFlatteningOptimizeFW, self).__init__(
            t,
            parents=parents,
            name=name,
            **kwargs)


@explicit_serialize
class RunQChemCustodian(FiretaskBase):
    """
    Run QChem using custodian "on rails", i.e. in a simple way that supports most common options.

    Required params:
        qchem_cmd (str): the name of the full executable for running QChem. Note that this is
                         explicitly different from qchem_cmd in RunQChemDirect because it does
                         not include any flags and should only be the call to the executable.
                         Supports env_chk.

    Optional params:
        multimode (str): Parallelization scheme, either openmp or mpi.
        input_file (str): Name of the QChem input file.
        output_file (str): Name of the QChem output file.
        max_cores (int): Maximum number of cores to parallelize over. Defaults to 32.
        qclog_file (str): Name of the file to redirect the standard output to. None means
                          not to record the standard output. Defaults to None.
        suffix (str): String to append to the file in postprocess.
        scratch_dir (str): QCSCRATCH directory. Defaults to "/dev/shm/qcscratch/".
        save_scratch (bool): Whether to save scratch directory contents. Defaults to False.
        save_name (str): Name of the saved scratch directory. Defaults to "default_save_name".
        max_errors (int): Maximum # of errors to fix before giving up (default=5)
        job_type (str): Choose from "normal" (default) and "opt_with_frequency_flattener"
        handler_group (str): Group of handlers to use. See handler_groups dict in the code
                             for the groups and complete list of handlers in each group.
        gzip_output (bool): gzip output (default=T)

        *** Just for opt_with_frequency_flattener ***
        max_iterations (int): Number of perturbation -> optimization -> frequency iterations
                              to perform. Defaults to 10.
        max_molecule_perturb_scale (float): The maximum scaled perturbation that can be
                                            applied to the molecule. Defaults to 0.3.
        reversed_direction (bool): Whether to reverse the direction of the vibrational
                                   frequency vectors. Defaults to False.

        *** Just for opt_with_freq_sp ***
        sp_params (dict): Describes parameters for single-point calculation, if
                          different from opt and freq.

    """
    required_params = ["qchem_cmd"]
    optional_params = [
        "multimode", "input_file", "output_file", "max_cores", "qclog_file",
        "suffix", "scratch_dir", "save_scratch", "save_name", "max_errors",
        "max_iterations", "max_molecule_perturb_scale", "reversed_direction",
        "job_type", "handler_group", "gzipped_output", "sp_params"
    ]

    def run_task(self, fw_spec):

        # initialize variables
        qchem_cmd = env_chk(self["qchem_cmd"], fw_spec)
        multimode = self.get("multimode", "openmp")
        input_file = self.get("input_file", "mol.qin")
        output_file = self.get("output_file", "mol.qout")
        max_cores = self.get("max_cores", 32)
        qclog_file = self.get("qclog_file", "mol.qclog")
        suffix = self.get("suffix", "")
        scratch_dir = env_chk(self.get("scratch_dir"), fw_spec)
        if scratch_dir == None:
            scratch_dir = "/dev/shm/qcscratch/"
        save_scratch = self.get("save_scratch", False)
        save_name = self.get("save_name", "default_save_name")
        max_errors = self.get("max_errors", 5)
        max_iterations = self.get("max_iterations", 10)
        max_molecule_perturb_scale = self.get("max_molecule_perturb_scale",
                                              0.3)
        job_type = self.get("job_type", "normal")
        gzipped_output = self.get("gzipped_output", True)
        sp_params = self.get("sp_params", None)

        handler_groups = {
            "default": [
                QChemErrorHandler(
                    input_file=input_file, output_file=output_file)
            ],
            "no_handler": []
        }

        # construct jobs
        if job_type == "normal":
            jobs = [
                QCJob(
                    qchem_command=qchem_cmd,
                    multimode=multimode,
                    input_file=input_file,
                    output_file=output_file,
                    max_cores=max_cores,
                    qclog_file=qclog_file,
                    suffix=suffix,
                    scratch_dir=scratch_dir,
                    save_scratch=save_scratch,
                    save_name=save_name)
            ]
        elif job_type == "opt_with_frequency_flattener":
            jobs = QCJob.opt_with_frequency_flattener(
                qchem_command=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                qclog_file=qclog_file,
                sp_params=sp_params,
                max_iterations=max_iterations,
                max_molecule_perturb_scale=max_molecule_perturb_scale,
                scratch_dir=scratch_dir,
                save_scratch=save_scratch,
                save_name=save_name,
                max_cores=max_cores)

        else:
            raise ValueError("Unsupported job type: {}".format(job_type))

        # construct handlers
        handlers = handler_groups[self.get("handler_group", "default")]

        c = Custodian(
            handlers,
            jobs,
            max_errors=max_errors,
            gzipped_output=gzipped_output)

        c.run()


@explicit_serialize
class WriteCustomInput(FiretaskBase):
    """
        Writes QChem Input files from custom input sets. This firetask gives the maximum flexibility when trying
        to define custom input parameters.

        required_params:
            qchem_input_custom (dict): Define custom input parameters to generate a qchem input file.
            This should be a dictionary of dictionaries (i.e. {{"rem": {"method": "b3lyp", basis": "6-31*G++", ...}
            Each QChem section should be a key with its own dictionary as the value. For more details on how
            the input should be structured look at pymatgen.io.qchem_io.inputs
            ***  ***

        optional_params:
            molecule (Molecule):
            input_file (str): Name of the QChem input file. Defaults to mol.qin
            write_to_dir (str): Path of the directory where the QChem input file will be written,
            the default is to write to the current working directory
        """

    required_params = ["rem"]
    # optional_params will need to be modified if more QChem sections are added QCInput
    optional_params = [
        "molecule", "opt", "pcm", "solvent", "smx", "input_file", "write_to_dir"
    ]

    def run_task(self, fw_spec):
        input_file = self.get("input_file", "mol.qin")
        # this adds the full path to the input_file
        if "write_to_dir" in self:
            input_file = os.path.join(self["write_to_dir"], input_file)
        # these if statements might need to be reordered at some point
        if "molecule" in self:
            molecule = self["molecule"]
        elif fw_spec.get("prev_calc_molecule"):
            molecule = fw_spec.get("prev_calc_molecule")
        else:
            raise KeyError(
                "No molecule present, add as an optional param or check fw_spec"
            )
        # in the current structure there needs to be a statement for every optional QChem section
        # the code below defaults the section to None if the variable is not passed
        opt = self.get("opt", None)
        pcm = self.get("pcm", None)
        solvent = self.get("solvent", None)
        smx = self.get("smx", None)

        qcin = QCInput(
            molecule=molecule,
            rem=self["rem"],
            opt=opt,
            pcm=pcm,
            solvent=solvent,
            smx=smx)
        qcin.write_file(input_file)

