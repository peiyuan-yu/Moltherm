from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.firetasks.parse_outputs import QChemToDb

from atomate.utils.utils import env_chk, get_logger
from fireworks import FiretaskBase, explicit_serialize
from fireworks import Firework

from custodian import Custodian
from custodian.qchem.new_handlers import QChemErrorHandler

from moltherm.compute.drones import MolThermDrone
from moltherm.compute.outputs import QCOutput
from moltherm.compute.jobs import QCJob


class OptFreqSPFW(Firework):
    def __init(self, molecule=None,
               name="opt_freq_sp",
               qchem_cmd="qchem",
               multimode="openmp",
               input_file="mol.qin",
               output_file="mol.qout",
               max_cores=64,
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
            max_cores (int): Maximum number of cores to parallelize over. Defaults to 32.
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
            WriteInputFromIOSet(
                molecule=molecule,
                qchem_input_set="OptSet",
                input_file=input_file,
                qchem_input_params=qchem_input_params["opt"]))
        t.append(
            RunQChemCustodian(
                qchem_cmd=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                max_cores=max_cores,
                sp_params=sp_params,
                job_type="opt_freq_sp",
                reversed_direction=reversed_direction,
            ))
        t.append(
            QChemToDb(
                db_file=db_file,
                input_file=input_file,
                output_file=output_file,
                additional_fields={
                    "task_label": name,
                    "special_run_type": "opt_freq_sp"
                }))
        super(OptFreqSPFW, self).__init__(
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
                max_iterations=max_iterations,
                max_molecule_perturb_scale=max_molecule_perturb_scale,
                scratch_dir=scratch_dir,
                save_scratch=save_scratch,
                save_name=save_name,
                max_cores=max_cores)

        elif job_type == "opt_freq_sp":
            jobs = QCJob.opt_with_freq_sp(
                qchem_command=qchem_cmd,
                multimode=multimode,
                input_file=input_file,
                output_file=output_file,
                qclog_file=qclog_file,
                sp_params=sp_params,
                max_cores=max_cores,
                scratch_dir=scratch_dir,
                save_scratch=save_scratch,
                save_name=save_name)


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
