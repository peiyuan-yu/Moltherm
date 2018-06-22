from os import listdir
from os.path import join, isfile, isdir
import operator

from pymatgen.io.babel import BabelMolAdaptor

from atomate.qchem.database import QChemCalcDb
from atomate.qchem.drones import QChemDrone
from atomate.qchem.fireworks.core import FrequencyFlatteningOptimizeFW

from pymatgen.core.structure import

from atomate.utils.utils import load_class
from fireworks import FiretaskBase, explicit_serialize
from pymatgen.io.qchem_io.inputs import QCInput
from pymatgen.io.qchem_io.outputs import QCOutputs
from pymatgen.io.qchem_io.sets import OptSet, SinglePointSet, FreqSet
from fireworks import Firework

from custodian.qchem.new_jobs import QCJob
from custodian.qchem.new_handlers import QChemErrorHandler


class OptFreqSPFW(Firework):
    def __init(self, molecule=None,
               name="opt-freq-sp",
               qchem_cmd="qchem",
               multimode="openmp",
               input_file_opt="mol_opt.qin",
               input_file_sp="mol_sp.qin",
               output_file_opt="mol_opt.qout",
               output_file_sp="mol_sp.qout",
               max_cores=32,
               qchem_input_params=None,
               max_iterations=1,
               max_molecule_perturb_scale=0.3,
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
                input_file=input_file_opt,
                qchem_input_params=qchem_input_params["opt"]))
        t.append(
            RunQChemCustodian(
                qchem_cmd=qchem_cmd,
                multimode=multimode,
                input_file=input_file_opt,
                output_file=output_file_opt,
                max_cores=max_cores,
                job_type="opt_with_frequency_flattener",
                max_iterations=max_iterations,
                max_molecule_perturb_scale=max_molecule_perturb_scale,
                reversed_direction=reversed_direction))
        t.append(
            WriteItemFromIOSet(
                qchem_input_set="SinglePointSet",
                input_file=input_file_sp,
                qchem_input_params=qchem_input_params["sp"]))
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
