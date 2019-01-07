# coding: utf-8

import os
from os.path import abspath, dirname, join
import shutil
import unittest
from unittest.mock import patch

from pymatgen.core.structure import Molecule
from pymatgen.io.qchem.inputs import QCInput
from pymatgen.io.qchem.outputs import QCOutput

from custodian.qchem.handlers import QChemErrorHandler

from atomate.utils.testing import AtomateTest
from atomate.qchem.firetasks.write_inputs import WriteInputFromIOSet
from atomate.qchem.firetasks.run_calc import RunQChemCustodian as RunQChemCustodianAtomate
from atomate.qchem.firetasks.parse_outputs import QChemToDb

from moltherm.compute.jobs import QCJob
from moltherm.compute.fireworks import OptFreqSPFW, FrequencyFlatteningOptimizeFW, RunQChemCustodian, WriteCustomInput

__author__ = "Evan Spotte-Smith, Samuel Blau, Brandon Wood"
__version__ = "0.2"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Beta"
__date__ = "Jan 2019"

module_dir = join(dirname(abspath(__file__)))
files_dir = join(module_dir, "..", "..", "..", "test_files")


class TestOptFreqSPFW(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_OptFreqSPFW_defaults(self):
        pass

    def test_OptFreqSPFw_not_defaults(self):
        pass


class TestFrequencyFlatteningOptimizeFW(AtomateTest):

    def setUp(self, lpad=False):
        out_file = os.path.join(files_dir, "qchem", "test.qout.opt_0")
        qc_out = QCOutput(filename=out_file)
        self.act_mol = qc_out.data["initial_molecule"]
        super(TestFrequencyFlatteningOptimizeFW, self).setUp(lpad=False)

        self.maxDiff = None

    def tearDown(self):
        pass

    def test_FrequencyFlatteningOptimizeFW_defaults(self):
        firework = FrequencyFlatteningOptimizeFW(molecule=self.act_mol)
        self.assertEqual(firework.tasks[0].as_dict(),
                         WriteInputFromIOSet(
                             molecule=self.act_mol,
                             qchem_input_set="OptSet",
                             input_file="mol.qin",
                             qchem_input_params={}).as_dict())
        self.assertEqual(firework.tasks[1].as_dict(),
                         RunQChemCustodianAtomate(
                             qchem_cmd=">>qchem_cmd<<",
                             multimode=">>multimode<<",
                             input_file="mol.qin",
                             output_file="mol.qout",
                             qclog_file="mol.qclog",
                             max_cores=">>max_cores<<",
                             job_type="opt_with_frequency_flattener",
                             max_iterations=10,
                             max_molecule_perturb_scale=0.3,
                             gzipped_output=False,
                             reversed_direction=False).as_dict())
        self.assertEqual(firework.tasks[2].as_dict(),
                         QChemToDb(
                             db_file=None,
                             input_file="mol.qin",
                             output_file="mol.qout",
                             additional_fields={"special_run_type":
                                                    "frequency_flattener",
                                                "task_label":
                                                    "frequency flattening structure optimization"}).as_dict())
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name,
                         "frequency flattening structure optimization")

    def test_FrequencyFlatteningOptimizeFW_not_defaults(self):
        firework = FrequencyFlatteningOptimizeFW(
            molecule=self.act_mol,
            name="special frequency flattening structure optimization",
            qchem_cmd="qchem -slurm",
            multimode="mpi",
            max_cores=12,
            qchem_input_params={"pcm_dielectric": 10.0},
            max_iterations=5,
            max_molecule_perturb_scale=0.2,
            db_file=os.path.join(files_dir, "db.json"),
            parents=None)
        self.assertEqual(firework.tasks[0].as_dict(),
                         WriteInputFromIOSet(
                             molecule=self.act_mol,
                             qchem_input_set="OptSet",
                             input_file="mol.qin",
                             qchem_input_params={
                                 "pcm_dielectric": 10.0
                             }).as_dict())
        self.assertEqual(firework.tasks[1].as_dict(),
                         RunQChemCustodianAtomate(
                             qchem_cmd="qchem -slurm",
                             multimode="mpi",
                             input_file="mol.qin",
                             output_file="mol.qout",
                             qclog_file="mol.qclog",
                             max_cores=12,
                             job_type="opt_with_frequency_flattener",
                             max_iterations=5,
                             max_molecule_perturb_scale=0.2,
                             gzipped_output=False,
                             reversed_direction=False).as_dict())
        self.assertEqual(
            firework.tasks[2].as_dict(),
            QChemToDb(
                db_file=os.path.join(files_dir, "db.json"),
                input_file="mol.qin",
                output_file="mol.qout",
                additional_fields={
                    "task_label":
                        "special frequency flattening structure optimization",
                    "special_run_type":
                        "frequency_flattener"
                }).as_dict())
        self.assertEqual(firework.parents, [])
        self.assertEqual(firework.name,
                         "special frequency flattening structure optimization")


class TestRunQChemCustodian(AtomateTest):

    def setUp(self, lpad=False):
        super(TestRunQChemCustodian, self).setUp(lpad=False)

    def tearDown(self):
        pass

    def test_RunQChemCustodian_basic_defaults(self):
        with patch("moltherm.compute.fireworks.Custodian"
                   ) as custodian_patch:
            firetask = RunQChemCustodian(
                qchem_cmd="qchem",
                input_file=os.path.join(files_dir, "qchem", "co_qc.in"),
                max_cores=32)
            firetask.run_task(fw_spec={})
            custodian_patch.assert_called_once()
            self.assertEqual(custodian_patch.call_args[0][0][0].as_dict(),
                             QChemErrorHandler(
                                 input_file=os.path.join(files_dir, "qchem",
                                                         "co_qc.in"),
                                 output_file="mol.qout").as_dict())
            self.assertEqual(custodian_patch.call_args[0][1][0].as_dict(),
                             QCJob(
                                 qchem_command="qchem",
                                 max_cores=32,
                                 multimode="openmp",
                                 input_file=os.path.join(files_dir, "qchem",
                                                         "co_qc.in"),
                                 output_file="mol.qout").as_dict())
            self.assertEqual(custodian_patch.call_args[1], {
                "max_errors": 5,
                "gzipped_output": True
            })

    def test_RunQChemCustodian_using_fw_spec_defaults(self):
        with patch("moltherm.compute.fireworks.Custodian"
                   ) as custodian_patch:
            firetask = RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                scratch_dir=">>scratch_dir<<",
                input_file=os.path.join(module_dir, "..", "..", "test_files",
                                        "co_qc.in"),
                max_cores=">>max_cores<<",
                multimode=">>multimode<<")
            firetask.run_task(
                fw_spec={
                    "_fw_env": {
                        "qchem_cmd": "qchem -slurm",
                        "scratch_dir": "/this/is/a/test",
                        "max_cores": 32,
                        "multimode": "openmp"
                    }
                })
            custodian_patch.assert_called_once()
            self.assertEqual(custodian_patch.call_args[0][0][0].as_dict(),
                             QChemErrorHandler(
                                 input_file=os.path.join(
                                     module_dir, "..", "..", "test_files",
                                     "co_qc.in"),
                                 output_file="mol.qout").as_dict())
            self.assertEqual(custodian_patch.call_args[0][1][0].as_dict(),
                             QCJob(
                                 qchem_command="qchem -slurm",
                                 max_cores=">>max_cores<<",
                                 multimode=">>multimode<<",
                                 input_file=os.path.join(
                                     module_dir, "..", "..", "test_files",
                                     "co_qc.in"),
                                 output_file="mol.qout",
                                 scratch_dir="/this/is/a/test").as_dict())
            self.assertEqual(custodian_patch.call_args[1], {
                "max_errors": 5,
                "gzipped_output": True
            })

    def test_RunQChemCustodian_basic_not_defaults(self):
        with patch("moltherm.compute.fireworks.Custodian"
                   ) as custodian_patch:
            firetask = RunQChemCustodian(
                qchem_cmd="qchem -slurm",
                multimode="mpi",
                input_file=os.path.join(module_dir, "..", "..", "test_files",
                                        "co_qc.in"),
                output_file="this_is_a_test.qout",
                max_cores=4,
                qclog_file="this_is_a_test.qclog",
                suffix="bad_idea",
                save_scratch=True,
                save_name="no_idea",
                max_errors=137,
                gzipped_output=False,
                handler_group="no_handler",
                scratch_dir="/this/is/a/test")
            firetask.run_task(fw_spec={})
            custodian_patch.assert_called_once()
            self.assertEqual(custodian_patch.call_args[0][0], [])
            self.assertEqual(custodian_patch.call_args[0][1][0].as_dict(),
                             QCJob(
                                 qchem_command="qchem -slurm",
                                 multimode="mpi",
                                 input_file=os.path.join(
                                     module_dir, "..", "..", "test_files",
                                     "co_qc.in"),
                                 output_file="this_is_a_test.qout",
                                 max_cores=4,
                                 qclog_file="this_is_a_test.qclog",
                                 suffix="bad_idea",
                                 save_scratch=True,
                                 save_name="no_idea",
                                 scratch_dir="/this/is/a/test").as_dict())
            self.assertEqual(custodian_patch.call_args[1], {
                "max_errors": 137,
                "gzipped_output": False
            })

    def test_RunQChemCustodian_using_fw_spec_not_defaults(self):
        with patch("moltherm.compute.fireworks.Custodian"
                   ) as custodian_patch:
            firetask = RunQChemCustodian(
                qchem_cmd=">>qchem_cmd<<",
                multimode=">>multimode<<",
                input_file=os.path.join(module_dir, "..", "..", "test_files",
                                        "co_qc.in"),
                output_file="this_is_a_test.qout",
                max_cores=4,
                qclog_file="this_is_a_test.qclog",
                suffix="bad_idea",
                save_scratch=True,
                save_name="no_idea",
                max_errors=137,
                gzipped_output=False,
                handler_group="no_handler",
                scratch_dir=">>scratch_dir<<",
            )
            firetask.run_task(
                fw_spec={
                    "_fw_env": {
                        "qchem_cmd": "qchem -slurm",
                        "scratch_dir": "/this/is/a/test",
                        "max_cores": 32,
                        "multimode": "mpi"
                    }
                })
            custodian_patch.assert_called_once()
            self.assertEqual(custodian_patch.call_args[0][0], [])
            self.assertEqual(custodian_patch.call_args[0][1][0].as_dict(),
                             QCJob(
                                 qchem_command="qchem -slurm",
                                 multimode=">>multimode<<",
                                 input_file=os.path.join(
                                     module_dir, "..", "..", "test_files",
                                     "co_qc.in"),
                                 output_file="this_is_a_test.qout",
                                 max_cores=4,
                                 qclog_file="this_is_a_test.qclog",
                                 suffix="bad_idea",
                                 save_scratch=True,
                                 save_name="no_idea",
                                 scratch_dir="/this/is/a/test").as_dict())
            self.assertEqual(custodian_patch.call_args[1], {
                "max_errors": 137,
                "gzipped_output": False
            })

    def test_RunQChemCustodian_FF_basic_defaults(self):
        self.maxDiff = None
        with patch("moltherm.compute.fireworks.Custodian"
                   ) as custodian_patch:
            with patch(
                    "moltherm.compute.fireworks.QCJob.opt_with_frequency_flattener"
            ) as FF_patch:
                firetask = RunQChemCustodian(
                    qchem_cmd="qchem",
                    max_cores=32,
                    input_file=join(files_dir, "qchem", "test_before_run.qin"),
                    output_file=join(files_dir, "qchem", "test_before_run.qout"),
                    job_type="opt_with_frequency_flattener"
                )
                firetask.run_task(fw_spec={})
                custodian_patch.assert_called_once()
                self.assertEqual(custodian_patch.call_args[0][0][0].as_dict(),
                                 QChemErrorHandler(
                                     input_file=join(files_dir, "qchem",
                                                     "test_before_run.qin"),
                                     output_file=join(files_dir, "qchem",
                                                      "test_before_run.qout")).as_dict())
                self.assertEqual(custodian_patch.call_args[1], {
                    "max_errors": 5,
                    "gzipped_output": True
                })
                self.assertEqual(
                    FF_patch.call_args[1], {
                        "qchem_command":
                        "qchem",
                        "multimode":
                        "openmp",
                        "input_file": os.path.join(files_dir, "qchem",
                                                   "test_before_run.qin"),
                        "output_file": os.path.join(files_dir, "qchem",
                                                    "test_before_run.qout"),
                        "qclog_file": "mol.qclog",
                        "max_iterations": 10,
                        "max_molecule_perturb_scale": 0.3,
                        "scratch_dir": "/dev/shm/qcscratch/",
                        "save_scratch": False,
                        "save_name": "default_save_name",
                        "max_cores": 32,
                        "sp_params": None
                    })

    def test_RunQChemCustodian_FF_using_fw_spec_defaults(self):
        with patch("moltherm.compute.fireworks.Custodian"
                   ) as custodian_patch:
            with patch(
                    "moltherm.compute.fireworks.QCJob.opt_with_frequency_flattener"
            ) as FF_patch:
                firetask = RunQChemCustodian(
                    qchem_cmd=">>qchem_cmd<<",
                    max_cores=">>max_cores<<",
                    multimode=">>multimode<<",
                    input_file=os.path.join(files_dir, "qchem",
                                            "test_before_run.qin"),
                    output_file=os.path.join(files_dir, "qchem",
                                             "test_before_run.qout"),
                    scratch_dir=">>scratch_dir<<",
                    job_type="opt_with_frequency_flattener")
                firetask.run_task(
                    fw_spec={
                        "_fw_env": {
                            "qchem_cmd": "qchem -slurm",
                            "scratch_dir": "/this/is/a/test",
                            "max_cores": 32,
                            "multimode": "openmp"
                        }
                    })
                custodian_patch.assert_called_once()

                self.assertEqual(custodian_patch.call_args[0][0][0].as_dict(),
                                 QChemErrorHandler(
                                     input_file=os.path.join(files_dir, "qchem",
                                                             "test_before_run.qin"),
                                     output_file=os.path.join(files_dir, "qchem",
                                                              "test_before_run.qout")).as_dict())
                self.assertEqual(custodian_patch.call_args[1], {
                    "max_errors": 5,
                    "gzipped_output": True
                })
                self.assertEqual(
                    FF_patch.call_args[1], {
                        "qchem_command": "qchem -slurm",
                        "multimode": ">>multimode<<",
                        "input_file": os.path.join(files_dir, "qchem",
                                                   "test_before_run.qin"),
                        "output_file": os.path.join(files_dir, "qchem",
                                                    "test_before_run.qout"),
                        "qclog_file": "mol.qclog",
                        "max_iterations": 10,
                        "max_molecule_perturb_scale": 0.3,
                        "scratch_dir": "/this/is/a/test",
                        "save_scratch": False,
                        "save_name": "default_save_name",
                        "max_cores": ">>max_cores<<",
                        "sp_params": None
                    })

    def test_RunQChemCustodian_FF_basic_not_defaults(self):
        with patch("moltherm.compute.fireworks.Custodian"
                   ) as custodian_patch:
            with patch(
                    "moltherm.compute.fireworks.QCJob.opt_with_frequency_flattener"
            ) as FF_patch:
                firetask = RunQChemCustodian(
                    qchem_cmd="qchem -slurm",
                    input_file=os.path.join(files_dir, "qchem",
                                            "test_before_run.qin"),
                    output_file=os.path.join(files_dir, "qchem",
                                             "test_before_run.qout"),
                    job_type="opt_with_frequency_flattener",
                    max_cores=4,
                    qclog_file="this_is_a_test.qclog",
                    suffix="bad_idea",
                    save_scratch=True,
                    save_name="no_idea",
                    max_errors=137,
                    gzipped_output=False,
                    handler_group="no_handler",
                    scratch_dir="/this/is/a/test",
                    max_iterations=1029,
                    max_molecule_perturb_scale=0.5,
                    multimode="mpi")
                firetask.run_task(fw_spec={})
                custodian_patch.assert_called_once()
                self.assertEqual(custodian_patch.call_args[0][0], [])
                self.assertEqual(custodian_patch.call_args[1], {
                    "max_errors": 137,
                    "gzipped_output": False
                })
                self.assertEqual(
                    FF_patch.call_args[1], {
                        "qchem_command": "qchem -slurm",
                        "multimode": "mpi",
                        "input_file": os.path.join(files_dir, "qchem",
                                                   "test_before_run.qin"),
                        "output_file": os.path.join(files_dir, "qchem",
                                                    "test_before_run.qout"),
                        "qclog_file": "this_is_a_test.qclog",
                        "max_iterations": 1029,
                        "max_molecule_perturb_scale": 0.5,
                        "scratch_dir": "/this/is/a/test",
                        "save_scratch": True,
                        "save_name": "no_idea",
                        "max_cores": 4,
                        "sp_params": None
                    })

    def test_RunQChemCustodian_FF_using_fw_spec_not_defaults(self):
        with patch("moltherm.compute.fireworks.Custodian"
                   ) as custodian_patch:
            with patch(
                    "moltherm.compute.fireworks.QCJob.opt_with_frequency_flattener"
            ) as FF_patch:
                firetask = RunQChemCustodian(
                    qchem_cmd=">>qchem_cmd<<",
                    input_file=os.path.join(files_dir, "qchem",
                                            "test_before_run.qin"),
                    output_file=os.path.join(files_dir, "qchem",
                                             "test_before_run.qout"),
                    job_type="opt_with_frequency_flattener",
                    max_cores=4,
                    qclog_file="this_is_a_test.qclog",
                    suffix="bad_idea",
                    save_scratch=True,
                    save_name="no_idea",
                    max_errors=137,
                    gzipped_output=False,
                    handler_group="no_handler",
                    scratch_dir=">>scratch_dir<<",
                    max_iterations=1029,
                    max_molecule_perturb_scale=0.5,
                    multimode=">>multimode<<")
                firetask.run_task(
                    fw_spec={
                        "_fw_env": {
                            "qchem_cmd": "qchem -slurm",
                            "scratch_dir": "/this/is/a/test",
                            "max_cores": 32,
                            "multimode": "mpi"
                        }
                    })
                custodian_patch.assert_called_once()
                self.assertEqual(custodian_patch.call_args[0][0], [])
                self.assertEqual(custodian_patch.call_args[1], {
                    "max_errors": 137,
                    "gzipped_output": False
                })

                self.assertEqual(
                    FF_patch.call_args[1], {
                        "qchem_command": "qchem -slurm",
                        "multimode": ">>multimode<<",
                        "input_file": join(files_dir, "qchem",
                                           "test_before_run.qin"),
                        "output_file":
                        join(files_dir, "qchem",
                             "test_before_run.qout"),
                        "qclog_file": "this_is_a_test.qclog",
                        "sp_params": None,
                        "max_iterations": 1029,
                        "max_molecule_perturb_scale": 0.5,
                        "scratch_dir": "/this/is/a/test",
                        "save_scratch": True,
                        "save_name": "no_idea",
                        "max_cores": 4
                    })

class TestWriteCustomInput(AtomateTest):
    @classmethod
    def setUpClass(cls):

        co_species = ["C", "O"]
        co_coords = [[0.0, 0.0, 0.0], [1.3, 0.0, 0.0]]
        cls.co_mol = Molecule(co_species, co_coords)
        cls.co_opt_ref_in = QCInput.from_file(
            os.path.join(files_dir, "qchem", "co_qc.in"))
        cls.opt_mol_ref_in = QCInput.from_file(
            os.path.join(files_dir, "qchem", "to_opt.qin"))
        cls.opt_mol = cls.opt_mol_ref_in.molecule
        cls.opt_mol_pcm_ref_in = QCInput.from_file(
            os.path.join(files_dir, "qchem", "to_opt_pcm.qin"))

    def setUp(self, lpad=False):
        super(TestWriteCustomInput, self).setUp(lpad=False)

    def tearDown(self):
        shutil.rmtree(self.scratch_dir)
        for x in ["mol.qin"]:
            if os.path.exists(os.path.join(module_dir, x)):
                os.remove(os.path.join(module_dir, x))

    def test_write_custom_input(self):
        mol = self.co_mol
        rem = {
            "job_type": "opt",
            "basis": "6-311++G*",
            "max_scf_cycles": 200,
            "method": "wB97xd",
            "geom_opt_max_cycles": 200,
            "gen_scfman": True,
            "scf_algorithm": "gdm"
        }
        ft = WriteCustomInput(molecule=mol, rem=rem)
        ft.run_task({})
        test_dict = QCInput.from_file("mol.qin").as_dict()
        for k, v in self.co_opt_ref_in.as_dict().items():
            self.assertEqual(v, test_dict[k])