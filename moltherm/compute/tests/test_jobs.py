# coding: utf-8

from __future__ import unicode_literals, division

import os
from os.path import join, dirname, abspath
import shutil
from unittest import TestCase
from unittest.mock import patch
import unittest

from moltherm.compute.jobs import QCJob
from pymatgen.io.qchem.inputs import QCInput

__author__ = "Samuel Blau, Evan Spotte-Smith"
__copyright__ = "Copyright 2018, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Evan Spotte-Smith"
__email__ = "espottesmith@gmail.com"
__status__ = "Beta"
__date__ = "Jan 2019"
__credits__ = "Shyam Dwaraknath"


test_dir = join(dirname(__file__), "..", "..", "..", "test_files", "qchem")
scr_dir = join(test_dir, "scr")
cwd = os.getcwd()

module_dir = join(dirname(abspath(__file__)))
files_dir = join(module_dir, "..", "..", "..", "test_files")


class QCJobTest(TestCase):

    def test_defaults(self):
        with patch("moltherm.compute.jobs.os.putenv") as putenv_patch:
            with patch("moltherm.compute.jobs.shutil.copy") as copy_patch:
                myjob = QCJob(qchem_command="qchem", max_cores=32)
                self.assertEqual(myjob.current_command, ["qchem", "-nt", "32", "mol.qin", "mol.qout"])
                myjob.setup()
                self.assertEqual(copy_patch.call_args_list[0][0][0], "mol.qin")
                self.assertEqual(copy_patch.call_args_list[0][0][1], "mol.qin.orig")
                self.assertEqual(putenv_patch.call_args_list[0][0][0], "QCSCRATCH")
                self.assertEqual(putenv_patch.call_args_list[0][0][1], "/dev/shm/qcscratch/")
                self.assertEqual(putenv_patch.call_args_list[1][0][0], "QCTHREADS")
                self.assertEqual(putenv_patch.call_args_list[1][0][1], "32")
                self.assertEqual(putenv_patch.call_args_list[2][0][0], "OMP_NUM_THREADS")
                self.assertEqual(putenv_patch.call_args_list[2][0][1], "32")

    def test_not_defaults(self):
        with patch("moltherm.compute.jobs.os.putenv") as putenv_patch:
            myjob = QCJob(qchem_command="qchem -slurm", multimode="mpi", input_file="different.qin", output_file="not_default.qout", max_cores=12, scratch_dir="/not/default/scratch/", backup=False)
            self.assertEqual(myjob.current_command, ["qchem", "-slurm", "-np", "12", "different.qin", "not_default.qout"])
            myjob.setup()
            self.assertEqual(putenv_patch.call_args_list[0][0][0], "QCSCRATCH")
            self.assertEqual(putenv_patch.call_args_list[0][0][1], "/not/default/scratch/")

    def test_save_scratch(self):
        with patch("moltherm.compute.jobs.os.putenv") as putenv_patch:
            with patch("moltherm.compute.jobs.shutil.copy") as copy_patch:
                myjob = QCJob(qchem_command="qchem -slurm", max_cores=32, scratch_dir=os.getcwd(), save_scratch=True, save_name="freq_scratch")
                self.assertEqual(myjob.current_command, ["qchem", "-slurm", "-save", "-nt", "32", "mol.qin", "mol.qout", "freq_scratch"])
                myjob.setup()
                self.assertEqual(copy_patch.call_args_list[0][0][0], "mol.qin")
                self.assertEqual(copy_patch.call_args_list[0][0][1], "mol.qin.orig")
                self.assertEqual(putenv_patch.call_args_list[0][0][0], "QCSCRATCH")
                self.assertEqual(putenv_patch.call_args_list[0][0][1], os.getcwd())
                self.assertEqual(putenv_patch.call_args_list[1][0][0], "QCTHREADS")
                self.assertEqual(putenv_patch.call_args_list[1][0][1], "32")
                self.assertEqual(putenv_patch.call_args_list[2][0][0], "OMP_NUM_THREADS")
                self.assertEqual(putenv_patch.call_args_list[2][0][1], "32")

    # def test_read_scratch(self):
    #     with patch("moltherm.compute.jobs.os.putenv") as putenv_patch:
    #         with patch("moltherm.compute.jobs.shutil.copy") as copy_patch:
    #             myjob = QCJob(qchem_command="qchem -slurm", max_cores=32, scratch_dir=os.getcwd(), read_scratch=True, save_name="freq_scratch")
    #             self.assertEqual(myjob.current_command, ["qchem", "-slurm", "-nt", "32", "mol.qin", "mol.qout", "freq_scratch"])
    #             myjob.setup()
    #             self.assertEqual(copy_patch.call_args_list[0][0][0], "mol.qin")
    #             self.assertEqual(copy_patch.call_args_list[0][0][1], "mol.qin.orig")
    #             self.assertEqual(putenv_patch.call_args_list[0][0][0], "QCSCRATCH")
    #             self.assertEqual(putenv_patch.call_args_list[0][0][1], os.getcwd())
    #             self.assertEqual(putenv_patch.call_args_list[1][0][0], "QCTHREADS")
    #             self.assertEqual(putenv_patch.call_args_list[1][0][1], "32")
    #             self.assertEqual(putenv_patch.call_args_list[2][0][0], "OMP_NUM_THREADS")
    #             self.assertEqual(putenv_patch.call_args_list[2][0][1], "32")


class OptFFTest(TestCase):
    def setUp(self):
        os.makedirs(scr_dir)
        shutil.copyfile(join(test_dir, "standard/rct_1_471171.in.opt_0"),
                        join(scr_dir, "test.qin"))
        shutil.copyfile(join(test_dir, "standard/rct_1_471171.out.opt_0"),
                        join(scr_dir, "test.qout.opt_0"))
        shutil.copyfile(join(test_dir, "standard/rct_1_471171.out.freq_0"),
                        join(scr_dir, "test.qout.freq_0"))
        shutil.copyfile(join(test_dir, "standard/rct_1_471171.out.sp"),
                        join(scr_dir, "test.qout.sp"))
        os.chdir(scr_dir)

        self.sp_params = {"rem": {"job_type": "sp",
                                  "method": "wb97x-d",
                                  "basis": "6-311++g(d,p)",
                                  "max_scf_cycles": 200,
                                  "gen_scfman": True,
                                  "scf_algorithm": "diis",
                                  "solvent_method": "pcm"},
                          "pcm": {
                              "theory": "iefpcm"},
                          "solvent": {
                              "dielectric": 80.4}
                          }

    def tearDown(self):
        os.chdir(cwd)
        shutil.rmtree(scr_dir)

    def test_OptFF(self):
        self.maxDiff = None
        myjob = QCJob.opt_with_frequency_flattener(qchem_command="qchem",
                                                   max_cores=32,
                                                   input_file="test.qin",
                                                   output_file="test.qout",
                                                   sp_params=self.sp_params)
        expected_next = QCJob(
                        qchem_command="qchem",
                        max_cores=32,
                        multimode="openmp",
                        input_file="test.qin",
                        output_file="test.qout",
                        suffix=".opt_0",
                        backup=True).as_dict()
        self.assertEqual(next(myjob).as_dict(),expected_next)
        expected_next = QCJob(
                        qchem_command="qchem",
                        max_cores=32,
                        multimode="openmp",
                        input_file="test.qin",
                        output_file="test.qout",
                        suffix=".freq_0",
                        backup=False).as_dict()
        self.assertEqual(next(myjob).as_dict(),expected_next)
        self.assertEqual(QCInput.from_file(join(test_dir,"standard/rct_1_471171.in.freq_0")).as_dict(),
                         QCInput.from_file(join(scr_dir,"test.qin")).as_dict())
        expected_next = QCJob(
                        qchem_command="qchem",
                        max_cores=32,
                        multimode="openmp",
                        input_file="test.qin",
                        output_file="test.qout",
                        suffix=".sp",
                        backup=True).as_dict()
        self.assertEqual(next(myjob).as_dict(), expected_next)
        self.assertEqual(QCInput.from_file(join(test_dir,"standard/rct_1_471171.in.sp")).as_dict(),
                         QCInput.from_file(join(scr_dir,"test.qin")).as_dict())
        self.assertRaises(StopIteration,myjob.__next__)


class TSFFTest(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_TSFF(self):
        pass


if __name__ == "__main__":
    unittest.main()
