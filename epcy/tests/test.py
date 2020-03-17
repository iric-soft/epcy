# cmd to run test: coverage run -m unittest discover ./epcy/tests

import unittest

import os
import sys
import io

from argparse import Namespace
from epcy.tools import pred as tp

from epcy.utils import Classifier as uc
from epcy.utils import readers as ur

from contextlib import contextmanager


@contextmanager
def captured_output():
    new_out, new_err = io.StringIO(), io.StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class epcyTest(unittest.TestCase):
    def test_pred(self):

        design = "./data/small_for_test/design.tsv"
        mat = "./data/small_for_test/exp_matrix.tsv"

        args = Namespace(
            C=1,
            DESIGN=design,
            EXP=0,
            L2FC=0.3,
            MATRIX=mat,
            THREAD=1,
            N_DRAW=100,
            N_BAGGING=1,
            BY=-1,
            BS=0,
            LOG=True,
            QUERY="Query",
            MIN_BW=0.0,
            CPM=False,
            ANNO=None,
            PATH_OUT=None,
            SUBGROUP="subgroup",
            UTEST=False,
            TTEST=False,
            FULL=False,
            AUC=False,
            NORMAL=False,
            RANDOM_SEED=42
        )

        with captured_output() as (out, err):
            tp.main_pred(args, None)

        output = out.getvalue()
        all_lines = output.split("\n")

        selected_line = all_lines[0].split("\t")
        self.assertEqual(selected_line[2],
                         "kernel_mcc",
                         "Test fail: test_pred -> header")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[1],
                         "2.8047473",
                         "Test fail: test_pred -> L2FC")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[2],
                         "1.0",
                         "Test fail: test_pred -> KERNEL_MCC")

        selected_line = all_lines[4].split("\t")
        self.assertEqual(selected_line[2],
                         "nan",
                         "Test fail: test_pred -> NaN")

#    def test_pred_pvalue(self):
#
#        design = "./data/small_for_test/design.tsv"
#        mat = "./data/small_for_test/exp_matrix.tsv"
#
#        args = Namespace(
#            C=1,
#            DESIGN=design,
#            EXP=0,
#            L2FC=0.3,
#            MATRIX=mat,
#            THREAD=1,
#            N_BAGGING=1,
#            BY=-1,
#            BS=0,
#            LOG=True,
#            QUERY="Query",
#            MIN_BW=0.0,
#            ANNO=None,
#            PATH_OUT=None,
#            SUBGROUP="subgroup",
#            UTEST=True,
#            TTEST=True,
#            FULL=True,
#            AUC=True,
#            NORMAL=False
#        )
#
#        with captured_output() as (out, err):
#            tp.main_pred(args, None)
#
#        output = out.getvalue()
#        all_lines = output.split("\n")
#
#        selected_line = all_lines[0].split("\t")
#        self.assertEqual(selected_line[5],
#                         "auc",
#                         "Test fail: test_pred_pvalue -> AUC")
#        self.assertEqual(selected_line[6],
#                         "u_pv",
#                         "Test fail: test_pred_pvalue -> UTEST")
#        self.assertEqual(selected_line[7],
#                         "t_pv",
#                         "Test fail: test_pred_pvalue -> TTEST")

    def test_pred_thread(self):
        design = "./data/small_for_test/design.tsv"
        mat = "./data/small_for_test/exp_matrix.tsv"

        args = Namespace(
            C=1,
            DESIGN=design,
            EXP=0,
            L2FC=0.3,
            MATRIX=mat,
            THREAD=2,
            N_DRAW=100,
            N_BAGGING=1,
            BY=-1,
            BS=0,
            LOG=True,
            QUERY="Query",
            MIN_BW=0.0,
            CPM=False,
            ANNO=None,
            PATH_OUT=None,
            SUBGROUP="subgroup",
            UTEST=False,
            TTEST=False,
            FULL=False,
            AUC=False,
            NORMAL=False,
            RANDOM_SEED=42
        )

        with captured_output() as (out, err):
            tp.main_pred(args, None)

        output = out.getvalue()
        all_lines = output.split("\n")

        selected_line = all_lines[0].split("\t")
        self.assertEqual(selected_line[2],
                         "kernel_mcc",
                         "Test fail: test_pred_thread normal -> header")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[1],
                         "2.8047473",
                         "Test fail: test_pred -> L2FC")


def runTests():
    unittest.main()


if __name__ == "__main__":
    runTests()
