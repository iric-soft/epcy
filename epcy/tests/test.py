# cmd to run test: coverage run -m unittest discover ./epcy/tests

import unittest

import os
import sys
import io

from argparse import Namespace
from epcy.tools import pred as tp
from epcy.tools import pred_rna as tpr

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

        selected_line = all_lines[8].split("\t")
        self.assertEqual(selected_line[2],
                         "nan",
                         "Test fail: test_pred -> NaN")

        selected_line = all_lines[3].split("\t")
        self.assertEqual(selected_line[2],
                         "0.7856584",
                         "Test fail: test_pred -> feature without missing value")

        selected_line = all_lines[4].split("\t")
        self.assertEqual(selected_line[2],
                         "0.5833333",
                         "Test fail: test_pred -> feature with missing value")

        selected_line = all_lines[5].split("\t")
        self.assertEqual(selected_line[2],
                         "nan",
                         "Test fail: test_pred -> feature too much missing value")

    def test_pred_rna_cpm(self):
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
            CPM=True,
            ANNO=None,
            PATH_OUT=None,
            SUBGROUP="subgroup",
            UTEST=False,
            TTEST=False,
            FULL=False,
            AUC=False,
            NORMAL=False,
            GENE=False,
            TPM=False,
            KAL=False,
            RANDOM_SEED=42
        )

        with captured_output() as (out, err):
            tpr.main_pred_rna(args, None)

        output = out.getvalue()
        all_lines = output.split("\n")

        selected_line = all_lines[0].split("\t")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[1],
                         "2.7649796",
                         "Test fail: test_pred_rna_cpm -> L2FC with CPM")

    def test_pred_rna_kall_gene_bagging(self):
        design = "./data/small_leucegene/5_inv16_vs_5/design.tsv"
        anno = "./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3"

        args = Namespace(
            C=1,
            DESIGN=design,
            EXP=0,
            L2FC=0,
            THREAD=1,
            N_DRAW=10,
            N_BAGGING=2,
            BY=-1,
            BS=0,
            LOG=True,
            QUERY="Query",
            MIN_BW=0.1,
            CPM=True,
            ANNO=anno,
            PATH_OUT=None,
            SUBGROUP="subgroup",
            UTEST=False,
            TTEST=False,
            FULL=False,
            AUC=False,
            NORMAL=True,
            GENE=True,
            TPM=False,
            KAL=True,
            RANDOM_SEED=42,

        )

        with captured_output() as (out, err):
            tpr.main_pred_rna(args, None)

        output = out.getvalue()
        all_lines = output.split("\n")

        selected_line = all_lines[0].split("\t")
        self.assertEqual(selected_line[2],
                         "kernel_mcc",
                         "Test fail: test_pred_rna_kall_gene_bagging -> header")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[1],
                         "1.6704607",
                         "Test fail: test_pred_rna_kall_gene_bagging -> L2FC")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[2],
                         "1.0",
                         "Test fail: test_pred_rna_kall_gene_bagging -> MCC")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[9],
                         "0.8",
                         "Test fail: test_pred_rna_kall_gene_bagging -> MCC")

    def test_pred_rna_kall_gene_bagging_tpm(self):
        design = "./data/small_leucegene/5_inv16_vs_5/design.tsv"
        anno = "./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3"

        args = Namespace(
            C=1,
            DESIGN=design,
            EXP=0,
            L2FC=0,
            THREAD=1,
            N_DRAW=10,
            N_BAGGING=2,
            BY=-1,
            BS=0,
            LOG=True,
            QUERY="Query",
            MIN_BW=0.1,
            CPM=False,
            ANNO=anno,
            PATH_OUT=None,
            SUBGROUP="subgroup",
            UTEST=False,
            TTEST=False,
            FULL=False,
            AUC=False,
            NORMAL=True,
            GENE=True,
            TPM=True,
            KAL=True,
            RANDOM_SEED=42,

        )

        with captured_output() as (out, err):
            tpr.main_pred_rna(args, None)

        output = out.getvalue()
        all_lines = output.split("\n")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[1],
                         "2.3213472",
                         "Test fail: test_pred_rna_kall_gene_bagging_tpm -> L2FC TPM")

    def test_pred_rna_kall_gene(self):
        design = "./data/small_leucegene/5_inv16_vs_5/design.tsv"
        anno = "./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3"

        args = Namespace(
            C=1,
            DESIGN=design,
            EXP=0,
            L2FC=0,
            THREAD=1,
            N_DRAW=100,
            N_BAGGING=1,
            BY=-1,
            BS=5,
            LOG=True,
            QUERY="Query",
            MIN_BW=0.1,
            CPM=True,
            ANNO=anno,
            PATH_OUT=None,
            SUBGROUP="subgroup",
            UTEST=False,
            TTEST=False,
            FULL=False,
            AUC=False,
            NORMAL=True,
            GENE=True,
            TPM=False,
            KAL=True,
            RANDOM_SEED=42,

        )

        with captured_output() as (out, err):
            tpr.main_pred_rna(args, None)

        output = out.getvalue()
        all_lines = output.split("\n")

        selected_line = all_lines[0].split("\t")
        self.assertEqual(selected_line[2],
                         "kernel_mcc",
                         "Test fail: test_pred_rna_kall_gene -> header")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[2],
                         "1.0",
                         "Test fail: test_pred_rna_kall_gene -> MCC")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[9],
                         "0.882861",
                         "Test fail: test_pred_rna_kall_gene -> MCC")

    def test_pred_rna_kall_trans(self):
        design = "./data/small_leucegene/5_inv16_vs_5/design.tsv"
        anno = "./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3"

        args = Namespace(
            C=1,
            DESIGN=design,
            EXP=0,
            L2FC=0,
            THREAD=1,
            N_DRAW=10,
            N_BAGGING=1,
            BY=-1,
            BS=0,
            LOG=True,
            QUERY="Query",
            MIN_BW=0.1,
            CPM=True,
            ANNO=anno,
            PATH_OUT=None,
            SUBGROUP="subgroup",
            UTEST=False,
            TTEST=False,
            FULL=False,
            AUC=False,
            NORMAL=False,
            GENE=False,
            TPM=False,
            KAL=True,
            RANDOM_SEED=42,

        )

        with captured_output() as (out, err):
            tpr.main_pred_rna(args, None)

        output = out.getvalue()
        all_lines = output.split("\n")

        selected_line = all_lines[0].split("\t")
        self.assertEqual(selected_line[2],
                         "kernel_mcc",
                         "Test fail: test_pred_rna_kall_trans -> header")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[2],
                         "0.51592886",
                         "Test fail: test_pred_rna_kall_trans -> MCC")

    def test_pred_rna_kall_miss_anno(self):
        design = "./data/small_leucegene/5_inv16_vs_5/design.tsv"
        anno = "./data/small_genome/Homo_sapiens.GRCh38.84.reduce2.gff3"

        args = Namespace(
            C=1,
            DESIGN=design,
            EXP=0,
            L2FC=0,
            THREAD=1,
            N_DRAW=10,
            N_BAGGING=1,
            BY=-1,
            BS=0,
            LOG=True,
            QUERY="Query",
            MIN_BW=0.1,
            CPM=True,
            ANNO=anno,
            PATH_OUT=None,
            SUBGROUP="subgroup",
            UTEST=False,
            TTEST=False,
            FULL=False,
            AUC=False,
            NORMAL=False,
            GENE=True,
            TPM=False,
            KAL=True,
            RANDOM_SEED=42,

        )

        with captured_output() as (out, err):
            tpr.main_pred_rna(args, None)

        output = out.getvalue()
        all_lines = output.split("\n")
        selected_line = all_lines[len(all_lines) - 2].split("\t")

        self.assertEqual(selected_line[0],
                         "ENST00000411957",
                         "Test fail: test_pred_rna_kall_miss_anno -> ENST")

    def test_pred_pvalue(self):

        design = "./data/small_for_test/design.tsv"
        mat = "./data/small_for_test/exp_matrix.tsv"

        args = Namespace(
            C=1,
            DESIGN=design,
            EXP=0,
            L2FC=0.3,
            MATRIX=mat,
            THREAD=1,
            N_BAGGING=1,
            BY=-1,
            BS=0,
            LOG=True,
            QUERY="Query",
            MIN_BW=0.0,
            ANNO=None,
            PATH_OUT=None,
            SUBGROUP="subgroup",
            UTEST=True,
            TTEST=True,
            FULL=True,
            AUC=True,
            NORMAL=False,
            RANDOM_SEED=42,
            N_DRAW=100
        )

        with captured_output() as (out, err):
            tp.main_pred(args, None)

        output = out.getvalue()
        all_lines = output.split("\n")

        selected_line = all_lines[0].split("\t")
        self.assertEqual(selected_line[9],
                         "auc",
                         "Test fail: test_pred_pvalue -> AUC")
        self.assertEqual(selected_line[10],
                         "u_pv",
                         "Test fail: test_pred_pvalue -> UTEST")
        self.assertEqual(selected_line[11],
                         "t_pv",
                         "Test fail: test_pred_pvalue -> TTEST")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[9],
                         "1.0",
                         "Test fail: test_pred_pvalue -> AUC value l1")
        self.assertEqual(selected_line[10],
                         "0.0047718217",
                         "Test fail: test_pred_pvalue -> UTEST value l1")
        self.assertEqual(selected_line[11],
                         "1.9629755e-05",
                         "Test fail: test_pred_pvalue -> TTEST value l1")

        selected_line = all_lines[8].split("\t")
        self.assertEqual(selected_line[9],
                         "nan",
                         "Test fail: test_pred_pvalue -> AUC value l4")
        self.assertEqual(selected_line[10],
                         "nan",
                         "Test fail: test_pred_pvalue -> UTEST value l4")
        self.assertEqual(selected_line[11],
                         "nan",
                         "Test fail: test_pred_pvalue -> TTEST value l4")

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
