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


default_args = Namespace(
    C=1,
    DESIGN=None,
    EXP=None,
    L2FC=0,
    MATRIX=None,
    THREAD=1,
    N_DRAW=100,
    N_BAGGING=1,
    BY=-1,
    BS=0,
    LOG=False,
    QUERY="Query",
    REF=None,
    MIN_BW=0.1,
    CPM=False,
    CPMED=False,
    ANNO=None,
    N_FOLD=None,
    PATH_OUT=None,
    CONDITION="condition",
    SHUFFLE=False,
    PPV=False,
    NPV=False,
    TPR=False,
    TNR=False,
    FNR=False,
    FPR=False,
    FDR=False,
    FOR=False,
    TS=False,
    ACC=False,
    F1=False,
    NORMAL=False,
    UTEST=False,
    TTEST=False,
    FULL=False,
    AUC=False,
    GENE=False,
    TPM=False,
    KAL=False,
    RANDOM_SEED=42
)


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

        args = Namespace(**vars(default_args))
        args.DESIGN = design
        args.MATRIX = mat
        args.L2FC = 0.3
        args.LOG = True

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
                         "2.8047472709397807",
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
                         "0.7856583920412332",
                         "Test fail: test_pred -> feature without missing value")

        selected_line = all_lines[4].split("\t")
        self.assertEqual(selected_line[2],
                         "0.5833333333333334",
                         "Test fail: test_pred -> feature with missing value")

        selected_line = all_lines[5].split("\t")
        self.assertEqual(selected_line[2],
                         "nan",
                         "Test fail: test_pred -> feature too much missing value")

    def test_pred_rna_cpm(self):
        design = "./data/small_for_test/design.tsv"
        mat = "./data/small_for_test/exp_matrix.tsv"

        args = Namespace(**vars(default_args))
        args.DESIGN = design
        args.MATRIX = mat
        args.LOG = True
        args.L2FC = 0.3
        args.CPM = True

        with captured_output() as (out, err):
            tpr.main_pred_rna(args, None)

        output = out.getvalue()
        all_lines = output.split("\n")

        selected_line = all_lines[0].split("\t")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[1],
                         "2.764979565971643",
                         "Test fail: test_pred_rna_cpm -> L2FC with CPM")

    def test_pred_rna_kall_gene_bagging(self):
        design = "./data/small_leucegene/5_inv16_vs_5/design.tsv"
        anno = "./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3"

        args = Namespace(**vars(default_args))
        args.DESIGN = design
        args.BS = 5
        args.LOG = True
        args.N_DRAW = 10
        args.N_BAGGING = 10
        args.CPM = True
        args.ANNO = anno
        args.GENE = True
        args.KAL = True

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
                         "1.6635006178253366",
                         "Test fail: test_pred_rna_kall_gene_bagging -> L2FC")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[2],
                         "1.0",
                         "Test fail: test_pred_rna_kall_gene_bagging -> MCC")

    def test_pred_rna_kall_gene_bagging_tpm(self):
        design = "./data/small_leucegene/5_inv16_vs_5/design.tsv"
        anno = "./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3"

        args = Namespace(**vars(default_args))
        args.DESIGN = design
        args.BS = 5
        args.LOG = True
        args.N_DRAW = 10
        args.N_BAGGING = 10
        args.TPM = True
        args.ANNO = anno
        args.GENE = True
        args.KAL = True

        with captured_output() as (out, err):
            tpr.main_pred_rna(args, None)

        output = out.getvalue()
        all_lines = output.split("\n")

        selected_line = all_lines[1].split("\t")
        self.assertEqual(selected_line[1],
                         "2.320041910153135",
                         "Test fail: test_pred_rna_kall_gene_bagging_tpm -> L2FC TPM")

    def test_pred_rna_kall_gene(self):
        design = "./data/small_leucegene/5_inv16_vs_5/design.tsv"
        anno = "./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3"

        args = Namespace(**vars(default_args))
        args.DESIGN = design
        args.BS = 5
        args.LOG = True
        args.CPM = True
        args.ANNO = anno
        args.GENE = True
        args.KAL = True

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

    def test_pred_rna_kall_trans(self):
        design = "./data/small_leucegene/5_inv16_vs_5/design.tsv"
        anno = "./data/small_genome/Homo_sapiens.GRCh38.84.reduce.gff3"

        args = Namespace(**vars(default_args))
        args.DESIGN = design
        args.N_DRAW = 10
        args.LOG = True
        args.CPM = True
        args.ANNO = anno
        args.KAL = True

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
                         "0.5827131783528432",
                         "Test fail: test_pred_rna_kall_trans -> MCC")

    def test_pred_rna_kall_miss_anno(self):
        design = "./data/small_leucegene/5_inv16_vs_5/design.tsv"
        anno = "./data/small_genome/Homo_sapiens.GRCh38.84.reduce2.gff3"

        args = Namespace(**vars(default_args))
        args.DESIGN = design
        args.N_DRAW = 10
        args.LOG = True
        args.CPM = True
        args.ANNO = anno
        args.GENE = True
        args.KAL = True

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

        args = Namespace(**vars(default_args))
        args.DESIGN = design
        args.MATRIX = mat
        args.L2FC = 0.3
        args.LOG = True
        args.UTEST = True
        args.TTEST = True
        args.FULL = True
        args.AUC = True

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
                         "0.004771821713797114",
                         "Test fail: test_pred_pvalue -> UTEST value l1")
        self.assertEqual(selected_line[11],
                         "1.9629755742218275e-05",
                         "Test fail: test_pred_pvalue -> TTEST value l1")

        selected_line = all_lines[8].split("\t")
        self.assertEqual(selected_line[9],
                         "nan",
                         "Test fail: test_pred_pvalue -> AUC value l8")
        self.assertEqual(selected_line[10],
                         "nan",
                         "Test fail: test_pred_pvalue -> UTEST value l8")
        self.assertEqual(selected_line[11],
                         "nan",
                         "Test fail: test_pred_pvalue -> TTEST value l8")

    def test_pred_all_predictive_score(self):

        design = "./data/small_for_test/design.tsv"
        mat = "./data/small_for_test/exp_matrix.tsv"

        args = Namespace(**vars(default_args))
        args.DESIGN = design
        args.MATRIX = mat
        args.L2FC = 0.3
        args.LOG = True
        args.PPV = True
        args.NPV = True
        args.TPR = True
        args.TNR = True
        args.FNR = True
        args.FPR = True
        args.FDR = True
        args.FOR = True
        args.TS = True
        args.ACC = True
        args.F1 = True

        with captured_output() as (out, err):
            tp.main_pred(args, None)

        output = out.getvalue()
        all_lines = output.split("\n")

        selected_line = all_lines[0].split("\t")
        self.assertEqual(selected_line[5],
                         "kernel_ppv",
                         "Test fail: test_pred_all_predictive_score -> PPV")
        self.assertEqual(selected_line[8],
                         "kernel_npv",
                         "Test fail: test_pred_all_predictive_score -> NPV")
        self.assertEqual(selected_line[11],
                         "kernel_tpr",
                         "Test fail: test_pred_all_predictive_score -> TPR")
        self.assertEqual(selected_line[14],
                         "kernel_tnr",
                         "Test fail: test_pred_all_predictive_score -> TNR")
        self.assertEqual(selected_line[17],
                         "kernel_fnr",
                         "Test fail: test_pred_all_predictive_score -> FNR")
        self.assertEqual(selected_line[20],
                         "kernel_fpr",
                         "Test fail: test_pred_all_predictive_score -> FPR")
        self.assertEqual(selected_line[23],
                         "kernel_fdr",
                         "Test fail: test_pred_all_predictive_score -> FDR")
        self.assertEqual(selected_line[26],
                         "kernel_for",
                         "Test fail: test_pred_all_predictive_score -> FOR")
        self.assertEqual(selected_line[29],
                         "kernel_ts",
                         "Test fail: test_pred_all_predictive_score -> TS")
        self.assertEqual(selected_line[32],
                         "kernel_acc",
                         "Test fail: test_pred_all_predictive_score -> ACC")
        self.assertEqual(selected_line[35],
                         "kernel_f1",
                         "Test fail: test_pred_all_predictive_score -> F1")

        selected_line = all_lines[3].split("\t")
        self.assertEqual(selected_line[5],
                         "0.8492063492063492",
                         "Test fail: test_pred_all_predictive_score -> PPV value l3")
        self.assertEqual(selected_line[8],
                         "0.9444444444444444",
                         "Test fail: test_pred_all_predictive_score -> NPV value l3")
        self.assertEqual(selected_line[11],
                         "0.9444444444444444",
                         "Test fail: test_pred_all_predictive_score -> TPR value l3")
        self.assertEqual(selected_line[14],
                         "0.8333333333333331",
                         "Test fail: test_pred_all_predictive_score -> TNR value l3")
        self.assertEqual(selected_line[17],
                         "0.05555555555555555",
                         "Test fail: test_pred_all_predictive_score -> FNR value l3")
        self.assertEqual(selected_line[20],
                         "0.16666666666666666",
                         "Test fail: test_pred_all_predictive_score -> FPR value l3")
        self.assertEqual(selected_line[23],
                         "0.15079365079365076",
                         "Test fail: test_pred_all_predictive_score -> FDR value l3")
        self.assertEqual(selected_line[26],
                         "0.05555555555555555",
                         "Test fail: test_pred_all_predictive_score -> FOR value l3")
        self.assertEqual(selected_line[29],
                         "0.8095238095238095",
                         "Test fail: test_pred_all_predictive_score -> TS value l3")
        self.assertEqual(selected_line[32],
                         "0.8888888888888888",
                         "Test fail: test_pred_all_predictive_score -> ACC value l3")
        self.assertEqual(selected_line[35],
                         "0.8931623931623932",
                         "Test fail: test_pred_all_predictive_score -> F1 value l3")

        selected_line = all_lines[8].split("\t")
        self.assertEqual(selected_line[5],
                         "nan",
                         "Test fail: test_pred_all_predictive_score -> PPV value l8")
        self.assertEqual(selected_line[8],
                         "nan",
                         "Test fail: test_pred_all_predictive_score -> NPV value l8")
        self.assertEqual(selected_line[11],
                         "nan",
                         "Test fail: test_pred_all_predictive_score -> TPR value l8")
        self.assertEqual(selected_line[14],
                         "nan",
                         "Test fail: test_pred_all_predictive_score -> TNR value l8")
        self.assertEqual(selected_line[17],
                         "nan",
                         "Test fail: test_pred_all_predictive_score -> FNR value l8")
        self.assertEqual(selected_line[20],
                         "nan",
                         "Test fail: test_pred_all_predictive_score -> FPR value l8")
        self.assertEqual(selected_line[23],
                         "nan",
                         "Test fail: test_pred_all_predictive_score -> FDR value l8")
        self.assertEqual(selected_line[26],
                         "nan",
                         "Test fail: test_pred_all_predictive_score -> FOR value l8")
        self.assertEqual(selected_line[29],
                         "nan",
                         "Test fail: test_pred_all_predictive_score -> TS value l8")
        self.assertEqual(selected_line[32],
                         "nan",
                         "Test fail: test_pred_all_predictive_score -> ACC value l8")
        self.assertEqual(selected_line[35],
                         "nan",
                         "Test fail: test_pred_all_predictive_score -> F1 value l8")

    def test_pred_thread(self):
        design = "./data/small_for_test/design.tsv"
        mat = "./data/small_for_test/exp_matrix.tsv"

        args = Namespace(**vars(default_args))
        args.DESIGN = design
        args.MATRIX = mat
        args.L2FC = 0.3
        args.LOG = True
        args.THREAD = 2

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
                         "2.8047472709397807",
                         "Test fail: test_pred -> L2FC")


def runTests():
    unittest.main()


if __name__ == "__main__":
    runTests()
