import os
import time
import sys
import h5py

import numpy as np
import pandas as pd

from .. utils import readers as ur
from .. utils import other as uo


def main_pred_rna(args, argparser):
    if args.KAL:
        if args.GENE:
            sys.stderr.write(time.strftime('%X') + ": Run EPCY on kallisto " +
                             "output on gene\n")
        else:
            sys.stderr.write(time.strftime('%X') + ": Run EPCY on kallisto " +
                             "output on transcript!!!\n")
            sys.stderr.write(time.strftime('%X') + ":\t add --gene to run " +
                             "on gene level\n")

    df_anno = ur.read_anno(args)

    sys.stderr.write(time.strftime('%X') + ": Read design and matrix " +
                     "features\n")
    (design, data, list_ids) = ur.read_design_matrix_rna(args, df_anno)

    if design is None:
        exit()

    num_pred = data.shape[0]

    uo.compute_pred(args, num_pred, list_ids, data, design)
    #uo.save_pred_res(args, all_classifier)
