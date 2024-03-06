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

    if args.PATH_OUT is not None:
        if not os.path.exists(args.PATH_OUT):
            os.makedirs(args.PATH_OUT)

        file_basename = ""
        if args.MATRIX is not None:
            file_basename = os.path.basename(args.MATRIX)
        
        file_out = os.path.join(args.PATH_OUT, file_basename)
        if args.CPM:
            if args.LOG:
                file_out = os.path.join(args.PATH_OUT, "log_cpm.tsv")
            else:
                file_out = os.path.join(args.PATH_OUT, "cpm.tsv")
        elif args.CPMED:
            if args.LOG:
                file_out = os.path.join(args.PATH_OUT, "log_cpmed.tsv")
            else:
                file_out = os.path.join(args.PATH_OUT, "cpmed.tsv")
        elif args.TPM:
            if args.LOG:
                file_out = os.path.join(args.PATH_OUT, "log_tpm.tsv")
            else:
                file_out = os.path.join(args.PATH_OUT, "tpm.tsv")
        elif args.LOG:
            file_out = os.path.join(args.PATH_OUT, "log_" + file_basename)

        df_data = pd.DataFrame(data=data, columns=design['sample'])
        df_data.insert(loc=0, column='ID', value=list_ids)

        df_data.to_csv(file_out, index=False, sep="\t")


    if design is None:
        exit()

    num_pred = data.shape[0]

    uo.compute_pred(args, num_pred, list_ids, data, design)
    #uo.save_pred_res(args, all_classifier)
