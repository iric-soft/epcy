import os
import time
import sys
import h5py

import numpy as np
import pandas as pd

from .. utils import readers as ur
from .. utils import other as uo


def main_kal2mat(args, argparser):
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

    df_data = pd.DataFrame(data=data, columns=design['sample'])
    df_data.insert(loc=0, column='ID', value=list_ids)

    if not os.path.exists(args.PATH_OUT):
        os.makedirs(args.PATH_OUT)

    file_out = os.path.join(args.PATH_OUT, "readcounts.xls")
    if args.CPM:
        if args.LOG:
            file_out = os.path.join(args.PATH_OUT, "readcounts_cpm_log.xls")
        else:
            file_out = os.path.join(args.PATH_OUT, "readcounts_cpm.xls")
    elif args.TPM:
        if args.LOG:
            file_out = os.path.join(args.PATH_OUT, "readcounts_tpm_log.xls")
        else:
            file_out = os.path.join(args.PATH_OUT, "readcounts_tpm_.xls")

    df_data.to_csv(file_out, index=False, sep="\t")
