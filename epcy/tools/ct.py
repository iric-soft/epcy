import time
import sys
import os

import numpy as np
import pandas as pd

from .. utils import readers as ur


def main_ct(args, argparser):
    sys.stderr.write(time.strftime('%X') + ": Read input files\n")
    design = ur.get_design(args)
    num_query = len(np.where(design[args.CONDITION] == 1)[0])
    num_ref = len(np.where(design[args.CONDITION] == 0)[0])

    df_subg = pd.read_csv(args.SUBG, sep="\t", index_col=0)
    df_subg.rename(str.upper, axis='columns', inplace=True)
    df_ct = df_subg.abs() > 0.5
    tp = df_ct.iloc[:, :num_query].sum(axis=1)
    tn = df_ct.iloc[:, num_query:].sum(axis=1)

    df_ct = {'TP': tp, 'FN': num_query-tp, 'FP': num_ref-tn, 'TN': tn}
    df_ct = pd.DataFrame(df_ct)
    df_ct.reset_index(drop=False, inplace=True)

    sys.stderr.write(time.strftime('%X') + ": Write output file\n")
    if not os.path.exists(args.PATH_OUT):
        os.makedirs(args.PATH_OUT)

    ct_out = os.path.join(args.PATH_OUT, "contingency_table.xls")
    df_ct.to_csv(ct_out, index=False, sep="\t")
    sys.stderr.write(time.strftime('%X') + ": End\n")
