import time
import sys

import numpy as np
import pandas as pd

from .. utils import other as uo
from .. utils import readers as ur
from .. utils import plot as up


def main_explore(args, argparser):
    sys.stderr.write(time.strftime('%X') + ": Read input files\n")
    design = ur.get_design(args)
    num_query = len(np.where(design[args.CONDITION] == 1)[0])

    df_pred = pd.read_csv(args.PRED, sep="\t")
    df_pred.rename(str.upper, axis='columns', inplace=True)
    df_pred = df_pred.reindex(
                df_pred.KERNEL_MCC.sort_values(ascending=False).index)

    df_subg = pd.read_csv(args.SUBG, sep="\t")
    df_subg.rename(str.upper, axis='columns', inplace=True)
    df_subg = df_subg.reindex(df_pred.index)

    df_heatmap = df_subg.iloc[:args.TOP]
    df_pred = df_pred.iloc[:args.TOP]

    sys.stderr.write(time.strftime('%X') + ": Plot fig\n")
    up.plot_explore_heatmap(df_heatmap, df_pred, args)
