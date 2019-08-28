import time
import sys

import numpy as np
import pandas as pd

from .. utils import other as uo
from .. utils import readers as ur

def read_design_matrix(args):
    design = ur.get_design(args)
    data = pd.io.parsers.read_csv(args.MATRIX, sep="\t", index_col=0)
    if sum(~design["sample"].isin(data.columns)) > 0:
        sys.stderr.write("WARNING: Some samples are present in the design, but not in the quantification matrix\n")
        sys.stderr.write("\t the analysis will be done without these samples:\n")
        sys.stderr.write(str(design[~design["sample"].isin(data.columns)]) + "\n")
        design = design[design["sample"].isin(data.columns)]

    data = data.reindex(design["sample"], axis=1)

    data = data[(data.T != 0).any()]

    list_ids = list(data.index)

    data = data.values
    if not args.LOG:
        data = np.log2(data + args.C)

    return(design, data, list_ids)

def main_pred(args, argparser):

    sys.stderr.write(time.strftime('%X') + ": Read design and matrix features\n")
    (design, data, list_ids) = read_design_matrix(args)

    num_pred = data.shape[0]

    all_classifier = uo.compute_pred(args, num_pred, list_ids, data, design)

    uo.save_pred_res(args, all_classifier)
