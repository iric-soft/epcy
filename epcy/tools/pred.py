import time
import sys
import os

import numpy as np
import pandas as pd

from .. utils import other as uo
from .. utils import readers as ur


def main_pred(args, argparser):

    sys.stderr.write(time.strftime('%X') + ": Read design and matrix " +
                     "features\n")
    (design, data, list_ids) = ur.read_design_matrix(args)

    if args.PATH_OUT is not None:
        if not os.path.exists(args.PATH_OUT):
            os.makedirs(args.PATH_OUT)

        file_basename = os.path.basename(args.MATRIX)
        file_out = os.path.join(args.PATH_OUT, file_basename)
        if args.NORM:
            if args.LOG:
                file_out = os.path.join(args.PATH_OUT, "log_norm_" + file_basename)
            else:
                file_out = os.path.join(args.PATH_OUT, "norm_" + file_basename)
        elif args.LOG:
            file_out = os.path.join(args.PATH_OUT, "log_" + file_basename)

        df_data = pd.DataFrame(data=data, columns=design['sample'])
        df_data.insert(loc=0, column='ID', value=list_ids)

        df_data.to_csv(file_out, index=False, sep="\t")

    num_pred = data.shape[0]

    all_classifier = uo.compute_pred(args, num_pred, list_ids, data, design)

    #uo.save_pred_res(args, all_classifier)
