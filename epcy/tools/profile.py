import time
import sys

import numpy as np

from .. utils import other as uo
from .. utils import readers as ur
from .. utils import plot as up
from .. utils import Classifier as uc

def main_profile(args, argparser):
    sys.stderr.write(time.strftime('%X') + ": Read design and matrix features\n")
    (design, data, list_ids) = ur.read_design_matrix(args)
    num_query = len(np.where(design[args.SUBGROUP] == 1)[0])

    sys.stderr.write(time.strftime('%X') + ": Start profiling:\n")
    for id in args.IDS:
        sys.stderr.write(time.strftime('%X') + ":\t" + id + "\n")
        pos = np.where(list_ids == id)[0]
        if pos.shape[0] == 1:
            row_data, row_num_query = uc.Classifier.rm_missing(data[pos,:][0], num_query)
            bw = uc.Classifier.bw_nrd0(row_data)
            query_exp = row_data[:row_num_query]
            ref_exp = row_data[row_num_query:]
            up.plot_profile(id, query_exp, ref_exp, bw, args)
        else:
            if pos.shape[0] == 0:
                sys.stderr.write("\t\t -> WARNING: This feasture is not found in your matrix (-m).\n")
            else:
                sys.stderr.write("\t\t -> WARNING: This feasture is find more than one time, check your matrix (-m).\n")


    sys.stderr.write(time.strftime('%X') + ": End\n")
