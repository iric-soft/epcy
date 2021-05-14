import time
import sys

import numpy as np

from .. utils import readers as ur
from .. utils import plot as up
from .. utils import Classifier as uc


def main_profile_rna(args, argparser):
    if args.KAL:
        if args.GENE:
            sys.stderr.write(time.strftime('%X') + ": Run EPCY on kallisto " +
                             "output on gene\n")
        else:
            if not hasattr(args, 'MATRIX') or args.MATRIX is None:
                sys.stderr.write(time.strftime('%X') + ": Run EPCY on " +
                                 "kallisto output on transcript!!!\n")
                sys.stderr.write(time.strftime('%X') + ":\t add --gene to " +
                                 "run on gene level\n")

    df_anno = None
    if args.GENE:
        if hasattr(args, 'ANNO') and args.ANNO is not None:
            df_anno = ur.read_anno(args)
        else:
            sys.stderr.write(time.strftime('%X') + ": An annotation file is " +
                             "need to switch transcripts quantification " +
                             "into gene quantification (see --anno)\n")
            exit()

    sys.stderr.write(time.strftime('%X') + ": Read design and matrix " +
                     "features\n")
    (design, data, list_ids) = ur.read_design_matrix_rna(args, df_anno)

    if design is None or data is None or list_ids is None:
        exit()

    num_query = len(np.where(design[args.CONDITION] == 1)[0])

    num_bs = 0
    if hasattr(args, 'BS') and args.BS is not None:
        num_bs = args.BS

    sys.stderr.write(time.strftime('%X') + ": Start profiling:\n")
    for id in args.IDS:
        sys.stderr.write(time.strftime('%X') + ":\t" + id + "\n")
        pos = np.where(list_ids == id)[0]
        if pos.shape[0] == 1:
            row_data, row_num_query, ids_na = uc.rm_missing(
                                                               data[pos, :][0],
                                                               num_query)
            row_query = row_data[:row_num_query]
            row_ref = row_data[row_num_query:]

            bw_query = uc.bw_nrd(row_query, num_bs)
            bw_ref = uc.bw_nrd(row_ref, num_bs)

            if bw_query < args.MIN_BW:
                bw_query = args.MIN_BW
            if bw_ref < args.MIN_BW:
                bw_ref = args.MIN_BW

            up.plot_profile(id, row_query, row_ref, bw_query, bw_ref, args)
        else:
            if pos.shape[0] == 0:
                sys.stderr.write("\t\t -> WARNING: This feasture is not " +
                                 "found  in your expression matrix.\n")
            else:
                sys.stderr.write("\t\t -> WARNING: This feasture is find " +
                                 "more than one time, check your expression " +
                                 "matrix.\n")

    sys.stderr.write(time.strftime('%X') + ": End\n")
