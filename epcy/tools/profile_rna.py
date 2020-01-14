import time
import sys

import numpy as np

from .. utils import readers as ur
from .. utils import plot as up
from .. utils import Classifier as uc

def main_profile_rna(args, argparser):
    if args.KAL:
        if args.GENE:
            sys.stderr.write(time.strftime('%X') + ": Run EPCY on kallisto output on gene\n")
        else:
            if not hasattr(args, 'MATRIX') or args.MATRIX == None:
                sys.stderr.write(time.strftime('%X') + ": Run EPCY on kallisto output on transcript!!!\n")
                sys.stderr.write(time.strftime('%X') + ":\t add --gene to run on gene level\n")

    df_anno = None
    if args.GENE:
        if hasattr(args, 'ANNO') and args.ANNO is not None:
            df_anno = ur.read_anno(args)
        else:
            sys.stderr.write(time.strftime('%X') + ": An annotation file is need to switch transcripts quantification into gene quantification (see --anno)\n")
            exit()

    sys.stderr.write(time.strftime('%X') + ": Read design and matrix features\n")
    (design, data, list_ids) = ur.read_design_matrix_rna(args, df_anno)

    if design is None or data is None or list_ids is None:
        exit()

    num_query = len(np.where(design[args.SUBGROUP] == 1)[0])

    sys.stderr.write(time.strftime('%X') + ": Start profiling:\n")
    for id in args.IDS:
        sys.stderr.write(time.strftime('%X') + ":\t" + id + "\n")
        pos = np.where(list_ids == id)[0]
        if pos.shape[0] == 1:
            bw = uc.Classifier.bw_nrd0(data[pos,:])
            query_exp = data[pos,:num_query][0]
            ref_exp = data[pos,num_query:][0]
            up.plot_profile(id, query_exp, ref_exp, bw, args)
        else:
            if pos.shape[0] == 0:
                sys.stderr.write("\t\t -> WARNING: This feasture is not found in your expression matrix.\n")
            else:
                sys.stderr.write("\t\t -> WARNING: This feasture is find more than one time, check your expression matrix.\n")


    sys.stderr.write(time.strftime('%X') + ": End\n")
