import os
import time
import sys

import numpy as np
import pandas as pd

from multiprocessing import Pool

from .. utils import Classifier as uc
from .. utils import readers as ur

def run_classifier(args, design, data, list_ids, start):
    data_selected = data[start:(start + args.BY)]
    list_ids_selected = list_ids[start:(start + args.BY)]
    c = uc.Classifier(args, design, data_selected, list_ids_selected, start)
    c.run()

    return(c)

def read_design_matrix(args):
    design = ur.get_design(args)
    data = pd.io.parsers.read_csv(args.MATRIX, sep="\t", index_col=0)
    if sum(~design["sample"].isin(data.columns)) > 0:
        sys.stderr.write("WARNING: Some samples are present in the design, but not in the quantification matrix\n")
        sys.stderr.write("\t the analysis will be done without these samples:\n")
        sys.stderr.write(str(design[~design["sample"].isin(data.columns)]) + "\n")
        design = design[design["sample"].isin(data.columns)]

    data = data.reindex(design["sample"], axis=1)

    if args.CPM:
        f_norm = 1e6 /  data.iloc[:,1:].sum()

    #data = data[start:(start+args.BY)]
    data = data[(data.T != 0).any()]

    if args.CPM:
        data.iloc[:,1:] = data.iloc[:,1:] * f_norm

    list_ids = list(data.index)

    data = data.values
    if not args.LOG:
        data = np.log2(data + args.C)

    return(design, data, list_ids)

def main_pred(args, argparser):

    if args.ANNO is not None:
        sys.stderr.write(time.strftime('%X') + ": Read annotation\n")
        df_anno = ur.gff_2_df(args.ANNO)
        df_anno.dropna(axis=0, subset=["ID"], inplace=True)
        df_anno["ID"] = df_anno["ID"].str.replace("gene:", "")
        df_anno["ID"] = df_anno["ID"].str.replace("transcript:", "")
        df_anno["Parent"] = df_anno["Parent"].str.replace("gene:", "")
        df_anno["Parent"] = df_anno["Parent"].str.replace("transcript:", "")
        df_anno.set_index("ID", inplace=True)

    num_pred = sum(1 for line in open(args.MATRIX, 'r')) - 1

    sys.stderr.write(time.strftime('%X') + ": Start epcy analysis\n")
    if args.BY == -1:
        args.BY = int(num_pred / args.THREAD) + 1

    sys.stderr.write(time.strftime('%X') + ": Read design and matrix features\n")
    (design, data, list_ids) = read_design_matrix(args)

    all_classifier = []
    if args.THREAD <= 1:
        tasks = [(args, design, data, list_ids, start) for start in np.arange(0, num_pred, args.BY)]
        for task in tasks:
            all_classifier.append(run_classifier(task[0], task[1], task[2], task[3], task[4]))

    else:
        with Pool(processes=args.THREAD) as p:
            tasks = [(args, design, data, start) for start in np.arange(0, num_pred, args.BY)]

            jobs = [p.apply_async(run_classifier, t) for t in tasks]
            all_classifier = [job.get(None) for job in jobs]
            p.close()
            p.join()


    sys.stderr.write(time.strftime('%X') + ": Save epcy results\n")
    if args.PATH_OUT is not None:
        if not os.path.exists(args.PATH_OUT):
            os.makedirs(args.PATH_OUT)

        file_out = args.PATH_OUT + "/prediction_capability.xls"
        file_pred_out = args.PATH_OUT + "/subgroup_predicted.xls"

        with open(file_out, 'w') as w_csv:
            all_classifier[0].print_feature_header(w_csv)
            for classifier in all_classifier:
                classifier.print_feature_pred(w_csv)

        with open(file_pred_out, 'w') as w_csv:
            all_classifier[0].print_subgroup_header(w_csv)
            for classifier in all_classifier:
                classifier.print_subgroup_predicted(w_csv, "kernel")
    else:
        all_classifier[0].print_feature_header(sys.stdout)
        for classifier in all_classifier:
            classifier.print_feature_pred(sys.stdout)

    sys.stderr.write(time.strftime('%X') + ": End\n")
