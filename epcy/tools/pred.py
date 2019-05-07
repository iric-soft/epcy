import os
import time
import sys

import numpy as np

from multiprocessing import Pool

from .. utils import Classifier as uc
from .. utils import readers as ur

def run_classifier(args, design, start):
    c = uc.Classifier(args, design, start)
    c.run()

    return(c)

def main_pred(args, argparser):

    if not os.path.exists(args.PATH_OUT):
        os.makedirs(args.PATH_OUT)

    file_out = args.PATH_OUT + "/prediction_capability.xls"
    file_pred_out = args.PATH_OUT + "/subgroup_predicted.xls"

    if args.ANNO is not None:
        sys.stderr.write(time.strftime('%X') + ": Read annotation\n")
        df_anno = ur.gff_2_df(args.ANNO)
        df_anno.dropna(axis=0, subset=["ID"], inplace=True)
        df_anno["ID"] = df_anno["ID"].str.replace("gene:", "")
        df_anno["ID"] = df_anno["ID"].str.replace("transcript:", "")
        df_anno["Parent"] = df_anno["Parent"].str.replace("gene:", "")
        df_anno["Parent"] = df_anno["Parent"].str.replace("transcript:", "")
        df_anno.set_index("ID", inplace=True)

    sys.stderr.write(time.strftime('%X') + ": Read design\n")
    design = ur.get_design(args)

    num_pred = sum(1 for line in open(args.MATRIX, 'r')) - 1
    # print(num_pred)

    sys.stderr.write(time.strftime('%X') + ": Start epcy analysis\n")
    if args.BY == -1:
        args.BY = int(num_pred / args.THREAD) + 1

    all_classifier = []
    if args.THREAD <= 1:
        tasks = [(args, design, start) for start in np.arange(0, num_pred, args.BY)]
        for task in tasks:
            all_classifier.append(run_classifier(task[0], task[1], task[2]))

    else:
        with Pool(processes=args.THREAD) as p:
            tasks = [(args, design, start) for start in np.arange(0, num_pred, args.BY)]

            jobs = [p.apply_async(run_classifier, t) for t in tasks]
            all_classifier = [job.get(None) for job in jobs]
            p.close()
            p.join()

    sys.stderr.write(time.strftime('%X') + ": Save epcy results\n")
    with open(file_out, 'w') as w_csv:
        all_classifier[0].print_feature_header(w_csv)
        for classifier in all_classifier:
            classifier.print_feature_pred(w_csv)

    with open(file_pred_out, 'w') as w_csv:
        all_classifier[0].print_subgroup_header(w_csv)
        for classifier in all_classifier:
            classifier.print_subgroup_predicted(w_csv, "kernel")

    sys.stderr.write(time.strftime('%X') + ": End\n")
