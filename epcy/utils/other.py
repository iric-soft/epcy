import sys
import os
import time

import numpy as np

from multiprocessing import Pool
from .. utils import Classifier as uc


def run_classifier(args, design, data, list_ids, start):
    """Run a classifier.
    """
    data_selected = data[start:(start + args.BY)]
    list_ids_selected = list_ids[start:(start + args.BY)]
    c = uc.Classifier(args, design, data_selected, list_ids_selected, start)
    c.run()

    return(c)


def save_pred_res(args, all_classifier):
    """Save results compute by all classifiers.
    """
    with_na = False
    for classifier in all_classifier:
        if classifier.with_na > 0:
            with_na = True
            break

    sys.stderr.write(time.strftime('%X') + ": Save epcy results\n")
    if args.PATH_OUT is not None:
        if not os.path.exists(args.PATH_OUT):
            os.makedirs(args.PATH_OUT)

        file_out = args.PATH_OUT + "/predictive_capability.xls"
        file_pred_out = args.PATH_OUT + "/subgroup_predicted.xls"

        with open(file_out, 'w') as w_csv:
            all_classifier[0].print_feature_header(w_csv, args, with_na)
            for classifier in all_classifier:
                classifier.print_feature_pred(w_csv, with_na)

        if args.FULL:
            with open(file_pred_out, 'w') as w_csv:
                all_classifier[0].print_subgroup_header(w_csv)
                for classifier in all_classifier:
                    classifier.print_subgroup_predicted(w_csv, "kernel")
    else:
        all_classifier[0].print_feature_header(sys.stdout, args)
        for classifier in all_classifier:
            classifier.print_feature_pred(sys.stdout)

    sys.stderr.write(time.strftime('%X') + ": End\n")


def cut_version(id_str):
    """Remove the version ensembl gene/transcript ids
    """
    index_version = id_str.find('.')
    if index_version != -1:
        id_str = id_str[0:index_version]

    return(id_str)


def compute_pred(args, num_pred, list_ids, data, design):
    sys.stderr.write(time.strftime('%X') + ": Start epcy analysis of " +
                     str(num_pred) + " features\n")
    if args.BY == -1:
        args.BY = int(num_pred / args.THREAD) + 1

    all_classifier = []
    if args.THREAD <= 1:
        tasks = [(args, design, data, list_ids, start)
                 for start in np.arange(0, num_pred, args.BY)]
        for task in tasks:
            all_classifier.append(run_classifier(task[0], task[1], task[2],
                                                 task[3], task[4]))

    else:
        with Pool(processes=args.THREAD) as p:
            tasks = [(args, design, data, list_ids, start)
                     for start in np.arange(0, num_pred, args.BY)]

            jobs = [p.apply_async(run_classifier, t) for t in tasks]
            all_classifier = [job.get(None) for job in jobs]
            p.close()
            p.join()

    return(all_classifier)
