import sys
import os
import time

import numpy as np

from .. utils import Classifier as uc

def cut_version(id_str):
    """Remove the version ensembl gene/transcript ids
    """
    index_version = id_str.find('.')
    if index_version != -1:
        id_str = id_str[0:index_version]

    return(id_str)


def get_back_bs_sample(folds, num_bs):
    folds_w_bs = []
    for fold in folds:
        all_ids = []
        for id in fold:
            id_bs = id * num_bs
            for x in range(num_bs):
                all_ids.append(id_bs + x)
        folds_w_bs.append(all_ids)

    return(folds_w_bs)


def compute_pred(args, num_pred, list_ids, data, design):
    sys.stderr.write(time.strftime('%X') + ": Start epcy analysis of " +
                     str(num_pred) + " features\n")

    num_bs = 0
    if hasattr(args, 'BS') and args.BS is not None:
        num_bs = args.BS

    random_state = np.random.RandomState()
    if args.RANDOM_SEED is not None:
        random_state = np.random.RandomState(args.RANDOM_SEED)

    # Prepare n-fold (or LOO) used for each gene
    folds = []
    folds_reorder = None
    N = design.shape[0]
    if args.N_FOLD is None:
        if num_bs > 0:
            folds = np.array_split(np.arange(N), N/num_bs)
        else:
            folds = np.array_split(np.arange(N), N)
    else:
        if num_bs > 0:
            nfold_decomposition = np.arange(N/num_bs)
            random_state.shuffle(nfold_decomposition)
            length = int(len(nfold_decomposition)/args.N_FOLD) #length of each fold
            for i in range(args.N_FOLD-1):
                folds += [nfold_decomposition[i*length:(i+1)*length]]
            folds += [nfold_decomposition[(args.N_FOLD-1)*length:len(nfold_decomposition)]]
            folds = get_back_bs_sample(folds, num_bs)

        else:
            nfold_decomposition = np.arange(N)
            random_state.shuffle(nfold_decomposition)
            length = int(len(nfold_decomposition)/args.N_FOLD) #length of each fold
            for i in range(args.N_FOLD-1):
                folds += [nfold_decomposition[i*length:(i+1)*length]]
            folds += [nfold_decomposition[(args.N_FOLD-1)*length:len(nfold_decomposition)]]

        if len(folds) != N:
            n_folds = np.ndarray(N, dtype=np.int)
            cpt_sample = 0
            for i, fold in enumerate(folds):
                for j, id in enumerate(fold):
                    n_folds[cpt_sample] = id
                    cpt_sample += 1
            folds_reorder = np.argsort(n_folds)

    # Prepare draws used for each gene
    draws = random_state.random(args.N_DRAW)

    c = uc.Classifier(args, design, data, list_ids, folds, draws, folds_reorder)
    del design
    del data
    c.run()
    c.pred2csv()

    sys.stderr.write(time.strftime('%X') + ": End\n")
    return(None)
