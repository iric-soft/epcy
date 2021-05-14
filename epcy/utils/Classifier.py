import sys
import math
import time
import os

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu, ttest_ind

from multiprocessing import Pool, RawArray
from collections import defaultdict

from itertools import product

import numexpr as ne
ne.set_num_threads(1)


def print_memory(fn):
    def wrapper(*args, **kwargs):
        print(psutil.virtual_memory())
        try:
            return fn(*args, **kwargs)
        finally:
            print(psutil.virtual_memory())
    return wrapper


def rm_missing(feature_data, num_query):
    ids_na = np.isnan(feature_data)
    if sum(ids_na) > 0:
        feature_data = feature_data[~np.isnan(feature_data)]
        num_query = num_query - sum(ids_na[:num_query])

    return(feature_data, num_query, np.where(ids_na)[0])


def get_foldchange(feature_data, num_query):
    mean_query = np.mean(feature_data[:num_query])
    mean_ref = np.mean(feature_data[num_query:])
    log2fc = mean_query - mean_ref

    return(log2fc, mean_query, mean_ref)


def auc_u_test(feature_data, num_query, num_ref):
    # print(feature_data)
    (u_value, p_value) = mannwhitneyu(feature_data[:num_query],
                                      feature_data[num_query:],
                                      alternative="two-sided")
    auc = u_value / (num_query * num_ref)

    return(auc, p_value)


def t_test_welch(feature_data, num_query):
    (t_value, p_value) = ttest_ind(feature_data[:num_query],
                                   feature_data[num_query:],
                                   equal_var=False)

    return(p_value)


def pred_fill_cont_table_normal(feature_data, num_query, n_folds, draws,
                                n_bagging=1, num_bs=0, random_seed=None):
    # Compute sample assignation using normal dist
    N = len(feature_data)
    return(get_ct_using_fx_normal(
        feature_data, num_query, N, n_folds, draws,
        n_bagging=n_bagging, num_bs=num_bs, random_seed=random_seed)
    )


def compute_kernel_fx(all_k, num_query):
    ''' return (1/(n*bw)) * sum k((x-xi) / bw) '''
    k_query = all_k[:num_query]
    k_ref = all_k[num_query:]

    sum_llh = k_query.sum() + k_ref.sum()
    if sum_llh == 0:
        res = 0.5
    else:
        res = k_query.sum() / sum_llh

    return(res)


def get_class_using_fx_kernel(feature_data, num_query, min_bw,
                              n_folds, folds_reorder, draws,
                              bw=None, n_bagging=1, num_bs=0,
                              random_seed=None):

    # get class [fold, sample_in_fold, bagging, draw]
    lst_sample_class = [get_class_gaussian_kernel_for_all_x(
                            feature_data, x_ids, num_query, draws, min_bw,
                            bw=bw, n_bagging=n_bagging, num_bs=num_bs,
                            random_seed=random_seed)
                    for x_ids in n_folds]

    n_sample = feature_data.shape[0]
    n_draw = draws.shape[0]
    # [sample, bagging, draw]
    sample_class = np.ndarray((n_sample, n_bagging, n_draw), dtype=np.bool)

    cpt_sample = 0
    for i, x_ids in enumerate(n_folds):
        for j, id in enumerate(x_ids):
            sample_class[cpt_sample, :, :] = np.asarray(
                lst_sample_class[i][j], dtype=np.bool)
            cpt_sample += 1

    if folds_reorder is not None:
        sample_class = sample_class[folds_reorder, :, :]

    return(sample_class)


def get_class_using_fx_normal(feature_data, num_query, n_folds, folds_reorder, draws,
                           n_bagging=1, num_bs=0, random_seed=None):

    # get class [fold, sample_in_fold, bagging]
    lst_sample_class = [get_class_normal_for_all_x(
                        feature_data, x_ids, num_query, draws,
                        n_bagging=n_bagging, num_bs=num_bs,
                        random_seed=random_seed)
                    for x_ids in n_folds]

    n_sample = feature_data.shape[0]
    n_draw = draws.shape[0]
    # [sample, bagging, draw]
    sample_class = np.ndarray((n_sample, n_bagging, n_draw), dtype=np.bool)

    cpt_sample = 0
    for i, x_ids in enumerate(n_folds):
        for j, id in enumerate(x_ids):
            sample_class[cpt_sample, :, :] = np.asarray(
                lst_sample_class[i][j], dtype=np.bool)
            cpt_sample += 1

    if folds_reorder is not None:
        sample_class = sample_class[folds_reorder, :, :]

    return(sample_class)


def get_classification(fx, draws):
    """return: draw < fx"""

    return([draw < fx for draw in draws])


def bw_var(x):
    return(np.var(x))


def bw_nrd(x, num_bs=0):
    # TODO need to improve speed of this part
    if num_bs != 0:
        x = x[range(1, len(x), num_bs)]

    hi = np.std(x)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    lo = min(hi, iqr/1.34)
    if (lo == 0):
        lo = hi
        if (lo == 0):
            lo = abs(x[1])
            if (lo == 0):
                lo = 1

    # this function can be run by ne.evaluate, with all lo pre-computed
    return(1.06 * lo * len(x)**(-0.2))


def bw_nrd0(x, num_bs=0):
    # TODO need to improve speed of this par
    if num_bs != 0:
        x = x[range(1, len(x), num_bs)]

    hi = np.std(x)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    lo = min(hi, iqr/1.34)
    if (lo == 0):
        lo = hi
        if (lo == 0):
            lo = abs(x[1])
            if (lo == 0):
                lo = 1

    # this function can be run by ne.evaluate, with all lo pre-computed
    return(0.9 * lo * len(x)**(-0.2))


def get_bagging_other(other, num_query, random_state=None):
    if random_state is None:
        random_state = np.random.RandomState()

    while True:
        # TODO fix randon seed to replicate same results
        ids = np.sort(random_state.choice(len(other), len(other)))
        # print(ids)
        bag_num_query = np.where(ids < num_query)[0].shape[0]
        if bag_num_query >= 2 and bag_num_query <= len(other) - 2:
            break

    return(other[ids], bag_num_query)


def get_class_gaussian_kernel_for_all_x(feature_data, x_ids, num_query, draws,
                                        min_bw, bw=None, n_bagging=1, num_bs=0,
                                        random_seed=None):

    x_values = feature_data[x_ids]
    other = np.delete(feature_data, x_ids)
    o_num_query = num_query - np.sum(x_ids < num_query)

    return([get_class_gaussian_kernel(
                    x, other, o_num_query, draws, min_bw,
                    bw, n_bagging, num_bs=num_bs, random_seed=random_seed)
            for x in x_values])


def get_class_gaussian_kernel(x, other, num_query, draws,
                              min_bw, bw, n_bagging,
                              num_bs=0, random_seed=None):
    if n_bagging > 1:
        random_state = np.random.RandomState()
        if random_state is not None:
            random_state = np.random.RandomState(random_seed)

        bag_others = [get_bagging_other(other, num_query, random_state)
                      for j in range(n_bagging)]

        return([get_classification(
                    compute_kernel_fx(
                        k_gaussian_kernel(
                            x, other, min_bw,
                            b_num_query, bw,
                            num_bs=num_bs),
                        b_num_query),
                    draws)
                for other, b_num_query in bag_others])
    else:
        return([get_classification(
                    compute_kernel_fx(
                        k_gaussian_kernel(
                            x, other, min_bw,
                            num_query, bw, num_bs=num_bs),
                        num_query),
                    draws)
                for j in range(1)])


def k_gaussian_kernel(x, other, min_bw, ids_split, bw, num_bs=0):
    other_query = other[:ids_split]
    other_ref = other[ids_split:]

    if bw is None:
        # bw_query = bw_nrd0(other_query)
        # bw_ref = bw_nrd0(other_ref)
        bw_query = bw_nrd(other_query, num_bs)
        bw_ref = bw_nrd(other_ref, num_bs)

        if bw_query < min_bw:
            bw_query = min_bw
        if bw_ref < min_bw:
            bw_ref = min_bw
    else:
        bw_query = bw
        bw_ref = bw

    norm_query = other_query.size * bw_query
    norm_ref = other_ref.size * bw_ref

    res_query = ne.evaluate('(0.3989423 * exp(-1/2*(((x - other_query) / bw_query)**2)))/norm_query')
    res_ref = ne.evaluate('(0.3989423 * exp(-1/2*(((x - other_ref) / bw_ref)**2)))/norm_ref')

    return(np.concatenate((res_query, res_ref)))


def get_class_normal_for_all_x(feature_data, x_ids, num_query, draws,
                                epsilon=0.001, n_bagging=1, num_bs=0,
                                random_seed=None):
    x_values = feature_data[x_ids]
    other = np.delete(feature_data, x_ids)
    o_num_query = num_query - np.sum(x_ids < num_query)

    return([get_class_normal(
                x, other, o_num_query, draws,
                epsilon, n_bagging, num_bs=num_bs, random_seed=random_seed)
            for x in x_values])


def get_class_normal(x, other, num_query, draws, epsilon, n_bagging,
                      num_bs=0, random_seed=None):
    if n_bagging > 1:
        random_state = np.random.RandomState()
        if random_state is not None:
            random_state = np.random.RandomState(random_seed)

        bag_others = (get_bagging_other(other, num_query, random_state)
                      for j in range(n_bagging))

        return([get_classification(
                    fx_normal(
                        x, other, b_num_query,
                        epsilon=epsilon, num_bs=num_bs),
                    draws)
                for other, b_num_query in bag_others])
    else:
        return([get_classification(
                    fx_normal(
                        x, other, num_query,
                        epsilon=epsilon, num_bs=num_bs),
                    draws)
                for j in range(1)])


def fx_normal(x, other, id_split, epsilon, num_bs=0):
    mu = np.mean(other[:id_split])
    var = np.var(other[:id_split]) + epsilon
    first_part = 1 / math.sqrt(2 * np.pi * var)
    fx_query = first_part * np.exp(-((x-mu)**2)/(2*var))

    mu = np.mean(other[id_split:])
    var = np.var(other[id_split:]) + epsilon
    first_part = 1 / math.sqrt(2 * np.pi * var)
    fx_ref = first_part * np.exp(-((x-mu)**2)/(2*var))

    if fx_query + fx_ref == 0:
        return(0.5)

    return(fx_query / (fx_query + fx_ref))


def get_mcc_pred(sample_class, num_query, args):
    # sample_class [sample_class, bagging, draw]

    n_sample = sample_class.shape[0]
    n_bag = sample_class.shape[1]
    n_draw = sample_class.shape[2]

    cont_tables = []
    pclass_by_sample = np.zeros(
        shape=(n_sample),
        dtype=np.float64)

    for i in range(n_bag):
        for j in range(n_draw):
            tp_fn = np.add.reduceat(sample_class[:, i, j], [0, num_query])

            cont_table = [tp_fn[0], num_query-tp_fn[0],
                          tp_fn[1], n_sample - num_query - tp_fn[1]]
            cont_tables.append(cont_table)

            pclass_by_sample[np.where(
                                sample_class[:num_query, i, j] == 1
                            )] += 1
            pclass_by_sample[np.array(
                                np.where(sample_class[num_query:, i, j] == 0)
                            ) + num_query] -= 1

    pclass_by_sample[:num_query] = pclass_by_sample[:num_query] / (n_draw * n_bag)
    pclass_by_sample[num_query:] = pclass_by_sample[num_query:] / (n_draw * n_bag)

    all_score = defaultdict(list)
    for ct in cont_tables:
        all_score['mcc'].append(get_mcc(ct))
        if args.PPV:
            all_score['ppv'].append(get_ppv(ct))
        if args.NPV:
            all_score['npv'].append(get_npv(ct))
        if args.TPR:
            all_score['tpr'].append(get_tpr(ct))
        if args.TNR:
            all_score['tnr'].append(get_tnr(ct))
        if args.FNR:
            all_score['fnr'].append(get_fnr(ct))
        if args.FPR:
            all_score['fpr'].append(get_fpr(ct))
        if args.FDR:
            all_score['fdr'].append(get_fdr(ct))
        if args.FOR:
            all_score['for'].append(get_for(ct))
        if args.TS:
            all_score['ts'].append(get_ts(ct))
        if args.ACC:
            all_score['acc'].append(get_acc(ct))
        if args.F1:
            all_score['f1'].append(get_f1(ct))

    dict_res = defaultdict(list)
    all_score['mcc'] = np.sort(all_score['mcc'])
    num_value = all_score['mcc'].size
    first_quantile = int(num_value * 0.05)
    last_quantile = int(num_value * 0.95)
    if last_quantile != 0:
        last_quantile = last_quantile - 1

    # set score name prefix
    score_name = "kernel_"
    if args.NORMAL:
        score_name = "normal_"

    for key, value in all_score.items():
        tmp = np.sort(all_score[key])
        mean_score = np.mean(tmp[first_quantile:(last_quantile+1)])
        dict_res[score_name + key] = [
            tmp[first_quantile],
            mean_score, tmp[last_quantile]
        ]

    return(dict_res, pclass_by_sample)


def get_mcc(ct):
    # 0: TP, 1:FN, 2:FP, 3:TN
    d1 = ct[0] + ct[1]
    d2 = ct[0] + ct[2]
    d3 = ct[3] + ct[1]
    d4 = ct[3] + ct[2]

    n1 = ct[0] * ct[3]
    n2 = ct[1] * ct[2]

    return((n1 - n2) / math.sqrt(d1 * d2 * d3 * d4)
           if d1 != 0 and d2 != 0 and d3 != 0 and d4 != 0 else 0)


def get_ppv(ct):
    # 0: TP, 1:FN, 2:FP, 3:TN
    # PPV = TP / (TP + FP)
    deno = ct[0] + ct[2]
    return(ct[0] / deno if deno != 0 else 0)


def get_npv(ct):
    # 0: TP, 1:FN, 2:FP, 3:TN
    # NPV = TN / (TN + FN)
    deno = ct[3] + ct[1]
    return(ct[3] / deno if deno != 0 else 0)


def get_tpr(ct):
    # 0: TP, 1:FN, 2:FP, 3:TN
    # TPR = TP / (TP + FN)
    deno = ct[0] + ct[1]
    return(ct[0] / deno if deno != 0 else 0)


def get_tnr(ct):
    # 0: TP, 1:FN, 2:FP, 3:TN
    # TNR = TN / (TN + FP)
    deno = ct[3] + ct[2]
    return(ct[3] / deno if deno != 0 else 0)


def get_fnr(ct):
    # 0: TP, 1:FN, 2:FP, 3:TN
    # FNR = FN / (FN + TP)
    deno = ct[1] + ct[0]
    return(ct[1] / deno if deno != 0 else 0)


def get_fpr(ct):
    # 0: TP, 1:FN, 2:FP, 3:TN
    # FPR = FP / (FP + TN)
    deno = (ct[2] + ct[3])
    return(ct[2] / deno if deno != 0 else 0)


def get_fdr(ct):
    # 0: TP, 1:FN, 2:FP, 3:TN
    # FDR = FP / (FP + TP)
    deno = ct[2] + ct[0]
    return(ct[2] / deno if deno != 0 else 0)


def get_for(ct):
    # 0: TP, 1:FN, 2:FP, 3:TN
    # FOR = FN / (FN + TN)
    deno = ct[1] + ct[3]
    return(ct[1] / deno if deno != 0 else 0)


def get_ts(ct):
    # 0: TP, 1:FN, 2:FP, 3:TN
    # TS = TP / (TP + FN + FP)
    deno = ct[0] + ct[1] + ct[2]
    return(ct[0] / deno if deno != 0 else 0)


def get_acc(ct):
    # 0: TP, 1:FN, 2:FP, 3:TN
    # ACC = (TP + TN) / (TP + TN + FP + FN)
    deno = ct[0] + ct[1] + ct[2] + ct[3]
    return((ct[0] + ct[3]) / deno if deno != 0 else 0)


def get_f1(ct):
    # 0: TP, 1:FN, 2:FP, 3:TN
    # F1 = 2TP / (2TP + FP + FN)
    deno = (2 * ct[0] + ct[2] + ct[1])
    return((2 * ct[0]) / deno if deno != 0 else 0)


def pred_feature(feature_data, num_query, num_ref,
                 n_folds, folds_reorder, draws, num_bs, args):

    dict_res = defaultdict(list)
    ids_na = None
    if np.isnan(feature_data).sum() > 0:
        num_ids = len(feature_data)
        # remove missing values
        feature_data, num_query, ids_na = rm_missing(feature_data,
                                                          num_query)
        # keep number of sample query and ref for the output
        # (after removing sample with missing value)
        num_ref = feature_data.size - num_query
        dict_res['sample_query'].append(num_query)
        dict_res['sample_ref'].append(num_ref)

        if num_query <= 2 or num_ref <= 2:
            return(dict_res)

        # update nfold decomposition
        ids2del = range(num_ids-len(ids_na), num_ids, 1)
        n_folds_saved = n_folds
        n_folds = []
        for fold in n_folds_saved:
            new_fold = np.setdiff1d(fold, ids2del)
            if len(new_fold) > 0:
                n_folds.append(new_fold)

        # update folds_reorder
        if folds_reorder is not None:
            folds_reorder = np.setdiff1d(folds_reorder, ids2del)

    sum_row = sum(feature_data)
    res = get_foldchange(feature_data, num_query)
    l2fc = res[0]
    dict_res['l2fc'].append(l2fc)
    dict_res['mean_query'].append(res[1])
    dict_res['mean_ref'].append(res[2])
    dict_res['bw_query'].append(
        bw_nrd(
            feature_data[:num_query],
            num_bs=num_bs
        )
    )
    dict_res['bw_ref'].append(
        bw_nrd(
            feature_data[num_query:],
            num_bs=num_bs
        )
    )

    if np.unique(feature_data).size == 1:
        return(dict_res)

    if args.EXP is not None and sum_row < args.EXP:
        return(dict_res)

    if abs(l2fc) < args.L2FC:
        return(dict_res)

    if args.AUC:
        res = auc_u_test(feature_data, num_query, num_ref)
        dict_res['auc'].append(res[0])

        if args.UTEST:
            dict_res['utest_pv'].append(res[1])

    if args.TTEST:
        dict_res['ttest_pv'].append(
            t_test_welch(
                feature_data,
                num_query
            )
        )

    if not args.NORMAL:
        sample_class = get_class_using_fx_kernel(
            feature_data, num_query, args.MIN_BW,
            n_folds, folds_reorder, draws,
            n_bagging=args.N_BAGGING, num_bs=num_bs,
            random_seed=args.RANDOM_SEED
        )
    else:
        sample_class = get_class_using_fx_normal(
            feature_data, num_query, n_folds, folds_reorder, draws,
            n_bagging=args.N_BAGGING, num_bs=num_bs,
            random_seed=args.RANDOM_SEED
        )

    dict_score, pred_by_sample = get_mcc_pred(sample_class, num_query, args)

    for key, value in dict_score.items():
        dict_res[key] = dict_score[key]

    if args.FULL:
        if ids_na is not None and len(ids_na) > 0:
            pred_by_sample = np.insert(pred_by_sample, ids_na, np.nan)
        dict_res['pred_by_sample'] = pred_by_sample

    return(dict_res)


def init_worker(raw_array, shape, dtype,
                num_query, num_ref, n_folds, folds_reorder,
                draws, num_bs, args):

    # The shared array pointer is a global variable so that it can be accessed by the
    # child processes. It is a tuple (pointer, dtype, shape).
    global shared_arr
    shared_arr = {}
    shared_arr['array'] = raw_array
    shared_arr['shape'] = shape
    shared_arr['dtype'] = dtype
    shared_arr['num_query'] = num_query
    shared_arr['num_ref'] = num_ref
    shared_arr['n_folds'] = n_folds
    shared_arr['folds_reorder'] = folds_reorder
    shared_arr['draws'] = draws
    shared_arr['num_bs'] = num_bs
    shared_arr['args'] = args


def worker_func(i):
    # Simply computes the sum of the i-th row of the input matrix X
    feature_data = np.frombuffer(
        shared_arr['array'], dtype=shared_arr['dtype'],
        offset=i * shared_arr['shape'][1] * 8,
        count=shared_arr['shape'][1]
    )

    return(
        pred_feature(
            feature_data,
            shared_arr['num_query'], shared_arr['num_ref'],
            shared_arr['n_folds'], shared_arr['folds_reorder'],
            shared_arr['draws'], shared_arr['num_bs'],
            shared_arr['args']
        )
    )


#def worker_shared_func(i, shm_name, num_query, num_ref, n_folds, folds_reorder,
#                       draws, num_bs, args):

#    existing_shm = shared_memory.SharedMemory(name=shm_name)
#    feature_data = np.ndarray(
#        (num_query+num_ref,),
#        dtype=np.float64,
#        buffer=existing_shm.buf,
#        offset=i * (num_query+num_ref) * 8
#    )

#    res = pred_feature(
#        feature_data, num_query, num_ref,
#        n_folds, folds_reorder, draws, num_bs, args
#    )

#    del feature_data
#    existing_shm.close()

#    return(res)


class Classifier:
    def __init__(self, args, design, data, list_ids, n_folds, draws, folds_reorder):
        self.args = args
        self.design = design
        self.data = data
        self.list_ids = list_ids
        self.n_folds = n_folds
        self.draws = draws
        self.folds_reorder = folds_reorder

        self.num_query = len(np.where(design[self.args.CONDITION] == 1)[0])
        self.num_ref = len(np.where(design[self.args.CONDITION] == 0)[0])

        self.with_na = np.isnan(data).sum()
        self.done = False

        self.result = []

    def run(self):
        self.__pred()

    def pred2csv(self):

        sys.stderr.write(time.strftime('%X') + ": Save epcy results\n")
        if self.args.PATH_OUT is not None:
            if not os.path.exists(self.args.PATH_OUT):
                os.makedirs(self.args.PATH_OUT)

            file_out = self.args.PATH_OUT + "/predictive_capability.xls"
            file_pred_out = self.args.PATH_OUT + "/condition_predicted.xls"
            file_pred_normal_out = self.args.PATH_OUT + "/condition_predicted_normal.xls"

            with open(file_out, 'w') as w_csv:
                self.print_feature_header(w_csv, self.args, self.with_na > 0)
                self.print_feature_pred(
                    self.result, self.list_ids,
                    self.num_query, self.num_ref,
                    w_csv, self.args, self.with_na > 0
                )

            if self.args.FULL:
                with open(file_pred_out, 'w') as w_csv:
                    self.print_condition_header(w_csv, self.design)
                    self.print_condition_predicted(
                        self.result, self.list_ids, w_csv, "kernel")
                if self.args.NORMAL:
                    with open(file_pred_normal_out, 'w') as w_csv:
                        self.print_condition_header(w_csv,  self.design)
                        self.print_condition_predicted(
                            self.result, self.list_ids, w_csv, "normal")

        else:
            self.print_feature_header(sys.stdout, self.args, self.with_na > 0)
            self.print_feature_pred(
                self.result, self.list_ids,
                self.num_query, self.num_ref,
                sys.stdout, self.args, self.with_na > 0
            )

    @staticmethod
    def print_feature_header(w_csv, args, with_na=False):
        score_prefix = "kernel_"
        if args.NORMAL:
            score_prefix = "normal_"

        header = "id\tl2fc"
        header = header + "\t" + score_prefix + "mcc"
        header = header + "\t" + score_prefix + "mcc_low"
        header = header + "\t" + score_prefix + "mcc_high"
        if args.PPV:
            header = header + "\t" + score_prefix + "ppv"
            header = header + "\t" + score_prefix + "ppv_low"
            header = header + "\t" + score_prefix + "ppv_high"
        if args.NPV:
            header = header + "\t" + score_prefix + "npv"
            header = header + "\t" + score_prefix + "npv_low"
            header = header + "\t" + score_prefix + "npv_high"
        if args.TPR:
            header = header + "\t" + score_prefix + "tpr"
            header = header + "\t" + score_prefix + "tpr_low"
            header = header + "\t" + score_prefix + "tpr_high"
        if args.TNR:
            header = header + "\t" + score_prefix + "tnr"
            header = header + "\t" + score_prefix + "tnr_low"
            header = header + "\t" + score_prefix + "tnr_high"
        if args.FNR:
            header = header + "\t" + score_prefix + "fnr"
            header = header + "\t" + score_prefix + "fnr_low"
            header = header + "\t" + score_prefix + "fnr_high"
        if args.FPR:
            header = header + "\t" + score_prefix + "fpr"
            header = header + "\t" + score_prefix + "fpr_low"
            header = header + "\t" + score_prefix + "fpr_high"
        if args.FDR:
            header = header + "\t" + score_prefix + "fdr"
            header = header + "\t" + score_prefix + "fdr_low"
            header = header + "\t" + score_prefix + "fdr_high"
        if args.FOR:
            header = header + "\t" + score_prefix + "for"
            header = header + "\t" + score_prefix + "for_low"
            header = header + "\t" + score_prefix + "for_high"
        if args.TS:
            header = header + "\t" + score_prefix + "ts"
            header = header + "\t" + score_prefix + "ts_low"
            header = header + "\t" + score_prefix + "ts_high"
        if args.ACC:
            header = header + "\t" + score_prefix + "acc"
            header = header + "\t" + score_prefix + "acc_low"
            header = header + "\t" + score_prefix + "acc_high"
        if args.F1:
            header = header + "\t" + score_prefix + "f1"
            header = header + "\t" + score_prefix + "f1_low"
            header = header + "\t" + score_prefix + "f1_high"
        header = header + "\tmean_query\tmean_ref"
        header = header + "\tbw_query\tbw_ref"

        if args.AUC:
            header = header + "\tauc"
            if args.UTEST:
                header = header + "\tu_pv"
        if args.TTEST:
            header = header + "\tt_pv"
        if with_na:
            header = header + "\tsample_query\tsample_ref"
        header = header + "\n"

        w_csv.write(header)

    @staticmethod
    def print_feature_pred(results, list_ids, num_query,
                           num_ref, w_csv, args, with_na=False):
        cpt_id = 0

        def print_score(score_name, res):
            if score_name in res:
                scores = res[score_name]
                line = str(scores[1]) + "\t"
                line = line + str(scores[0]) + "\t"
                line = line + str(scores[2]) + "\t"
            else:
                line = "nan\tnan\tnan\t"

            return(line)

        score_prefix = "kernel_"
        if args.NORMAL:
            score_prefix = "normal_"

        for res in results:
            line = str(list_ids[cpt_id]) + "\t"
            if "l2fc" in res:
                line = line + str(res['l2fc'][0]) + "\t"
            else:
                line = line + "nan\t"

            score_name = score_prefix + "mcc"
            line = line + print_score(score_name, res)

            if args.PPV:
                score_name = score_prefix + "ppv"
                line = line + print_score(score_name, res)
            if args.NPV:
                score_name = score_prefix + "npv"
                line = line + print_score(score_name, res)
            if args.TPR:
                score_name = score_prefix + "tpr"
                line = line + print_score(score_name, res)
            if args.TNR:
                score_name = score_prefix + "tnr"
                line = line + print_score(score_name, res)
            if args.FNR:
                score_name = score_prefix + "fnr"
                line = line + print_score(score_name, res)
            if args.FPR:
                score_name = score_prefix + "fpr"
                line = line + print_score(score_name, res)
            if args.FDR:
                score_name = score_prefix + "fdr"
                line = line + print_score(score_name, res)
            if args.FOR:
                score_name = score_prefix + "for"
                line = line + print_score(score_name, res)
            if args.TS:
                score_name = score_prefix + "ts"
                line = line + print_score(score_name, res)
            if args.ACC:
                score_name = score_prefix + "acc"
                line = line + print_score(score_name, res)
            if args.F1:
                score_name = score_prefix + "f1"
                line = line + print_score(score_name, res)

            if 'mean_query' in res:
                line = line + str(res['mean_query'][0]) + "\t"
                line = line + str(res['mean_ref'][0]) + "\t"
                line = line + str(res['bw_query'][0]) + "\t"
                line = line + str(res['bw_ref'][0])
            else:
                line = line + "nan\tnan\tnan\tnan"

            if args.AUC:
                if 'auc' in res:
                    auc = res['auc'][0] if res['auc'][0] >= 0.5 else 1 - res['auc'][0]
                    line = line + "\t" + str(auc)
                    if args.UTEST:
                        line = line + "\t" + str(res['utest_pv'][0])
                else:
                    line = line + "\tnan"
                    if args.UTEST:
                        line = line + "\tnan"
            if args.TTEST:
                if 'ttest_pv' in res:
                    line = line + "\t" + str(res['ttest_pv'][0])
                else:
                    line = line + "\tnan"
            if with_na:
                if 'sample_query' in res:
                    line = line + "\t" + str(res['sample_query'][0])
                    line = line + "\t" + str(res['sample_ref'][0])
                else:
                    line = line + "\t" + str(num_query)
                    line = line + "\t" + str(num_ref)

            line = line + "\n"

            w_csv.write(line)
            cpt_id += 1

    @staticmethod
    def print_condition_header(w_csv, design):
        line = "id" + '\t'
        line = line + '\t'.join(design['sample'])
        line = line + "\n"
        w_csv.write(line)

    @staticmethod
    def print_condition_predicted(results, list_ids, w_csv, pred):
        key = 'kernel_pred'
        if pred == "normal":
            key = 'normal_pred'

        cpt_id = 0
        for res in results:
            line = str(list_ids[cpt_id]) + "\t"
            if key in res:
                line = line + '\t'.join([str(x) for x in res[key][0]])
            line = line + "\n"

            w_csv.write(line)
            cpt_id += 1

    def __pred(self):
        num_bs = 0
        if hasattr(self.args, 'BS') and self.args.BS is not None:
            num_bs = self.args.BS

        result = []
        if self.args.THREAD <= 1:
            for i in range(self.data.shape[0]):
                feature_data = self.data[i,:]
                self.result.append(
                    pred_feature(
                        feature_data, self.num_query, self.num_ref,
                        self.n_folds, self.folds_reorder, self.draws,
                        num_bs, self.args
                    )
                )
        else:
            # Use shared_memory
            #shm = shared_memory.SharedMemory(create=True, size=self.data.nbytes)
            #data_shared = np.ndarray(self.data.shape, dtype=self.data.dtype, buffer=shm.buf)
            #np.copyto(data_shared, self.data)

            #params = [(
            #        x, shm.name,
            #        self.num_query, self.num_ref,
            #        self.n_folds, self.draws, num_bs,
            #        self.args
            #    )
            #    for x in range(len(self.list_ids))
            #]
            #with Pool(processes=self.args.THREAD) as p:
            #    self.result = p.starmap(
            #        worker_shared_func,
            #        params
            #    )

            #del data_shared
            #shm.close()
            #shm.unlink()

            # Use RawArray
            # https://research.wmz.ninja/articles/2018/03/on-sharing-large-arrays-when-using-pythons-multiprocessing.html
            dtype = np.float64
            cdtype = np.ctypeslib.as_ctypes_type(dtype)
            data_shape = self.data.shape
            raw_array = RawArray(cdtype, range(data_shape[0] * data_shape[1]))
            raw_array_np = np.frombuffer(raw_array, dtype=dtype).reshape(data_shape)
            np.copyto(raw_array_np, self.data)
            del self.data
            with Pool(
                processes=self.args.THREAD, initializer=init_worker,
                initargs=(
                    raw_array, data_shape, dtype,
                    self.num_query, self.num_ref,
                    self.n_folds, self.folds_reorder, self.draws,
                    num_bs, self.args
                )
            ) as p:
                self.result = p.map(worker_func, range(data_shape[0]))
