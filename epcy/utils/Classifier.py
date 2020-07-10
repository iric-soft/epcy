import sys
import math
import time
import os

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu, ttest_ind
from statistics import median

from multiprocessing import Pool, RawArray
from collections import defaultdict

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


def pred_fill_cont_table_kernel(feature_data, num_query, min_bw, bw=None,
                                n_bagging=1, num_bs=0, num_draw=100,
                                random_state=np.random.RandomState()):
    # Compute sample assignation using kernel
    N = len(feature_data)
    return(get_ct_using_fx_kernel(
        feature_data, num_query, N, min_bw, bw=bw, n_bagging=n_bagging,
        num_bs=num_bs, num_draw=num_draw, random_state=random_state)
    )


def pred_fill_cont_table_normal(feature_data, num_query, n_bagging=1,
                                num_draw=100, num_bs=0,
                                random_state=np.random.RandomState()):
    # Compute sample assignation using normal dist
    N = len(feature_data)
    return(get_ct_using_fx_normal(
        feature_data, num_query, N, n_bagging=n_bagging, num_draw=num_draw,
        random_state=random_state, num_bs=num_bs)
    )


def compute_kernel_fx(all_k, i, num_bs):
    ''' return (1/(n*bw)) * sum k((x-xi) / bw) '''
    # TODO: clean this
    all_k, num_query = all_k

    k_query = all_k[:num_query]
    k_ref = all_k[num_query:]

    sum_llh = k_query.sum() + k_ref.sum()
    if sum_llh == 0:
        res = 0.5
    else:
        res = k_query.sum() / sum_llh

    return(res)


def compute_kernel_fx_and_ct(k_bw_nq_gen, num_query, N, num_bs=0,
                             num_draw=100,
                             random_state=np.random.RandomState()):
    fx_by_sample = [compute_kernel_fx(k_bw_nq, i, num_bs=num_bs)
                    for i, k_bw_nq in enumerate(k_bw_nq_gen)]

    return(
        fx_to_tables(
            fx_by_sample, num_query, N,
            num_draw=num_draw, random_state=random_state
        )
    )


def get_ct_using_fx_kernel(feature_data, num_query, N, min_bw, bw=None,
                           n_bagging=1, num_bs=0, num_draw=100,
                           random_state=np.random.RandomState()):
    # Create leave one out index and manage kallisto bootstrap
    # TODO: implement n_folds
    n_folds = None
    if num_bs > 0:
        n_folds = np.array_split(np.arange(N), N/num_bs)
    else:
        n_folds = np.array_split(np.arange(N), N)

    # Compute k((x-xi) / bw) for each leave one out
    k_bw_gen_by_fold = [get_k_gaussian_kernel_for_all_x(
                            feature_data, x_ids, num_query, min_bw,
                            bw=bw, n_bagging=n_bagging, num_bs=num_bs,
                            random_state=random_state)
                        for x_ids in n_folds]

    k_bw_gen_by_bag = np.transpose(
                        np.asarray(k_bw_gen_by_fold),
                        (2, 0, 1, 3))
    k_bw_gen_by_bag = np.reshape(k_bw_gen_by_bag, (n_bagging, N, 2))

    # (1/(n*bw)) * sum k((x-xi) / bw)
    ct_by_bag = np.asarray([compute_kernel_fx_and_ct(
                                k_bw_nq_gen, num_query, N,
                                num_bs, num_draw=num_draw,
                                random_state=random_state)
                            for k_bw_nq_gen in k_bw_gen_by_bag])

    return(ct_by_bag)


def get_ct_using_fx_normal(feature_data, num_query, N, n_bagging=1, num_bs=0,
                           num_draw=100,
                           random_state=np.random.RandomState()):
    # Create leave one out index
    n_folds = None
    if num_bs > 0:
        n_folds = np.array_split(np.arange(N), N/num_bs)
    else:
        n_folds = np.array_split(np.arange(N), N)

    fx_by_fold = [compute_normal_fx_for_all_x(
                    feature_data, x_ids, num_query,
                    n_bagging=n_bagging, num_bs=num_bs,
                    random_state=random_state)
                  for x_ids in n_folds]

    fx_by_fold = np.asarray(fx_by_fold)
    fx_by_bag = np.transpose(fx_by_fold, (2, 0, 1))
    fx_by_bag = np.reshape(fx_by_bag, (n_bagging, N))

    ct_by_bag = (fx_to_tables(
                    fx_by_sample, num_query, N,
                    num_draw=num_draw,
                    random_state=random_state)
                 for fx_by_sample in fx_by_bag)

    return(ct_by_bag)


def fx_to_tables(fx_by_sample, num_query, N, num_draw=100,
                 random_state=np.random.RandomState()):
    """return:
       cont_table[tp, fp, fn, tn]
       pred_by_sample: 1=tp, 2=fn, 3=fp, tn=4"""

    cont_tables = []
    pclass_by_sample = np.zeros(
        shape=(len(fx_by_sample)),
        dtype=np.float16)
    for i in range(num_draw):
        pred_by_sample = np.array([random_state.random() < fx
                                  for fx in fx_by_sample])
        tp_fn = np.add.reduceat(pred_by_sample, [0, num_query])

        cont_table = [tp_fn[0], num_query-tp_fn[0],
                      tp_fn[1], N - num_query - tp_fn[1]]
        cont_tables.append(cont_table)

        pclass_by_sample[np.where(pred_by_sample[:num_query] == 1)] += 1
        pclass_by_sample[np.array(
                        np.where(pred_by_sample[num_query:] == 0)
                     ) + num_query] -= 1

    pclass_by_sample[:num_query] = pclass_by_sample[:num_query] / num_draw
    pclass_by_sample[num_query:] = pclass_by_sample[num_query:] / num_draw
    return(cont_tables, pclass_by_sample)


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


def get_bagging_other(other, num_query,
                      random_state=np.random.RandomState()):
    while True:
        # TODO fix randon seed to replicate same results
        ids = np.sort(random_state.choice(len(other), len(other)))
        # print(ids)
        bag_num_query = np.where(ids < num_query)[0].shape[0]
        if bag_num_query >= 2 and bag_num_query <= len(other) - 2:
            break

    return(other[ids], bag_num_query)


def get_k_gaussian_kernel_for_all_x(feature_data, x_ids, num_query, min_bw,
                                    bw=None, n_bagging=1, num_bs=0,
                                    random_state=np.random.RandomState()):
    x_values = feature_data[x_ids]
    other = np.delete(feature_data, x_ids)
    o_num_query = num_query - np.sum(x_ids < num_query)

    return([get_k_gaussian_kernel(
                x, other, o_num_query, min_bw,
                bw, n_bagging, random_state,
                num_bs=num_bs)
            for x in x_values])


def get_k_gaussian_kernel(x, other, num_query, min_bw, bw, n_bagging,
                          random_state, num_bs=0):
    if n_bagging > 1:
        bag_others = [get_bagging_other(
                        other, num_query,
                        random_state=random_state)
                      for j in range(n_bagging)]

        return([k_gaussian_kernel(
                    x, other, min_bw,
                    b_num_query, bw,
                    num_bs=num_bs)
                for other, b_num_query in bag_others])
    else:
        return([k_gaussian_kernel(
                    x, other, min_bw,
                    num_query, bw, num_bs=num_bs)
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

    return(np.concatenate((res_query, res_ref)), ids_split)


def compute_normal_fx_for_all_x(feature_data, x_ids, num_query, epsilon=0.001,
                                n_bagging=1, num_bs=0,
                                random_state=np.random.RandomState()):
    x_values = feature_data[x_ids]
    other = np.delete(feature_data, x_ids)
    o_num_query = num_query - np.sum(x_ids < num_query)

    return([compute_normal_fx(
                x, other, o_num_query,
                epsilon, n_bagging, num_bs=num_bs,
                random_state=random_state)
            for x in x_values])


def compute_normal_fx(x, other, num_query, epsilon, n_bagging,
                      random_state, num_bs=0):
    if n_bagging > 1:
        bag_others = (get_bagging_other(
                        other, num_query,
                        random_state=random_state)
                      for j in range(n_bagging))

        return([fx_normal(
                    x, other, b_num_query,
                    epsilon=epsilon, num_bs=num_bs)
                for other, b_num_query in bag_others])
    else:
        return([fx_normal(
                    x, other, num_query,
                    epsilon=epsilon, num_bs=num_bs)
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


def get_mcc_ped_by_bagging(ct_by_bagging):
    all_mcc = []
    all_pred = []
    for cts, pred in ct_by_bagging:
        for ct in cts:
            all_mcc.append(get_mcc(ct))
        all_pred.append(pred)
    pred_by_sample = np.mean(np.asarray(all_pred), axis=0)

    all_mcc = np.sort(all_mcc)
    num_value = all_mcc.size
    first_quantile = int(num_value * 0.05)
    last_quantile = int(num_value * 0.95)

    if last_quantile != 0:
        last_quantile = last_quantile - 1
    mean_mcc = np.mean(all_mcc[first_quantile:(last_quantile+1)])

    return([all_mcc[first_quantile], mean_mcc,
            all_mcc[last_quantile]], pred_by_sample)


def get_mcc(ct):
    d1 = ct[0] + ct[1]
    d2 = ct[0] + ct[2]
    d3 = ct[3] + ct[1]
    d4 = ct[3] + ct[2]

    n1 = ct[0] * ct[3]
    n2 = ct[1] * ct[2]

    return((n1 - n2) / math.sqrt(d1 * d2 * d3 * d4)
           if d1 != 0 and d2 != 0 and d3 != 0 and d4 != 0 else 0)


def pred_feature(feature_data, num_query, num_ref, num_bs, args, random_state):
    dict_res = defaultdict(list)
    ids_na = None
    if np.isnan(feature_data).sum() > 0:
        feature_data, num_query, ids_na = rm_missing(feature_data,
                                                          num_query)

        num_ref = feature_data.size - num_query
        dict_res['sample_query'].append(num_query)
        dict_res['sample_ref'].append(num_ref)

        if num_query <= 2 or num_ref <= 2:
            return(dict_res)

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

    if args.EXP is not None and sum_row >= args.EXP:
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

    ct_by_bagging_and_pred_by_sample = pred_fill_cont_table_kernel(
        feature_data, num_query, args.MIN_BW,
        n_bagging=args.N_BAGGING, num_bs=num_bs,
        num_draw=args.N_DRAW,
        random_state=random_state
    )

    mcc, pred_by_sample = get_mcc_ped_by_bagging(
        ct_by_bagging_and_pred_by_sample
    )

    dict_res['kernel_mcc'].append(mcc)

    if args.FULL:
        if ids_na is not None and len(ids_na) > 0:
            pred_by_sample = np.insert(pred_by_sample, ids_na, np.nan)
        dict_res['kernel_pred'].append(pred_by_sample)

    if args.NORMAL:
        ct_by_bagging_and_pred_by_sample = pred_fill_cont_table_normal(
            feature_data, num_query, n_bagging=args.N_BAGGING,
            num_draw=args.N_DRAW,
            random_state=random_state,
            num_bs=num_bs
        )

        mcc, pred_by_sample = get_mcc_ped_by_bagging(
            ct_by_bagging_and_pred_by_sample
        )

        dict_res['normal_mcc'].append(mcc)
        if args.FULL:
            if ids_na is not None and len(ids_na) > 0:
                pred_by_sample = np.insert(pred_by_sample, ids_na, np.nan)
            dict_res['normal_pred'].append(pred_by_sample)

    return(dict_res)


def init_worker(raw_array, shape, dtype,
                num_query, num_ref, num_bs, args, random_state):

    # The shared array pointer is a global variable so that it can be accessed by the
    # child processes. It is a tuple (pointer, dtype, shape).
    global shared_arr
    shared_arr = {}
    shared_arr['array'] = raw_array
    shared_arr['shape'] = shape
    shared_arr['dtype'] = dtype
    shared_arr['num_query'] = num_query
    shared_arr['num_ref'] = num_ref
    shared_arr['num_bs'] = num_bs
    shared_arr['args'] = args
    shared_arr['random_state'] = random_state


def worker_func(i):
    # Simply computes the sum of the i-th row of the input matrix X
    raw_array_np = np.frombuffer(
        shared_arr['array'], dtype=shared_arr['dtype']
    ).reshape(shared_arr['shape'])

    feature_data = raw_array_np[i,:]
    return(
        pred_feature(
            feature_data,
            shared_arr['num_query'], shared_arr['num_ref'],
            shared_arr['num_bs'], shared_arr['args'],
            shared_arr['random_state']
        )
    )


class Classifier:
    def __init__(self, args, design, data, list_ids):
        self.args = args
        self.design = design
        self.data = data
        self.list_ids = list_ids

        self.num_query = len(np.where(design[self.args.SUBGROUP] == 1)[0])
        self.num_ref = len(np.where(design[self.args.SUBGROUP] == 0)[0])

        self.with_na = np.isnan(data).sum()
        self.done = False

        if args.RANDOM_SEED is not None:
            self.random_state = np.random.RandomState(args.RANDOM_SEED)
        else:
            self.random_state = np.random.RandomState()

        self.result = []

    def run(self):
        self.__pred()

    def pred2csv(self):

        sys.stderr.write(time.strftime('%X') + ": Save epcy results\n")
        if self.args.PATH_OUT is not None:
            if not os.path.exists(self.args.PATH_OUT):
                os.makedirs(self.args.PATH_OUT)

            file_out = self.args.PATH_OUT + "/predictive_capability.xls"
            file_pred_out = self.args.PATH_OUT + "/subgroup_predicted.xls"
            file_pred_normal_out = self.args.PATH_OUT + "/subgroup_predicted_normal.xls"

            with open(file_out, 'w') as w_csv:
                self.print_feature_header(w_csv, self.args, self.with_na > 0)
                self.print_feature_pred(
                    self.result, self.list_ids,
                    self.num_query, self.num_ref,
                    w_csv, self.args, self.with_na > 0
                )

            if self.args.FULL:
                with open(file_pred_out, 'w') as w_csv:
                    self.print_subgroup_header(w_csv, self.design)
                    self.print_subgroup_predicted(
                        self.result, self.list_ids, w_csv, "kernel")
                if self.args.NORMAL:
                    with open(file_pred_normal_out, 'w') as w_csv:
                        self.print_subgroup_header(w_csv,  self.design)
                        self.print_subgroup_predicted(
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
        header = "id\tl2fc\tkernel_mcc"
        header = header + "\tkernel_mcc_low"
        header = header + "\tkernel_mcc_high"
        header = header + "\tmean_query\tmean_ref"
        header = header + "\tbw_query\tbw_ref"
        if args.NORMAL:
            header = header + "\tnormal_mcc"
            if args.N_BAGGING > 1:
                header = header + "\tnormal_mcc_low"
                header = header + "\tnormal_mcc_high"
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
        for res in results:
            line = str(list_ids[cpt_id]) + "\t"
            if "l2fc" in res:
                line = line + str(res['l2fc'][0]) + "\t"
            else:
                line = line + "nan\t"

            if 'kernel_mcc' in res:
                k_mcc = res['kernel_mcc'][0]
                line = line + str(k_mcc[1]) + "\t"
                line = line + str(k_mcc[0]) + "\t"
                line = line + str(k_mcc[2]) + "\t"
            else:
                line = line + "nan\tnan\tnan\t"

            if 'mean_query' in res:
                line = line + str(res['mean_query'][0]) + "\t"
                line = line + str(res['mean_ref'][0]) + "\t"
                line = line + str(res['bw_query'][0]) + "\t"
                line = line + str(res['bw_ref'][0])
            else:
                line = line + "nan\tnan\tnan\tnan"
            if args.NORMAL:
                if 'normal_mcc' in res:
                    n_mcc = res['normal_mcc'][0]
                    line = line + "\t" + str(n_mcc[1])
                    if args.N_BAGGING > 1:
                        line = line + "\t" + str(n_mcc[0])
                        line = line + "\t" + str(n_mcc[2])
                else:
                    line = line + "nan\t"
                    if args.N_BAGGING > 1:
                        line = line + "nan\tnan\t"

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
    def print_subgroup_header(w_csv, design):
        line = "id" + '\t'
        line = line + '\t'.join(design['sample'])
        line = line + "\n"
        w_csv.write(line)

    @staticmethod
    def print_subgroup_predicted(results, list_ids, w_csv, pred):
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
                        num_bs, self.args, self.random_state
                    )
                )
        else:
            dtype = np.float64
            cdtype = np.ctypeslib.as_ctypes_type(dtype)
            data_shape = self.data.shape
            raw_array = RawArray(cdtype, data_shape[0] * data_shape[1])
            raw_array_np = np.frombuffer(raw_array, dtype=dtype).reshape(data_shape)
            np.copyto(raw_array_np, self.data)

            with Pool(
                processes=self.args.THREAD, initializer=init_worker,
                initargs=(
                    raw_array, data_shape, dtype,
                    self.num_query, self.num_ref, num_bs,
                    self.args, self.random_state
                )
            ) as p:
                self.result = p.map(worker_func, range(data_shape[0]))
