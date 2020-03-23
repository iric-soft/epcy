import sys
import math

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu, ttest_ind
from statistics import median

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


class Classifier:
    def __init__(self, args, design, data, list_ids, start):
        self.args = args
        self.start = start
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

    def run(self):
        self.__pred()

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

    def print_feature_pred(self, w_csv, with_na=False):
        if self.done:
            cpt_id = 0
            for select_id in self.list_ids:
                line = str(self.list_ids[cpt_id]) + "\t"
                line = line + str(self.l2fc[cpt_id]) + "\t"
                line = line + str(self.kernel_mcc[cpt_id][1]) + "\t"
                line = line + str(self.kernel_mcc[cpt_id][0]) + "\t"
                line = line + str(self.kernel_mcc[cpt_id][2]) + "\t"
                line = line + str(self.mean_query[cpt_id]) + "\t"
                line = line + str(self.mean_ref[cpt_id]) + "\t"
                line = line + str(self.bw_query[cpt_id]) + "\t"
                line = line + str(self.bw_ref[cpt_id])
                if self.args.NORMAL:
                    line = line + "\t" + str(self.normal_mcc[cpt_id][1])
                    if self.args.N_BAGGING > 1:
                        line = line + "\t" + str(self.normal_mcc[cpt_id][0])
                        line = line + "\t" + str(self.normal_mcc[cpt_id][2])

                if self.args.AUC:
                    auc = self.auc[cpt_id] if self.auc[cpt_id] >= 0.5 else 1 - self.auc[cpt_id]
                    line = line + "\t" + str(auc)
                    if self.args.UTEST:
                        line = line + "\t" + str(self.utest_pv[cpt_id])
                if self.args.TTEST:
                    line = line + "\t" + str(self.ttest_pv[cpt_id])
                if with_na:
                    line = line + "\t" + str(self.sample_query[cpt_id])
                    line = line + "\t" + str(self.sample_ref[cpt_id])

                line = line + "\n"

                w_csv.write(line)
                cpt_id += 1

    def print_subgroup_header(self, w_csv):
        line = "id" + '\t'
        line = line + '\t'.join(self.design['sample'])
        line = line + "\n"
        w_csv.write(line)

    def print_subgroup_predicted(self, w_csv, pred):
        if pred == "kernel":
            mat_pred = self.kernel_pred
        if pred == "normal":
            mat_pred = self.normal_pred

        if self.done:
            cpt_id = 0
            for select_id in self.list_ids:
                line = str(self.list_ids[cpt_id]) + "\t"
                line = line + '\t'.join([str(x) for x in mat_pred[cpt_id, :]])
                line = line + "\n"

                w_csv.write(line)
                cpt_id += 1

    def __create_empty_res(self):
        self.l2fc = np.empty(shape=(len(self.list_ids)), dtype=np.float32)
        self.l2fc.fill(np.nan)
        self.mean_query = np.empty(shape=(len(self.list_ids)),
                                   dtype=np.float32)
        self.mean_query.fill(np.nan)
        self.mean_ref = np.empty(shape=(len(self.list_ids)), dtype=np.float32)
        self.mean_ref.fill(np.nan)
        self.bw_query = np.empty(shape=(len(self.list_ids)), dtype=np.float32)
        self.bw_query.fill(np.nan)
        self.bw_ref = np.empty(shape=(len(self.list_ids)), dtype=np.float32)
        self.bw_ref.fill(np.nan)

        if self.with_na > 0:
            self.sample_query = np.empty(shape=(len(self.list_ids)),
                                         dtype=np.int32)
            self.sample_ref = np.empty(shape=(len(self.list_ids)),
                                       dtype=np.int32)
        if self.args.AUC:
            self.auc = np.empty(shape=(len(self.list_ids)), dtype=np.float32)
            self.auc.fill(np.nan)
            if self.args.UTEST:
                self.utest_pv = np.empty(shape=(len(self.list_ids)),
                                         dtype=np.float32)
                self.utest_pv.fill(np.nan)
        if self.args.TTEST:
            self.ttest_pv = np.empty(shape=(len(self.list_ids)),
                                     dtype=np.float32)
            self.ttest_pv.fill(np.nan)
        if self.args.NORMAL:
            self.normal_mcc = np.empty(shape=(len(self.list_ids), 3),
                                       dtype=np.float32)
            self.normal_mcc.fill(np.nan)
        self.kernel_mcc = np.empty(shape=(len(self.list_ids), 3),
                                   dtype=np.float32)
        self.kernel_mcc.fill(np.nan)

        if self.args.FULL:
            if self.args.NORMAL:
                self.normal_pred = np.empty(shape=(len(self.list_ids),
                                                   len(self.design["sample"])),
                                            dtype=np.float16)
                self.normal_pred.fill(np.nan)
            self.kernel_pred = np.empty(shape=(len(self.list_ids),
                                               len(self.design["sample"])),
                                        dtype=np.float16)
            self.kernel_pred.fill(np.nan)

    def __pred(self):
        num_bs = 0
        if hasattr(self.args, 'BS') and self.args.BS is not None:
            num_bs = self.args.BS

        self.__create_empty_res()

        cpt_id = 0
        for select_id in self.list_ids:
            row_data = self.data[cpt_id, :]

            num_query = self.num_query
            if self.with_na > 0:
                row_data, num_query = self.rm_missing(row_data, num_query)

                self.sample_query[cpt_id] = num_query
                self.sample_ref[cpt_id] = row_data.size - num_query

                if self.sample_query[cpt_id] <= 2 or self.sample_ref[cpt_id] <= 2:
                    cpt_id += 1
                    continue

            sum_row = sum(row_data)
            res = self.get_foldchange(row_data, num_query)
            self.l2fc[cpt_id] = res[0]
            self.mean_query[cpt_id] = res[1]
            self.mean_ref[cpt_id] = res[2]
            self.bw_query[cpt_id] = Classifier.bw_nrd(row_data[:num_query],
                                                      num_bs=num_bs)
            self.bw_ref[cpt_id] = Classifier.bw_nrd(row_data[num_query:],
                                                    num_bs=num_bs)
            # if select_id == "ENSG00000261541.1":
            #    print(select_id)
            #    print(row_data)
            #    print(self.l2fc[cpt_id])
            #    print(self.mean_query[cpt_id])
            #    print(self.mean_ref[cpt_id])

            if np.unique(row_data).size != 1:
                if sum_row >= self.args.EXP and abs(self.l2fc[cpt_id]) >= self.args.L2FC:
                    if self.args.AUC:
                        res = self.auc_u_test(row_data, num_query,
                                              self.num_ref)
                        self.auc[cpt_id] = res[0]

                        if self.args.UTEST:
                            self.utest_pv[cpt_id] = res[1]

                    if self.args.TTEST:
                        self.ttest_pv[cpt_id] = self.t_test_welch(row_data,
                                                                  num_query)

                    ct_by_bagging_and_pred_by_sample = self.pred_fill_cont_table_kernel(
                        row_data, num_query, self.args.MIN_BW,
                        n_bagging=self.args.N_BAGGING, num_bs=num_bs,
                        num_draw=self.args.N_DRAW,
                        random_state=self.random_state)
                    mcc, pred_by_sample = self.get_mcc_ped_by_bagging(ct_by_bagging_and_pred_by_sample)

                    self.kernel_mcc[cpt_id] = mcc

                    if self.args.FULL:
                        self.kernel_pred[cpt_id, :] = pred_by_sample

                    if self.args.NORMAL:
                        ct_by_bagging_and_pred_by_sample = self.pred_fill_cont_table_normal(
                            row_data, num_query, n_bagging=self.args.N_BAGGING,
                            num_draw=self.args.N_DRAW,
                            random_state=self.random_state,
                            num_bs=num_bs)
                        mcc, pred_by_sample = self.get_mcc_ped_by_bagging(ct_by_bagging_and_pred_by_sample)

                        self.normal_mcc[cpt_id] = mcc
                        if self.args.FULL:
                            self.normal_pred[cpt_id, :] = pred_by_sample

            cpt_id += 1

        del self.data
        if not self.args.FULL:
            del self.design

        self.done = True

    @staticmethod
    def rm_missing(row_data, num_query):
        ids_na = np.isnan(row_data)
        if sum(ids_na) > 0:
            row_data = row_data[~np.isnan(row_data)]
            num_query = num_query - sum(ids_na[:num_query])

        return(row_data, num_query)

    @staticmethod
    def get_foldchange(row_data, num_query):
        mean_query = np.mean(row_data[:num_query])
        mean_ref = np.mean(row_data[num_query:])
        log2fc = mean_query - mean_ref

        return(log2fc, mean_query, mean_ref)

    @staticmethod
    def auc_u_test(row_data, num_query, num_ref):
        # print(row_data)
        (u_value, p_value) = mannwhitneyu(row_data[:num_query],
                                          row_data[num_query:],
                                          alternative="two-sided")
        auc = u_value / (num_query * num_ref)

        return(auc, p_value)

    @staticmethod
    def t_test_welch(row_data, num_query):
        (t_value, p_value) = ttest_ind(row_data[:num_query],
                                       row_data[num_query:],
                                       equal_var=False)

        return(p_value)

    @staticmethod
    def pred_fill_cont_table_kernel(row_data, num_query, min_bw, bw=None,
                                    n_bagging=1, num_bs=0, num_draw=100,
                                    random_state=np.random.RandomState()):
        # Compute sample assignation using kernel
        N = len(row_data)
        return(Classifier.get_ct_using_fx_kernel(
            row_data, num_query, N, min_bw, bw=bw, n_bagging=n_bagging,
            num_bs=num_bs, num_draw=num_draw, random_state=random_state)
        )

    @staticmethod
    def pred_fill_cont_table_normal(row_data, num_query, n_bagging=1,
                                    num_draw=100, num_bs=0,
                                    random_state=np.random.RandomState()):
        # Compute sample assignation using normal dist
        N = len(row_data)
        return(Classifier.get_ct_using_fx_normal(
            row_data, num_query, N, n_bagging=n_bagging, num_draw=num_draw,
            random_state=random_state, num_bs=num_bs)
        )

    @staticmethod
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

    @staticmethod
    def compute_kernel_fx_and_ct(k_bw_nq_gen, num_query, N, num_bs,
                                 num_draw=100,
                                 random_state=np.random.RandomState()):
        fx_by_sample = [Classifier.compute_kernel_fx(k_bw_nq, i, num_bs=num_bs)
                        for i, k_bw_nq in enumerate(k_bw_nq_gen)]

        return(Classifier.fx_to_tables(fx_by_sample, num_query, N,
                                       num_draw=num_draw,
                                       random_state=random_state))

    @staticmethod
    def get_ct_using_fx_kernel(row_data, num_query, N, min_bw, bw=None,
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
        k_bw_gen_by_fold = [Classifier.get_k_gaussian_kernel_for_all_x(
                                row_data, x_ids, num_query, min_bw,
                                bw=bw, n_bagging=n_bagging, num_bs=num_bs,
                                random_state=random_state)
                            for x_ids in n_folds]

        k_bw_gen_by_bag = np.transpose(
                            np.asarray(k_bw_gen_by_fold),
                            (2, 0, 1, 3))
        k_bw_gen_by_bag = np.reshape(k_bw_gen_by_bag, (n_bagging, N, 2))

        # print("BAG")
        # print(k_bw_gen_by_bag)

        # (1/(n*bw)) * sum k((x-xi) / bw)
        ct_by_bag = np.asarray([Classifier.compute_kernel_fx_and_ct(
                                                k_bw_nq_gen, num_query, N,
                                                num_bs, num_draw=num_draw,
                                                random_state=random_state)
                                for k_bw_nq_gen in k_bw_gen_by_bag])

        return(ct_by_bag)

    @staticmethod
    def get_ct_using_fx_normal(row_data, num_query, N, n_bagging=1, num_bs=0,
                               num_draw=100,
                               random_state=np.random.RandomState()):
        # Create leave one out index
        n_folds = None
        if num_bs > 0:
            n_folds = np.array_split(np.arange(N), N/num_bs)
        else:
            n_folds = np.array_split(np.arange(N), N)

        fx_by_fold = [Classifier.compute_normal_fx_for_all_x(
                                    row_data, x_ids, num_query,
                                    n_bagging=n_bagging, num_bs=num_bs,
                                    random_state=random_state)
                      for x_ids in n_folds]
        fx_by_fold = np.asarray(fx_by_fold)
        # print(fx_by_fold.shape)
        # print(fx_by_fold)
        fx_by_bag = np.transpose(fx_by_fold, (2, 0, 1))
        fx_by_bag = np.reshape(fx_by_bag, (n_bagging, N, N - num_bs))

        ct_by_bag = (Classifier.fx_to_tables(
                                    fx_by_sample, num_query, N,
                                    num_draw=num_draw,
                                    random_state=random_state)
                     for fx_by_sample in fx_by_bag)

        return(ct_by_bag)

    @staticmethod
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

    @staticmethod
    def bw_var(x):
        return(np.var(x))

    @staticmethod
    def bw_nrd(x, num_bs=0):
        # TODO need to improve speed of this part
        if num_bs != 0:
            x = x[range(1, len(x), 4)]

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

    @staticmethod
    def bw_nrd0(x, num_bs=0):
        # TODO need to improve speed of this part
        if num_bs != 0:
            x = x[range(1, len(x), 4)]

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

    @staticmethod
    def get_bagging_other(other, num_query,
                          random_state=np.random.RandomState()):
        while True:
            ids = np.sort(random_state.choice(len(other), len(other)))
            bag_num_query = np.where(ids < num_query)[0].shape[0]
            if bag_num_query >= 2 and bag_num_query <= len(other) - 2:
                break

        return(other[ids], bag_num_query)

    @staticmethod
    def get_k_gaussian_kernel_for_all_x(row_data, x_ids, num_query, min_bw,
                                        bw=None, n_bagging=1, num_bs=0,
                                        random_state=np.random.RandomState()):
        x_values = row_data[x_ids]
        other = np.delete(row_data, x_ids)
        o_num_query = num_query - np.sum(x_ids < num_query)

        return([Classifier.get_k_gaussian_kernel(x, other, o_num_query, min_bw,
                                                 bw, n_bagging, random_state,
                                                 num_bs=num_bs)
                for x in x_values])

    @staticmethod
    def get_k_gaussian_kernel(x, other, num_query, min_bw, bw, n_bagging,
                              random_state, num_bs=0):
        if n_bagging > 1:
            bag_others = [Classifier.get_bagging_other(other, num_query,
                                                       random_state=random_state)
                          for j in range(n_bagging)]
            return([Classifier.k_gaussian_kernel(x, other, min_bw,
                                                 b_num_query, bw,
                                                 num_bs=num_bs)
                    for other, b_num_query in bag_others])
        else:
            return([Classifier.k_gaussian_kernel(x, other, min_bw,
                                                 num_query, bw, num_bs=num_bs)
                    for j in range(1)])

    @staticmethod
    def k_gaussian_kernel(x, other, min_bw, ids_split, bw, num_bs=0):
        other_query = other[:ids_split]
        other_ref = other[ids_split:]

        if bw is None:
            # bw_query = Classifier.bw_nrd0(other_query)
            # bw_ref = Classifier.bw_nrd0(other_ref)
            bw_query = Classifier.bw_nrd(other_query, num_bs)
            bw_ref = Classifier.bw_nrd(other_ref, num_bs)

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

    @staticmethod
    def compute_normal_fx_for_all_x(row_data, x_ids, num_query, epsilon=0.001,
                                    n_bagging=1, num_bs=0,
                                    random_state=np.random.RandomState()):
        x_values = row_data[x_ids]
        other = np.delete(row_data, x_ids)
        o_num_query = num_query - np.sum(x_ids < num_query)

        return([Classifier.compute_normal_fx(x, other, o_num_query,
                                             epsilon, n_bagging, num_bs=num_bs,
                                             random_state=random_state)
                for x in x_values])

    @staticmethod
    def compute_normal_fx(x, other, num_query, epsilon, n_bagging,
                          random_state, num_bs=0):
        if n_bagging > 1:
            bag_others = (Classifier.get_bagging_other(other, num_query,
                                                       random_state=random_state)
                          for j in range(n_bagging))
            return([Classifier.fx_normal(x, other, b_num_query,
                                         epsilon=epsilon, num_bs=num_bs)
                    for other, b_num_query in bag_others])
        else:
            return([Classifier.fx_normal(x, other, num_query, epsilon=epsilon,
                                         num_bs=num_bs)
                    for j in range(1)])

    @staticmethod
    def fx_normal(x, other, id_split, epsilon, num_bs=0):
        if num_bs != 0:
            x = x[range(1, len(x), 4)]
            other = other[range(1, len(x), 4)]

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

    @staticmethod
    def get_mcc_ped_by_bagging(ct_by_bagging):
        all_mcc = []
        all_pred = []
        for cts, pred in ct_by_bagging:
            for ct in cts:
                all_mcc.append(Classifier.get_mcc(ct))
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

    @staticmethod
    def get_mcc(ct):
        d1 = ct[0] + ct[1]
        d2 = ct[0] + ct[2]
        d3 = ct[3] + ct[1]
        d4 = ct[3] + ct[2]

        n1 = ct[0] * ct[3]
        n2 = ct[1] * ct[2]

        return((n1 - n2) / math.sqrt(d1 * d2 * d3 * d4)
               if d1 != 0 and d2 != 0 and d3 != 0 and d4 != 0 else 0)
