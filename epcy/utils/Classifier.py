import sys
import math

import numpy as np
import pandas as pd

import numexpr as ne
ne.set_num_threads(1)

from scipy.stats import mannwhitneyu, ttest_ind
from statistics import median
#from random import random

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

        self.done = False

    def run(self):
        self.__pred()

    @staticmethod
    def print_feature_header(w_csv, args):
        header = "id\tl2fc\tkernel_mcc"
        if args.N_BAGGING > 1:
            header = header + "\tkernel_mcc_low"
            header = header + "\tkernel_mcc_high"
        header = header + "\tmean_query\tmean_ref"
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
        header = header + "\n"

        w_csv.write(header)

    def print_feature_pred(self, w_csv):
        if self.done:
            cpt_id = 0
            for select_id in self.list_ids:
                line = str(self.list_ids[cpt_id]) + "\t"
                line = line + str(self.l2fc[cpt_id]) + "\t"
                line = line + str(self.kernel_mcc[cpt_id][1]) + "\t"
                if self.args.N_BAGGING > 1:
                    line = line + str(self.kernel_mcc[cpt_id][0]) + "\t"
                    line = line + str(self.kernel_mcc[cpt_id][2]) + "\t"
                line = line + str(self.mean_query[cpt_id]) + "\t"
                line = line + str(self.mean_ref[cpt_id])
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
        self.mean_query = np.empty(shape=(len(self.list_ids)), dtype=np.float32)
        self.mean_query.fill(np.nan)
        self.mean_ref = np.empty(shape=(len(self.list_ids)), dtype=np.float32)
        self.mean_ref.fill(np.nan)


        if self.args.AUC:
            self.auc = np.empty(shape=(len(self.list_ids)), dtype=np.float32)
            self.auc.fill(np.nan)
            if self.args.UTEST:
                self.utest_pv = np.empty(shape=(len(self.list_ids)), dtype=np.float32)
                self.utest_pv.fill(np.nan)
        if self.args.TTEST:
            self.ttest_pv = np.empty(shape=(len(self.list_ids)), dtype=np.float32)
            self.ttest_pv.fill(np.nan)
        if self.args.NORMAL:
            self.normal_mcc = np.empty(shape=(len(self.list_ids), 3), dtype=np.float32)
            self.normal_mcc.fill(np.nan)
        self.kernel_mcc = np.empty(shape=(len(self.list_ids), 3), dtype=np.float32)
        self.kernel_mcc.fill(np.nan)

        if self.args.FULL:
            if self.args.NORMAL:
                self.normal_pred = np.empty(shape=(len(self.list_ids), len(self.design["sample"])), dtype=np.int8)
                self.normal_pred.fill(np.nan)
            self.kernel_pred = np.empty(shape=(len(self.list_ids), len(self.design["sample"])), dtype=np.int8)
            self.kernel_pred.fill(np.nan)

    def __pred(self):
        self.__create_empty_res()

        cpt_id = 0
        for select_id in self.list_ids:
            row_data = self.data[cpt_id,:]

            sum_row = sum(row_data)
            res = self.get_foldchange(row_data, self.num_query)
            self.l2fc[cpt_id] = res[0]
            self.mean_query[cpt_id] = res[1]
            self.mean_ref[cpt_id] = res[2]
            #if select_id == "ENSG00000261541.1":
            #    print(select_id)
            #    print(row_data)
            #    print(self.l2fc[cpt_id])
            #    print(self.mean_query[cpt_id])
            #    print(self.mean_ref[cpt_id])

            if np.unique(row_data).size != 1:
                if sum_row >= self.args.EXP and abs(self.l2fc[cpt_id]) >= self.args.L2FC:
                    if self.args.AUC:
                        res = self.auc_u_test(row_data, self.num_query, self.num_ref)
                        self.auc[cpt_id] = res[0]

                        if self.args.UTEST:
                            self.utest_pv[cpt_id] = res[1]

                    if self.args.TTEST:
                        self.ttest_pv[cpt_id] = self.t_test_welch(row_data, self.num_query)

                    ct_by_bagging_and_pred_by_sample = self.pred_fill_cont_table_kernel(row_data, self.num_query, self.args.MIN_BW, n_bagging=self.args.N_BAGGING)
                    mcc, pred_by_sample = self.get_mcc_ped_by_bagging(ct_by_bagging_and_pred_by_sample)

                    self.kernel_mcc[cpt_id] = mcc

                    if self.args.FULL:
                        self.kernel_pred[cpt_id, :] = pred_by_sample
                    #if self.kernel_mcc[cpt_id] < 0:
                    #    print(select_id)
                    #    print(self.l2fc[cpt_id])
                    #    # print(row_data)
                    #    print(ct)
                    #    print(self.kernel_mcc[cpt_id])
                    #########
                    if self.args.NORMAL:
                        ct_by_bagging_and_pred_by_sample = self.pred_fill_cont_table_normal(row_data, self.num_query, n_bagging=self.args.N_BAGGING)
                        mcc, pred_by_sample = self.get_mcc_ped_by_bagging(ct_by_bagging_and_pred_by_sample)

                        self.normal_mcc[cpt_id] = mcc
                        if self.args.FULL:
                            self.normal_pred[cpt_id, :] = pred_by_sample

            #else:
            #    print(select_id)
            #    print(row_data)
            #    print(sum_row)
            #    print(abs(self.l2fc[cpt_id]))

            cpt_id += 1

        del self.data
        if not self.args.FULL:
            del self.design

        self.done = True

    @staticmethod
    def get_foldchange(row_data, num_query):
        mean_query = np.mean(row_data[:num_query])
        mean_ref = np.mean(row_data[num_query:])
        log2fc = mean_query - mean_ref

        return(log2fc, mean_query, mean_ref)

    @staticmethod
    def auc_u_test(row_data, num_query, num_ref):
        #print(row_data)
        (u_value, p_value) = mannwhitneyu(row_data[:num_query], row_data[num_query:], alternative="two-sided")
        auc = u_value / (num_query * num_ref)

        return(auc, p_value)

    @staticmethod
    def t_test_welch(row_data, num_query):
        (t_value, p_value) = ttest_ind(row_data[:num_query], row_data[num_query:], equal_var=False)

        return(p_value)

    @staticmethod
    def pred_fill_cont_table_kernel(row_data, num_query, min_bw, bw=None, n_bagging=1):
        # Compute sample assignation using kernel
        N = len(row_data)
        return(Classifier.get_ct_using_fx_kernel(row_data, num_query, N, min_bw, bw=bw, n_bagging=n_bagging))

    @staticmethod
    def pred_fill_cont_table_normal(row_data, num_query, n_bagging=1):
        # Compute sample assignation using normal dist
        N = len(row_data)
        return(Classifier.get_ct_using_fx_normal(row_data, num_query, N, n_bagging=n_bagging))

    @staticmethod
    def compute_kernel_fx(k_query, k_ref):
        ''' return (1/(n*bw)) * sum k((x-xi) / bw) '''
        sum_llh = k_query.sum() + k_ref.sum()
        if sum_llh == 0:
            res = 0.5
        else:
            res = k_query.sum() / sum_llh

        return(res)

    @staticmethod
    def compute_kernel_fx_and_ct(k_bw_nq_gen, num_query, N):
        #print("compute_kernel_fx_and_ct")
        fx_by_sample = [Classifier.compute_kernel_fx(k_bw_nq[0], k_bw_nq[1]) for k_bw_nq in k_bw_nq_gen]
        #print("#1")
        #print(fx_by_sample)
        #print("compute_kernel_fx end")
        return(Classifier.fx_to_tables(fx_by_sample, num_query, N))

    @staticmethod
    def get_ct_using_fx_kernel(row_data, num_query, N, min_bw, bw=None, n_bagging=1):
        # Create leave one out index
        idx_gen = (np.arange(1, N) - ([1]*i + [0]*(N-i-1)) for i in range(N))

        # Compute k((x-xi) / bw) for each leave one out
        k_bw_gen_by_fold = [Classifier.get_k_gausian_kernel(row_data, i, idx, num_query, min_bw, bw=bw, n_bagging=n_bagging)
                            for i, idx in enumerate(idx_gen)]
        #print(k_bw_gen_by_fold)
        k_bw_gen_by_bag = np.transpose(np.asarray(k_bw_gen_by_fold),(1, 0, 2))
        #print("BAG")
        #print(k_bw_gen_by_bag)

        # (1/(n*bw)) * sum k((x-xi) / bw)
        ct_by_bag = [Classifier.compute_kernel_fx_and_ct(k_bw_nq_gen, num_query, N)
                        for k_bw_nq_gen in k_bw_gen_by_bag]

        return(ct_by_bag)

    @staticmethod
    def get_ct_using_fx_normal(row_data, num_query, N, n_bagging=1):
        # Create leave one out index
        idx_gen = (np.arange(1, N) - ([1]*i + [0]*(N-i-1)) for i in range(N))

        fx_by_fold = [Classifier.compute_normal_fx(row_data, i, idx, num_query, n_bagging=n_bagging)
                            for i, idx in enumerate(idx_gen)]
        fx_by_bag = np.transpose(np.asarray(fx_by_fold),(1, 0))

        #print("####################################### normal")
        ct_by_bag = (Classifier.fx_to_tables(fx_by_sample, num_query, N)
                        for fx_by_sample in fx_by_bag)

        return(ct_by_bag)

    @staticmethod
    def fx_to_tables(fx_by_sample, num_query, N):
        # print(fx_by_sample)
        # return:
        #  cont_table[tp, fp, fn, tn]
        #  pred_by_sample: 1=tp, 2=fn, 3=fp, tn=4

        cont_tables = []
        for i in range(5):
            pred_by_sample = [np.random.random() < fx for fx in fx_by_sample]
            #pred_by_sample = [True if fx > 0 else False for fx in fx_by_sample]
            tp_fn = np.add.reduceat(pred_by_sample,[0,num_query])

            cont_table = [tp_fn[0], num_query-tp_fn[0], tp_fn[1], N - num_query - tp_fn[1]]
            cont_tables.append(cont_table)

        pred_by_sample = np.fromiter((1 if pred else 2 for pred in pred_by_sample), np.int8, len(pred_by_sample))
        pred_by_sample[num_query:] += 2

        # switch 2 and 3 for dendrogram distance
        id2 = np.where(pred_by_sample == 2)
        id3 = np.where(pred_by_sample == 3)
        pred_by_sample[id2] = 3
        pred_by_sample[id3] = 2

        return(cont_tables, pred_by_sample)

    @staticmethod
    def bw_nrd0(x):
        #TODO need to improve speed of this part
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
    def get_bagging_other(other, num_query):
        while True:
            ids = np.sort(np.random.choice(len(other), len(other)))
            bag_num_query = np.where(ids<num_query)[0].shape[0]
            if bag_num_query >= 2 and bag_num_query <= len(other) - 2:
                break

        return(other[ids], bag_num_query)


    @staticmethod
    def get_k_gausian_kernel(row_data, i, idx, num_query, min_bw, bw=None, n_bagging=1):
        x = row_data[i]
        other = row_data[idx]
        o_num_query = num_query - 1 if i < num_query else num_query
        if n_bagging > 1:
            bag_others = [Classifier.get_bagging_other(other, o_num_query) for j in range(n_bagging)]
            return( [Classifier.k_gausian_kernel(x, other, i, min_bw, b_num_query, bw=bw) for other, b_num_query in bag_others] )
        else:
            return( [Classifier.k_gausian_kernel(x, other, i, min_bw, o_num_query, bw=bw) for j in range(1)] )

    def k_gausian_kernel(x, other, i, min_bw, ids_split, bw=None):
        other_query = other[:ids_split]
        other_ref = other[ids_split:]
        
        if bw is None:
            bw = Classifier.bw_nrd0(other)
            if bw < min_bw:
                bw = min_bw

        norm_query = other_query.shape[0] * bw
        norm_ref = other_ref.shape[0] * bw

        return([
            ne.evaluate('(0.3989423 * exp(-1/2*(((x - other_query) / bw)**2)))/norm_query'),
            ne.evaluate('(0.3989423 * exp(-1/2*(((x - other_ref) / bw)**2)))/norm_ref'),
        ])

    @staticmethod
    def compute_normal_fx(row_data, i, idx, num_query, epsilon=0.001, n_bagging=1):
        o_num_query = num_query - 1 if i < num_query else num_query
        if n_bagging > 1:
            other = row_data[idx]
            bag_others = (Classifier.get_bagging_other(other, o_num_query) for j in range(n_bagging))
            return( [Classifier.fx_normal(row_data[i], other, i, b_num_query, epsilon=epsilon) for other, b_num_query in bag_others] )
        else:
            return( [Classifier.fx_normal(row_data[i], row_data[idx], i, o_num_query, epsilon=epsilon) for j in range(1)] )

    def fx_normal(x, other, i, id_split, epsilon=0.001):
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
            mccs = []
            for ct in cts:
                mccs.append(Classifier.get_mcc(ct))
            all_mcc.append(median(mccs))
            all_pred.append(pred)
        pred_by_sample = np.median(np.asarray(all_pred),axis=0)

        return(np.quantile(all_mcc, [0.05, 0.5, 0.95]), pred_by_sample)

    @staticmethod
    def get_mcc(ct):
        d1 = ct[0] + ct[1]
        d2 = ct[0] + ct[2]
        d3 = ct[3] + ct[1]
        d4 = ct[3] + ct[2]

        n1 = ct[0] * ct[3]
        n2 = ct[1] * ct[2]

        return((n1 - n2) / math.sqrt(d1 * d2 * d3 * d4) if d1!=0 and d2!=0 and d3!=0 and d4!=0 else 0)
