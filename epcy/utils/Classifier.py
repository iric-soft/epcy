import math

import numpy as np
import pandas as pd

import numexpr as ne
ne.set_num_threads(1)

from scipy.stats import mannwhitneyu, ttest_ind

class Classifier:
    def __init__(self, args, design, start):
        self.args = args
        self.start = start
        self.design = design

        self.num_query = len(np.where(design[self.args.SUBGROUP] == 1)[0])
        self.num_ref = len(np.where(design[self.args.SUBGROUP] == 0)[0])

        self.done = False

    def run(self):
        self.__load_data()
        self.__pred()

    @staticmethod
    def print_feature_header(w_csv):
        w_csv.write("ID\tL2FC\tKERNEL_MCC\tNORMAL_MCC\tAUC\tMEAN_QUERY\tMEAN_REF\tU_PV\tT_PV\n")

    def print_feature_pred(self, w_csv):
        if self.done:
            cpt_id = 0
            for select_id in self.list_ids:
                line = str(self.list_ids[cpt_id]) + "\t"
                line = line + str(self.l2fc[cpt_id]) + "\t"
                line = line + str(abs(self.kernel_mcc[cpt_id])) + "\t"
                line = line + str(abs(self.normal_mcc[cpt_id])) + "\t"
                auc = self.auc[cpt_id] if self.auc[cpt_id] >= 0.5 else 1 - self.auc[cpt_id]
                line = line + str(auc) + "\t"
                line = line + str(self.mean_query[cpt_id]) + "\t"
                line = line + str(self.mean_ref[cpt_id]) + "\t"
                line = line + str(self.utest_pv[cpt_id]) + "\t"
                line = line + str(self.ttest_pv[cpt_id]) + "\n"

                w_csv.write(line)
                cpt_id += 1

    def print_subgroup_header(self, w_csv):
        line = "ID" + '\t'
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


    def __load_data(self):
        self.data = pd.io.parsers.read_csv(self.args.MATRIX, sep="\t", index_col=0)
        self.data = self.data[self.start:(self.start+self.args.BY)]
        self.data = self.data.reindex(self.design["sample"], axis=1)
        self.list_ids = list(self.data.index)

        self.data = self.data.values
        if not self.args.LOG:
            self.data = np.log2(self.data + self.args.C)

    def __create_empty_res(self):
        self.l2fc = np.empty(shape=(len(self.list_ids)), dtype=np.float64)
        self.l2fc.fill(np.nan)
        self.mean_query = np.empty(shape=(len(self.list_ids)), dtype=np.float64)
        self.mean_query.fill(np.nan)
        self.mean_ref = np.empty(shape=(len(self.list_ids)), dtype=np.float64)
        self.mean_ref.fill(np.nan)

        self.auc = np.empty(shape=(len(self.list_ids)), dtype=np.float64)
        self.auc.fill(np.nan)
        self.utest_pv = np.empty(shape=(len(self.list_ids)), dtype=np.float64)
        self.utest_pv.fill(np.nan)
        self.ttest_pv = np.empty(shape=(len(self.list_ids)), dtype=np.float64)
        self.ttest_pv.fill(np.nan)
        self.normal_mcc = np.empty(shape=(len(self.list_ids)), dtype=np.float64)
        self.normal_mcc.fill(np.nan)
        self.kernel_mcc = np.empty(shape=(len(self.list_ids)), dtype=np.float64)
        self.kernel_mcc.fill(np.nan)

        self.normal_pred = np.empty(shape=(len(self.list_ids), len(self.design["sample"])), dtype=np.float64)
        self.normal_pred.fill(np.nan)
        self.kernel_pred = np.empty(shape=(len(self.list_ids), len(self.design["sample"])), dtype=np.float64)
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

            if sum_row >= self.args.EXP and abs(self.l2fc[cpt_id]) >= self.args.L2FC:
                res = self.auc_u_test(row_data, self.num_query, self.num_ref)
                self.auc[cpt_id] = res[0]
                self.utest_pv[cpt_id] = res[1]

                self.ttest_pv[cpt_id] = self.t_test_welch(row_data, self.num_query)
                #########
                #TODO (linked with bw_nrd0 function)
                #This version compute one bw (using all samples) which will be use for all leave-one-out.
                #It's fastest but not clean.
                #(ct_kernel, pred_by_sample) = pred_fill_cont_table_kernel(scores_tpm, num_query, bw_nrd0(scores_tpm))
                (ct, pred_by_sample) = self.pred_fill_cont_table_kernel(row_data, self.num_query)
                self.kernel_mcc[cpt_id] = self.get_mcc(ct)
                self.kernel_pred[cpt_id, :] = pred_by_sample
                #########
                (ct, pred_by_sample) = self.pred_fill_cont_table_normal(row_data, self.num_query)
                self.normal_mcc[cpt_id] = self.get_mcc(ct)
                self.normal_pred[cpt_id, :] = pred_by_sample
            else:
                print(select_id)
                print(row_data)
                print(sum_row)
                print(abs(self.l2fc[cpt_id]))

            cpt_id += 1

        self.done = True


    @staticmethod
    def get_foldchange(row_data, num_query):
        mean_query = np.mean(row_data[:num_query])
        mean_ref = np.mean(row_data[num_query:])
        log2fc = mean_query - mean_ref

        return(log2fc, mean_query, mean_ref)

    @staticmethod
    def auc_u_test(row_data, num_query, num_ref):
        (u_value, p_value) = mannwhitneyu(row_data[:num_query], row_data[num_query:], alternative="two-sided")
        auc = u_value / (num_query * num_ref)

        return(auc, p_value)

    @staticmethod
    def t_test_welch(row_data, num_query):
        (t_value, p_value) = ttest_ind(row_data[:num_query], row_data[num_query:], equal_var=False)

        return(p_value)

    @staticmethod
    def pred_fill_cont_table_kernel(row_data, num_query, bw=None):
        # Compute sample assignation using kernel
        N = len(row_data)
        fx_by_sample = Classifier.get_fx_kernel(row_data, num_query, N, bw)

        return(Classifier.fx_to_tables(fx_by_sample, num_query, N))

    @staticmethod
    def pred_fill_cont_table_normal(row_data, num_query):
        # Compute sample assignation using normal dist
        N = len(row_data)
        fx_by_sample = Classifier.get_fx_normal(row_data, num_query, N)

        return(Classifier.fx_to_tables(fx_by_sample, num_query, N))

    @staticmethod
    def get_fx_kernel(row_data, num_query, N, bw = None):
        # Create leave one out index
        idx = np.arange(1, N) - np.tri(N, N-1, k=-1, dtype=bool)
        # Compute k((x-xi) / h) for each leave one out
        k_bw = [Classifier.k_gausian_kernel(row_data[i], row_data[idx[i]], bw)
                    for i in range(N)]

        # (1/(n*bw)) * sum k((x-xi) / h)
        sum_k = [np.add.reduceat(sample_out,[0,num_query])
                            for sample_out, bw in k_bw]
        fx_by_sample = [sum_k[i] / (np.array([num_query,N-num_query]) * k_bw[i][1])
                            for i in range(N)]

        return(fx_by_sample)

    @staticmethod
    def get_fx_normal(row_data, num_query, N):
        # Create leave one out index
        idx = np.arange(1, N) - np.tri(N, N-1, k=-1, dtype=bool)
        fx_by_sample = [Classifier.exp_normal(row_data[i], row_data[idx[i]], i, num_query)
                            for i in range(N)]

        return(fx_by_sample)

    @staticmethod
    def fx_to_tables(fx_by_sample, num_query, N):
        # return:
        #  cont_table[tp, fp, fn, tn]
        #  pred_by_sample: 1=tp, 2=fn, 3=fp, tn=4
        pred_by_sample = [True if fx[0] > fx[1] else False for fx in fx_by_sample]
        tp_fn = np.add.reduceat(pred_by_sample,[0,num_query])
        cont_table = [tp_fn[0], num_query-tp_fn[0], tp_fn[1], N-num_query-tp_fn[1]]

        pred_by_sample = np.fromiter((1 if pred else 2 for pred in pred_by_sample), np.int32, len(pred_by_sample))
        pred_by_sample[num_query:] += 2

        # switch 2 and 3 for dendrogram distance
        id2 = np.where(pred_by_sample == 2)
        id3 = np.where(pred_by_sample == 3)
        pred_by_sample[id2] = 3
        pred_by_sample[id3] = 2

        return(cont_table, pred_by_sample)

    @staticmethod
    def bw_nrd0(x):
        #TODO need to improve speed of this part
        hi = np.std(x)
        iqr = np.subtract(*np.percentile(x, [75, 25]))
        lo = min(hi, iqr/1.34)
        if (lo == 0):
            lo = hi
            if (lo == 0):
                lo = math.abs(x[1])
                if (lo == 0):
                    lo = 1

        # this function can be run by ne.evaluate, with all lo pre-computed
        return(0.9 * lo * len(x)**(-0.2))

    @staticmethod
    def k_gausian_kernel(x, other, bw = None):
        if bw is None:
            bw = Classifier.bw_nrd0(other)
        pi = np.pi
        return(ne.evaluate('(1/sqrt(2 * pi)) * exp(-1/2*(((x - other) / bw)**2))'), bw)

    @staticmethod
    def gausian_kernel(u, pi):
        pi = np.pi
        return(ne.evaluate('(1/sqrt(2 * pi)) * exp(-1/2*(u**2))'))

    @staticmethod
    def exp_normal(x, other, i, num_query, epsilon=0.0000000000001):
        id_split = num_query
        if i < num_query:
            sep = num_query-1

        mu = np.mean(other[:id_split])
        var = np.var(other[:id_split]) + epsilon
        first_part = 1 / math.sqrt(2 * np.pi * var)
        fx_query = first_part * np.exp(-((x-mu)**2)/(2*var))

        mu = np.mean(other[id_split:])
        var = np.var(other[id_split:]) + epsilon
        first_part = 1 / math.sqrt(2 * np.pi * var)
        fx_ref = first_part * np.exp(-((x-mu)**2)/(2*var))

        #return(ne.evaluate('exp(-((x-mu)**2)/(2*var))'))
        return(np.array([fx_query,fx_ref]))

    @staticmethod
    def get_mcc(ct):
        d1 = ct[0] + ct[1]
        d2 = ct[0] + ct[2]
        d3 = ct[3] + ct[1]
        d4 = ct[3] + ct[2]

        n1 = ct[0] * ct[3]
        n2 = ct[1] * ct[2]

        return((n1 - n2) / math.sqrt(d1 * d2 * d3 * d4) if d1!=0 and d2!=0 and d3!=0 and d4!=0 else 0)



class Classifier_rsem(Classifier):
    def __load_data():
        self.data = np.zeros((self.design.shape[0], self.args.BY), dtype=np.float32)
        cpt = 0
        for index, row in self.design.iterrows():
            sample = row['sample']
            file_name = str(row['rsem'])


class Classifier_kallisto(Classifier):
    def __init__(self, args, design, start):
        Classifier.__init__(self, args, design, start)

    def __load_data():
        self.data = np.zeros((self.design.shape[0],
                              self.args.SAMPLE_BS,
                              self.args.BY
                             ),
                             dtype=np.float32
                            )

        cpt = 0
        for index, row in self.design.iterrows():
            sample = row['sample']
            file_name = os.path.join(str(row['kallisto']), "abundance.h5")
            tmp = add_kall_sample(file_name, kallisto, sample, num_kal_bs,
                                  transcripts, transcripts_len, args,
                                  ids_genes, uniq_genes, start)
            #!!! Don't use index (which is the num of sample line in project file) !!!
            self.data[cpt,:,:] = tmp[0]
            trans_used = tmp[1]
            genes_used = tmp[2]
            ids_genes_used = tmp[3]

            cpt += 1
