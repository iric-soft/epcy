import os
import time
import sys
import h5py

import numpy as np
import pandas as pd

from .. utils import readers as ur
from .. utils import other as uo

def create_kal_mat(args, design, df_anno):
    """ Create features matrix from kallisto output using kallisto
        column in the dsign file.
    """
    file_name = os.path.join(str(design["kallisto"][0]), "abundance.h5")
    sys.stderr.write(time.strftime('%X') + ":\tRead kallisto info\n")
    f = h5py.File(file_name, 'r')

    if args.BS == -1:
        args.BS = f["aux/num_bootstrap"][0]

    transcripts = f["aux/ids"][:]
    transcripts = [uo.cut_version(trans.decode('utf-8')) for trans in transcripts]
    transcripts = np.array(transcripts)
    uniq_genes = None

    ids_genes = None
    if args.GENE:
        sys.stderr.write(time.strftime('%X') + ":\tTranscripts to genes\n")
        parents = df_anno.reindex(transcripts)["Parent"].values
        #TODO find an other way to works with other ids not ENSG
        uniq_genes = df_anno.index[df_anno.index.str.contains("ENSG")]
        uniq_genes = np.array(uniq_genes)
        sys.stderr.write(time.strftime('%X') + ":\t\tkeep " + str(uniq_genes.size) + " ENSG ids\n")
        sys.stderr.write(time.strftime('%X') + ":\t\tfound genes ids for each trans\n")
        ids_genes = [np.where(parents == gene)[0] for gene in uniq_genes]
        ids_genes = np.array(ids_genes)

    transcripts_len = f["aux/eff_lengths"][:]
    f.close()

    rows_id = transcripts
    num_features = transcripts.size
    if args.GENE:
        rows_id = uniq_genes
        num_features = uniq_genes.size

    if args.BS == 0:
        kallisto = np.zeros((num_features, design.shape[0]), dtype=np.float32)
    else:
        kallisto = np.zeros((num_features, design.shape[0]*args.BS), dtype=np.float32)

    sys.stderr.write(time.strftime('%X') + ":\tRead samples kallisto quantification\n")
    cpt = 0
    for index, row in design.iterrows():
        sample = row['sample']
        #print(str(index) + " " + sample)

        file_h5 = os.path.join(str(row['kallisto']), "abundance.h5")
        sample_data = ur.read_kall_project(file_h5, num_features,
                                           transcripts_len, args, ids_genes)
        if args.BS != 0:
            kallisto[:,cpt : (cpt + args.BS)] = sample_data
            cpt += args.BS
        else:
            kallisto[:,cpt : (cpt + 1)] = sample_data
            cpt += 1


    if args.BS == 0:
        header = design['sample']
    else:
        header = np.repeat(design['sample'].values, args.BS)
        header_id = np.tile(["_" + str(x) for x in range(0,args.BS)], design.shape[0])
        header = np.core.defchararray.add(header.astype(str), header_id.astype(str))
        design = design.loc[design.index.repeat(args.BS)]
        design.is_copy = None
        design['sample'] = header

    data = pd.DataFrame(data=kallisto, index=rows_id, columns=header)
    data = data.loc[(data.sum(axis=1) != 0), :]
    data.is_copy = None

    return(data, design)

def read_design_matrix(args, df_anno=None):
    design = ur.get_design(args)
    if args.KAL:
        data, design = create_kal_mat(args, design, df_anno)
    else:
        data = pd.io.parsers.read_csv(args.MATRIX, sep="\t", index_col=0)
        if args.GENE:
            #TODO
            sys.stderr.write("Not implemented!!! the analysis will be made on transcript.\n")

    if sum(~design["sample"].isin(data.columns)) > 0:
        sys.stderr.write("WARNING: Some samples are present in the design, but not in the quantification matrix\n")
        sys.stderr.write("\t the analysis will be done without these samples:\n")
        sys.stderr.write(str(design[~design["sample"].isin(data.columns)]) + "\n")
        design = design[design["sample"].isin(data.columns)]

    data = data.reindex(design["sample"], axis=1)

    if args.CPM:
        f_norm = 1e6 /  data.iloc[:,1:].sum()

    data = data[(data.T != 0).any()]

    if args.CPM:
        data.iloc[:,1:] = data.iloc[:,1:] * f_norm

    list_ids = list(data.index)

    data = data.values
    if not args.LOG:
        data = np.log2(data + args.C)

    return(design, data, list_ids)

def main_pred_rna(args, argparser):

    df_anno = None
    if args.ANNO is not None:
        sys.stderr.write(time.strftime('%X') + ": Read annotation\n")
        df_anno = ur.gff_2_df(args.ANNO)
        df_anno.dropna(axis=0, subset=["ID"], inplace=True)
        df_anno["ID"] = df_anno["ID"].str.replace("gene:", "")
        df_anno["ID"] = df_anno["ID"].str.replace("transcript:", "")
        df_anno["Parent"] = df_anno["Parent"].str.replace("gene:", "")
        df_anno["Parent"] = df_anno["Parent"].str.replace("transcript:", "")
        df_anno.set_index("ID", inplace=True)
        print(df_anno)

    sys.stderr.write(time.strftime('%X') + ": Read design and matrix features\n")
    (design, data, list_ids) = read_design_matrix(args, df_anno)

    num_pred = data.shape[0]

    all_classifier = uo.compute_pred(args, num_pred, list_ids, data, design)
    uo.save_pred_res(args, all_classifier)
