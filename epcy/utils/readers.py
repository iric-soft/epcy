import os
import re
import sys
import time

import gzip
import h5py

import pandas as pd
import numpy as np

from random import gauss
from collections import defaultdict

from . import other as uo

GTF_HEADER = ['seqname', 'source', 'feature', 'start', 'end', 'score',
              'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA = re.compile(r'\s*,\s*')
R_KEYVALUE = re.compile(r'(\s+|\s*=\s*)')

# GTF/GFF reader FROM: https://gist.github.com/slowkow/8101481
def gff_2_df(filename):
    """Open an optionally gzipped GTF (or GFF) file and return a pandas.DataFrame.
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)

    filename = str(filename)
    for i, line in enumerate(lines(filename)):
        # Ensure this row has some value for each column.
        for key in ["seqname", "feature", "ID", "Parent", "Name", "biotype", "description"]:
            result[key].append(line.get(key, None))

    return(pd.DataFrame(result))


def lines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open
    filtred_fields = ["exon", "chromosome", "CDS"]

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue

            fields = line.split('\t')
            if fields[2] in filtred_fields:
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single GTF line and return a dict.
    """
    result = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)

    return(result)


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return(value)


def get_trans_by_gene(transcripts, ids_genes):
    trans_by_gene = []
    for ids in ids_genes:
        trans_by_gene.append(','.join(transcripts[ids]))

    return(trans_by_gene)


def trans_to_gene(values, ids_genes):
    counts = np.zeros(len(ids_genes))

    cpt = 0
    for ids in ids_genes:
        counts[cpt] = sum(values[ids])
        cpt += 1

    return(counts)

def counts2tpm(counts, args):
    """ Transform counts into tpm.
    """
    total_mass = np.sum(counts)
    tpm_value = (counts / total_mass) * 1e6

    return(np.log2(tpm_value + args.C))

def read_kall_project(file_h5, num_features, transcripts_len, args,
                      ids_genes):
    """ Open an kallisto h5 file and return tpm as a dict.
    """
    f = h5py.File(file_h5, 'r', libver='latest')

    if args.BS != 0:
        sample_data = np.zeros([num_features, args.BS], dtype=np.float32)
        for cpt_bs in range(0, args.BS):
            counts = f["bootstrap/bs" + str(cpt_bs)][:] / transcripts_len
            if args.GENE:
                counts = trans_to_gene(counts, ids_genes)

            if args.CPM:
                sample_data[:, cpt_bs] = counts
            else:
                sample_data[:, cpt_bs] = counts2tpm(counts, args)
    else:
        sample_data = np.zeros([num_features, 1], dtype=np.float32)
        counts = f["est_counts"][:] / transcripts_len
        if args.GENE:
            counts = trans_to_gene(counts, ids_genes)

        if args.CPM:
            sample_data[:, 0] = counts
        else:
            sample_data[:, 0] = counts2tpm(counts, args)

    f.close()

    return(sample_data)


def get_design(args):
    design = pd.read_csv(args.DESIGN, sep="\t")
    drop_ids = design[ design[args.SUBGROUP] == 'None' ].index
    design.drop(drop_ids , inplace=True)
    design[args.SUBGROUP] = [1 if condition == args.QUERY else 0 for condition in design[args.SUBGROUP]]
    design = design.sort_values(by=[args.SUBGROUP, 'sample'], ascending=[False, True])

    return(design)


def read_design_matrix(args):
    design = get_design(args)
    data = pd.io.parsers.read_csv(args.MATRIX, sep="\t", index_col=0)
    if sum(~design["sample"].isin(data.columns)) > 0:
        sys.stderr.write("WARNING: Some samples are present in the design, but not in the quantification matrix\n")
        sys.stderr.write("\t the analysis will be done without these samples:\n")
        sys.stderr.write(str(design[~design["sample"].isin(data.columns)]) + "\n")
        design = design[design["sample"].isin(data.columns)]

    data = data.reindex(design["sample"], axis=1)

    data = data[(data.T != 0).any()]

    list_ids = np.array(data.index)

    data = data.values
    if not args.LOG:
        data = np.log2(data + args.C)

    return(design, data, list_ids)


def create_kal_mat_rna(args, design, df_anno):
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
        sample_data = read_kall_project(file_h5, num_features,
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


def read_design_matrix_rna(args, df_anno=None):
    design = get_design(args)
    if args.KAL:
        data, design = create_kal_mat_rna(args, design, df_anno)
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

    list_ids = np.array(data.index)

    data = data.values
    if not args.LOG:
        data = np.log2(data + args.C)

    return(design, data, list_ids)

def read_anno(args):
    df_anno = None
    if args.ANNO is not None:
        sys.stderr.write(time.strftime('%X') + ": Read annotation\n")
        df_anno = gff_2_df(args.ANNO)
        df_anno.dropna(axis=0, subset=["ID"], inplace=True)
        df_anno["ID"] = df_anno["ID"].str.replace("gene:", "")
        df_anno["ID"] = df_anno["ID"].str.replace("transcript:", "")
        df_anno["Parent"] = df_anno["Parent"].str.replace("gene:", "")
        df_anno["Parent"] = df_anno["Parent"].str.replace("transcript:", "")
        df_anno.set_index("ID", inplace=True)
        #print(df_anno)

    return(df_anno)
