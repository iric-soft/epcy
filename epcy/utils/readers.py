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
    """Open an optionally gzipped GTF (or GFF) file and return a
       pandas.DataFrame.
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)

    filename = str(filename)
    for i, line in enumerate(lines(filename)):
        # Ensure this row has some value for each column.
        keys = ["seqname", "feature", "ID", "Parent", "Name", "biotype",
                "description"]
        for key in keys:
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
    """ Merge transcript quantification to work on gene level.
    """
    counts = np.zeros(len(ids_genes))

    cpt = 0
    for ids in ids_genes:
        counts[cpt] = sum(values[ids])
        cpt += 1

    return(counts)


def counts2cpmed(counts):
    """ Transform counts matrix into cpmed.
    """
    total_mass = np.nansum(counts, axis=0)
    med_count = np.median(total_mass)
    counts = (counts / total_mass) * med_count

    return(counts)


def counts2cpm(counts):
    """ Transform counts matrix into cpm.
    """
    total_mass = np.nansum(counts, axis=0)
    counts = (counts / total_mass) * 1e6

    return(counts)


def counts2tpm(counts, transcripts_len):
    """ Transform counts arrays into tpm.
    """
    counts = counts / transcripts_len
    total_mass = np.sum(counts)
    tpm_value = (counts / total_mass) * 1e6

    return(tpm_value)


def read_kall_project(file_h5, num_features, transcripts_len, args,
                      ids_genes):
    """ Open an kallisto h5 file and return tpm as a dict.
    """
    f = h5py.File(file_h5, 'r', libver='latest')

    if args.BS != 0:
        sample_data = np.zeros([num_features, args.BS], dtype=np.float64)
        for cpt_bs in range(0, args.BS):
            counts = f["bootstrap/bs" + str(cpt_bs)][:]

            if args.TPM:
                counts = counts2tpm(counts, transcripts_len)

            if args.GENE:
                counts = trans_to_gene(counts, ids_genes)

            sample_data[:, cpt_bs] = counts
    else:
        sample_data = np.zeros([num_features, 1], dtype=np.float64)
        counts = f["est_counts"][:]

        if args.TPM:
            counts = counts2tpm(counts, transcripts_len)

        if args.GENE:
            counts = trans_to_gene(counts, ids_genes)

        sample_data[:, 0] = counts

    f.close()

    return(sample_data)


def get_design(args):
    """ Read and format design.
    """
    design = pd.read_csv(args.DESIGN, sep="\t")

    if "sample" not in design.columns.values:
        sys.stderr.write("ERROR: The design file need to have a " +
                         "column 'sample'.\n")
        exit(-1)

    if args.CONDITION not in design.columns.values:
        sys.stderr.write("ERROR: The design file need to have a " +
                         "column '" + args.CONDITION + "' or specify a " +
                         "column name using --condition.\n")
        exit(-1)

    if args.QUERY not in design[args.CONDITION].values:
        sys.stderr.write("ERROR: " + args.QUERY + " is not present in the  " +
                         "column '" + args.CONDITION + "' of the design " +
                         "file.\n")
        exit(-1)

    if args.REF is not None and args.REF not in design[args.CONDITION].values:
        sys.stderr.write("ERROR: " + args.REF + " is not present in the  " +
                         "column '" + args.CONDITION + "' of the design " +
                         "file.\n")
        exit(-1)

    design['sample'] = design['sample'].apply(str)
    drop_ids = design[design[args.CONDITION] == 'None'].index
    design.drop(drop_ids, inplace=True)
    if args.REF is not None:
        design.drop(
            design[~design[args.CONDITION].isin([args.REF, args.QUERY])].index,
            inplace=True
        )

    if hasattr(args, 'SHUFFLE') and args.SHUFFLE:
        design[args.CONDITION] = np.random.permutation(
            design[args.CONDITION].values
        )
        sys.stderr.write(time.strftime('%X') + ": Save shuffled design\n")
        if args.PATH_OUT is not None:
            if not os.path.exists(args.PATH_OUT):
                os.makedirs(args.PATH_OUT)

            file_out = args.PATH_OUT + "/shuffled_design.tsv"
            design.to_csv(file_out, index=False, sep="\t")

    design[args.CONDITION] = [1 if condition == args.QUERY else 0
                              for condition in design[args.CONDITION]]

    design = design.sort_values(by=[args.CONDITION, 'sample'],
                                ascending=[False, True])

    return(design)


def read_design_matrix(args):
    design = get_design(args)

    f = open(args.MATRIX, "r")
    lines = f.readlines()
    matrix_samples = lines[0]
    matrix_samples = matrix_samples.rstrip().split('\t')
    matrix_samples = matrix_samples[1:]
    matrix_num_samples = len(matrix_samples)
    if sum(~design["sample"].isin(matrix_samples)) > 0:
        sys.stderr.write("WARNING: Some samples are present in the design, " +
                         "but not in the quantification matrix\n")
        sys.stderr.write("\t the analysis will be made without these " +
                         "samples:\n")
        sys.stderr.write(str(design[~design["sample"].isin(matrix_samples)]) +
                         "\n")
        design = design[design["sample"].isin(matrix_samples)]

    list_ids = [x.split('\t')[0] for x in lines[1:]]
    f.close()

    data = np.loadtxt(
        args.MATRIX,
        usecols=range(1, matrix_num_samples+1, 1),
        skiprows=1, dtype=np.float64
    )
    if hasattr(args, 'REPLACE_NA') and args.REPLACE_NA is not None:
        data = np.nan_to_num(data, nan=args.REPLACE_NA)

    row_ids_0 = ~np.all(data == 0, axis=1)
    num_removed = len(row_ids_0) - sum(row_ids_0)
    if num_removed != 0:
        sys.stderr.write(time.strftime('%X') + ": " + str(num_removed) +
                         " features with sum==0 have been removed.\n")

    list_ids = [list_ids[x] for x in np.where(row_ids_0)[0]]
    list_ids = np.asarray(list_ids)
    data = data[row_ids_0]

    ids_sorted = [i for x in design["sample"] for i,y in enumerate(matrix_samples) if y == x]
    data = data[:, ids_sorted]

    #data = pd.io.parsers.read_csv(args.MATRIX, sep="\t", index_col=0)
    #if sum(~design["sample"].isin(data.columns)) > 0:
    #    sys.stderr.write("WARNING: Some samples are present in the design, " +
    #                     "but not in the quantification matrix\n")
    #    sys.stderr.write("\t the analysis will be made without these " +
    #                     "samples:\n")
    #    sys.stderr.write(str(design[~design["sample"].isin(data.columns)]) +
    #                     "\n")
    #    design = design[design["sample"].isin(data.columns)]
    #data = data.reindex(design["sample"], axis=1)
    #data = data[(data.T != 0).any()]
    #list_ids = np.array(data.index)
    #data = data.values

    if args.LOG:
        data = np.log2(data + args.C)

    return(design, data, list_ids)


def create_kal_mat(args, design, design_bootstrapped, df_anno):
    """ Create features matrix from kallisto h5 files.
        This function use kallisto column in the design file.
    """

    file_name = os.path.join(str(design["kallisto"][0]), "abundance.h5")
    sys.stderr.write(time.strftime('%X') + ":\tRead kallisto info\n")
    f = h5py.File(file_name, 'r')

    if args.BS == -1:
        args.BS = f["aux/num_bootstrap"][0]

    transcripts = f["aux/ids"][:]
    transcripts = [uo.cut_version(trans.decode('utf-8'))
                   for trans in transcripts]
    transcripts = np.array(transcripts)
    uniq_genes = None

    ids_genes = None
    if args.GENE:
        sys.stderr.write(time.strftime('%X') + ":\tTranscripts to genes\n")
        parents = df_anno.reindex(transcripts)["Parent"].values
        ids_notin_anno = np.argwhere(pd.isnull(parents))
        parents[ids_notin_anno] = transcripts[ids_notin_anno]
        # TODO find an other way to works with other ids not ENSG
        uniq_genes = df_anno.index[df_anno.index.str.contains("ENSG")]
        uniq_genes = np.array(uniq_genes)
        # To add transcropts not in annotation file
        uniq_genes = np.append(uniq_genes, transcripts[ids_notin_anno])
        sys.stderr.write(time.strftime('%X') + ":\t\tkeep " +
                         str(uniq_genes.size) + " ENSG ids\n")
        sys.stderr.write(time.strftime('%X') + ":\t\tfound genes ids for " +
                         "each trans\n")
        ids_genes = [np.where(parents == gene)[0] for gene in uniq_genes]
        ids_genes = np.array(ids_genes, dtype=object)

    transcripts_len = f["aux/eff_lengths"][:]
    f.close()

    rows_id = transcripts
    num_features = transcripts.size
    if args.GENE:
        rows_id = uniq_genes
        num_features = uniq_genes.size

    sys.stderr.write(time.strftime('%X') + ":\tRead samples kallisto " +
                     "quantification from h5\n")
    if args.BS == 0:
        kallisto = np.zeros((num_features, design.shape[0]), dtype=np.float64)
    else:
        kallisto = np.zeros((num_features, design.shape[0]*args.BS),
                            dtype=np.float64)

    cpt = 0
    for index, row in design.iterrows():
        sample = row['sample']
        # print(str(index) + " " + sample)

        file_h5 = os.path.join(str(row['kallisto']), "abundance.h5")
        sample_data = read_kall_project(file_h5, num_features,
                                        transcripts_len, args, ids_genes)
        if args.BS != 0:
            kallisto[:, cpt:(cpt + args.BS)] = sample_data
            cpt += args.BS
        else:
            kallisto[:, cpt:(cpt + 1)] = sample_data
            cpt += 1

    data = pd.DataFrame(data=kallisto, index=rows_id,
                        columns=design_bootstrapped['sample'])

    return(data)


def bootstrapped_design(design, args):
    """ Update generic design to match with the number of botstrap.
    """
    if args.BS != 0:
        header = np.repeat(design['sample'].values, args.BS)
        header_id = np.tile(["_" + str(x)
                             for x in range(0, args.BS)], design.shape[0])
        header = np.core.defchararray.add(header.astype(str),
                                          header_id.astype(str))
        design_bs = design.loc[design.index.repeat(args.BS)]
        design_bs['sample'] = header

        return(design_bs)

    return(design)


def read_design_matrix_rna(args, df_anno=None):
    design = get_design(args)
    design_bootstrapped = design
    if args.BS is not None and args.BS > 0:
        design_bootstrapped = bootstrapped_design(design, args)

    num_query = len(np.where(design[args.CONDITION] == 1)[0])
    if num_query == 0:
        sys.stderr.write("ERROR: EPCY havn't found Query samples in your " +
                         "design!\n\tCheck your design file and --condition, " +
                         "--query options\n")
        return(None, None, None)

    num_ref = len(np.where(design[args.CONDITION] == 0)[0])
    if num_ref == 0:
        sys.stderr.write("ERROR: EPCY havn't found Query samples in your " +
                         "design!\n\tCheck your design file and --condition, " +
                         "--ref options\n")
        return(None, None, None)

    if hasattr(args, 'MATRIX') and args.MATRIX is not None:
        if args.BS is not None and args.BS > 0:
            design = design_bootstrapped

        if args.GENE:
            # TODO
            sys.stderr.write("ERROR: Sorry, switch transcript to gene " +
                             "quantification from a matrix file is not " +
                             "implemented!!!\n")
            return(None, None, None)
        f = open(args.MATRIX, "r")
        lines = f.readlines()
        matrix_samples = lines[0]
        matrix_samples = matrix_samples.rstrip().split('\t')
        matrix_samples = matrix_samples[1:]
        matrix_num_samples = len(matrix_samples)
        if sum(~design["sample"].isin(matrix_samples)) > 0:
            sys.stderr.write("WARNING: Some samples are present in the design, " +
                             "but not in the quantification matrix\n")
            sys.stderr.write("\t the analysis will be made without these " +
                             "samples:\n")
            sys.stderr.write(str(design[~design["sample"].isin(matrix_samples)]) +
                             "\n")
            design = design[design["sample"].isin(matrix_samples)]

        list_ids = [x.split('\t')[0] for x in lines[1:]]
        f.close()

        data = np.loadtxt(
            args.MATRIX,
            usecols=range(1, matrix_num_samples+1, 1),
            skiprows=1, dtype=np.float64
        )
        if hasattr(args, 'REPLACE_NA') and args.REPLACE_NA is not None:
            data = np.nan_to_num(data, nan=args.REPLACE_NA)

        row_ids_0 = ~np.all(data == 0, axis=1)
        num_removed = len(row_ids_0) - sum(row_ids_0)
        if num_removed != 0:
            sys.stderr.write(time.strftime('%X') + ": " + str(num_removed) +
                             " features with sum==0 have been removed.\n")

        list_ids = [list_ids[x] for x in np.where(row_ids_0)[0]]
        list_ids = np.asarray(list_ids)
        data = data[row_ids_0]

        ids_sorted = [i for x in design["sample"] for i,y in enumerate(matrix_samples) if y == x]
        data = data[:, ids_sorted]
    else:
        if args.KAL:
            data = create_kal_mat(args, design, design_bootstrapped, df_anno)

            if hasattr(args, 'REPLACE_NA') and args.REPLACE_NA is not None:
                data = data.fillna(args.REPLACE_NA)

            if args.KAL:
                design = design_bootstrapped

            if sum(~design["sample"].isin(data.columns)) > 0:
                sys.stderr.write("WARNING: Some samples are present in the design, " +
                                 "but not in the quantification matrix\n")
                sys.stderr.write("\t the analysis will be made without these \
                                 samples:\n")
                sys.stderr.write(str(design[~design["sample"].isin(data.columns)]) +
                                 "\n")
                sys.stderr.write("Samples in matrix:\n")
                sys.stderr.write(str(data.columns) + "\n")
                design = design[design["sample"].isin(data.columns)]

            # Select and order sample column in fonction of design
            data = data.reindex(design["sample"], axis=1)

            # delete rows with only 0 into it
            data = data[(data.T != 0).any()]

            list_ids = np.array(data.index)

            data = data.values
        else:
            sys.stderr.write("ERROR: No quantification matrix can be find!\n")
            return(None, None, None)

    if not args.TPM:
        if args.CPM:
            data = counts2cpm(data)
        if args.CPMED:
            data = counts2cpmed(data)

    if args.LOG:
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
        # print(df_anno)

    return(df_anno)
