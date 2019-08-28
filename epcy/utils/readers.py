from collections import defaultdict
import gzip
import re
import pandas as pd
import numpy as np
from random import gauss

import h5py


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
    design[args.SUBGROUP] = [1 if condition == args.QUERY else 0 for condition in design[args.SUBGROUP]]
    design = design.sort_values(by=[args.SUBGROUP, 'sample'], ascending=[False, True])

    return(design)
