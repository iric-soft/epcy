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


def read_kall_project(file_h5, num_kal_bs, transcripts_len, args,
                      ids_genes, start=-1, end=-1):
    """ Open an kallisto h5 file and return tpm as a dict.
    """
    f = h5py.File(file_h5, 'r', libver='latest')

    all_tpm = np.zeros([num_kal_bs, args.BY], dtype=np.float32)
    total_mass = 0.0

    if args.KAL and args.BS > 0:
        for cpt_bs in range(0, num_kal_bs):
            counts = f["bootstrap/bs" + str(cpt_bs)][:] / transcripts_len
            if args.GENE:
                counts = trans_to_gene(counts, ids_genes)

            total_mass = np.sum(counts)
            if start != -1 and end != -1:
                counts = counts[start:end:1]

            for cpt_elt in range(0, counts.size):
                epsilon = abs(gauss(0, 0.0000000000001))
                tpm_value = (counts[cpt_elt] / total_mass) * 1e6
                all_tpm[cpt_bs][cpt_elt] = np.log2(tpm_value + args.C) + epsilon
    else:
        counts = f["est_counts"][:] / transcripts_len
        if args.GENE:
            counts = trans_to_gene(counts, ids_genes)
        total_mass = np.sum(counts)
        if start != -1 and end != -1:
            counts = counts[start:end:1]

        for cpt_elt in range(0, counts.size):
            epsilon = abs(gauss(0, 0.0000000000001))
            tpm_value = (counts[cpt_elt] / total_mass) * 1e6
            all_tpm[0][cpt_elt] = np.log2(tpm_value + args.C) + epsilon


    f.close()

    return(all_tpm)


def get_design(args):
    design = pd.read_csv(args.DESIGN, sep="\t")
    design[args.SUBGROUP] = [1 if condition == args.QUERY else 0 for condition in design[args.SUBGROUP]]
    design = design.sort_values(by=[args.SUBGROUP, 'sample'], ascending=[False, True])

    return(design)
