from .common import *
from .pred import *

def get_argparser_pred_rna(parser):
    get_argparser_pred(parser)

    parser.add_argument("--anno",
        dest="ANNO",
        help="(Optional) gff3 file of the feautres annotation.",
        type=lambda x: is_valid_file(parser, x))

    parser.add_argument("--cpm",
        dest="CPM",
        help="To apply a Count Par Million (CPM) normalization to the matrix given in input (-m)",
        action='store_true')

    parser.add_argument("--kal",
                        dest="KAL",
                        help='To work with kallisto quantification output. A collumn "kallisto" which contains the path of kallisto output folder for each sample, need to be added to the design file. By default EPCY will works on TPM unless --cpm is specified.',
                        action='store_true')

    parser.add_argument("--bs",
                        dest="BS",
                        help="Number of bootstrap (BS) used for each sample (Default: No BS). Use -1 to let EPCY find the number of BS in kallisto output.",
                        type=int,
                        default=0)

    parser.add_argument("--gene",
                        dest="GENE",
                        help="If the quantification is compute on transcripts, this option allow to calculate predictive capability on genes, using annotation file (--anno).",
                        action='store_true')

    parser.set_defaults(CPM=False)
    parser.set_defaults(GENE=False)
