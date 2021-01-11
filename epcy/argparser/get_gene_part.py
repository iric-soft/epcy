from .common import *


def get_argparser_gene_part(parser):

    parser.add_argument(
        "--gene",
        dest="GENE",
        help="If the quantification is compute on " +
             "transcripts, this option allow to calculate " +
             "predictive capability on genes, using " +
             "annotation file (--anno).",
        action='store_true'
    )

    parser.add_argument(
        "--anno",
        dest="ANNO",
        help="gff3 file of the feautres's annotation.",
        type=lambda x: is_valid_file(parser, x)
    )

    parser.set_defaults(GENE=False)
