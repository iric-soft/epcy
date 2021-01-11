from .common import *

from .get_rna_norm_part import *
from .get_kallisto_part import *
from .get_gene_part import *
from .get_common_pred import *


def get_argparser_pred_rna(parser):

    requiredNamed = parser.add_argument_group('required arguments')

    parser.add_argument(
        "-m",
        dest="MATRIX",
        help="tsv file of features quantification matrix.",
        type=lambda x: is_valid_file(parser, x)
    )

    get_argparser_common_pred_part(parser, requiredNamed)

    get_argparser_gene_part(parser)
    get_argparser_kallisto_part(parser)
    get_argparser_rna_norm_part(parser)
