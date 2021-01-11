from .common import *

from .get_common_profile import *
from .get_rna_norm_part import *
from .get_kallisto_part import *
from .get_gene_part import *


def get_argparser_profile_rna(parser):

    requiredNamed = parser.add_argument_group('required arguments')

    parser.add_argument(
        "-m",
        dest="MATRIX",
        help="tsv file of features matrix quantification.",
        type=lambda x: is_valid_file(parser, x)
    )

    get_argparser_common_profile_part(parser, requiredNamed)
    get_argparser_rna_norm_part(parser)
    get_argparser_kallisto_part(parser)
    get_argparser_gene_part(parser)
