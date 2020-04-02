from .common import *

from .get_design_part import *
from .get_gene_part import *
from .get_rna_norm_part import *
from .get_log_part import *
from .get_bootstrap_part import *
from .get_output_part import *


def get_argparser_kal2mat(parser):

    requiredNamed = parser.add_argument_group('required arguments')

    get_argparser_design_part(parser, requiredNamed)
    get_argparser_rna_norm_part(parser)
    get_argparser_log_part(parser)
    get_argparser_gene_part(parser)
    get_argparser_bootstrap_part(parser)
    get_argparser_output_part(parser)

    parser.set_defaults(KAL=True)
