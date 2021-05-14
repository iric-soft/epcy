from .common import *

from .get_output_part import *
from .get_design_part import *


def get_argparser_ct(parser):

    requiredNamed = parser.add_argument_group('required arguments')

    get_argparser_design_part(parser, requiredNamed)

    parser.add_argument(
        "-s",
        dest="SUBG",
        help="Path to EPCY condition_predicted output file.",
        type=lambda x: is_valid_file(parser, x)
    )

    get_argparser_output_part(parser)
