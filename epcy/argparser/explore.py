from .common import *

from .get_output_part import *
from .get_design_part import *


def get_argparser_explore(parser):

    requiredNamed = parser.add_argument_group('required arguments')

    get_argparser_design_part(parser, requiredNamed)
    parser.add_argument(
        "-p",
        dest="PRED",
        required=True,
        help="Path to EPCY predictive_capability output file.",
        type=lambda x: is_valid_file(parser, x)
    )

    parser.add_argument(
        "-s",
        dest="SUBG",
        help="Path to EPCY condition_predicted output file.",
        type=lambda x: is_valid_file(parser, x)
    )

    requiredNamed.add_argument(
        "--top",
        dest="TOP",
        required=True,
        help="Top x features used to explore. (Default: 10)",
        type=int,
        default=10
    )

    get_argparser_output_part(parser)
