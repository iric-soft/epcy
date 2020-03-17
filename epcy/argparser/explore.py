from .common import *

from .get_output_part import *
from .get_design_part import *


def get_argparser_explore(parser):

    get_argparser_design_part(parser)
    parser.add_argument("-p",
                        dest="PRED",
                        help="Path to EPCY predictive_capability output file.",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument("-s",
                        dest="SUBG",
                        help="Path to EPCY subgroup_predicted output file.",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument("--top",
                        dest="TOP",
                        help="Top x features used to explore. (Default: 10)",
                        type=int,
                        default=10)

    get_argparser_output_part(parser)
