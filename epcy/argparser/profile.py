from .common import *

from .get_log_part import *
from .get_design_part import *
from .get_bandwidth_part import *
from .get_output_part import *
from .get_matrix_part import *


def get_argparser_profile(parser):

    get_argparser_matrix(parser)
    get_argparser_log_part(parser)
    get_argparser_design_part(parser)
    get_argparser_bandwidth_part(parser)
    get_argparser_output_part(parser)

    parser.add_argument("--ids",
                        dest="IDS",
                        help='List of features id to plot',
                        type=str,
                        nargs='+',
                        default=[])
