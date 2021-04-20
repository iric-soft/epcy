from .common import *

from .get_design_part import *
from .get_log_part import *
from .get_bandwidth_part import *
from .get_output_part import *
from .get_matrix_part import *


def get_argparser_common_profile_part(parser, requiredNamed):

    get_argparser_matrix(parser)
    get_argparser_log_part(parser)
    get_argparser_design_part(parser, requiredNamed)
    get_argparser_bandwidth_part(parser)
    get_argparser_output_part(parser)

    parser.add_argument("--strip",
                        dest="STRIP",
                        help='To create a strip plot.',
                        action='store_true')

    parser.add_argument("--violin",
                        dest="VIOLIN",
                        help='To create a violin plot.',
                        action='store_true')

    parser.add_argument("--no_density",
                        dest="NO_DENSITY",
                        help='remove density plot.',
                        action='store_true')

    parser.add_argument("--size",
                        dest="SIZE",
                        help="Radius of a dot, in points. (Default:5.0)",
                        type=float,
                        default=5.0)

    requiredNamed.add_argument("--ids",
                               dest="IDS",
                               required=True,
                               help='List of features id to plot',
                               type=str,
                               nargs='+',
                               default=[])

    parser.set_defaults(STRIP=False)
    parser.set_defaults(NO_DENSITY=False)
