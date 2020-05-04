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

    parser.add_argument("--size",
                        dest="SIZE",
                        help="Radius of a dot, in points. (Default:6.0)",
                        type=float,
                        default=6.0)

    parser.add_argument("--alpha",
                        dest="ALPHA",
                        help="Transparency of dot. (Default:1.0, no transparency)",
                        type=float,
                        default=1.0)

    requiredNamed.add_argument("--ids",
                               dest="IDS",
                               required=True,
                               help='List of features id to plot',
                               type=str,
                               nargs='+',
                               default=[])

    parser.set_defaults(STRIP=False)
