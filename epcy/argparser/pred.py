from .common import *

from .get_design_part import *
from .get_log_part import *
from .get_filter_part import *
from .get_other_pred_part import *
from .get_bandwidth_part import *
from .get_output_part import *
from .get_matrix_part import *


def get_argparser_pred(parser):

    get_argparser_matrix(parser)
    get_argparser_design_part(parser)
    get_argparser_log_part(parser)
    get_argparser_filter_part(parser)
    get_argparser_other_pred_part(parser)
    get_argparser_bandwidth_part(parser)
    get_argparser_output_part(parser)

    parser.add_argument("-b",
                        dest="N_BAGGING",
                        help="Number of bagging performed for each " +
                             "prediction (Default: 1 = no bagging).",
                        type=int,
                        default=1)

    parser.add_argument("-t",
                        dest="THREAD",
                        help="Number of thread.",
                        type=int,
                        default=1)
    parser.add_argument("--by",
                        dest="BY",
                        help="Number of feature by thread. By default this " +
                             "number is automaticaly set " +
                             "(#features/#thread). If you encounter memory " +
                             "issues, you can try using lower values.",
                        type=int,
                        default=-1)

    parser.add_argument("--full",
                        dest="FULL",
                        help="enable full output files.",
                        action='store_true')

    parser.add_argument("--ndraw",
                        dest="N_DRAW",
                        help="Number of time that EPCY will draw a " +
                             "predicted class (default 100).",
                        type=int,
                        default=100)

    parser.add_argument("--randomseed",
                        dest="RANDOM_SEED",
                        help="To specify a random seed (Int). If None, the " +
                             "random number generator is the RandomState " +
                             "instance used by np.random.",
                        type=int,
                        default=None)

    parser.set_defaults(FULL=False)
