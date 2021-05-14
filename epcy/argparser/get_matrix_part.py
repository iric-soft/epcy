from .common import *


def get_argparser_matrix(parser):

    parser.add_argument(
        "--replacena",
        dest="REPLACE_NA",
        help="To specify a constant to replace missing " +
             "value. Without (and if missing value is " +
             "present) predictive capabilty will be compute " +
             "without samples with missing value, for each " +
             "feature. In that case, EPCY output will " +
             "report the number of samples used in each " +
             "condition and each features.",
        type=float,
        default=None
    )
