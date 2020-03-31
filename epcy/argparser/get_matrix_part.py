from .common import *


def get_argparser_matrix(parser):

    parser.add_argument("-m",
                        dest="MATRIX",
                        help="tsv file of features matrix quantification.",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument("--replacena",
                        dest="REPLACE_NA",
                        help="To specify a value to replace missing value. " +
                             "Without (and if missing value is present) " +
                             "predictive capabilty will be compute without " +
                             "samples with missing value into the feature. " +
                             "In that case, the number of samples used for " +
                             "each group will be reported for each features.",
                        type=float,
                        default=None)
