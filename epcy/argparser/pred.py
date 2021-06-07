from .common import *

from .get_common_pred import *


def get_argparser_pred(parser):

    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument(
        "-m",
        dest="MATRIX",
        required=True,
        help="tsv file of features quantification matrix.",
        type=lambda x: is_valid_file(parser, x)
    )

    parser.add_argument(
        "--norm",
        dest="NORM",
        help="To apply a depth normalization on quantification matrix.",
        action='store_true'
    )

    parser.add_argument(
        "--nolog",
        dest="NOLOG",
        help="To specify that quantification matrix do not need to be log " +
             "transformed (like with percentage).",
        action='store_true'
    )

    parser.set_defaults(NORM=False)
    parser.set_defaults(NOLOG=False)

    get_argparser_common_pred_part(parser, requiredNamed)
