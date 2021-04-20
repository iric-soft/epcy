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

    get_argparser_common_pred_part(parser, requiredNamed)
