from .common import *

from .get_common_profile import *


def get_argparser_profile(parser):

    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument(
        "-m",
        dest="MATRIX",
        required=True,
        help="tsv file of features matrix quantification.",
        type=lambda x: is_valid_file(parser, x)
    )

    get_argparser_common_profile_part(parser, requiredNamed)
