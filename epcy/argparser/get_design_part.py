from .common import *


def get_argparser_design_part(parser, requiredNamed):

    requiredNamed.add_argument(
        "-d",
        dest="DESIGN",
        required=True,
        help="A path to tabulated file who discribe your design.",
        type=lambda x: is_valid_file(parser, x)
    )

    parser.add_argument(
        "--subgroup",
        dest="SUBGROUP",
        help="In design file, the column name use for this analysis " +
             "(Default: subgroup).",
        type=str,
        default="subgroup"
    )

    parser.add_argument(
        "--query",
        dest="QUERY",
        help="Value of queried subgroup to compare to all other samples " +
             "(Default: Query).",
        type=str,
        default="Query"
    )
