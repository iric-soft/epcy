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
        "--condition",
        dest="CONDITION",
        help="In design file, the column name use for this analysis " +
             "(Default: condition).",
        type=str,
        default="condition"
    )

    parser.add_argument(
        "--query",
        dest="QUERY",
        help="To specify a query condition to compare to a reference " +
             "condition (Default: Query).",
        type=str,
        default="Query"
    )

    parser.add_argument(
        "--ref",
        dest="REF",
        help="To specify a reference condition to compare to a query " +
             "condition. Without all other samples will be use as reference " +
             "condition. (Default: None).",
        type=str,
        default=None
    )
