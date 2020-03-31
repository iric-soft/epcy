from .common import *


def get_argparser_design_part(parser):

    parser.add_argument("-d",
                        dest="DESIGN",
                        help="Tabulated file who discribe your design.",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument("--query",
                        dest="QUERY",
                        help="Query value specify in the subgroup column " +
                             "for samples queried (Default: Query).",
                        type=str,
                        default="Query")

    parser.add_argument("--subgroup",
                        dest="SUBGROUP",
                        help="Header name of the subgroup column in your " +
                             "design file (Default: subgroup).",
                        type=str,
                        default="subgroup")
