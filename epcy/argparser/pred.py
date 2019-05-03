from .common import *


def get_argparser_diff(parser):
    parser.add_argument("-c",
                        dest="C",
                        help="Constant value used during log transformation, log2(x+C) (Default: C=1).",
                        type=float,
                        default=1.0)
    parser.add_argument("-d",
                        dest="DESIGN",
                        help="Tabulated file who discribe your design.",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-e",
                        dest="EXP",
                        help="filter features with sum(values(features)) < (Default:0).",
                        type=float,
                        default=0)
    parser.add_argument("-l",
                        dest="L2FC",
                        help="filter trans/genes with log2FC < (Default:0.3).",
                        type=float,
                        default=0.3)
    parser.add_argument("-m",
                        dest="MATRIX",
                        help="tsv file of the annotation used.",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-o",
                        dest="PATH_OUT",
                        help="Path to the directory output.")
    parser.add_argument("-t",
                        dest="THREAD",
                        help="Number of thread.",
                        type=int,
                        default=1)
    parser.add_argument("--anno",
                        dest="ANNO",
                        help="(Optional) gff3 file of the annotation used.",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("--by",
                        dest="BY",
                        help="Number of feature by thread. By default this number automaticaly compute in function of number of thread to optimise execution time, but you can specify a specific value if you encounter memory issue.",
                        type=int,
                        default=-1)
    parser.add_argument("--log",
                        dest="LOG",
                        help="Add this parameter, if your data are already log transformed.",
                        action='store_true')
    parser.add_argument("--query",
                        dest="QUERY",
                        help="Query value specify in the subgroup column for samples queried (Default: Query).",
                        type=str,
                        default="Query")
    parser.add_argument("--subgroup",
                        dest="SUBGROUP",
                        help="Header name of the subgroup column in your design file (Default: subgroup).",
                        type=str,
                        default="subgroup")

    parser.set_defaults(LOG=False)
