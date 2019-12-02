from .common import *


def get_argparser_qc(parser):

    parser.add_argument("-p",
                        dest="PRED",
                        help="EPCY predictive_capability output file.",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument("-o",
                        dest="PATH_OUT",
                        help="Path to the directory output.",
                        type=str,
                        default=None)

    parser.add_argument("--l2fc",
                        dest="L2FC",
                        help="Filled in foction of log2(FoldChange). (Default: max of mean expression)",
                        action='store_true')

    parser.add_argument("--ylog",
                        dest="YLOG",
                        help="To set yscale to log",
                        action='store_true')

    parser.set_defaults(L2FC=False)
    parser.set_defaults(YLOG=False)
