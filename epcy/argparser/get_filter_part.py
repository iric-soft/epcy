
def get_argparser_filter_part(parser):

    parser.add_argument(
        "-e",
        dest="EXP",
        help="filter features with " +
             "sum(values(features)) < (Default:None, disable).",
        type=float,
        default=None
    )

    parser.add_argument(
        "-l",
        dest="L2FC",
        help="filter trans/genes with abs(log2FC) < (Default:0).",
        type=float,
        default=0
    )
