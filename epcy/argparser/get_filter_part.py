
def get_argparser_filter_part(parser):

    parser.add_argument("-e",
                        dest="EXP",
                        help="filter features with " +
                             "sum(values(features)) < (Default:0).",
                        type=float,
                        default=0)
    parser.add_argument("-l",
                        dest="L2FC",
                        help="filter trans/genes with log2FC < (Default:0.3).",
                        type=float,
                        default=0.3)
