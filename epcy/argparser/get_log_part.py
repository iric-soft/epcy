
def get_argparser_log_part(parser):
    parser.add_argument(
        "--log",
        dest="LOG",
        help="To apply a log2 transformation log2(x + C) " +
             "(see -c) ).",
        action='store_true'
    )

    parser.add_argument(
        "-c",
        dest="C",
        help="Constant value used during log " +
             "transformation, log2(x+C) (Default: C=1).",
        type=float,
        default=1.0
    )

    parser.set_defaults(LOG=False)
