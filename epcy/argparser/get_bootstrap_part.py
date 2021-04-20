
def get_argparser_bootstrap_part(parser):

    parser.add_argument(
        "--bs",
        dest="BS",
        help="Number of bootstrap (BS) used for each sample " +
             "(Default: No BS).",
        type=int,
        default=0
    )
