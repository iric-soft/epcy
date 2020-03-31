
def get_argparser_bootstrap_part(parser):

    parser.add_argument("--bs",
                        dest="BS",
                        help="Number of bootstrap (BS) used for each sample " +
                             "(Default: No BS). Use -1 to let EPCY find the " +
                             "number of BS in kallisto output.",
                        type=int,
                        default=0)
