
def get_argparser_bandwidth_part(parser):

    parser.add_argument("--min_bw",
                        dest="MIN_BW",
                        help="To compute KDE MCC a bandwidth need to estimate from data using bw_nrd0. To avoid very small bw you can use this parameter to set a minimum (Default:0.0).",
                        type=float,
                        default=0.0)
