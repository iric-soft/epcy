
def get_argparser_output_part(parser):

    parser.add_argument(
        "-o",
        dest="PATH_OUT",
        help="Path to the directory output.",
        type=str,
        default=None
    )
