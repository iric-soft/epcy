from .get_bootstrap_part import *


def get_argparser_kallisto_part(parser):

    parser.add_argument(
        "--kal",
        dest="KAL",
        help='To work with kallisto quantification output. ' +
             'A collumn "kallisto" which contains the path ' +
             'of kallisto output folder for each sample, ' +
             'need to be added to the design file.',
        action='store_true')

    get_argparser_bootstrap_part(parser)

    parser.set_defaults(KAL=False)
