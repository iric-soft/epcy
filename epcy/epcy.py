import argparse

from .argparser.pred import *
from .argparser.pred_rna import *

from .tools.pred import main_pred
from .tools.pred_rna import main_pred_rna


# ###########################################################################
# Main function
def main():

    argparser = argparse.ArgumentParser(prog='PROG')
    subparsers = argparser.add_subparsers(help='sub-command help')

    # create the argparser for the "pred" command
    pred = subparsers.add_parser(
        'pred',
        help='Compute predictive capability of each normalized features (Generic case).'
    )
    pred.set_defaults(func=main_pred)
    get_argparser_pred(pred)

    # create the argparser for the "pred_rna" command
    pred_rna = subparsers.add_parser(
        'pred_rna',
        help='Compute predictive capability of each genes/transcripts.'
    )
    pred_rna.set_defaults(func=main_pred_rna)
    get_argparser_pred_rna(pred_rna)


    # recover arguments
    args = argparser.parse_args()

    # execute the command
    args.func(args, argparser)
