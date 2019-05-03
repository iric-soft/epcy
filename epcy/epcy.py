import argparse

from .argparser.pred import *

from .tools.pred import main_pred


# ###########################################################################
# Main function
def main():

    argparser = argparse.ArgumentParser(prog='PROG')
    subparsers = argparser.add_subparsers(help='sub-command help')

    # create the argparser for the "diff" command
    pred = subparsers.add_parser(
        'pred',
        help='Compute predictive capability of each features.'
    )
    pred.set_defaults(func=main_pred)
    get_argparser_diff(pred)

    # recover arguments
    args = argparser.parse_args()

    # execute the command
    args.func(args, argparser)
