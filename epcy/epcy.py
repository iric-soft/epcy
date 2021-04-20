import argparse

from .argparser.pred import *
from .argparser.pred_rna import *
from .argparser.profile import *
from .argparser.profile_rna import *
from .argparser.qc import *
from .argparser.kal2mat import *
from .argparser.explore import *
from .argparser.ct import *

from .tools.pred import main_pred
from .tools.pred_rna import main_pred_rna
from .tools.profile import main_profile
from .tools.profile_rna import main_profile_rna
from .tools.qc import main_qc
from .tools.kal2mat import main_kal2mat
from .tools.explore import main_explore
from .tools.ct import main_ct

import sys

# ###########################################################################
# Main function
def main():

    argparser = argparse.ArgumentParser(prog='PROG')
    subparsers = argparser.add_subparsers(help='sub-command help')

    # create the argparser for the "pred" command
    pred = subparsers.add_parser(
        'pred',
        help="Compute predictive capability of each normalized features " +
             "(Generic case)."
    )
    pred.set_defaults(func=main_pred)
    get_argparser_pred(pred)

    # create the argparser for the "pred_rna" command
    pred_rna = subparsers.add_parser(
        'pred_rna',
        help="Compute predictive capability of each genes/transcripts " +
             "expression."
    )
    pred_rna.set_defaults(func=main_pred_rna)
    get_argparser_pred_rna(pred_rna)

    # create the argparser for the "profile" command
    profile = subparsers.add_parser(
        'profile',
        help='Plot profile of a list of features.'
    )
    profile.set_defaults(func=main_profile)
    get_argparser_profile(profile)

    # create the argparser for the "profile" command
    profile_rna = subparsers.add_parser(
        'profile_rna',
        help='Plot profile of a list of genes/transcipts.'
    )
    profile_rna.set_defaults(func=main_profile_rna)
    get_argparser_profile_rna(profile_rna)

    # create the argparser for the "qc" command
    qc = subparsers.add_parser(
        'qc',
        help='Plot quality conrol gaph.'
    )
    qc.set_defaults(func=main_qc)
    get_argparser_qc(qc)

    # create the argparser for the "kal2mat" command
    kal2mat = subparsers.add_parser(
        'kal2mat',
        help="Build and save matrix expression from kallisto quantification " +
             "h5 files."
    )
    kal2mat.set_defaults(func=main_kal2mat)
    get_argparser_kal2mat(kal2mat)

    # create the argparser for the "explore" command
    explore = subparsers.add_parser(
        'explore',
        help='Create figures to explore subgroup_predicted.xls.'
    )
    explore.set_defaults(func=main_explore)
    get_argparser_explore(explore)

    # create the argparser for the "ct" command
    ct = subparsers.add_parser(
        'ct',
        help="Return a contingency table by feature, " +
             "using subgroup_predicted.xls."
    )
    ct.set_defaults(func=main_ct)
    get_argparser_ct(ct)

    # recover arguments
    args = argparser.parse_args()

    if not len(sys.argv) > 1:
        sys.stderr.write("WARNING: use epcy -h, if you need help\n")
        exit(0)

    # execute the command
    args.func(args, argparser)
