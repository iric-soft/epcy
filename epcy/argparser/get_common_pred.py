from .common import *

from .get_design_part import *
from .get_log_part import *
from .get_filter_part import *
from .get_bandwidth_part import *
from .get_output_part import *
from .get_matrix_part import *


def get_argparser_common_pred_part(parser, requiredNamed):

    get_argparser_matrix(parser)
    get_argparser_design_part(parser, requiredNamed)
    get_argparser_log_part(parser)
    get_argparser_filter_part(parser)
    get_argparser_bandwidth_part(parser)
    get_argparser_output_part(parser)

    parser.add_argument(
        "-b",
        dest="N_BAGGING",
        help="Number of bagging performed for each " +
             "prediction (Default: 1, for no bagging).",
        type=int,
        default=1
    )

    parser.add_argument(
        "-t",
        dest="THREAD",
        help="Number of thread.",
        type=int,
        default=1
    )

    parser.add_argument(
        "--npv",
        dest="NPV",
        help="Compute NPV (Negative Predictive Value) for each features " +
             "(e.g. genes).",
        action='store_true'
    )

    parser.add_argument(
        "--ppv",
        dest="PPV",
        help="Compute PPV (Positive Predictive Value, or precision) for " +
             "each features (e.g. genes).",
        action='store_true'
    )

    parser.add_argument(
        "--tpr",
        dest="TPR",
        help="Compute TPR (True Positive Rate or sensitivity) for " +
             "each features (e.g. genes).",
        action='store_true'
    )

    parser.add_argument(
        "--tnr",
        dest="TNR",
        help="Compute TNR (True Negative Rate or specificity) for " +
             "each features (e.g. genes).",
        action='store_true'
    )

    parser.add_argument(
        "--fnr",
        dest="FNR",
        help="Compute FNR (False Negative Rate or miss rate) for " +
             "each features (e.g. genes).",
        action='store_true'
    )

    parser.add_argument(
        "--fpr",
        dest="FPR",
        help="Compute FPR (False Positive Rate, or fall-out) for " +
             "each features (e.g. genes).",
        action='store_true'
    )

    parser.add_argument(
        "--fdr",
        dest="FDR",
        help="Compute FDR (False Discovery Rate) for " +
             "each features (e.g. genes).",
        action='store_true'
    )

    parser.add_argument(
        "--for",
        dest="FOR",
        help="Compute FOR (False Omission Rate) for " +
             "each features (e.g. genes).",
        action='store_true'
    )

    parser.add_argument(
        "--ts",
        dest="TS",
        help="Compute TS (Treat Score or critical success index) for " +
             "each features (e.g. genes).",
        action='store_true'
    )

    parser.add_argument(
        "--acc",
        dest="ACC",
        help="Compute ACC scores (Accuracy) for each features (e.g. genes).",
        action='store_true'
    )

    parser.add_argument(
        "--f1",
        dest="F1",
        help="Compute F1 scores for each features (e.g. genes).",
        action='store_true'
    )

    parser.add_argument(
        "--auc",
        dest="AUC",
        help="Compute AUC for each features (e.g. genes).",
        action='store_true'
    )

    parser.add_argument(
        "--full",
        dest="FULL",
        help="Enable full output files.",
        action='store_true'
    )

    parser.add_argument(
        "--nfold",
        dest="N_FOLD",
        help="Number of fold (default is None, to perform a " +
             "Leave-one-out cross validation).",
        type=int,
        default=None
    )

    parser.add_argument(
        "--ndraw",
        dest="N_DRAW",
        help="Number of time that EPCY will draw a " +
             "predicted class (Default: 100).",
        type=int,
        default=100
    )

    parser.add_argument(
        "--normal",
        dest="NORMAL",
        help="Compute sample assignation using normal distributions " +
             "predictor (to replace KDE).",
        action='store_true'
    )

    parser.add_argument(
        "--randomseed",
        dest="RANDOM_SEED",
        help="To specify a random seed (Int). If None, the " +
             "random number generator is the RandomState " +
             "instance used by numpyp.random. (Default: None)",
        type=int,
        default=None
    )

    parser.add_argument(
        "--ttest",
        dest="TTEST",
        help="Compute a p-value using ttest_ind from scipy.stats.",
        action='store_true'
    )

    parser.add_argument(
        "--utest",
        dest="UTEST",
        help="Compute a p-value using Mann-Whitney from " +
             "scipy.stats. (NEED --auc)",
        action='store_true'
    )

    parser.add_argument(
        "--shuffle",
        dest="SHUFFLE",
        help="Perform an analysis on shuffled design. Useful to create a " +
             "null distribution and find a cutoff.",
        action='store_true'
    )

    parser.set_defaults(AUC=False)
    parser.set_defaults(NORMAL=False)
    parser.set_defaults(UTEST=False)
    parser.set_defaults(TTEST=False)
    parser.set_defaults(FULL=False)
    parser.set_defaults(SHUFFLE=False)
