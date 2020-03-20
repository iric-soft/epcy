
def get_argparser_other_pred_part(parser):

    parser.add_argument("--auc",
                        dest="AUC",
                        help="Compute sample assignation using normal dist.",
                        action='store_true')
    parser.add_argument("--normal",
                        dest="NORMAL",
                        help="Compute sample assignation using normal dist.",
                        action='store_true')
    parser.add_argument("--ttest",
                        dest="TTEST",
                        help="Compute a p-value using ttest_ind " +
                             "from scipy.stats.",
                        action='store_true')
    parser.add_argument("--utest",
                        dest="UTEST",
                        help="Compute a p-value using Mann-Whitney from " +
                             "scipy.stats. (NEED --auc)",
                        action='store_true')

    parser.set_defaults(AUC=False)
    parser.set_defaults(NORMAL=False)
    parser.set_defaults(UTEST=False)
    parser.set_defaults(TTEST=False)
