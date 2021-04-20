
def get_argparser_rna_norm_part(parser):

    parser.add_argument(
        "--cpm",
        dest="CPM",
        help="To apply a Count Par Million (CPM) " +
             "normalization to the matrix given in input (-m)",
        action='store_true'
    )

    parser.add_argument(
        "--cpmed",
        dest="CPMED",
        help="To apply a Count Par Median " +
             "normalization to the matrix given in input (-m)",
        action='store_true'
    )

    parser.add_argument(
        "--tpm",
        dest="TPM",
        help="Compute TPM from readcounts. (Need --anno)",
        action='store_true'
    )

    parser.set_defaults(TPM=False)
    parser.set_defaults(CPM=False)
    parser.set_defaults(CPMED=False)
