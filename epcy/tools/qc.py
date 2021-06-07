import sys

import pandas as pd
import numpy as np
import math

from .. utils import plot as up


def set_color_legend():
    df_pred['abs_l2fc'] = np.abs(df_pred['l2fc'])
    quantiles = [math.trunc(np.quantile(df_pred['abs_l2fc'], x) * 10000) / 10000
                 for x in [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99]]
    legend_quantile = [
        "abs(l2fc) <" + str(quantiles[0]),
        str(quantiles[0]) + "<= abs(l2fc) <=" + str(quantiles[1]),
        str(quantiles[1]) + "<= abs(l2fc) <" + str(quantiles[2]),
        str(quantiles[2]) + "<= abs(l2fc) <" + str(quantiles[3]),
        str(quantiles[3]) + "<= abs(l2fc) <" + str(quantiles[4]),
        str(quantiles[4]) + "<= abs(l2fc) <" + str(quantiles[5]),
        str(quantiles[6]) + "<= abs(l2fc)"
    ]

    df_pred = df_pred.sort_values(['abs_l2fc'], ascending=True)
    df_pred['color_legend'] = legend_quantile[0]
    df_pred.loc[(df_pred['abs_l2fc'] >= quantiles[0]) &
                (df_pred['abs_l2fc'] < quantiles[1]), 'color_legend'] = legend_quantile[1]
    df_pred.loc[(df_pred['abs_l2fc'] >= quantiles[1]) &
                (df_pred['abs_l2fc'] < quantiles[2]), 'color_legend'] = legend_quantile[2]
    df_pred.loc[(df_pred['abs_l2fc'] >= quantiles[2]) &
                (df_pred['abs_l2fc'] < quantiles[3]), 'color_legend'] = legend_quantile[3]
    df_pred.loc[(df_pred['abs_l2fc'] >= quantiles[3]) &
                (df_pred['abs_l2fc'] < quantiles[4]), 'color_legend'] = legend_quantile[4]
    df_pred.loc[(df_pred['abs_l2fc'] >= quantiles[4]) &
                (df_pred['abs_l2fc'] < quantiles[5]), 'color_legend'] = legend_quantile[5]
    df_pred.loc[(df_pred['abs_l2fc'] >= quantiles[5]), 'color_legend'] = legend_quantile[6]


def main_qc(args, argparser):
    #Check param
    if args.EXP and args.L2FC:
        sys.stderr.write("We can't fill histogram on EXP and L2FC.\n" +
                         "You need choose one of these options.\n")
        exit()


    # Import Data
    df_pred = pd.read_csv(args.PRED, sep="\t")
    df_pred = df_pred.dropna()

    # Prepare data
    min_mcc = df_pred['kernel_mcc'].min()
    start_bin = math.trunc((round(min_mcc, 1) - 0.05) * 100)
    mcc_bins = [x/100 for x in range(start_bin, 105, 5)]

    if (args.L2FC):
        df_pred['abs_l2fc'] = np.abs(df_pred['l2fc'])
        quantiles = [math.trunc(np.quantile(df_pred['abs_l2fc'], x) * 10000) / 10000
                     for x in [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99]]
        legend_quantile = [
            "abs(l2fc) <" + str(quantiles[0]),
            str(quantiles[0]) + "<= abs(l2fc) <" + str(quantiles[1]),
            str(quantiles[1]) + "<= abs(l2fc) <" + str(quantiles[2]),
            str(quantiles[2]) + "<= abs(l2fc) <" + str(quantiles[3]),
            str(quantiles[3]) + "<= abs(l2fc) <" + str(quantiles[4]),
            str(quantiles[4]) + "<= abs(l2fc) <" + str(quantiles[5]),
            str(quantiles[6]) + "<= abs(l2fc)"
        ]

        df_pred = df_pred.sort_values(['abs_l2fc'], ascending=True)
        df_pred['color_legend'] = legend_quantile[0]
        df_pred.loc[(df_pred['abs_l2fc'] >= quantiles[0]) &
                    (df_pred['abs_l2fc']< quantiles[1]), 'color_legend'] = legend_quantile[1]
        df_pred.loc[(df_pred['abs_l2fc'] >= quantiles[1]) &
                    (df_pred['abs_l2fc'] < quantiles[2]), 'color_legend'] = legend_quantile[2]
        df_pred.loc[(df_pred['abs_l2fc'] >= quantiles[2]) &
                    (df_pred['abs_l2fc'] < quantiles[3]), 'color_legend'] = legend_quantile[3]
        df_pred.loc[(df_pred['abs_l2fc'] >= quantiles[3]) &
                    (df_pred['abs_l2fc'] < quantiles[4]), 'color_legend'] = legend_quantile[4]
        df_pred.loc[(df_pred['abs_l2fc'] >= quantiles[4]) &
                    (df_pred['abs_l2fc'] < quantiles[5]), 'color_legend'] = legend_quantile[5]
        df_pred.loc[(df_pred['abs_l2fc'] >= quantiles[5]), 'color_legend'] = legend_quantile[6]
    else:
        if 'mean_log2_query' in df_pred:
            df_pred['max(query, ref)'] = df_pred[
                ['mean_log2_query', 'mean_log2_ref']
            ].max(axis=1)
        else:
            df_pred['max(query, ref)'] = df_pred[
                ['mean_query', 'mean_ref']
            ].max(axis=1)

        quantiles = [math.trunc(np.quantile(df_pred['max(query, ref)'], x) * 10000) / 10000
                     for x in [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99]]
        legend_quantile = [
            "max(mean_query, mean_ref) <" + str(quantiles[0]),
            str(quantiles[0]) + "<= max(mean_query, mean_ref) <" + str(quantiles[1]),
            str(quantiles[1]) + "<= max(mean_query, mean_ref) <" + str(quantiles[2]),
            str(quantiles[2]) + "<= max(mean_query, mean_ref) <" + str(quantiles[3]),
            str(quantiles[3]) + "<= max(mean_query, mean_ref) <" + str(quantiles[4]),
            str(quantiles[4]) + "<= max(mean_query, mean_ref) <" + str(quantiles[5]),
            str(quantiles[6]) + "<= max(mean_query, mean_ref)"
        ]

        df_pred = df_pred.sort_values(['max(query, ref)'], ascending=True)
        df_pred['color_legend'] = legend_quantile[0]
        df_pred.loc[(df_pred['max(query, ref)'] >= quantiles[0]) &
                    (df_pred['max(query, ref)'] < quantiles[1]), 'color_legend'] = legend_quantile[1]
        df_pred.loc[(df_pred['max(query, ref)'] >= quantiles[1]) &
                    (df_pred['max(query, ref)'] < quantiles[2]), 'color_legend'] = legend_quantile[2]
        df_pred.loc[(df_pred['max(query, ref)'] >= quantiles[2]) &
                    (df_pred['max(query, ref)'] < quantiles[3]), 'color_legend'] = legend_quantile[3]
        df_pred.loc[(df_pred['max(query, ref)'] >= quantiles[3]) &
                    (df_pred['max(query, ref)'] < quantiles[4]), 'color_legend'] = legend_quantile[4]
        df_pred.loc[(df_pred['max(query, ref)'] >= quantiles[4]) &
                    (df_pred['max(query, ref)'] < quantiles[5]), 'color_legend'] = legend_quantile[5]
        df_pred.loc[(df_pred['max(query, ref)'] >= quantiles[5]), 'color_legend'] = legend_quantile[6]

    up.plot_qc_histo(df_pred, quantiles, legend_quantile, mcc_bins, args)
