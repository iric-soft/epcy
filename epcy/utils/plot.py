import os

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42

col_pal = [
    mpl.colors.hex2color('#E62528'),
    mpl.colors.hex2color('#0D6BAC')
]


def plot_qc_histo(df_pred, quantiles, legend_quantile, mcc_bins, args):
    if not os.path.exists(args.PATH_OUT):
        os.makedirs(args.PATH_OUT)


    x_var = 'kernel_mcc'
    groupby_var = 'color_legend'
    df_pred_agg = df_pred.loc[:, [x_var, groupby_var]].groupby(groupby_var)
    vals = [df[x_var].values.tolist() for i, df in df_pred_agg]

    # Draw
    plt.figure(figsize=(16,9), dpi=150)
    colors = [plt.cm.copper(i/float(len(vals)-1)) for i in range(len(vals))]
    n, bins, patches = plt.hist(vals, mcc_bins, stacked=True, density=False, color=colors[:len(vals)])

    # Decoration
    plt.legend({group:col for group, col in zip(legend_quantile, colors[:len(vals)])})
    if args.L2FC:
        plt.title("QC histogram of kernel_mcc colored by abs_l2fc", fontsize=22)
    else:
        plt.title("QC histogram of kernel_mcc colored by max of mean expression", fontsize=22)
    plt.xlabel(x_var)
    plt.ylabel("# features")

    # add xlimit
    plt.axvline(x=0, color='r', linestyle='--')

    file_out = os.path.join(args.PATH_OUT, "qc_histogram_exp.pdf")
    if args.YLOG:
        plt.yscale("log")
        if args.L2FC:
            file_out = os.path.join(args.PATH_OUT, "qc_histogram_l2fc_ylog.pdf")
        else:
            file_out = os.path.join(args.PATH_OUT, "qc_histogram_exp_ylog.pdf")
    else:
        if args.L2FC:
            file_out = os.path.join(args.PATH_OUT, "qc_histogram_l2fc.pdf")

    plt.savefig(file_out)


def plot_profile(id, query_exp, ref_exp, bw, args):
    df_swarn = pd.DataFrame(
        data={
            'x': np.append(query_exp, ref_exp),
            'subgroup': np.append(
                np.repeat(args.QUERY, len(query_exp)),
                np.repeat("Other", len(ref_exp))
            )
        }
    )

    fig = plt.figure(figsize=(5, 5))
    gs = plt.GridSpec(4, 1)

    ax_kde = fig.add_subplot(gs[0, 0])
    ax_swarm = fig.add_subplot(gs[1:, 0], sharex=ax_kde)

    fig.subplots_adjust(hspace=0)
    sns.despine(ax=ax_kde, top=True, right=True, left=True, bottom=True)
    sns.despine(ax=ax_swarm, top=True, right=True, left=True, bottom=False)

    # Turn off kde axis visibility
    plt.setp(ax_kde.get_xticklabels(), visible=False)
    plt.setp(ax_kde.get_yticklabels(), visible=False)
    plt.setp(ax_kde.yaxis.get_majorticklines(), visible=False)
    plt.setp(ax_kde.yaxis.get_minorticklines(), visible=False)
    plt.setp(ax_kde.xaxis.get_majorticklines(), visible=False)
    plt.setp(ax_kde.xaxis.get_minorticklines(), visible=False)
    ax_kde.yaxis.grid(False)
    ax_kde.xaxis.grid(False)

    # Turn off swarm axis visibility
    plt.setp(ax_swarm.get_yticklabels(), visible=False)
    plt.setp(ax_swarm.yaxis.get_majorticklines(), visible=False)
    plt.setp(ax_swarm.yaxis.get_minorticklines(), visible=False)
    ax_swarm.yaxis.grid(False)

    sns_plot = sns.kdeplot(query_exp, shade=True, bw=bw, color = col_pal[0], label=args.QUERY, ax=ax_kde)
    #sns_plot = sns.rugplot(query_exp, color = "r")
    sns_plot = sns.kdeplot(ref_exp, shade=True, bw=bw, color = col_pal[1], label="Other", ax=ax_kde)
    #sns_plot = sns.rugplot(ref_exp, color = "b")
    sns_plot.set_title(str(id) + " " + args.QUERY + "\nbw=" + str(bw))

    sns_plot = sns.swarmplot(
        x="x", y="subgroup", data=df_swarn, ax=ax_swarm,
        palette=sns.color_palette([col_pal[0], col_pal[1]])
    )

    x_label = "x"
    if hasattr(args, 'CPM') and args.CPM:
        x_label = "cpm(x)"
    if hasattr(args, 'LOG') and args.LOG:
        x_label = 'log2(' + x_label + "+"+ str(args.C) + ')'

    sns_plot.set(xlabel=x_label)

    fig_dir = os.path.join(args.PATH_OUT)
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    fig_file = os.path.join(fig_dir, id + "_density.pdf")
    sns_plot.figure.savefig(fig_file)
    plt.close()
