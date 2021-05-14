import os

import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18

col_pal = [
    mpl.colors.hex2color('#D21417'),
    mpl.colors.hex2color('#1483D2')
]

def plot_explore_heatmap(df_heatmap, df_pred, args):
    if not os.path.exists(args.PATH_OUT):
        os.makedirs(args.PATH_OUT)

    df_heatmap = df_heatmap.set_index('ID')

    mcc_colors_palette = sns.light_palette("green", reverse=True, n_colors=3)


    mcc_colors = []
    for x in df_pred.KERNEL_MCC.values:
        if x >= 0.9:
            mcc_colors.append(mcc_colors_palette[0])
        elif x >=0.5:
            mcc_colors.append(mcc_colors_palette[1])
        else:
            mcc_colors.append(mcc_colors_palette[2])

    #DOESN'T WORKS
    #mcc = df_pred.KERNEL_MCC.values.copy()
    #mcc_colors = df_pred.KERNEL_MCC.values.copy()
    #mcc_colors[np.where(mcc >= 0.9)] = [mcc_colors_palette[0]] * len(mcc_colors[np.where(mcc >= 0.9)])
    #mcc_colors[np.where((mcc < 0.9) & (mcc >= 0.5))] = [mcc_colors_palette[1]] * len(mcc_colors[np.where((mcc < 0.9) & (mcc >= 0.5))])
    #mcc_colors[np.where(mcc < 0.5)] = [mcc_colors_palette[2]] * len(mcc_colors[np.where(mcc < 0.5)])


    sns_plot = sns.clustermap(
        df_heatmap, linewidths=0, metric='euclidean',
        xticklabels=True, yticklabels=True,
        vmin=-1, vmax=1, row_cluster=False,
        row_colors=mcc_colors,
        cmap="vlag"
    )
    sns_plot.fig.suptitle("Heatmap of predicted condition using top " + str(args.TOP) + " features")
    sns_plot.ax_heatmap.set_xticklabels(sns_plot.ax_heatmap.get_xmajorticklabels(), fontsize = 4)
    sns_plot.ax_heatmap.set_yticklabels(sns_plot.ax_heatmap.get_ymajorticklabels(), fontsize = 4)

    file_out = os.path.join(args.PATH_OUT, "explore_heatmap.pdf")
    sns_plot.fig.savefig(file_out)

def plot_explore_heatmap(df_heatmap, df_pred, args):
    if not os.path.exists(args.PATH_OUT):
        os.makedirs(args.PATH_OUT)

    df_heatmap = df_heatmap.set_index('ID')

    mcc_colors_palette = sns.light_palette("green", reverse=True, n_colors=3)

    mcc_colors = []
    for x in df_pred.KERNEL_MCC.values:
        if x >= 0.9:
            mcc_colors.append(mcc_colors_palette[0])
        elif x >= 0.5:
            mcc_colors.append(mcc_colors_palette[1])
        else:
            mcc_colors.append(mcc_colors_palette[2])

    sns_plot = sns.clustermap(
        df_heatmap, linewidths=0, metric='euclidean',
        xticklabels=True, yticklabels=True,
        vmin=-1, vmax=1, row_cluster=False,
        row_colors=mcc_colors,
        cmap="vlag"
    )
    sns_plot.fig.suptitle("Heatmap of predicted condition using top " +
                          str(args.TOP) + " features")
    sns_plot.ax_heatmap.set_xticklabels(
        sns_plot.ax_heatmap.get_xmajorticklabels(),
        fontsize=4)
    sns_plot.ax_heatmap.set_yticklabels(
        sns_plot.ax_heatmap.get_ymajorticklabels(),
        fontsize=4)

    file_out = os.path.join(args.PATH_OUT, "explore_heatmap.pdf")
    sns_plot.fig.savefig(file_out)


def plot_qc_histo(df_pred, quantiles, legend_quantile, mcc_bins, args):
    if not os.path.exists(args.PATH_OUT):
        os.makedirs(args.PATH_OUT)

    x_var = 'kernel_mcc'
    groupby_var = 'color_legend'

    if args.EXP or args.L2FC:
        df_pred_agg = df_pred.loc[:, [x_var, groupby_var]].groupby(groupby_var)
        vals = [df[x_var].values.tolist() for i, df in df_pred_agg]
    else:
        vals = df_pred[x_var].values.tolist()

    # Draw
    plt.figure(figsize=(16, 9), dpi=150)
    colors = [plt.cm.copper(i/float(len(vals)-1)) for i in range(len(vals))]

    if args.EXP or args.L2FC:
        n, bins, patches = plt.hist(vals, mcc_bins, stacked=True, density=False,
                                    color=colors[:len(vals)])
        # Decoration
        plt.legend({group: col
                    for group, col in zip(legend_quantile, colors[:len(vals)])})
    else:
        n, bins, patches = plt.hist(vals, mcc_bins, stacked=True, density=False)

    plt.title("QC histogram of kernel_mcc",
              fontsize=28)

    if args.L2FC:
        plt.title("QC histogram of kernel_mcc colored by abs_l2fc",
                  fontsize=28)
    if args.EXP:
        plt.title(
            "QC histogram of kernel_mcc colored by max of mean expression",
            fontsize=28
        )
    plt.xlabel(x_var, fontsize=20)
    plt.ylabel("# features", fontsize=20)

    # add xlimit
    plt.axvline(x=0, color='r', linestyle='--')

    file_out = os.path.join(args.PATH_OUT, "qc_mcc.pdf")

    if args.EXP:
        file_out = os.path.join(args.PATH_OUT, "qc_mcc_exp.pdf")

    if args.YLOG:
        plt.yscale("log")
        if args.L2FC:
            file_out = os.path.join(args.PATH_OUT,
                                    "qc_mcc_l2fc_ylog.pdf")
        else:
            file_out = os.path.join(args.PATH_OUT,
                                    "qc_mcc_exp_ylog.pdf")
    else:
        if args.L2FC:
            file_out = os.path.join(args.PATH_OUT, "qc_mcc_l2fc.pdf")

    plt.savefig(file_out)

    x_var = 'kernel_mcc'
    # Draw
    plt.figure(figsize=(16, 9), dpi=150)
    n, bins, patches = plt.hist([df_pred['bw_query'], df_pred['bw_ref']],
                                bins=100, color=['r', 'b'],
                                label=['Query', 'Ref'])
    plt.title("QC histogram of bandwidth", fontsize=28)
    plt.xlabel("bandwidth", fontsize=20)
    plt.ylabel("# features", fontsize=20)

    # add xlimit
    plt.axvline(x=args.MIN_BW, color='black', linestyle='--', label="min_bw")

    plt.legend()
    file_out = os.path.join(args.PATH_OUT, "qc_bw.pdf")
    plt.savefig(file_out)


def plot_profile(id, query_exp, ref_exp, bw_query, bw_ref, args):
    df_swarn = pd.DataFrame(
        data={
            'x': np.append(query_exp, ref_exp),
            'condition': np.append(
                np.repeat(args.QUERY, len(query_exp)),
                np.repeat("Other", len(ref_exp))
            )
        }
    )

    df_swarn.x = df_swarn.x + np.random.normal(0, 0.01, df_swarn.shape[0])
    # dummy plots, just to get the Path objects
    fig, ax = plt.subplots(1, 1)
    a = ax.scatter([1, 2], [3, 4], marker='s')
    b = ax.scatter([1, 2], [3, 4])
    square_mk, = a.get_paths()
    circle_mk, = b.get_paths()

    fig = plt.figure(figsize=(5, 5.9))
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

    if not args.NO_DENSITY:
        sns_plot = sns.kdeplot(query_exp, shade=True, bw=bw_query,
                               color=col_pal[0], label=args.QUERY,
                               ax=ax_kde)
        sns_plot = sns.kdeplot(ref_exp, shade=True, bw=bw_ref,
                               color=col_pal[1], label="Other",
                               ax=ax_kde)
        sns_plot.set_title(str(id) + " " + args.QUERY + "\nbw_query=" +
                           str(round(bw_query,2)) + ", bw_ref=" +
                           str(round(bw_ref,2)))

    if args.STRIP:
        sns_plot = sns.stripplot(
            x="x", y="condition", data=df_swarn, ax=ax_swarm,
            size=args.SIZE, jitter=0.4, hue="condition",
            palette=sns.color_palette([col_pal[0], col_pal[1]]),
            edgecolor="gray"
        )
    elif args.VIOLIN:
        bw = (bw_query + bw_ref) / 2
        sns_plot = sns.violinplot(
            x="x", y="condition", data=df_swarn, ax=ax_swarm,
            inner=None, linewidth=None, bw=bw,
            palette=sns.color_palette([col_pal[0], col_pal[1]])
        )
    else:
        sns_plot = sns.swarmplot(
            x="x", y="condition", data=df_swarn, ax=ax_swarm, size=args.SIZE,
            palette=sns.color_palette([col_pal[0], col_pal[1]])
        )
        ax_swarm.set_ylabel('')

    if args.NO_DENSITY:
        sns_plot.set_title(str(id) + " " + args.QUERY + "\nbw_query=" +
                           str(round(bw_query,2)) + ", bw_ref=" +
                           str(round(bw_ref,2)))

    # Change shape in function of condition
    if not args.VIOLIN:
        collections = sns_plot.collections
        unique_colors = [list(col_pal[0]) + [1], list(col_pal[1]) + [1]]
        markers = [circle_mk, square_mk]
        for collection in collections:
            paths = []
            for current_color in collection.get_facecolors():
                for possible_marker, possible_color in zip(markers, unique_colors):
                    if np.array_equal(current_color, possible_color):
                        paths.append(possible_marker)
                        break
            collection.set_paths(paths)
    # sns_plot.legend(collections[-2:], pd.unique(df_swarn.condition))

    x_label = "x"
    if hasattr(args, 'CPM') and args.CPM:
        x_label = "cpm(x)"
    if hasattr(args, 'LOG') and args.LOG:
        x_label = 'log2(' + x_label + "+" + str(args.C) + ')'

    sns_plot.set_xlabel(x_label, fontsize=18)

    fig_dir = os.path.join(args.PATH_OUT)
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    fig_file = os.path.join(fig_dir, id + "_density.pdf")
    sns_plot.figure.savefig(fig_file)
    plt.close('all')
