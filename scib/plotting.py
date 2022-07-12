import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def metrics(
    metrics_df,
    method_column="method",
    metric_column="metric",
    value_column="value",
    batch_metrics=None,
    bio_metrics=None,
    palette=None,
    overall=True,
    return_fig=False,
):
    """
    :param metrics_df: dataframe with columns for methods, metrics and metric values
    :param method_column: column in ``metrics_df`` of methods
    :param metric_column: column in ``metrics_df`` of metrics
    :param value_column: column in ``metrics_df`` with metric values
    :param batch_metrics: list of batch correction metrics in metrics column for annotating metric type
    :param bio_metrics: list of biological conservation metrics in the metrics column for annotating metric type
    :param palette: color map as input for ``seaborn.scatterplot``
    :param overall: whether to include a column for the overall score
    :param return_fig: whether to return a fig object
    """
    sns.set_context("paper")
    sns.set_style("white")

    if palette is None:
        palette = "viridis_r"
        # sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True)

    if batch_metrics is None:
        batch_metrics = ["ASW_batch", "PCR_batch", "graph_conn", "kBET", "iLISI"]

    if bio_metrics is None:
        bio_metrics = [
            "NMI_cluster",
            "ARI_cluster",
            "ASW_label",
            "cell_cycle_conservation",
            "isolated_label_F1",
            "isolated_label_silhouette",
            "cLISI",
            "hvg_overlap",
            "trajectory",
        ]

    df = metrics_df.copy()

    conditions = [
        (df[metric_column].isin(batch_metrics)),
        (df[metric_column].isin(bio_metrics)),
    ]
    metric_type = ["Batch Correction", "Biological Conservation"]
    df["metric_type"] = np.select(conditions, metric_type)
    df[metric_column] = df[metric_column].str.replace("_", " ")

    # overall score
    df_list = [
        df,
        df.groupby([method_column, "metric_type"])[value_column]
        .mean()
        .reset_index()
        .assign(metric="Overall"),
        df.groupby(method_column)[value_column]
        .mean()
        .reset_index()
        .assign(metric_type="Overall", metric="Overall"),
    ]
    df = pd.concat(df_list)

    # rank metrics
    df["rank"] = (
        df.groupby([metric_column, "metric_type"])[value_column]
        .rank(
            method="min",
            ascending=False,
            na_option="bottom",
        )
        .astype(int)
    )
    method_rank = df.query('metric_type == "Overall"').sort_values(
        "rank", ascending=True
    )[method_column]
    df[method_column] = pd.Categorical(df[method_column], categories=method_rank)

    # get plot dimensions
    dims = (
        df[["metric_type", metric_column]]
        .drop_duplicates()["metric_type"]
        .value_counts()
    )
    n_metric_types = dims.shape[0]
    n_metrics = dims.sum()
    n_methods = df[method_column].nunique()
    metric_len = df[metric_column].str.len().max()
    dim_x = np.max([4, (n_metrics + n_metric_types) * 0.4 + (metric_len / 10)])
    dim_y = np.max([2.5, n_methods * 0.9])

    # Build plot
    fig, axs = plt.subplots(
        nrows=1,
        ncols=n_metric_types,
        figsize=(dim_x, dim_y),
        sharey=True,
        gridspec_kw=dict(width_ratios=list(dims)),
    )

    for i, metric_type in enumerate(dims.index):
        legend = "brief" if i == 0 else None
        df_sub = df.query(f'metric_type == "{metric_type}"')
        ax = axs if n_metric_types == 1 else axs[i]
        sns.scatterplot(
            data=df_sub,
            x=metric_column,
            y=method_column,
            hue="rank",
            palette=palette,
            size=value_column,
            sizes=(df_sub["value"].min() * 100, df_sub["value"].max() * 100),
            # sizes={x: int(x * 200) for x in df_sub['value'].dropna().unique()},
            edgecolor="black",
            legend=legend,
            ax=ax,
        )
        ax.set(title=metric_type, xlabel=None, ylabel=None)
        ax.tick_params(axis="x", rotation=90)
        if legend is not None:
            ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
            for t in ax.legend_.texts:
                t.set_text(t.get_text()[:5])
        sns.despine(bottom=True, left=True)

    fig.tight_layout()

    if return_fig:
        return fig
