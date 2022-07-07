import numpy as np
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
):
    """
    :param metrics_df: dataframe with columns for methods, metrics and metric values
    :param metric_column: column in ``metrics_df`` of metrics
    :param method_column: column in ``metrics_df`` of methods
    :param batch_metrics: list of batch correction metrics in metrics column for annotating metric type
    :param bio_metrics: list of biological conservation metrics in the metrics column for annotating metric type
    :param palette: color map as input for ``seaborn.scatterplot``
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
    df["rank"] = df.groupby(metric_column)[value_column].rank(ascending=False)

    dims = (
        df[["metric_type", metric_column]]
        .drop_duplicates()["metric_type"]
        .value_counts()
    )
    n_metric_types = dims.shape[0]
    n_metrics = dims.sum()
    n_methods = df[method_column].nunique()
    dim_x = np.max([4, (n_metrics + n_metric_types) * 0.4])
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
            legend="brief",
            ax=ax,
        )
        ax.set(title=metric_type, xlabel=None, ylabel=None)
        ax.tick_params(axis="x", rotation=90)
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
        for t in ax.legend_.texts:
            t.set_text(t.get_text()[:5])
        sns.despine(bottom=True, left=True)

    fig.tight_layout()
