from pygam import LinearGAM, s as spline_term
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap as lsc
from typing import List, Tuple, Union, Optional
import scipy as sp


def get_gam_trend(
    x: np.ndarray,
    y: np.ndarray,
    n_splines: int = 8, #8
    spline_order: int = 4, #3
    x_res: int = 500,
    weights: Optional[np.ndarray] = None,
):
    x_pred = np.linspace(x.min(), x.max(), x_res)
    spline = spline_term(0, n_splines=n_splines, spline_order=spline_order)
    gam = LinearGAM(spline).fit(x, y, weights)
    y_pred = gam.predict(x_pred)
    p = gam.predict(x)
    n = len(x)
    mu = np.mean(x)
    sigma = np.sqrt(((y - p) ** 2).sum() / (n - 2))
    std = (
        np.sqrt(
            1 + 1 / n + (x_pred - mu) ** 2 / ((x - mu) ** 2).sum()
        ) * sigma / 2
    )
    return pd.Series(y_pred, index=x_pred), pd.Series(std, index=x_pred)


def plot_palantir_trends(
    adata: sc.AnnData,
    features: List[str],
    ps_column: str,
    branch_column: str,
    fig: plt.Figure,
    feature_colors: Optional[dict] = None,
    gs: Optional[gridspec.GridSpec] = None,
    **kwargs,
):
    # Colors for each trend/row
    if feature_colors is None:
        feature_colors = dict()
        cycler = plt.rcParams['axes.prop_cycle']
        for f, props in zip(features, cycler):
            feature_colors[f] = props['color']

    # Plot organization
    if gs is None:
        gs = fig.add_gridspec(1, 1)
    gs_0 = gs.subgridspec(2, 1, hspace=0.5)
    ax_1 = fig.add_subplot(gs_0[0])
    ax_2 = fig.add_subplot(gs_0[1], sharex=ax_1)

    # Calculate and plot feature trends for this branch
    ps_max = adata.obs.loc[adata.obs[branch_column] > 0.9, ps_column].max()
    mask = adata.obs[ps_column].lt(ps_max)
    x = adata.obs.loc[mask, ps_column]
    weights = adata.obs.loc[mask, branch_column]
    # Strengthen weights to account for instability near w=0
    eps = 1e-5
    weights = weights.clip(eps, 1-eps) ** 2
    for f in features:
        if f in adata.obs:
            y = adata.obs.loc[mask, f]
        else:
            y = adata[mask, f].layers['palantir_imputed'].flatten()
        # Plot trend
        trend, std = get_gam_trend(x, y, weights=weights, **kwargs)
        trend = sp.stats.zscore(trend)
        color = feature_colors[f]
        ax_1.plot(trend, color=color)
        # Plot first derivative
        diffs = pd.Series(np.diff(trend.values, n=1), trend.index[:-1])
        diffs = sp.stats.zscore(diffs)
        ax_2.plot(diffs, color=color)
        # Mark maxima
        signs = np.sign(np.diff(trend.values, n=2))
        try:
            max_idx = np.argwhere(signs[:-1] > signs[1:])[0]
        except:
            max_idx = np.argmax(diffs)
        styles = dict(lw=0.5, linestyle='--', color=color)
        for ax in [ax_1, ax_2]:
            t = ax.get_xaxis_transform()
            ax.vlines(trend.index[max_idx], 0, 1, transform=t, **styles)
    
    # Formatting
    ax_1.xaxis.set_visible(False)
    for ax in [ax_1, ax_2]:
        ax.tick_params(direction='in', length=2, labelsize=5, pad=2)


def plot_pseudotime_trends(
    adata: sc.AnnData,
    ps_column: str,
    feature_columns: Union[str, List[str]],
    ax,
    feature_colors: Optional[dict] = None,
    clip: Tuple[float] = (0.2, 1.0),
    labels: Optional[List[str]] = None,
    **kwargs,
):
    # Colors for each trend/row
    if feature_colors is None:
        feature_colors = dict()
        cycler = plt.rcParams['axes.prop_cycle']
        for f, props in zip(feature_columns, cycler):
            feature_colors[f] = props['color']

    # Calculate feature trends and reformat them as rows in a heatmap
    # 2) Find positions where trends increase (maxima in 2nd derivative)
    heatmap_trends = []
    tick_positions = []
    x = adata.obs[ps_column]
    for f in feature_columns:
        y = adata.obs[f]
        trend, _ = get_gam_trend(x, y, **kwargs)
        norm = plt.Normalize(*trend.quantile(clip))
        trend_cmap = lsc.from_list("", [[1,1,1], feature_colors[f]])
        row = trend.apply(norm).apply(trend_cmap).values.tolist()
        heatmap_trends.append(row)

    # Plotting
    ax.imshow(heatmap_trends, aspect="auto", interpolation='none',)
    borders = np.arange(-0.5, len(feature_columns)+0.5, 1)
    for y_pos, x_pos in enumerate(tick_positions):
        ax.plot([x_pos]*2, [y_pos-0.5, y_pos+0.5], c='k')

    # Formatting
    ax.hlines(borders, -1, 501, color="w", lw=1, zorder=2, clip_on=False)
    ax.set_xticks([0, 500], [0., 1.0])
    ax.set_xlim(-5, 505)
    ax.set_ylim(len(feature_columns)-0.5+0.2, -0.5-0.2)
    # Hide y-axis
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)

    # Row colors
    for i, f in enumerate(feature_columns):
        color = feature_colors[f]
        kwargs = dict(
            fill=True, facecolor=color, lw=1.5, edgecolor='w',
            clip_on=False, zorder=0
        )
        row_color = plt.Rectangle(
            (-0.05, i-0.5), 0.04, 1,
            transform=ax.get_yaxis_transform(), **kwargs
        )
        ax.add_patch(row_color)


def plot_module_progressions(
    adata: sc.AnnData,
    named_colors: dict,
    features: List = None,
    peaks: List = None,
    figsize: Tuple = (2, 1.2),
):
    # Setup figure
    fig, axes = plt.subplots(
        2,1,
        figsize=figsize,
        gridspec_kw=dict(height_ratios=[1, 0.33], hspace=0.33/(1.33/2)),
    )

    # Standard columns to plot
    if features is None:
        features = [
            'Module Absorptive Intestine Score',
            'Module Secretory Intestine Score',
            'Module Intestine Score',
            'Module Tumor ISC-like Score',
            'Module Endoderm Development Score',
            f"Module {adata.uns['DC Terminal State']} Score",
            'Fetal, Conserved'
        ]

    # Plot DC trends for all modules
    ax = axes[0]
    ps_column = adata.uns['DC']
    plot_pseudotime_trends(adata, ps_column, features, ax, named_colors)

    # Mark positions in fetal and terminal states where trend crosses 0.8
    if peaks is None: peaks = features[-2:]
    for feature in peaks:
        # Calculate trends
        trend, _ = get_gam_trend(adata.obs[ps_column], adata.obs[feature])
        trend -= trend.min()
        trend /= trend.max()
        # Mark positions
        signs = np.sign(trend.values - 0.75)
        idx = np.argwhere(signs[:-1] < signs[1:]).flatten()
        t = ax.get_xaxis_transform()
        styles = dict(lw=0.5, linestyle='--', color=named_colors[feature])
        ax.vlines(idx, 0, 1, transform=t, **styles)

    # Plot sample type positions
    ax = axes[1]
    x = adata.obs[ps_column]
    y = adata.obs['Sample Type'].str.contains('Primary')
    c = adata.obs['Sample Type'].map(named_colors).tolist()
    styles = dict(s=8, lw=0, alpha=0.25, rasterized=True, clip_on=False)
    ax.scatter(x, y, c=c, **styles)

    # Formatting
    offset = (x.max() - x.min()) / 100
    ax.set_xlim(x.min() - offset, x.max() + offset)
    ax.set_xticks([x.min(), x.max()], [0., 1.])
    ax.set_ylim(-0.67, 1.67)
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)

    return fig
