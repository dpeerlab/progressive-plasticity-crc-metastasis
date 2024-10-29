import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from typing import List, Optional


def plot_overlaps(
    overlaps: pd.DataFrame,
    feature_colors: dict,
    labels: List[str] = None,
):

    # Plot overlaps
    fig, ax = plt.subplots(1, 1, figsize=(2, 2))
    sns.heatmap(
        overlaps,
        cmap='Purples',
        cbar_kws={'shrink': 0.2, 'aspect': 6},
        vmin=0,
        vmax=1,
        ax=ax,
    )
    n = overlaps.shape[0]
    ax.hlines(range(n+1), -1, n+0.1, color='w', lw=3, clip_on=False)

    # Adjust tick labels
    ax.tick_params(which='major', length=0, labelsize=6.5)
    ax.tick_params(axis='x', pad=7)
    ax.tick_params(axis='y', pad=8)
    if labels is None:
        labels = overlaps.columns.tolist()
    ax.set_yticklabels(labels, rotation=0)
    ax.set_xticks(np.arange(len(labels))+0.75)  # offset to align with center
    ax.set_xticklabels(labels, rotation=45, ha='right', va='top')

    # Add row and column color annotations
    for i, f in enumerate(overlaps):
        color = feature_colors[f]
        kwargs = dict(
            fill=True, facecolor=color, lw=1.5, edgecolor='w',
            clip_on=False, zorder=0
        )
        p = 0.075
        row_color = plt.Rectangle(
            (-p, i), p*0.75, 1, transform=ax.get_yaxis_transform(), **kwargs
        )
        ax.add_patch(row_color)
        p = 0.055
        col_color = plt.Rectangle(
            (i, -p), 1, p, transform=ax.get_xaxis_transform(), **kwargs
        )
        ax.add_patch(col_color)

    return fig


def plot_ratios(
    ratios: pd.Series,
    cmap,
    ax,
    row_colors: Optional[pd.Series] = None,
):
    # Cap infinity and negative infinity values
    mask = ratios.ne(0) & ratios.ne(np.inf)
    log_ratios = ratios.copy()
    log_ratios.loc[mask] = np.log(ratios.loc[mask])
    offset = log_ratios[mask].abs().max() * 0.05
    vmin = log_ratios[mask].min()
    vmax = log_ratios[mask].max()
    log_ratios.replace({0: vmin - offset, np.inf: vmax + offset}, inplace=True)

    # Setup colors
    colors = log_ratios.copy()
    colors[colors > 0] += vmax / 10
    colors[colors < 0] -= vmin / 10
    absmax = colors.abs().max()
    colors = (colors + absmax) / (absmax * 2)
    colors = colors.apply(cmap)

    # Plot log-ratios
    x = np.power(np.e, log_ratios) - 1
    n = log_ratios.shape[0]
    y = np.arange(n)
    left = np.repeat(1, n)
    ax.barh(y, x, height=0.75, lw=0, left=left, color=colors)

    # Formatting
    ax.spines[['top', 'left', 'right']].set_visible(False)
    ax.set_ylim(-1.25, n+0.25)
    ax.set_yticks(y, log_ratios.index)
    ax.tick_params(axis='y', length=0, labelsize=6, pad=6)
    ax.tick_params(axis='x', direction='in', labelsize=6)
    ax.set_xscale('log')

    # Add row and column color annotations
    for i, color in enumerate(row_colors.loc[log_ratios.index]):
        kwargs = dict(
            fill=True, facecolor=color, lw=1.5, edgecolor='w',
            clip_on=False, zorder=0
        )
        row_color = plt.Rectangle(
            (-0.08, i-0.5), 
            0.06, 1, transform=ax.get_yaxis_transform(), **kwargs
        )
        ax.add_patch(row_color)


def plot_fractions(
    fractions: pd.DataFrame,
    cmap: dict,
    ax,
    row_colors: pd.Series,
):
    # Plot cumulative fractions
    lefts = np.zeros(fractions.shape[0])
    yvals = np.arange(fractions.shape[0])
    for col in fractions.columns:
        ax.barh(
            y=yvals,
            width=fractions[col],
            height=0.75,
            left=lefts,
            color=cmap[col],
            lw=0.,
            edgecolor='w',
        )
        lefts += fractions[col]
        if 'Osteoblast' in col:
            styles = dict(lw=0.5, color='k', clip_on=False)
            for x, y in zip(lefts, yvals):
                ax.plot([x]*2, [y-0.5, y+0.5], **styles)

    # Formatting
    ax.spines[['top', 'left', 'right']].set_visible(False)
    ax.set_ylim(-1.25, fractions.shape[0] + 0.25)
    ax.set_yticks(yvals, fractions.index)
    ax.tick_params(axis='y', length=0, labelsize=6, pad=6)
    ax.tick_params(axis='x', direction='in', labelsize=6)

    # Add row and column color annotations
    for i, color in enumerate(row_colors.loc[fractions.index]):
        kwargs = dict(
            fill=True, facecolor=color, lw=1.5, edgecolor='w',
            clip_on=False, zorder=0
        )
        row_color = plt.Rectangle(
            (-0.08, i-0.5), 
            0.06, 1, transform=ax.get_yaxis_transform(), **kwargs
        )
        ax.add_patch(row_color)