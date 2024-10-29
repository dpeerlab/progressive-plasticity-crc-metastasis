import scipy as sp
import numpy as np
import pandas as pd
import ternary
from ternary.helpers import (
    project_point,
    planar_to_coordinates,
    simplex_iterator
)
from matplotlib.colors import cnames, to_rgb
from matplotlib import collections, lines, pyplot as plt
import seaborn as sns
import colorsys
from typing import List


def ternary_kde(
    points,
    tax,
    n_levels=10,
    cmap="viridis",
    outline=False,
    bw_method='scott',
):
    # Project to 2D simplex
    simplex_points = np.apply_along_axis(project_point, 0, points)

    # Fit density model to 2D projection
    kde = sp.stats.gaussian_kde(simplex_points, bw_method=bw_method)

    # Evaluate density on triangular grid
    n = 100
    tri_grid = np.array(list(simplex_iterator(n))).T
    simplex_grid = np.apply_along_axis(project_point, 0, tri_grid)/n
    densities = kde(simplex_grid)
    levels = np.linspace(
        np.percentile(densities, 5), 
        np.percentile(densities, 95), 
        n_levels,
    ) 
    tax.ax.tricontourf(
        simplex_grid[0], simplex_grid[1], 
        densities, levels=levels, cmap=cmap,
        extend="both"
    )
    if outline:
        tax.ax.tricontour(
            simplex_grid[0], simplex_grid[1], 
            densities, levels=levels,
            colors=[[0,0,0,0.25]]+[[0,0,0,0.25]]*4,
            linewidths=[0.25]+[0.25]*4
        )


def format_tax(
    tax, 
    labels,
    fontsize,
    tick_width,
    boundary_width,
    pad,
):
    tax.gridlines(
        color="gray", lw=tick_width, linestyle='--', alpha=0.5, multiple=0.5
    )
    tax.ticks(
        axis='lbr', lw=tick_width, fontsize=fontsize, tick_formats='%.1f',
        offset=0.05, multiple=1.0,
    )
    tax.boundary(linewidth=boundary_width, zorder=4)
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')


def plot_ternary(
    points,
    cmap,
    titles,
    ax,
    n_pts: int = 1000,
):
    _, tax = ternary.figure(ax=ax)
    ternary_kde(
        points.T, 
        tax, 
        n_levels=9,
        cmap=cmap, 
        outline=True,
        bw_method=0.3,
    )
    idx = np.random.choice(points.shape[0], n_pts)
    tax.scatter(points[idx], s=0.5, lw=0, color="gray", alpha=0.33)
    format_tax(
        tax=tax, 
        labels=titles,
        fontsize=8,
        tick_width=0.5,
        boundary_width=1,
        pad=2,
    )


def lighten_color(color, amount: float):
    # Lookup color in matplotlib named colors
    try:
        color = cnames[color]
    except:
        pass
    h, l, s = colorsys.rgb_to_hls(*to_rgb(color))
    color = colorsys.hls_to_rgb(h, 1 - amount * (1 - l), s)
    return color


def patch_violinplot(ax):
    children = ax.get_children()
    i = 0
    for n in range(0, len(children), 4):
        art = children[n: n+4]
        is_violin = len(art) == 4
        is_violin &= isinstance(art[0], collections.PolyCollection)
        is_violin &= all([isinstance(a, lines.Line2D) for a in art[1:]])
        if is_violin:
            violin, q1, q2, q3 = art
            c = violin.get_facecolor()
            if i%2==1: c = lighten_color(c, 0.5)
            violin.set_facecolor(c)
            violin.set_edgecolor(c)
            violin.set_linewidth(0.1)
            q2.set_linestyle('solid')
            q2.set_linewidth(0.5)
            q2.set_solid_capstyle('butt')
            for q in [q1, q3]:
                q.set_alpha(0)
            i += 1


def plot_kde(
    data: pd.DataFrame,
    row: str,
    row_order: List,
    **kwargs,
):
    # Update with default styles
    styles = dict(
        cut=0,
        common_norm=False, 
        density_norm='width',
        width=0.75,
        gap=0,
        split=True,
        fill=True,
        linewidth=0.5,
        legend=False,
        inner='quart',
    )
    kwargs.update(styles)

    # Plot KDEs
    grid = sns.FacetGrid(
        data,
        row=row,
        row_order=row_order,
        height=1.5,
        aspect=1.5,
        sharex=False,
        gridspec_kws=dict(hspace=0.5),
        despine=False,
    )
    fig = grid.map_dataframe(sns.violinplot, **kwargs)

    # Formatting
    for ax in fig.axes.flat:
        patch_violinplot(ax)
        ax.set_title('')
        ax.set_xlim(0, 1)
        ax.yaxis.set_visible(False)
        ax.spines[['top', 'left', 'right']].set_visible(False)
        ax.set_xlabel(kwargs['x'])