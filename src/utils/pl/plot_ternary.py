import scipy as sp
import numpy as np
import ternary
from ternary.helpers import (
    project_point,
    planar_to_coordinates,
    simplex_iterator
)


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