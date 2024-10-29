import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from sklearn.preprocessing import MinMaxScaler
import sys
import os


# Class to suppress output from PhenoGraph
class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def quartile_to_level(data, quantile):
    """Return data levels corresponding to quantile cuts of mass."""
    isoprop = np.asarray(quantile)
    values = np.ravel(data)
    sorted_values = np.sort(values)[::-1]
    normalized_values = np.cumsum(sorted_values) / values.sum()
    idx = np.searchsorted(normalized_values, 1 - isoprop)
    levels = np.take(sorted_values, idx, mode="clip")
    return levels


def get_kde(
    data,
    grid_size=500,
    min_q=0.0,
    **kwargs
):
    kernel = gaussian_kde(data, **kwargs)
    positions = np.linspace(data.min(), data.max(), grid_size)
    estimate = kernel(positions)
    level = quartile_to_level(estimate, min_q)
    mask = estimate>=level
    return positions[mask], estimate[mask]-level


def plot_violin(
    data, x, y,
    palette,
    ax,
    bw_adjust=1,
    grid_size=100,
    h=2,
    norm=True,
    min_q=0.01,
    lw=0.5,
):
    n_groups = data[x].nunique()
    x_ticks = np.linspace(0, -h * (n_groups - 1), n_groups)

    for (name, group), xpos in zip(data.sort_values(x).groupby(x), x_ticks):
        
        # Violin Plot
        n = group.shape[0]
        bw_method = group.shape[0] ** (-1/(5)) * bw_adjust
        pos, est = get_kde(group[y], min_q=min_q, bw_method=bw_method)
        if norm:
            est = MinMaxScaler().fit_transform(est.reshape(-1,1)).flatten()
        ax.fill_betweenx(
            pos, xpos - est, xpos + est,
            facecolor=palette[name], alpha=0.75, lw=0,
            zorder=1,
        )
        ax.plot(xpos + est, pos, lw=lw, color="w", zorder=1)
        ax.plot(xpos - est, pos, lw=lw, color="w", zorder=1)
        # IQR Plot
        ymin, ymid, ymax = group[y].quantile([0.25, 0.5, 0.75])
        ax.plot([xpos, xpos], [ymin, ymax], lw=0.5, color="black", alpha=0.75)
        ax.scatter([xpos], [ymid], s=2, color="black", alpha=0.9)
    
    ax.set_xticks(x_ticks)


def plot_heatmap(
    data: pd.DataFrame,
    k: int = 35,
 ):
    pass
