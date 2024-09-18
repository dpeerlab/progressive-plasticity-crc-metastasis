import scanpy as sc
from collections.abc import Iterable
from typing import Union, List
import itertools
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import pyplot as plt
import math
import os

# Set style sheet on import
cwd = os.path.dirname(os.path.realpath(__file__))
plt.style.use(f'{cwd}/assets/default.mplstyle')

def format_ax(
    fig, ax,
    style="umap",
    title="",
    cbar=True,
    dim_label="UMAP",
    fs=12,
    lw=1.5,
    arrow_len=0.2,
    draw_arrows=True,
):
    ax.set_facecolor('white')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.grid(False)
    ax.spines[list(ax.spines)].set_visible(False)

    if style == "umap":
        change_aspect(ax)
        if draw_arrows:
            arrowed_spines(ax, arrow_len, text=dim_label, fs=fs, lw=lw)
        ax.set_title(title, weight="bold")
    if cbar:
        format_cbar(fig, ax)


def format_cbar(fig, ax):

    cbar = ax.get_children()[0].colorbar
    if cbar: 
        cbar.remove()
    data = ax.get_children()[0]

    # Create colorbar ax
    bbox = ax.get_position()
    cax = fig.add_axes([
        bbox.x1+bbox.width*0.025, #min x
        bbox.y0+bbox.height*0.25, #min y
        bbox.width*0.03, #width
        bbox.height*0.5 #height
    ])
    cax.grid(False)
    new_cbar = fig.colorbar(
        data, ax=ax, cax=cax,
    )
    new_cbar.outline.set_visible(False)
    if not cbar:
        bbox = ax.get_position()
        ax.get_children()[0].colorbar.remove()
        ax.set_position(bbox)


def change_aspect(ax):

    # Reset x and y limits for square plotting
    xmin, xmax = ax.get_xlim()
    xrange = xmax - xmin
    xcenter = (xrange/2) + xmin

    ymin, ymax = ax.get_ylim()
    yrange = ymax - ymin
    ycenter = (yrange/2) + ymin

    axrange = max(xrange, yrange)/2

    xmin = xcenter - (axrange)
    xmax = xcenter + (axrange)
    ax.set_xlim(xmin, xmax)

    ymin = ycenter - (axrange)
    ymax = ycenter + (axrange)
    ax.set_ylim(ymin, ymax)

    ax.set_aspect('equal', adjustable = 'box')


def arrowed_spines(
    ax,
    length = 0.2,
    text = None,
    fs = None,
    lw = 1.5,
):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    hw = 1./30.*(ymax-ymin)
    hl = 1./30.*(xmax-xmin)
    lw = lw # axis line width
    ohg = 0.0 # arrow overhang

    ax.spines[list(ax.spines)].set_visible(False)
    ax.arrow(
        xmin, ymin, (xmax-xmin)*length, 0, fc='k', ec='k', lw = lw, 
        head_width=hw, head_length=hl, overhang = ohg, 
        length_includes_head= True, clip_on = False
    )
    ax.arrow(
        xmin, ymin, 0, (ymax-ymin)*length, fc='k', ec='k', lw = lw, 
        head_width=hw, head_length=hl, overhang = ohg, 
        length_includes_head= True, clip_on = False
    )
    if fs == None:
        fs = plt.rcParams["xtick.labelsize"]
    ax.text(
        s=f"{text}1", 
        y=ymin-(ymax-ymin)*0.05, x=xmin+(xmax-xmin)*length/2,
        ha="center", va="top",
        fontsize = fs
    )
    ax.text(
        s=f"{text}2", 
        x=xmin-(xmax-xmin)*0.05, y=ymin+(ymax-ymin)*length/2,
        ha="right", va="center", rotation=90,
        fontsize = fs
    )
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)


def plot_embedding(
    adata: sc.AnnData,
    features: Union[str, List[str]],
    basis: str = 'X_umap',
    palette: str = "tab20",
    cmap: str = "plasma",
    titles: Union[str, List[str]] = None,
    ncols: int = 5,
    dim: int = 5,
    layer: str = "imputed",
    dim_label = "UMAP",
    ax = None,
    fs: int = 12,
    lw: float = 1.5,
    arrow_len: float = 0.2,
    draw_arrows=False,
    rasterized=False,
    **kwargs,
):
    iterify = lambda x: x if isinstance(x, Iterable) and not isinstance(x, str) else [x]
    features = iterify(features)
    titles = iterify(titles)
    if not ax:
        nrows = math.ceil(len(features))
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(dim*ncols,dim*nrows),
        )
        fig.tight_layout(pad=dim*0.75)
        axes = axes.flat if isinstance(axes, Iterable) else [axes]
    else:
        assert (len(features)==1) and (len(titles)==1)
        fig = ax.get_figure()
        axes = [ax]
    for ax, feature, title in itertools.zip_longest(axes, features, titles):
        if not title: title = feature
        if feature:
            sc.pl.embedding(
                adata,
                basis=basis,
                color=feature,
                ax=ax,
                show=False,
                palette=palette, cmap=cmap,
                layer=layer,
                **kwargs,
            )
            if rasterized:
                ax.get_children()[0].set_rasterized(True)
            format_ax(
                fig, ax, style="umap",
                title=title, dim_label=dim_label, fs=fs,
                arrow_len=arrow_len, lw=lw, draw_arrows=draw_arrows,
            )
        else:
            ax.set_visible(False)
    return fig


def lighten_color(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def saturate(c, s=1.0):
    from matplotlib import colors
    from collections.abc import Iterable
    # Assumed to be hex if string
    if isinstance(c, str):
        rgb = colors.to_rgb(c)
        hex = True
    # Assumed to be RGB if iterable
    elif isinstance(c, Iterable):
        rgb = c
        hex = False
    hsv = colors.rgb_to_hsv(rgb)
    hsv[1] *= s
    c = colors.hsv_to_rgb(hsv)
    if hex:
        return colors.to_hex(c)
    else:
        return c