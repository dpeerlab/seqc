import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.colors import hex2color
from matplotlib import font_manager
from scipy.stats import gaussian_kde
from cycler import cycler
from mpl_toolkits.axes_grid1 import make_axes_locatable

# make matplotlib logger less verbose
import logging

logging.getLogger("matplotlib").setLevel(logging.WARNING)

try:
    os.environ["DISPLAY"]
except KeyError:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

with warnings.catch_warnings():
    warnings.simplefilter("ignore")  # catch warnings that system can't find fonts
    fm = font_manager.fontManager
    fm.findfont("Raleway")
    fm.findfont("Lato")

warnings.filterwarnings(action="ignore", module="matplotlib", message="^tight_layout")

dark_gray = ".15"

_colors = ["#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD"]

style_dictionary = {
    "figure.figsize": (3, 3),
    "figure.facecolor": "white",
    "figure.dpi": 200,
    "savefig.dpi": 200,
    "text.color": "k",
    "legend.frameon": False,
    "legend.numpoints": 1,
    "legend.scatterpoints": 1,
    "font.family": ["sans-serif"],
    "font.serif": ["Computer Modern Roman", "serif"],
    "font.monospace": ["Inconsolata", "Computer Modern Typewriter", "Monaco"],
    "font.sans-serif": ["Helvetica", "Lato", "sans-serif"],
    "patch.facecolor": _colors[0],
    "patch.edgecolor": "none",
    "grid.linestyle": "-",
    "axes.labelcolor": dark_gray,
    "axes.facecolor": "white",
    "axes.linewidth": 1.0,
    "axes.grid": False,
    "axes.axisbelow": False,
    "axes.edgecolor": dark_gray,
    "axes.prop_cycle": cycler("color", _colors),
    "lines.solid_capstyle": "round",
    "lines.color": _colors[0],
    "lines.markersize": 4,
    "image.cmap": "viridis",
    "image.interpolation": "none",
    "xtick.direction": "in",
    "xtick.major.size": 4,
    "xtick.minor.size": 2,
    "xtick.color": dark_gray,
    "ytick.direction": "in",
    "ytick.major.size": 4,
    "ytick.minor.size": 2,
    "ytick.color": dark_gray,
}

matplotlib.rcParams.update(style_dictionary)


def refresh_rc():
    matplotlib.rcParams.update(style_dictionary)
    print("rcParams updated")


class FigureGrid:
    """
    Generates a grid of axes for plotting

    axes can be iterated over or selected by number. e.g.:

    >>> # iterate over axes and plot some nonsense
    >>> fig = FigureGrid(4, max_cols=2)
    >>> for i, ax in enumerate(fig):
    >>>     plt.plot(np.arange(10) * i)

    >>> # select axis using indexing
    >>> ax3 = fig[3]
    >>> ax3.set_title("I'm axis 3")
    """

    def __init__(self, n: int, max_cols=3, scale=3):
        """
        :param n: number of axes to generate
        :param max_cols: maximum number of axes in a given row
        """

        self.n = n
        self.nrows = int(np.ceil(n / max_cols))
        self.ncols = int(min((max_cols, n)))
        figsize = self.ncols * scale, self.nrows * scale

        # create figure
        self.gs = plt.GridSpec(nrows=self.nrows, ncols=self.ncols)
        self.figure = plt.figure(figsize=figsize)

        # create axes
        self.axes = {}
        for i in range(n):
            row = int(i // self.ncols)
            col = int(i % self.ncols)
            self.axes[i] = plt.subplot(self.gs[row, col])

    def __getitem__(self, item):
        return self.axes[item]

    def __iter__(self):
        for i in range(self.n):
            yield self[i]

    def tight_layout(self, **kwargs):
        """wrapper for plt.tight_layout"""
        self.gs.tight_layout(self.figure, **kwargs)

    def despine(self, top=True, right=True, bottom=False, left=False):
        """removes axis spines (default=remove top and right)"""
        despine(ax=self, top=top, right=right, bottom=bottom, left=left)

    def detick(self, x=True, y=True):
        """
        removes tick labels

        :param x: bool, if True, remove tick labels from x-axis
        :param y: bool, if True, remove tick labels from y-axis
        """

        for ax in self:
            detick(ax, x=x, y=y)

    def savefig(self, filename, pad_inches=0.1, bbox_inches="tight", *args, **kwargs):
        """
        wrapper for savefig, including necessary paramters to avoid cut-off

        :param filename: str, name of output file
        :param pad_inches: float, number of inches to pad
        :param bbox_inches: str, method to use when considering bbox inches
        :param args: additional args for plt.savefig()
        :param kwargs: additional kwargs for plt.savefig()
        :return:
        """
        self.figure.savefig(
            filename, pad_inches=pad_inches, bbox_inches=bbox_inches, *args, **kwargs
        )


def detick(ax=None, x=True, y=True):
    """helper function for removing tick labels from an axis"""
    if not ax:
        ax = plt.gca()
    if x:
        ax.xaxis.set_major_locator(plt.NullLocator())
    if y:
        ax.yaxis.set_major_locator(plt.NullLocator())


def despine(ax=None, top=True, right=True, bottom=False, left=False) -> None:
    """helper function for removing axis spines"""
    if not ax:
        ax = plt.gca()

    # set spines
    if top:
        ax.spines["top"].set_visible(False)
    if right:
        ax.spines["right"].set_visible(False)
    if bottom:
        ax.spines["bottom"].set_visible(False)
    if left:
        ax.spines["left"].set_visible(False)

    # set ticks
    if top and bottom:
        ax.xaxis.set_ticks_position("none")
    elif top:
        ax.xaxis.set_ticks_position("bottom")
    elif bottom:
        ax.xaxis.set_ticks_position("top")
    if left and right:
        ax.yaxis.set_ticks_position("none")
    elif left:
        ax.yaxis.set_ticks_position("right")
    elif right:
        ax.yaxis.set_ticks_position("left")


def xtick_vertical(ax=None):
    """set xticklabels on ax to vertical instead of the horizontal default orientation"""
    if ax is None:
        ax = plt.gca()
    xt = ax.get_xticks()
    if np.all(xt.astype(int) == xt):  # ax.get_xticks() returns floats
        xt = xt.astype(int)
    ax.set_xticklabels(xt, rotation="vertical")


def equalize_numerical_tick_number(ax=None):
    if ax is None:
        ax = plt.gca()
    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    nticks = min(len(xticks), len(yticks))
    ax.set_xticks(np.round(np.linspace(min(xticks), max(xticks), nticks), 1))
    ax.set_yticks(np.round(np.linspace(min(yticks), max(yticks), nticks), 1))


def equalize_axis_size(ax=None):
    if ax is None:
        ax = plt.gca()
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax_min = min(xlim[0], ylim[0])
    ax_max = max(xlim[1], ylim[1])
    ax.set_xlim((ax_min, ax_max))
    ax.set_ylim((ax_min, ax_max))


def map_categorical_to_cmap(data: np.ndarray, cmap=plt.get_cmap()):
    """
    create a discrete colormap from cmap appropriate for data

    :param data: categorical vector to map to colormap
    :param cmap: cmap to discretize, or 'random'
    :return np.ndarray, dict: vector of colors matching input data, dictionary of labels
      to their respective colors
    """
    categories = np.unique(data)
    n = len(categories)
    if isinstance(cmap, str) and "random" in cmap:
        colors = np.random.rand(n, 3)
    else:
        colors = cmap(np.linspace(0, 1, n))
    category_to_color = dict(zip(categories, colors))
    return np.array([category_to_color[i] for i in data]), category_to_color


def add_legend_to_categorical_vector(
    colors: np.ndarray,
    labels,
    ax,
    loc="best",  # bbox_to_anchor=(0.98, 0.5),
    markerscale=0.75,
    **kwargs
):
    """
    Add a legend to a plot where the color scale was set by discretizing a colormap.

    :param colors: np.ndarray, output of map_categorical_vector_to_cmap()
    :param labels: np.ndarray, category labels
    :param ax: axis on which the legend should be plotted
    :param kwargs: additional kwargs for legend
    :return: None
    """
    artists = []
    for c in colors:
        artists.append(plt.Line2D((0, 1), (0, 0), color=c, marker="o", linestyle=""))
    ax.legend(
        artists,
        labels,
        loc=loc,
        markerscale=markerscale,  # bbox_to_anchor=bbox_to_anchor,
        **kwargs
    )


class scatter:
    @staticmethod
    def categorical(
        x,
        y,
        c,
        ax=None,
        cmap=plt.get_cmap(),
        legend=True,
        legend_kwargs=None,
        randomize=True,
        remove_ticks=False,
        *args,
        **kwargs
    ):
        """
        wrapper for scatter wherein the output should be colored by a categorical vector
        c

        :param x, y: np.ndarray, coordinate data to be scattered
        :param c: categories for data
        :param ax: axis on which to scatter data
        :param cmap: color map
        :param legend: bool, if True, plot legend
        :param legend_kwargs: additional kwargs for legend
        :param randomize: if True, randomize order of plotting
        :param remove_ticks: if True, removes axes ticks and labels
        :param args: additional args for scatter
        :param kwargs: additional kwargs for scatter
        :return: ax
        """
        if not ax:  # todo replace with plt.gridspec() method
            ax = plt.gca()

        if legend_kwargs is None:
            legend_kwargs = dict()

        color_vector, category_to_color = map_categorical_to_cmap(c, cmap)

        if randomize:
            ind = np.random.permutation(len(x))
        else:
            ind = np.argsort(np.ravel(c))

        ax.scatter(
            np.ravel(x)[ind], np.ravel(y)[ind], c=color_vector[ind], *args, **kwargs
        )
        if remove_ticks:
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())

        labels, colors = zip(*sorted(category_to_color.items()))
        if legend:
            add_legend_to_categorical_vector(
                colors, labels, ax, markerscale=2, **legend_kwargs
            )
        return ax

    @staticmethod
    def continuous(
        x,
        y,
        c=None,
        ax=None,
        colorbar=True,
        randomize=True,
        remove_ticks=False,
        **kwargs
    ):
        """
        wrapper for scatter wherein the coordinates x and y are colored according to a
        continuous vector c
        :param x, y: np.ndarray, coordinate data
        :param c: np.ndarray, continuous vector by which to color data points
        :param remove_ticks: remove axis ticks and labels
        :param args: additional args for scatter
        :param kwargs: additional kwargs for scatter
        :return: ax
        """

        if ax is None:
            ax = plt.gca()

        if c is None:  # plot density if no color vector is provided
            x, y, c = scatter.density_2d(x, y)

        if randomize:
            ind = np.random.permutation(len(x))
        else:
            ind = np.argsort(c)

        sm = ax.scatter(x[ind], y[ind], c=c[ind], **kwargs)
        if remove_ticks:
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())
        if colorbar:
            cb = plt.colorbar(sm)
            cb.ax.xaxis.set_major_locator(plt.NullLocator())
            cb.ax.yaxis.set_major_locator(plt.NullLocator())
        return ax

    @staticmethod
    def density_2d(x, y):
        """return x and y and their density z, sorted by their density (smallest to largest)

        :param x, y: np.ndarray: coordinate data
        :return: sorted x, y, and density
        """
        xy = np.vstack([np.ravel(x), np.ravel(y)])
        z = gaussian_kde(xy)(xy)
        return np.ravel(x), np.ravel(y), np.arcsinh(z)


def tatarize(n):
    """
    Return n-by-3 RGB color matrix using the "tatarize" color alphabet (n <= 269)
    :param n:
    :return:
    """

    with open(os.path.expanduser("~/.seqc/tools/tatarize_269.txt")) as f:
        s = f.read().split('","')
    s[0] = s[0].replace('{"', "")
    s[-1] = s[-1].replace('"}', "")
    s = [hex2color(s) for s in s]
    return s[:n]


class Diagnostics:
    @staticmethod
    def mitochondrial_fraction(data: pd.DataFrame, ax=None):
        """plot the fraction of mRNA that are of mitochondrial origin for each cell.

        :param data: DataFrame of cells x genes containing gene expression information
        :param ax: matplotlib axis
        :return: ax
        """

        mt_genes = data.molecules.columns[data.molecules.columns.str.contains("MT-")]
        mt_counts = data.molecules[mt_genes].sum(axis=1)
        library_size = data.molecules.sum(axis=1)

        if ax is None:
            ax = plt.gca()

        scatter.continuous(library_size, mt_counts / library_size)
        ax.set_title("Mitochondrial Fraction")
        ax.set_xlabel("Total Gene Expression")
        ax.set_ylabel("Mitochondrial Gene Expression")
        _, xmax = ax.get_xlim()
        ax.set_xlim((None, xmax))
        _, ymax = ax.get_ylim()
        ax.set_ylim((None, ymax))
        despine(ax)
        return ax

    @staticmethod
    def pca_components(fig_name, variance_ratio, pca_comps):
        """
        :param fig_name:    name for the figure
        :param variance_ratio:    variance ratios of at least 20 pca components
        :param pca_comps:    pca components of cells
        """

        fig = FigureGrid(4, max_cols=2)
        ax_pca, ax_pca12, ax_pca13, ax_pca23 = iter(fig)
        ax_pca.plot(variance_ratio[0:20] * 100.0, c="#1f77b4")
        ax_pca.set_xlabel("pca components")
        ax_pca.set_ylabel("explained variance")
        ax_pca.set_xlim([0, 20.5])

        ax_pca12.scatter(pca_comps[:, 0], pca_comps[:, 1], s=3, c="#1f77b4")
        ax_pca12.set_xlabel("pca 1")
        ax_pca12.set_ylabel("pca 2")
        xtick_vertical(ax=ax_pca12)

        ax_pca13.scatter(pca_comps[:, 0], pca_comps[:, 2], s=3, c="#1f77b4")
        ax_pca13.set_xlabel("pca 1")
        ax_pca13.set_ylabel("pca 3")
        xtick_vertical(ax=ax_pca13)

        ax_pca23.scatter(pca_comps[:, 1], pca_comps[:, 2], s=3, c="#1f77b4")
        ax_pca23.set_xlabel("pca 2")
        ax_pca23.set_ylabel("pca 3")
        xtick_vertical(ax=ax_pca23)

        fig.tight_layout()
        fig.savefig(fig_name, dpi=300, transparent=True)

    @staticmethod
    def phenograph_clustering(fig_name, cell_sizes, clust_info, tsne_comps):
        # sketching tSNE and Phenograph figure
        fig = FigureGrid(2, max_cols=2)
        ax_tsne, ax_phenograph = iter(fig)

        cl = np.log10(cell_sizes)
        splot = ax_tsne.scatter(
            tsne_comps[:, 0],
            tsne_comps[:, 1],
            c=cl,
            s=3,
            cmap=plt.cm.coolwarm,
            vmin=np.min(cl),
            vmax=np.percentile(cl, 98),
        )

        ax_tsne.set_title("UMI Counts (log10)")
        ax_tsne.set_xticks([])
        ax_tsne.set_yticks([])
        divider = make_axes_locatable(ax_tsne)
        cax = divider.append_axes("right", size="3%", pad=0.04)
        fig.figure.colorbar(splot, cax=cax, orientation="vertical")

        # this is a list of contrast colors for clutering
        cmap = [
            "#010067",
            "#D5FF00",
            "#FF0056",
            "#9E008E",
            "#0E4CA1",
            "#FFE502",
            "#005F39",
            "#00FF00",
            "#95003A",
            "#FF937E",
            "#A42400",
            "#001544",
            "#91D0CB",
            "#620E00",
            "#6B6882",
            "#0000FF",
            "#007DB5",
            "#6A826C",
            "#00AE7E",
            "#C28C9F",
            "#BE9970",
            "#008F9C",
            "#5FAD4E",
            "#FF0000",
            "#FF00F6",
            "#FF029D",
            "#683D3B",
            "#FF74A3",
            "#968AE8",
            "#98FF52",
            "#A75740",
            "#01FFFE",
            "#FFEEE8",
            "#FE8900",
            "#BDC6FF",
            "#01D0FF",
            "#BB8800",
            "#7544B1",
            "#A5FFD2",
            "#FFA6FE",
            "#774D00",
            "#7A4782",
            "#263400",
            "#004754",
            "#43002C",
            "#B500FF",
            "#FFB167",
            "#FFDB66",
            "#90FB92",
            "#7E2DD2",
            "#BDD393",
            "#E56FFE",
            "#DEFF74",
            "#00FF78",
            "#009BFF",
            "#006401",
            "#0076FF",
            "#85A900",
            "#00B917",
            "#788231",
            "#00FFC6",
            "#FF6E41",
            "#E85EBE",
        ]

        colors = []
        for i in range(len(clust_info)):
            colors.append(cmap[clust_info[i]])

        for ci in range(np.min(clust_info), np.max(clust_info) + 1):
            x1 = []
            y1 = []
            for i in range(len(clust_info)):
                if clust_info[i] == ci:
                    x1.append(tsne_comps[i, 0])
                    y1.append(tsne_comps[i, 1])
                    cl = colors[i]
            ax_phenograph.scatter(x1, y1, c=cl, s=3, label="C" + str(ci + 1))
        ax_phenograph.set_title("Phenograph Clustering")
        ax_phenograph.set_xticks([])
        ax_phenograph.set_yticks([])
        ax_phenograph.legend(
            bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.0, markerscale=2
        )

        fig.tight_layout()
        fig.savefig(fig_name, dpi=300, transparent=True)

    @staticmethod
    def cell_size_histogram(data, f=None, ax=None, save=None):
        if ax is None:
            f, ax = plt.subplots(figsize=(3.5, 3.5))
        if f is None:
            f = plt.gcf()

        cell_size = data.sum(axis=1)

        plt.hist(np.log10(cell_size), bins=25, log=True)
        ax.set_xlabel("log10(cell size)")
        ax.set_ylabel("frequency")
        despine(ax)
        xtick_vertical(ax)

        if save is not None:
            if not isinstance(save, str):
                raise TypeError(
                    "save must be the string filename of the " "figure-to-be-saved"
                )
            plt.tight_layout()
            f.savefig(save, dpi=300)
