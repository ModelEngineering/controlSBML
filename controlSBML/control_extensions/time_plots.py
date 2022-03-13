"""Plots for DataFrames where the index is time."""


import controlSBML.constants as cn
from controlSBML.option_management.option_manager import OptionManager

from docstring_expander.expander import Expander
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@Expander(cn.KWARGS, cn.PLOT_KWARGS)
def plotOneDF(df, **kwargs):
    """
    Plots a dataframe as multiple lines on a single plot. The
    columns are legends.

    Parameters
    ----------
    df: pd.DataFrame
    #@expand
    """
    mgr = OptionManager(kwargs)
    mgr.plot_opts.set(cn.O_XLABEL, default="time")
    plt.plot(df.index, df)
    legend_spec = cn.LegendSpec(df.columns, crd=mgr.plot_opts[cn.O_LEGEND_CRD])
    mgr.plot_opts.set(cn.O_LEGEND_SPEC, default=legend_spec)
    mgr.doPlotOpts()
    mgr.doFigOpts()

@Expander(cn.KWARGS, cn.PLOT_KWARGS)
def plotManyDF(*dfs, ncol=1, names=None, **kwargs):
    """
    Plots multiple dataframes with the same columns so that each column
    is a different plot.

    Parameters
    ----------
    dfs: sequence of pd.DataFrame
    ncol: int (number of columns)
    names: list-str
        Names of the dataframes
    #@expand
    """
    mgr = OptionManager(kwargs)
    mgr.plot_opts.set(cn.O_XLABEL, default="time")
    nrow = int(np.ceil(len(dfs)/ncol))
    fig, axes = plt.subplots(nrow, ncol)
    axes = np.reshape(axes, (nrow, ncol))
    columns = dfs[0].columns
    for idx, col in enumerate(columns):
        irow = int(np.floor(idx/ncol))
        icol = idx - irow*ncol
        ax = axes[irow, icol]
        mgr.plot_opts[cn.O_AX] = ax
        for df_idx, df in enumerate(dfs):
            try:
                values = df[col]
            except KeyError:
                raise ValueError("DataFrame %d lacks column %s"
                      % (df_idx, col))
            ax.plot(df.index, df[col])
        if names is not None:
            legend_spec = cn.LegendSpec(names,
                   crd=mgr.plot_opts[cn.O_LEGEND_CRD])
            mgr.plot_opts.set(cn.O_LEGEND_SPEC, default=legend_spec)
        mgr.doPlotOpts()
    mgr.doFigOpts()

