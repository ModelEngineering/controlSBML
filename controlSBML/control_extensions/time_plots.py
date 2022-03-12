"""Plots for DataFrames where the index is time."""


import controlSBML.constants as cn
from controlSBML.option_management.option_manager import OptionManager

from docstring_expander.expander import Expander
import matplotlib.pyplot as plt
import pandas as pd


@Expander(cn.KWARGS, cn.PLOT_KWARGS)
def plot1DF(df, **kwargs):
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

