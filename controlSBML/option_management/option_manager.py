"""Manages options used in ControlSBML."""

import controlSBML.constants as cn
from controlSBML.option_management.options import Options

from docstring_expander.expander import Expander
import matplotlib.pyplot as plt


class OptionManager(object):

    def __init__(self, kwargs, default_dcts):
        """
        Parameters
        ----------
        kwargs: dict
            key word arguments
        default_dcts: list-dict
            dictionaries to parse
        """
        if len(kwargs) > 0:
            self.options = Options(kwargs, default_dcts)
            self.plot_opts, self.fig_opts, self.sim_opts = self.options.parse()

    def copy(self):
        new_mgr = self.__class__({}, [{}])
        new_mgr.options = self.options
        new_mgr.plot_opts = self.plot_opts
        new_mgr.fig_opts = self.fig_opts
        new_mgr.sim_opts = self.sim_opts
        return new_mgr
    
    @Expander(cn.KWARGS, cn.PLOT_KWARGS)
    def doPlotOpts(self):
        """
        Executes codes for the single plot options

        Parameters
        ----------
        #@expand

        Returns
        -------
        Axes
        """
        new_kwargs = {k: self.plot_opts[k] if k in self.plot_opts else v for k, v in
             cn.PLOT_DCT.items()}
        ax  = new_kwargs[cn.O_AX]
        if ax is None:
             _, ax  = plt.subplots(1)
             new_kwargs[cn.O_AX]  = ax
        if new_kwargs[cn.O_LEGEND_SPEC] is not None:
            legend_spec = new_kwargs[cn.O_LEGEND_SPEC]
            ax.legend(legend_spec.names,
                  bbox_to_anchor=legend_spec.crd,
                  loc=legend_spec.loc)
        if new_kwargs[cn.O_TITLE] != cn.PLOT_DCT[cn.O_TITLE]:
            ax.set_title(new_kwargs[cn.O_TITLE])
        if new_kwargs[cn.O_XLABEL] != cn.PLOT_DCT[cn.O_XLABEL]:
            ax.set_xlabel(new_kwargs[cn.O_XLABEL])
        if new_kwargs[cn.O_XLIM] is not None:
            ax.set_ylim(new_kwargs[cn.O_XLIM])
        if new_kwargs[cn.O_XTICKLABELS] is not None:
            ax.set_xticklabels(new_kwargs[cn.O_XTICKLABELS])
        if new_kwargs[cn.O_YLABEL] != cn.PLOT_DCT[cn.O_YLABEL]:
            ax.set_ylabel(new_kwargs[cn.O_YLABEL])
        if new_kwargs[cn.O_YLIM] is not None:
            ax.set_ylim(new_kwargs[cn.O_YLIM])
        if new_kwargs[cn.O_YTICKLABELS] is not None:
            ax.set_yticklabels(new_kwargs[cn.O_YTICKLABELS])
        return new_kwargs[cn.O_AX]

    @Expander(cn.KWARGS, cn.FIG_KWARGS)
    def doFigOpts(self):
        """
        Executes figure options.

        Parameters
        ----------
        #@expand
        """
        new_kwargs = {k: self.fig_opts[k] if k in self.fig_opts else v for k, v in
             cn.FIG_DCT.items()}
        plt.suptitle(new_kwargs[cn.O_SUPTITLE])
        if new_kwargs[cn.O_IS_PLOT]:
            plt.show()
