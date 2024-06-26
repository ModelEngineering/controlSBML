"""Manages options used in ControlSBML."""

import controlSBML.constants as cn
from controlSBML.option_management.options import Options

from docstring_expander.expander import Expander # type: ignore
import matplotlib.pyplot as plt
import numpy as np
import os


class OptionManager(object):
    figure_idx = 0  # Index for automated figure generation

    def __init__(self, kwargs):
        """
        Parameters
        ----------
        kwargs: dict
            key word arguments
        default_dcts: list-dict
            dictionaries to parse
        """
        options = Options(kwargs, cn.DEFAULT_DCTS)
        self.plot_opts, self.fig_opts, self.sim_opts = options.parse()

    @property
    def options(self):
        options = dict(self.plot_opts)
        options.update(self.fig_opts)
        options.update(self.sim_opts)
        return options
    
    @property
    def plot_fig_options(self):
        options = dict(self.plot_opts)
        options.update(self.fig_opts)
        return options

    def copy(self):
        new_mgr = self.__class__({})
        new_mgr.plot_opts = Options(self.plot_opts, [cn.PLOT_DCT])
        new_mgr.fig_opts = Options(self.fig_opts, [cn.FIG_DCT])
        new_mgr.sim_opts = Options(self.sim_opts, [cn.SIM_DCT])
        return new_mgr

    def setYlim(self, values, is_override=False):
        """
        Sets the range for y based on the values provided.

        Parameters
        ----------
        values: list-float
        """
        arr = np.array(values)
        y_abs_max = max(max(arr), min(np.abs(arr)))
        ylim = [-y_abs_max, y_abs_max]
        override = None
        default = None
        if is_override:
            override = ylim
        else:
            default = ylim
        self.plot_opts.set(cn.O_YLIM, default=default, override=override)

    def setFigure(self, is_override=False):
        """
        Sets the current figure.
        
        Parameters
        ----------
        is_override: bool
            Set the current figure even if there is an existing setting
        """
        fig = plt.gcf()
        override = None
        default = None
        if is_override:
            override = fig
        else:
            default = fig
        self.fig_opts.set(cn.O_FIGURE, override=override, default=default)

    def setAx(self, is_override=False):
        """
        Sets the current axis.

        Parameters
        ----------
        is_override: bool
            Set the current axis even if there is an existing setting
        """
        ax = plt.gca()
        override = None
        default = None
        if is_override:
            override = ax
        else:
            default = ax
        self.plot_opts.set(cn.O_AX, override=override, default=default)

    def getAx(self, is_override=False):
        """
        Sets the current axis.

        Returns
        -------
        matplotlib.Axes
        """
        ax = plt.gca()
        self.plot_opts.set(cn.O_AX, default=ax)
        return self.plot_opts.get(cn.O_AX)
        
    @Expander(cn.KWARGS, cn.PLOT_KWARGS)
    def doPlotOpts(self):
        """
        Executes codes for the single plot options

        Parameters
        ----------
        #@expand
        """
        self.setAx()
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
        if new_kwargs[cn.O_LEGEND] is not None:
            o_legend = new_kwargs[cn.O_LEGEND]
            if isinstance(o_legend, list):
                ax.legend(new_kwargs[cn.O_LEGEND])
            elif not o_legend:
                ax.legend([])
        if new_kwargs[cn.O_TITLE] != cn.PLOT_DCT[cn.O_TITLE]:
            ax.set_title(new_kwargs[cn.O_TITLE])
        if new_kwargs[cn.O_XLABEL] != cn.PLOT_DCT[cn.O_XLABEL]:
            ax.set_xlabel(new_kwargs[cn.O_XLABEL])
        if new_kwargs[cn.O_XLABEL_ANGLE] is not None:
            labels = ax.get_xticklabels()
            xticks = ax.get_xticks()
            ax.set_xticks(xticks)
            ax.set_xticklabels(labels, rotation=new_kwargs[cn.O_XLABEL_ANGLE], ha='right')
        if new_kwargs[cn.O_XLIM] is not None:
            ax.set_xlim(new_kwargs[cn.O_XLIM])
        if new_kwargs[cn.O_XTICKLABELS] is not None:
            ax.set_xticklabels(new_kwargs[cn.O_XTICKLABELS])
        if new_kwargs[cn.O_YLABEL] != cn.PLOT_DCT[cn.O_YLABEL]:
            ax.set_ylabel(new_kwargs[cn.O_YLABEL])
        if new_kwargs[cn.O_YLIM] is not None:
            ax.set_ylim(new_kwargs[cn.O_YLIM])
        if new_kwargs[cn.O_YTICKLABELS] is not None:
            ax.set_yticklabels(new_kwargs[cn.O_YTICKLABELS])

    @Expander(cn.KWARGS, cn.FIG_KWARGS)
    def doFigOpts(self):
        """
        Executes figure options.

        Parameters
        ----------
        #@expand
        """
        self.setFigure()
        fig = self.fig_opts[cn.O_FIGURE]
        new_kwargs = {k: self.fig_opts[k] if k in self.fig_opts else v for k, v in
             cn.FIG_DCT.items()}
        plt.suptitle(new_kwargs[cn.O_SUPTITLE])
        fig_width, fig_length = self.fig_opts[cn.O_FIGSIZE]
        fig.set_size_inches(fig_width, fig_length)
        if new_kwargs[cn.O_IS_PLOT]:
            plt.show()
            self.writeFigure(new_kwargs[cn.O_WRITEFIG])

    @classmethod
    def writeFigure(cls, writefig):
        """
        Writes a figure to a path specified by the option.
        """
        is_write = False
        file_path = None
        if isinstance(writefig, str):
            filepath = writefig
            is_write = True
        elif isinstance(writefig, bool):
            if writefig:
                filepath = "figure_%d.pdf" % cls.figure_idx
                cls.figure_idx  += 1
                is_write = True
        if is_write:
            if file_path is None:
                file_path = os.path.join(cn.PLOT_DIR, filepath)
            plt.savefig(file_path)
