"""LTI Control for SBML models"""

"""
TO DO:

1. Plot difference between time jacoabian at reference vs. Current.
2. Plot TOTAL residual SSQ vs. jacobian difference
"""

import controlSBML.constants as cn
from controlSBML.control_analysis import ControlAnalysis
from controlSBML.options import Options

from docstring_expander.expander import Expander
import matplotlib.pyplot as plt
import numpy as np


class ControlPlot(ControlAnalysis):

    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotTrueModel(self, **kwargs):
        """
        Plots the underlying SBML model.

        Parameters
        ----------
        #@expand
        """
        # Parse the options
        options = Options(kwargs, cn.DEFAULT_DCTS)
        plot_opts, fig_opts, sim_opts = options.parse()
        # Run the simulation
        df = self.simulateRoadrunner(**sim_opts)
        # Adjust the option values
        plot_opts.set(cn.O_XLABEL, default=cn.TIME)
        y_max = df.max().max()
        plot_opts.set(cn.O_YLIM, default=[0, y_max])
        if cn.O_LEGEND_CRD in plot_opts.keys():
            legend_spec =cn.LegendSpec(df.columns, crd=plot_opts[cn.O_LEGEND_CRD])
        else:
            legend_spec =cn.LegendSpec(df.columns)
        plot_opts.set(cn.O_LEGEND_SPEC, default=legend_spec)
        ax = self._doPlotOpts(**plot_opts)
        # Do the plot
        for col in df.columns:
            ax.plot(df.index, df[col])
        plot_opts.set(cn.O_AX, ax)
        _ = self._doPlotOpts(**plot_opts)  # Recover lost plot options
        # Finalize the figure
        self._doFigOpts(**fig_opts)

    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotLinearApproximation(self, A_mat=None, is_reduced=True, **kwargs):
        """
        Creates a plot that compares the linear approximation with the true model.

        Parameters
        ----------
        A_mat: A matrix of approximation model
            default is Jacobian at current time
        is_reduced: bool
            only available if A_mat is None. Construct reduced order model.
        #@expand
        """
        options = Options(kwargs, cn.DEFAULT_DCTS)
        plot_opts, fig_opts, sim_opts = options.parse()
        start_time = sim_opts[cn.O_START_TIME]
        rr_df = self.simulateRoadrunner(**sim_opts)
        nrow = 1
        ncol = len(rr_df.columns)
        _, axes = plt.subplots(nrow, ncol, figsize=fig_opts[cn.O_FIGSIZE])
        axes = np.reshape(axes, (nrow, ncol))
        linear_df = self.simulateLinearSystem(timepoint=start_time,
              A_df=A_mat, is_reduced=is_reduced, **sim_opts)
        y_min = min(linear_df.min().min(), rr_df.min().min())
        y_max = max(linear_df.max().max(), rr_df.max().max())
        plot_opts[cn.O_YLIM] = [y_min, y_max]
        irow = 0
        base_plot_opts = Options(plot_opts, plot_opts.default_dcts)
        for icol, column in enumerate(linear_df.columns):
            plot_opts = Options(base_plot_opts, base_plot_opts.default_dcts)
            ax = axes[irow, icol]
            ax.plot(linear_df.index, linear_df[column], color="red")
            ax.plot(rr_df.index, rr_df[column], color="blue")
            if irow < nrow - 1:
                plot_opts[cn.O_XTICKLABELS] = []
            if irow == 0:
                ax.set_title(column, rotation=45)
                if icol == 0:
                    names = ["approximation", "true"]
                    if cn.O_LEGEND_CRD in plot_opts.keys():
                        legend_spec = cn.LegendSpec(
                              names, crd=plot_opts[cn.O_LEGEND_CRD])
                    else:
                        legend_spec =cn.LegendSpec(name)
                    plot_opts.set(cn.O_LEGEND_SPEC, default=legend_spec)
                else:
                    plot_opts[cn.O_LEGEND_SPEC] = None
            if icol > 0:
                ax.set_yticklabels([])
            plot_opts[cn.O_AX] = ax
            _ = self._doPlotOpts(**plot_opts)
        self._doFigOpts(**fig_opts)

    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotAccuracy(self, model_reference=None, timepoints=[0], **kwargs):
        """
        Creates a plot that evaluates the
        accouract of a linear model where the Jacobian is calculated
        at multiple timepoints. X0 is taken at the start of the simulation.

        Parameters
        ----------
        model_reference: as required by constructor
            default: current model
        timepoints: list-float
            Time at which Jacobian is calculated
            default: [0]
        #@expand
        """
        options = Options(kwargs, cn.DEFAULT_DCTS)
        plot_opts, fig_opts, sim_opts = options.parse()
        if isinstance(timepoints, float) or isinstance(timepoints, int):
            timepoints = [timepoints]
        if model_reference is not None:
            ctlsb = self.__class__(model_reference)
        else:
            ctlsb = self
        rr_df = ctlsb.simulateRoadrunner(**sim_opts)
        nrow = len(timepoints)
        ncol = len(rr_df.columns)
        _, axes = plt.subplots(nrow, ncol, figsize=fig_opts[cn.O_FIGSIZE])
        axes = np.reshape(axes, (nrow, ncol))
        for irow, timepoint in enumerate(timepoints):
            linear_df = ctlsb.simulateLinearSystem(timepoint=timepoint, **sim_opts)
            if plot_opts[cn.O_YLIM] is None:
                y_min = min(linear_df.min().min(), rr_df.min().min())
                y_max = max(linear_df.max().max(), rr_df.max().max())
                plot_opts[cn.O_YLIM] = [y_min, y_max]
            y_min, y_max = plot_opts[cn.O_YLIM]
            base_plot_opts = dict(plot_opts)
            for icol, column in enumerate(linear_df.columns):
                plot_opts = dict(base_plot_opts)
                ax = axes[irow, icol]
                plot_opts[cn.O_AX] = ax
                ax.plot(linear_df.index, linear_df[column], color="red")
                ax.plot(rr_df.index, rr_df[column], color="blue")
                ax.scatter(timepoint, y_min, s=40,
                       marker="o", color="g")
                if irow < nrow - 1:
                    plot_opts[cn.O_XTICKLABELS] = []
                if irow == 0:
                    ax.set_title(column, rotation=45)
                    if icol == 0:
                        ax.text(-3, 0.75*y_max, "Jacobian Time")
                        ax.legend(["linear", "nonlinear"])
                if icol > 0:
                    ax.set_yticklabels([])
                else:
                    ax.text(-2, y_max/2, "%2.1f" % timepoint)
                _ = self._doPlotOpts(**plot_opts)
        self._doFigOpts(**fig_opts)

    @classmethod
    @Expander(cn.KWARGS, cn.PLOT_KWARGS)
    def _doPlotOpts(cls, **kwargs):
        """
        Executes codes for the single plot options

        Parameters
        ----------
        #@expand

        Returns
        -------
        Axes
        """
        new_kwargs = {k: kwargs[k] if k in kwargs else v for k, v in
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

    @classmethod
    @Expander(cn.KWARGS, cn.FIG_KWARGS)
    def _doFigOpts(cls, **kwargs):
        """
        Executes figure options.

        Parameters
        ----------
        #@expand
        """
        new_kwargs = {k: kwargs[k] if k in kwargs else v for k, v in
             cn.FIG_DCT.items()}
        plt.suptitle(new_kwargs[cn.O_SUPTITLE])
        if new_kwargs[cn.O_IS_PLOT]:
            plt.show()
