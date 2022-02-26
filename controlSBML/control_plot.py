"""LTI Control for SBML models"""


import controlSBML.constants as cn
from controlSBML.control_analysis import ControlAnalysis
from controlSBML.option_management.option_manager import OptionManager

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
        mgr = OptionManager(kwargs)
        # Run the simulation
        df = self.simulateRoadrunner(**mgr.sim_opts)
        # Adjust the option values
        mgr.plot_opts.set(cn.O_XLABEL, default=cn.TIME)
        y_max = df.max().max()
        mgr.plot_opts.set(cn.O_YLIM, default=[0, y_max])
        if cn.O_LEGEND_CRD in mgr.plot_opts.keys():
            legend_spec =cn.LegendSpec(df.columns,
                  crd=mgr.plot_opts[cn.O_LEGEND_CRD])
        else:
            legend_spec =cn.LegendSpec(df.columns)
        mgr.plot_opts.set(cn.O_LEGEND_SPEC, default=legend_spec)
        ax = mgr.doPlotOpts()
        # Do the plot
        for col in df.columns:
            ax.plot(df.index, df[col])
        mgr.plot_opts.set(cn.O_AX, ax)
        _ = mgr.doPlotOpts()  # Recover lost plot options
        # Finalize the figure
        mgr.doFigOpts()

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
        mgr = OptionManager(kwargs)
        start_time = mgr.sim_opts[cn.O_START_TIME]
        rr_df = self.simulateRoadrunner(**mgr.sim_opts)
        nrow = 1
        ncol = len(rr_df.columns)
        _, axes = plt.subplots(nrow, ncol, figsize=mgr.fig_opts[cn.O_FIGSIZE])
        axes = np.reshape(axes, (nrow, ncol))
        linear_df = self.simulateLinearSystem(timepoint=start_time,
              A_df=A_mat, is_reduced=is_reduced, **mgr.sim_opts)
        y_min = min(linear_df.min().min(), rr_df.min().min())
        y_max = max(linear_df.max().max(), rr_df.max().max())
        mgr.plot_opts[cn.O_YLIM] = [y_min, y_max]
        irow = 0
        for icol, column in enumerate(linear_df.columns):
            new_mgr = mgr.copy()
            plot_opts = new_mgr.plot_opts
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
            _ = new_mgr.doPlotOpts()
        mgr.doFigOpts()

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
        mgr = OptionManager(kwargs)
        if isinstance(timepoints, float) or isinstance(timepoints, int):
            timepoints = [timepoints]
        if model_reference is not None:
            ctlsb = self.__class__(model_reference)
        else:
            ctlsb = self
        rr_df = ctlsb.simulateRoadrunner(**mgr.sim_opts)
        nrow = len(timepoints)
        ncol = len(rr_df.columns)
        _, axes = plt.subplots(nrow, ncol, figsize=mgr.fig_opts[cn.O_FIGSIZE])
        axes = np.reshape(axes, (nrow, ncol))
        for irow, timepoint in enumerate(timepoints):
            linear_df = ctlsb.simulateLinearSystem(timepoint=timepoint,
                  **mgr.sim_opts)
            if mgr.plot_opts[cn.O_YLIM] is None:
                y_min = min(linear_df.min().min(), rr_df.min().min())
                y_max = max(linear_df.max().max(), rr_df.max().max())
                mgr.plot_opts[cn.O_YLIM] = [y_min, y_max]
            y_min, y_max = mgr.plot_opts[cn.O_YLIM]
            for icol, column in enumerate(linear_df.columns):
                new_mgr = mgr.copy()
                ax = axes[irow, icol]
                new_mgr.plot_opts[cn.O_AX] = ax
                ax.plot(linear_df.index, linear_df[column], color="red")
                ax.plot(rr_df.index, rr_df[column], color="blue")
                ax.scatter(timepoint, y_min, s=40,
                       marker="o", color="g")
                if irow < nrow - 1:
                    new_mgr.plot_opts[cn.O_XTICKLABELS] = []
                if irow == 0:
                    ax.set_title(column, rotation=45)
                    if icol == 0:
                        ax.text(-3, 0.75*y_max, "Jacobian Time")
                        ax.legend(["linear", "nonlinear"])
                if icol > 0:
                    ax.set_yticklabels([])
                else:
                    ax.text(-2, y_max/2, "%2.1f" % timepoint)
                _ = new_mgr.doPlotOpts()
        mgr.doFigOpts()
