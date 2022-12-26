"""LTI Control for SBML models"""


import controlSBML.constants as cn
from controlSBML.control_analysis import ControlAnalysis
from controlSBML.control_extensions.state_space_tf import StateSpaceTF
from controlSBML.option_management.option_manager import OptionManager

from docstring_expander.expander import Expander
import matplotlib.pyplot as plt
import numpy as np


class ControlPlot(ControlAnalysis):

    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotTrueModel(self, names=None, **kwargs):
        """
        Plots the underlying SBML model.

        Parameters
        ----------
        names: list-str
            variables to to plot
        #@expand
        """
        # Parse the options
        mgr = OptionManager(kwargs)
        # Run the simulation
        ts = self.simulateRoadrunner(**mgr.sim_opts)
        if names is not None:
            ts = ts[names]
        # Adjust the option values
        mgr.plot_opts.set(cn.O_XLABEL, default=cn.TIME)
        y_max = ts.max().max()
        mgr.plot_opts.set(cn.O_YLIM, default=[0, y_max])
        if cn.O_LEGEND_CRD in mgr.plot_opts.keys():
            legend_spec =cn.LegendSpec(ts.columns,
                  crd=mgr.plot_opts[cn.O_LEGEND_CRD])
        else:
            legend_spec =cn.LegendSpec(ts.columns)
        mgr.plot_opts.set(cn.O_LEGEND_SPEC, default=legend_spec)
        ax = mgr.getAx()
        # Do the plot
        for col in ts.columns:
            ax.plot(ts.times, ts[col])
        mgr.plot_opts.set(cn.O_AX, ax)
        mgr.doPlotOpts()  # Recover lost plot options
        # Finalize the figure
        mgr.doFigOpts()

    # TODO: Deprecate plotLinearApproximation. Use plotAccuracy instead.
    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotLinearApproximation(self, A_mat=None, **kwargs):
        """
        Creates a plot that compares the linear approximation with the true model.

        Parameters
        ----------
        A_mat: A matrix of approximation model
            default is Jacobian at current time
        #@expand
        """
        mgr = OptionManager(kwargs)
        start_time = mgr.sim_opts[cn.O_START_TIME]
        rr_ts = self.simulateRoadrunner(**mgr.sim_opts)
        nrow = 1
        ncol = len(self.output_names)
        _, axes = plt.subplots(nrow, ncol, figsize=mgr.fig_opts[cn.O_FIGSIZE])
        axes = np.reshape(axes, (nrow, ncol))
        linear_ts = self.simulateLinearSystem(time=start_time,
              **mgr.sim_opts)
        y_min = min(linear_ts.min().min(), rr_ts.min().min())
        y_max = max(linear_ts.max().max(), rr_ts.max().max())
        mgr.plot_opts[cn.O_YLIM] = [y_min, y_max]
        irow = 0
        for icol, column in enumerate(linear_ts.columns):
            new_mgr = mgr.copy()
            plot_opts = new_mgr.plot_opts
            ax = axes[irow, icol]
            ax.plot(linear_ts.times, linear_ts[column], color="red")
            ax.plot(rr_ts.times, rr_ts[column], color="blue")
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
                        legend_spec =cn.LegendSpec(names)
                    plot_opts.set(cn.O_LEGEND_SPEC, default=legend_spec)
                else:
                    plot_opts[cn.O_LEGEND_SPEC] = None
            if icol > 0:
                ax.set_yticklabels([])
            new_mgr.doPlotOpts()
        mgr.doFigOpts()

    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotAccuracy(self, model_reference=None, times=0, **kwargs):
        """
        Creates a plot that evaluates the
        accouract of a linear model where the Jacobian is calculated
        at multiple times. X0 is taken at the start of the simulation.

        Parameters
        ----------
        model_reference: as required by constructor
            default: current model
        times: list-float
            Time at which Jacobian is calculated
            default: [0]
        #@expand
        """
        if times == 0:
            times = [0]
        mgr = OptionManager(kwargs)
        if isinstance(times, float) or isinstance(times, int):
            times = [times]
        if model_reference is not None:
            ctlsb = self.__class__(model_reference)
        else:
            ctlsb = self
        rr_ts = ctlsb.simulateRoadrunner(**mgr.sim_opts)
        nrow = len(times)
        ncol = len(self.output_names)
        _, axes = plt.subplots(nrow, ncol, figsize=mgr.fig_opts[cn.O_FIGSIZE])
        axes = np.reshape(axes, (nrow, ncol))
        for irow, time in enumerate(times):
            linear_ts = ctlsb.simulateLinearSystem(time=time,
                  **mgr.sim_opts)
            #
            for icol, column in enumerate(self.output_names):
                start_time = mgr.sim_opts[cn.O_START_TIME]
                y_min = min(linear_ts.min().min(), rr_ts.min().min())
                y_max = max(linear_ts.max().max(), rr_ts.max().max())
                mgr.plot_opts.set(cn.O_YLIM, default=[y_min, y_max])
                y_min, y_max = mgr.plot_opts[cn.O_YLIM]
                x_min = start_time
                x_max = mgr.sim_opts[cn.O_END_TIME]
                mgr.plot_opts.set(cn.O_XLIM, [x_min, x_max])
                x_min, x_max = mgr.plot_opts[cn.O_XLIM]
                new_mgr = mgr.copy()
                ax = axes[irow, icol]
                new_mgr.plot_opts[cn.O_AX] = ax
                ax.plot(linear_ts.times, linear_ts[column], color="red")
                ax.plot(rr_ts.times, rr_ts[column], color="blue")
                ax.scatter(time, y_min, s=40,
                       marker="o", color="g")
                if irow < nrow - 1:
                    new_mgr.plot_opts[cn.O_XTICKLABELS] = []
                if irow == 0:
                    ax.set_title(column, rotation=45)
                    if icol == 0:
                        #ax.text(start_time-3, 0.75*y_max, "Jacobian Time")
                        ax.legend(["linear", "nonlinear"])
                if icol > 0:
                    ax.set_yticklabels([])
                else:
                    pass
                    #ax.text(-2, y_max/2, "%2.1f" % time)
                new_mgr.doPlotOpts()
        mgr.doFigOpts()

    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotBode(self, is_magnitude=True, is_phase=True, is_plot=True,
          **kwargs):
        """
        Constructs bode plots for a State Spaesystem.
        This is done by constructing n*n
        SISO systems where there n states.
    
        Parameters
        ----------
        is_magnitude: bool
            Do magnitude plots
        is_phase: bool
            Do phase plots
        is_plot: bool
            Display plots
        #@expand
        """
        mimo_sys = self.makeStateSpace()
        tf = StateSpaceTF(mimo_sys, input_names=self.input_names,
              output_names=self.output_names)
        tf.plotBode(is_magnitude==is_magnitude, is_phase==is_phase,
              is_plot=is_plot, **kwargs)
