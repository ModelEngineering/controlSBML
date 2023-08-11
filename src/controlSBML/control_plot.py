"""LTI Control for SBML models"""


import controlSBML.constants as cn
from controlSBML.control_analysis import ControlAnalysis
from controlSBML.option_management.option_manager import OptionManager
from controlSBML.util import PlotResult
from controlSBML import util
from controlSBML.nonlinear_io_system import NonlinearIOSystem

import control
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

        Returns
        ------
        PlotResult
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
        return PlotResult(ax=ax)
    
    def makeNonlinearIOSystem(self, name="", **kwargs):
        """
        Creates an object that can be used in connections with the
        control package.

        Parameters
        ----------
        name: str (name of the system)
        kwargs: dict (additional arguments for NonlinearIOSystem)

        Returns
        -------
        controlSBML.NonelinearIOSystem
        """
        if "input_names" not in kwargs:
            kwargs["input_names"] = self.input_names
        if "output_names" not in kwargs:
            kwargs["output_names"] = self.output_names
        return NonlinearIOSystem(name, self, **kwargs)

    # TODO: Deprecate plotLinearApproximation. Use plotAccuracy instead.
    @Expander(cn.KWARGS, cn.PLOT_KWARGS)
    def plotLinearApproximation(self, transfer_function_df, **kwargs):
        """
        Creates a plot that compares the linear approximation with the true model.

        Parameters
        ----------
        transfer_function_df: pd.DataFrame
            Transfer functions
        #@expand

        Returns
        ------
        PlotResult
        """
        mgr = OptionManager(kwargs)
        start_time = mgr.sim_opts[cn.O_START_TIME]
        rr_ts = self.simulateRoadrunner(**mgr.sim_opts)
        nrow = 1
        ncol = len(self.output_names)
        fig, axes = plt.subplots(nrow, ncol, figsize=mgr.fig_opts[cn.O_FIGSIZE])
        axes = np.reshape(axes, (nrow, ncol))
        linear_ts = self.simulateLinearSystem(op_time=start_time,
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
        return PlotResult(ax=axes, fig=fig, time_series=linear_ts)

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

        Returns
        ------
        PlotResult
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
        fig, axes = plt.subplots(nrow, ncol, figsize=mgr.fig_opts[cn.O_FIGSIZE])
        axes = np.reshape(axes, (nrow, ncol))
        for irow, time in enumerate(times):
            linear_ts = ctlsb.simulateLinearSystem(op_time=time,
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
        return PlotResult(ax=axes, fig=fig, time_series=rr_ts)

    @staticmethod
    def _setYAxColor(ax, position, color):
        # Set the colors of the labels, axes, and spines
        ax.tick_params(axis='y', labelcolor=color)
        ax.spines[position].set_color(color)
        ax.yaxis.label.set_color(color)

    # TESTME
    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotBodeTF(self, transfer_function_df, is_magnitude=True, is_phase=True, **kwargs):
        """
        Constructs bode plots for the transfer functions provided.

        Parameters
        ----------
        transfer_function_df: pd.DataFrame
        is_magnitude: bool
            Do magnitude plots
        is_phase: bool
            Do phase plots
        is_plot: bool
            Display plots
        #@expand

        Returns
        ------
        PlotResult
        """
        mgr = OptionManager(kwargs)
        figure, axes = plt.subplots(len(transfer_function_df.index), 
                                    len(transfer_function_df.columns),
                                    figsize=mgr.fig_opts[cn.O_FIGSIZE])
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.4, hspace=0.2)
        # FIXME: Display fitted transfer function on Bode plot
        mgr.fig_opts[cn.O_FIGURE] = figure
        for input_name in transfer_function_df.index:
            for output_name in transfer_function_df.columns:
                icol = self.input_names.index(input_name)
                irow = self.output_names.index(output_name)
                if (icol ==0) and (irow == 0):
                    is_label = True
                else:
                    is_label = False
                ax = axes[irow, icol]
                mgr.plot_opts[cn.O_AX] = ax
                #
                tf = transfer_function_df.loc[input_name, output_name]
                if tf is None:
                    continue
                #
                magnitude_arr, phase_arr, frequency_arr  = control.bode(tf, plot=False, deg=False, dB=True)
                ax.set_title("%s->%s" % (input_name, output_name))
                if is_magnitude:
                    color="red"
                    ax.plot(frequency_arr, magnitude_arr, color=color)
                    ax.set_yscale("log")
                    ax.set_xscale("log")
                    self._setYAxColor(ax, "left", color)
                    ax.grid(axis="x")
                if is_phase:
                    ax2 = ax.twinx()
                    color="blue"
                    ax2.plot(frequency_arr, phase_arr, color=color)
                    ax2.set_xscale("log")
                    self._setYAxColor(ax2, "right", color)
                    ax2.grid(axis="x")
                latex = util.latexifyTransferFunction(tf)
                if len(mgr.plot_opts[cn.O_TITLE]) == 0:
                    title = latex
                else:
                    title = mgr.plot_opts[cn.O_TITLE]
                ax.set_title(title, y=0.2, pad=-14, fontsize=14, loc="right")
                ax.legend(["magnitude (dB)"], loc="upper left", fontsize=8)
                ax2.legend(["phase (rad)"], loc="upper right", fontsize=8)
                if irow == len(self.output_names) - 1:
                    ax.set_xlabel("frequency (rad)")
                mgr.doPlotOpts()
        mgr.doFigOpts()