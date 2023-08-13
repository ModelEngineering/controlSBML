"""LTI Control for SBML models"""


import controlSBML.constants as cn
from controlSBML.control_analysis import ControlAnalysis
from controlSBML.option_management.option_manager import OptionManager
from controlSBML.util import PlotResult
from controlSBML import util

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
        mgr.doPlotOpts()
        # Finalize the figure
        mgr.doFigOpts()
        return PlotResult(ax=ax)
    
    @Expander(cn.KWARGS, cn.ALL_KWARGS)
    def plotAccuracy(self, transfer_function_df, op_time=None, times=None, **kwargs):
        """
        Creates a plot that evaluates the
        accouract of a set of linear models described by transfer functions. The transfer function
        input should be a species with a fixed value.

        Parameters
        ----------
        transfer_function_df: pd.DataFrame / control.TransferFunction
            columns: output species
            index: input species
            values: control.TransferFunction
        op_time: float (time at which operating point values are obtained); if None, use steady state
        times: np.array (times for simulation)
        #@expand                     
        """
        mgr = OptionManager(kwargs)
        rr_ts = self.simulateRoadrunner(**mgr.sim_opts)
        times = np.array(rr_ts.index)/cn.MS_IN_SEC
        input_names = list(transfer_function_df.index)
        output_names = list(transfer_function_df.columns)
        nrow = len(input_names)
        ncol = len(output_names)
        fig, axes = plt.subplots(nrow, ncol, figsize=mgr.fig_opts[cn.O_FIGSIZE])
        axes = np.reshape(axes, (nrow, ncol))
        for irow, input_name in enumerate(input_names):
            for icol, output_name in enumerate(output_names):
                tf = transfer_function_df.loc[input_name, output_name]
                linear_ts = self.simulateLinearSystem(tf, input_name, output_name,
                    op_time=op_time, U=self.get(input_name), times=times,
                    **mgr.sim_opts)
                predicted_ts = linear_ts[output_name]
                simulated_ts = rr_ts[output_name]
                x_min = min(predicted_ts.min(), simulated_ts.min())
                x_max = max(predicted_ts.max(), simulated_ts.max())
                mgr.plot_opts.set(cn.O_XLIM, default=[x_min, x_max])
                y_min = min(predicted_ts.min(), simulated_ts.min())
                y_max = max(predicted_ts.max(), simulated_ts.max())
                mgr.plot_opts.set(cn.O_YLIM, default=[y_min, y_max])
                ax = axes[irow, icol]
                mgr.plot_opts[cn.O_AX] = ax
                ax.scatter(rr_ts[output_name], linear_ts[output_name])
                ax.plot([x_min, x_max], [y_min, y_max], color='red')
                if irow < nrow - 1:
                    ax.set_xticklabels([])
                else:
                    ax.set_xlabel("simulated")
                ax.set_title("%s->%s" % (input_name, output_name))
                if icol > 0:
                    ax.set_yticklabels([])
                else:
                    ax.set_ylabel("predicted")
                if (irow == 0) and (icol == 0):
                    ax.legend(["actual", "ideal"])
                mgr.doPlotOpts()
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