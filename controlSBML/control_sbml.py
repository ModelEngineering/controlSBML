"""LTI Control for SBML models"""

"""
TO DO:

1. Plot difference between time jacoabian at reference vs. Current.
2. Plot TOTAL residual SSQ vs. jacobian difference
"""

from controlSBML.make_roadrunner import makeRoadrunner
from controlSBML.options import Options

import control
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tellurium as te


TIME = "time"
START_TIME = 0  # Default start time
END_TIME = 5  # Default endtime
POINTS_PER_TIME = 10
TIME = "time"
# Options
O_AX = "ax"
O_END_TIME = "end_time"
O_FIGSIZE = "figsize"
O_IS_PLOT = "is_plot"
O_LEGEND_SPEC = "legend_spec"
O_POINTS_PER_TIME = "points_per_time"
O_START_TIME = "start_time"
O_SUPTITLE = "suptitle"
O_TITLE = "title"
O_XLABEL = "xlabel"
O_XTICKLABELS = "xticklabels"
O_YTICKLABELS = "yticklabels"
O_YLABEL = "ylabel"
O_YLIM = "ylim"

"""
        ax: Matplotlib.axes
        end_time: float (end time of simulation)
        figsize: float, float (widith, height)
        is_plot: bool (do the plot)
        legend_spec: LegendSpec (position of the legend)
"""

        

# Legend specification
class LegendSpec():

    def __init__(self, names, crd=(1.4, 1), loc="upper right"):
        """
        
        Parameters
        ----------
        names: list-srs - names of the legends
        crd: (float, float) - coordinate of the legend
        loc: str - position of the legend
        """
        self.names = list(names)
        self.crd = crd
        self.loc = loc


# OPtions for simulation methods
SIM_OPTS = Options(dict(
      start_time=START_TIME,  # Start time of the simulation
      end_time=END_TIME,      # End time of the simulation
      points_per_time=POINTS_PER_TIME,    # Number of points in the simulation
      ))
# Options for a single plot
PLOT_OPTS = Options(dict(
      ylim=None,           # maximum and minimum value of y
      xlabel="",           
      ylabel="",           
      title="",             # plot title
      legend_spec=None,     # LegendSpec
      ax=None,              # axis to plot
      xticklabels=None,
      yticklabels=None,
      ))
# Options for the full figure
FIG_OPTS = Options(dict(
      is_plot=True,         # Is a figure generated
      figsize=(10, 10),     # Size of the figure
      suptitle="",          # Title for the figure
      ))
OPTS_LST = [PLOT_OPTS, FIG_OPTS, SIM_OPTS]


class ControlSBML(object):

    def __init__(self, model_reference):
        """
        Initializes instance variables
        :param str model_reference: string or SBML file or Roadrunner object
        """
        ##### PUBLIC #####
        self.model_reference = model_reference
        self.roadrunner = makeRoadrunner(model_reference)
        self.antimony = self.roadrunner.getAntimony()
        # Do the initializations
        self.boundary_species = self.roadrunner.getBoundarySpeciesIds()
        self._species_names = None
        self.roadrunner.reset()

    def _mkBoundarySpeciesFloating(self):
        for name in self.boundary_species:
            self.roadrunner.setBoundary(name, False)

    def _unmkBoundarySpeciesFloating(self):
        for name in self.boundary_species:
            self.roadrunner.setBoundary(name, True)

    @property
    def species_names(self):
        if self._species_names is None:
            self._mkBoundarySpeciesFloating()
            try:
                mat = self.roadrunner.getFullJacobian()
                self._unmkBoundarySpeciesFloating()
            except Exception:
                self._unmkBoundarySpeciesFloating()
                mat = self.roadrunner.getFullJacobian()
            self._species_names = list(mat.colnames)
        return self._species_names

    @property
    def jacobian(self):
        """
        Returns
        -------
        NamedArray with names for rows (rownames) and columns (colnames)
        Handle boundary species.
        """
        self._mkBoundarySpeciesFloating()
        try:
            mat = self.roadrunner.getFullJacobian()
            self._unmkBoundarySpeciesFloating()
        except Exception:
            self._unmkBoundarySpeciesFloating()
            mat = self.roadrunner.getFullJacobian()
        for idx, name in enumerate(mat.rownames):
            if name in self.boundary_species:
                mat[idx, :] = 0
        df = pd.DataFrame(mat, columns=mat.colnames, index=self.species_names)
        return df

    @property
    def current_state(self):
        """
        Contructs vector of current state values (floating and boundary species)

        Returns
        -------
        Series
            index: str (state names)
        """
        values = list(self.get(self.species_names).values())
        return pd.Series(values, index=self.species_names)

    @staticmethod
    def isRoadrunnerKey(key):
        return not ((key[0] == "_") or ("(" in key) or (key[-1] == "'"))

    def setTime(self, time):
        self.roadrunner.reset()
        _ = self.roadrunner.simulate(0, time)

    def copy(self):
        """
        Creates a copy of the object.

        Returns
        -------
        controlSBML
        """
        ctlsb = ControlSBML(self.model_reference)
        # Update roadrunner
        for key, value in self.roadrunner.items():
            if self.isRoadrunnerKey(key):
                ctlsb.roadrunner[key] = value
        return ctlsb

    def equals(self, other):
        """
        Checks that they have the same information

        Parameters
        ----------
        other: ControlSBML

        Returns
        -------
        bool
        """
        bValue = self.antimony == other.antimony
        bValue = bValue and all([s1 == s2 for s1, s2
              in zip(self.species_names, self.jacobian.columns)])
        diff = set(self.roadrunner.keys()).symmetric_difference(
              other.roadrunner.keys())
        bValue = bValue and (len(diff) == 0)
        # Check the roadrunner state
        if bValue:
            for key, value in self.roadrunner.items():
                if self.isRoadrunnerKey(key):
                    bValue = bValue and (other.roadrunner[key] == value)
        return bValue

    def get(self, names=None):
        """
        Provides the roadrunner values for a name. If no name,
        then all values are given.

        Parameters
        ----------
        name: str/list-str

        Returns
        -------
        object/dict
        """
        if isinstance(names, str):
            return self.roadrunner[names]
        if names is None:
            names = self.roadrunner.keys()
        return {n: self.roadrunner[n] for n in names}

    def set(self, name_dct):
        """
        Sets the values of names and values.

        Parameters
        ----------
        name_dct: dict
            key: str
            value: value
        """
        for name, value in name_dct.items():
            self.roadrunner[name] = value

    def makeStateSpace(self, A=None, B=None, C=None, D=None):
        """
        Creates a control system object for the n X n jacobian.

        Parameters
        ----------
        A: np.array(n X n)
        B: np.array(n X p)
        C: np.array(q X n)
        D: np.array(q X p)

        Returns
        -------
        control.StateSpace
        """
        # Construct the matrices
        if A is None:
            A = self.jacobian.values
        if B is None:
            B = np.repeat(0, A.shape[0])
            B = np.reshape(B, (A.shape[0], 1))
        if C is None:
            C = np.identity(A.shape[0])
        if D is None:
            D = B
        return control.StateSpace(A, B, C, D)

    @staticmethod
    def _getSimulationParameters(sim_opts):
        start_time = sim_opts[O_START_TIME]
        end_time = sim_opts[O_END_TIME]
        points_per_time = sim_opts[O_POINTS_PER_TIME]
        num_points = int(points_per_time*(start_time - end_time))
        return start_time, end_time, num_points

    def simulateLinearSystem(self, A_mat=None, timepoint=0, **kwargs):
        """
        Creates an approximation of the SBML model based on the Jacobian, and
        constructs predictions based on this Jacobian and the values of
        floating species at the start_time.

        Parameters
        ----------
        kwargs: dict
            SIM_OPTS
        
        Returns
        -------
        pd.dataframe
            columns: floating species
            index: time
        """
        options = Options(kwargs)
        sim_opts = options.parse(SIM_OPTS)
        start_time, end_time, num_points = self._getSimulationParameters(sim_opts)
        cur_time = self.get(TIME)
        self.setTime(timepoint)
        sys = self.makeStateSpace(A=A_mat)
        self.setTime(start_time)
        x0 = self.current_state
        self.setTime(cur_time)  # Restore the time
        # Run the linear simulation
        dt = (end_time - start_time)/num_point
        times = [start_time + n*dt for n in range(num_point)]
        times, y_vals = control.forced_response(sys, T=times, X0=x0)
        df = pd.DataFrame(y_vals.transpose(), index=times)
        df.columns = self.species_names
        return df

    def simulateRoadrunner(self, **kwargs):
        """
        Runs a new roadrunner simulation.

        Parameters
        ----------
        kwargs: dict
            SIM_OPTS
        
        Returns
        -------
        pd.dataframe
            columns: floating species
            index: time
        """
        options = Options(kwargs)
        sim_opts = options.parse(SIM_OPTS)
        start_time, end_time, num_points = self._getSimulationParameters(sim_opts)
        #
        self.roadrunner.reset()
        data = self.roadrunner.simulate(start_time, end_time, num_point)
        columns = [c[1:-1] if c[0] =="[" else c for c in data.colnames]
        df = pd.DataFrame(data, columns=columns)
        df = df.set_index(TIME)
        return df

    def plotTrueModel(self, **kwargs):
        """
        Plots the underlying SBML model.

        Parameters
        ----------
        kwargs: dict
            PLOT_OPTS, FIG_OPTS
        """
        # Parse the options
        options = Options(kwargs)
        options = options.parse(options)
        plot_opts, fig_opts, sim_opts = options.parse(OPTS_LST)
        # Run the simulation
        df = self.simulateRoadrunner(**sim_opts)
        # Adjust the option values
        plot_opts.set(O_XLABEL, default=TIME)
        y_max = df.max().max()
        plot_opts.set(O_YLIM, default=[0, y_max])
        plot_opts.set(O_LEGEND_SPEC, default=LegendSpec(df.columns))
        ax = self._doPlotOpts(**plot_opts)
        # Do the plot
        for col in df.columns:
            ax.plot(df.index, df[col])
        # Finalize the figure
        self._doFigOpts(**fig_opts)

    def plotLinearApproximation(self, A_mat=None, **kwargs):
        """
        Creates a plot that compares the linear approximation with the true model.

        Parameters
        ----------
        A_mat: A matrix of approximation model
            default is Jacobian at current time
        kwargs: dict
            a combination of plot, figure, and simulation options
        """
        options = Options(kwargs)
        plot_opts, fig_opts, sim_opts = options.parse(OPTS_LST)
        start_time = sim_opts[O_START_TIME]
        rr_df = self.simulateRoadrunner(**sim_opts)
        nrow = 1
        ncol = len(rr_df.columns)
        fig, axes = plt.subplots(nrow, ncol, figsize=fig_opts[O_FIGSIZE])
        axes = np.reshape(axes, (nrow, ncol))
        linear_df = self.simulateLinearSystem(timepoint=start_time,
              A_mat=A_mat, **sim_opts)
        y_min = min(linear_df.min().min(), rr_df.min().min())
        y_max = max(linear_df.max().max(), rr_df.max().max())
        plot_opts[O_YLIM] = [y_min, y_max]
        irow = 0
        base_plot_opts = dict(plot_opts)
        for icol, column in enumerate(rr_df.columns):
            plot_opts = dict(base_plot_opts)
            ax = axes[irow, icol]
            ax.plot(linear_df.index, linear_df[column], color="red")
            ax.plot(rr_df.index, rr_df[column], color="blue")
            if irow < nrow - 1:
                plot_opts[O_XTICKLABELS] = []
            if irow == 0:
                ax.set_title(column, rotation=45)
                if icol == 0:
                    plot_opts[O_LEGEND_SPEC] = LegendSpec(
                          ["approximation", "true"])
                else:
                    plot_opts[O_LEGEND_SPEC] = None
            if icol > 0:
                ax.set_yticklabels([])
            plot_opts[O_AX] = ax
            _ = self._doPlotOpts(**plot_opts)
        self._doFigOpts(**fig_opts)

    @classmethod
    def evaluateAccuracy(cls, model_reference, timepoints, **kwargs):
        """
        Creates a plot that evaluates the
        accouract of a linear model where the Jacobian is calculated
        at multiple timepoints. X0 is taken at the start of the simulation.

        Parameters
        ----------
        model_reference: as required by constructor
        timepoints: list-float
            Time at which Jacobian is calculated
        kwargs: dict
            SIM_OPTS, PLOT_OPTS, FIG_OPTS
        """
        options = Options(kwargs)
        plot_opts, fig_opts, sim_opts = options.parse(OPTS_LST)
        if isinstance(timepoints, float) or isinstance(timepoints, int):
            timepoints = [timepoints]
        ctlsb = cls(model_reference)
        rr_df = ctlsb.simulateRoadrunner(**sim_opts)
        nrow = len(timepoints)
        ncol = len(rr_df.columns)
        fig, axes = plt.subplots(nrow, ncol, figsize=fig_opts[O_FIGSIZE])
        axes = np.reshape(axes, (nrow, ncol))
        for irow, timepoint in enumerate(timepoints):
            linear_df = ctlsb.simulateLinearSystem(timepoint=timepoint, **sim_opts)
            if plot_opts[O_YLIM] is None:
                y_min = min(linear_df.min().min(), rr_df.min().min())
                y_max = max(linear_df.max().max(), rr_df.max().max())
                plot_opts[O_YLIM] = [y_min, y_max]
            base_plot_opts = dict(plot_opts)
            for icol, column in enumerate(rr_df.columns):
                plot_opts = dict(base_plot_opts)
                ax = axes[irow, icol]
                plot_opts[O_AX] = ax
                ax.plot(linear_df.index, linear_df[column], color="red")
                ax.plot(rr_df.index, rr_df[column], color="blue")
                ax.scatter(timepoint, y_min, s=40, marker="o", color="g")
                if irow < nrow - 1:
                    plot_opts[O_XTICKLABELS] = []
                if irow == 0:
                    ax.set_title(column, rotation=45)
                    if icol == 0:
                        ax.text(-3, 0.75*y_max, "Jacobian Time")
                        ax.legend(["linear", "nonlinear"])
                if icol > 0:
                    ax.set_yticklabels([])
                else:
                    ax.text(-2, y_max/2, "%2.1f" % timepoint)
                _ = cls._doPlotOpts(**plot_opts)
        cls._doFigOpts(**fig_opts)

    @classmethod
    def _doPlotOpts(cls, **kwargs):
        """
        Executes codes for the single plot options

        Parameters
        ----------
        kwargs: dict
               see PLOT_OPTS

        Returns
        -------
        Axes
        """
        new_kwargs = {k: kwargs[k] if k in kwargs else v for k, v in
             PLOT_OPTS.items()}
        ax  = new_kwargs[O_AX]
        if ax is None:
             _, ax  = plt.subplots(1)
             new_kwargs[O_AX]  = ax
        if new_kwargs[O_YLIM] is not None:
            ax.set_ylim(new_kwargs[O_YLIM])
        if new_kwargs[O_XLABEL] is not None:
            ax.set_xlabel(new_kwargs[O_XLABEL])
        if new_kwargs[O_YLABEL] is not None:
            ax.set_ylabel(new_kwargs[O_YLABEL])
        if new_kwargs[O_TITLE] is not None:
            ax.set_title(new_kwargs[O_TITLE])
        if new_kwargs[O_XTICKLABELS] is not None:
            ax.set_xticklabels(new_kwargs[O_XTICKLABELS])
        if new_kwargs[O_YTICKLABELS] is not None:
            ax.set_yticklabels(new_kwargs[O_YTICKLABELS])
        if new_kwargs[O_LEGEND_SPEC] is not None:
            legend_spec = new_kwargs[O_LEGEND_SPEC]
            ax.legend(legend_spec.names,
                  bbox_to_anchor=legend_spec.crd,
                  loc=legend_spec.loc)
        return new_kwargs[O_AX]

    @classmethod
    def _doFigOpts(cls, **kwargs):
        """
        Executes figure options.

        Parameters
        ----------
        kwargs: dict
            see FIG_OPTS       
        """
        new_kwargs = {k: kwargs[k] if k in kwargs else v for k, v in
             FIG_OPTS.items()}
        plt.suptitle(new_kwargs[O_SUPTITLE])
        if new_kwargs[O_IS_PLOT]:
            plt.show()
