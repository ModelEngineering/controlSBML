"""LTI Control for SBML models"""

"""
TO DO:

1. Plot difference between time jacoabian at reference vs. Current.
2. Plot TOTAL residual SSQ vs. jacobian difference
"""

from controlSBML.make_roadrunner import makeRoadrunner

import control
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tellurium as te


TIME = "time"
START_TIME = 0  # Default start time
END_TIME = 5  # Default endtime
NUM_POINT = 101
TIME = "time"
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
SIM_OPTS = dict(
      start_time=START_TIME,  # Start time of the simulation
      end_time=END_TIME,      # End time of the simulation
      num_point=NUM_POINT,    # Number of points in the simulation
      )
# Options for a single plot
PLOT_OPTS = dict(
      ylim=None,           # maximum and minimum value of y
      xlabel="",           
      ylabel="",           
      title="",             # plot title
      legend_spec=None,     # LegendSpec
      ax=None,              # axis to plot
      xticklabels=None,
      yticklabels=None,
      )
# Options for the full figure
FIG_OPTS = dict(
      is_plot=True,         # Is a figure generated
      figsize=(10, 10),     # Size of the figure
      suptitle="",          # Title for the figure
      )
ALL_OPTS = list(SIM_OPTS.keys())
ALL_OPTS.extend(list(PLOT_OPTS.keys()))
ALL_OPTS.extend(list(FIG_OPTS.keys()))


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

    def simulateLinearSystem(self, A_mat=None, timepoint=0, start_time=0,
          end_time=END_TIME, num_point=NUM_POINT):
        """
        Creates an approximation of the SBML model based on the Jacobian, and
        constructs predictions based on this Jacobian and the values of
        floating species at the start_time.

        Parameters
        ----------
        timepoint: float
        start_time: float
        end_time: float
        
        Returns
        -------
        pd.dataframe
            columns: floating species
            index: time
        """
        cur_time = self.get("time")
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

    def simulateRoadrunner(self, start_time=0, end_time=END_TIME,
          num_point=NUM_POINT):
        """
        Runs a new roadrunner simulation.

        Parameters
        ----------
        start_time: float
        end_time: float
        num_point: int
        
        Returns
        -------
        pd.dataframe
            columns: floating species
            index: time
        """
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
        is_plot: bool
        start_time: float
        end_time: float
        num_point: int
        y_max: float
            max value of y
        legend_crd: (float, float)
            coordinates for the legend
        """
        # Parse the options
        plot_opts, fig_opts, sim_opts = self._parseOpts(**kwargs)
        # Run the simulation
        df = self.simulateRoadrunner(**sim_opts)
        # Adjust the option values
        plot_opts["xlabel"] = TIME
        if plot_opts["ylim"] is None:
            y_max = df.max().max()
            plot_opts["ylim"] = [0, y_max]
        if plot_opts["legend_spec"] is None:
            legend_spec = LegendSpec(df.columns)
        ax = self._doPlotOpts(**plot_opts)
        # Do the plot
        for col in df.columns:
            ax.plot(df.index, df[col])
        # Finalize the figure
        self._doFigOpts(**fig_opts)

    @classmethod
    def _parseOpts(cls, **kwargs):
        """
        Parses options into plot, figure, and simulation

        Parameters
        ----------
        kwargs: dict
        
        Returns
        -------
        plot_opts: dict
        fig_opts: dict
        sim_opts: dict
        """
        # Validate
        unknown_options = set(kwargs.keys()).difference(ALL_OPTS)
        if len(unknown_options) > 0:
            raise ValueError("Unknown options: %s" % str(unknown_options))
        #
        plot_opts = {k: kwargs[k] if k in kwargs.keys() else v
              for k, v in PLOT_OPTS.items()}
        fig_opts = {k: kwargs[k] if k in kwargs.keys() else v
              for k, v in FIG_OPTS.items()}
        sim_opts = {k: kwargs[k] if k in kwargs.keys() else v
              for k, v in SIM_OPTS.items()}
        return plot_opts, fig_opts, sim_opts

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
        plot_opts, fig_opts, sim_opts = self._parseOpts(**kwargs)
        start_time = sim_opts["start_time"]
        rr_df = self.simulateRoadrunner(**sim_opts)
        nrow = 1
        ncol = len(rr_df.columns)
        fig, axes = plt.subplots(nrow, ncol, figsize=fig_opts["figsize"])
        axes = np.reshape(axes, (nrow, ncol))
        linear_df = self.simulateLinearSystem(timepoint=start_time,
              A_mat=A_mat, **sim_opts)
        y_min = min(linear_df.min().min(), rr_df.min().min())
        y_max = max(linear_df.max().max(), rr_df.max().max())
        plot_opts["ylim"] = [y_min, y_max]
        irow = 0
        base_plot_opts = dict(plot_opts)
        for icol, column in enumerate(rr_df.columns):
            plot_opts = dict(base_plot_opts)
            ax = axes[irow, icol]
            ax.plot(linear_df.index, linear_df[column], color="red")
            ax.plot(rr_df.index, rr_df[column], color="blue")
            if irow < nrow - 1:
                plot_opts["xticklabels"] = []
            if irow == 0:
                ax.set_title(column, rotation=45)
                if icol == 0:
                    plot_opts["legend_spec"] = LegendSpec(
                          ["approximation", "true"])
                else:
                    plot_opts["legend_spec"] = None
            if icol > 0:
                ax.set_yticklabels([])
            plot_opts["ax"] = ax
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
        plot_opts, fig_opts, sim_opts = cls._parseOpts(**kwargs)
        if isinstance(timepoints, float) or isinstance(timepoints, int):
            timepoints = [timepoints]
        ctlsb = cls(model_reference)
        rr_df = ctlsb.simulateRoadrunner(**sim_opts)
        nrow = len(timepoints)
        ncol = len(rr_df.columns)
        fig, axes = plt.subplots(nrow, ncol, figsize=fig_opts["figsize"])
        axes = np.reshape(axes, (nrow, ncol))
        for irow, timepoint in enumerate(timepoints):
            linear_df = ctlsb.simulateLinearSystem(timepoint=timepoint, **sim_opts)
            if plot_opts["ylim"] is None:
                y_min = min(linear_df.min().min(), rr_df.min().min())
                y_max = max(linear_df.max().max(), rr_df.max().max())
                plot_opts["ylim"] = [y_min, y_max]
            base_plot_opts = dict(plot_opts)
            for icol, column in enumerate(rr_df.columns):
                plot_opts = dict(base_plot_opts)
                ax = axes[irow, icol]
                plot_opts["ax"] = ax
                ax.plot(linear_df.index, linear_df[column], color="red")
                ax.plot(rr_df.index, rr_df[column], color="blue")
                ax.scatter(timepoint, y_min, s=40, marker="o", color="g")
                if irow < nrow - 1:
                    plot_opts["xticklabels"] = []
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
        ax  = new_kwargs["ax"]
        if ax is None:
             _, ax  = plt.subplots(1)
             new_kwargs["ax"]  = ax
        if new_kwargs["ylim"] is not None:
            ax.set_ylim(new_kwargs["ylim"])
        if new_kwargs["xlabel"] is not None:
            ax.set_xlabel(new_kwargs["xlabel"])
        if new_kwargs["ylabel"] is not None:
            ax.set_ylabel(new_kwargs["ylabel"])
        if new_kwargs["title"] is not None:
            ax.set_title(new_kwargs["title"])
        if new_kwargs["xticklabels"] is not None:
            ax.set_xticklabels(new_kwargs["xticklabels"])
        if new_kwargs["yticklabels"] is not None:
            ax.set_yticklabels(new_kwargs["yticklabels"])
        if new_kwargs["legend_spec"] is not None:
            legend_spec = new_kwargs["legend_spec"]
            ax.legend(legend_spec.names,
                  bbox_to_anchor=legend_spec.crd,
                  loc=legend_spec.loc)
        return new_kwargs["ax"]

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
        plt.suptitle(new_kwargs["suptitle"])
        if new_kwargs["is_plot"]:
            plt.show()
