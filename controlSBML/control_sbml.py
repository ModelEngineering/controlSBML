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
        Creates an approximation of the SBML model based on the Jacobian, and
        constructs predictions based on this Jacobian and the values of
        floating species at the start_time.

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

    def plotTrueModel(self, is_plot=True, start_time=START_TIME,
          end_time=END_TIME, num_point=NUM_POINT, y_max=None,
          figsize=(10, 10), legend_crd=(1.4, 1),  ax=None):
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
        if ax is None:
            _, ax = plt.subplots(1, figsize=figsize)
        self.roadrunner.reset()
        data = self.roadrunner.simulate(start_time, end_time, num_point)
        for col in data.colnames[1:]:
            ax.plot(data[TIME], data[col])
        ax.set_xlabel("time")
        if y_max is None:
            y_max = max((data[:, 1:].flatten()))
        ax.set_ylim([0, y_max])
        ax.legend(data.colnames[1:], bbox_to_anchor=legend_crd,
               loc="upper right")
        if is_plot:
            plt.show()

    def plotLinearApproximation(self, A_mat=None, suptitle="",
          is_plot=True, start_time=START_TIME, end_time=END_TIME,
          figsize=(20, 10), num_point=NUM_POINT):
        """
        Creates a plot that compares the linear approximation with the true model.

        Parameters
        ----------
        A_mat: A matrix of approximation model
            default is Jacobian at current time
        suptitle: str
        is_plot: bool
        start_time: float
        end_time: float
        num_point: int
        """
        rr_df = self.simulateRoadrunner(start_time, end_time, num_point)
        nrow = 1
        ncol = len(rr_df.columns)
        fig, axes = plt.subplots(nrow, ncol, figsize=figsize)
        axes = np.reshape(axes, (nrow, ncol))
        linear_df = self.simulateLinearSystem(timepoint=start_time,
              A_mat=A_mat,
              start_time=start_time, end_time=end_time, num_point=num_point)
        y_min = min(linear_df.min().min(), rr_df.min().min())
        y_max = max(linear_df.max().max(), rr_df.max().max())
        irow = 0
        for icol, column in enumerate(rr_df.columns):
            ax = axes[irow, icol]
            ax.plot(linear_df.index, linear_df[column], color="red")
            ax.plot(rr_df.index, rr_df[column], color="blue")
            ax.set_ylim([y_min, y_max])
            if irow < nrow - 1:
                ax.set_xticklabels([])
            if irow == 0:
                ax.set_title(column, rotation=45)
                if icol == 0:
                    ax.legend(["approximation", "true"])
            if icol > 0:
                ax.set_yticklabels([])
        plt.suptitle(suptitle)
        if is_plot:
            plt.show()
        else:
            plt.close()

    @classmethod
    def evaluateAccuracy(cls, model_reference, timepoints, suptitle="",
         is_plot=True, y_max=None, **kwargs):
        """
        Creates a plot that evaluates the
        accouract of a linear model where the Jacobian is calculated
        at multiple timepoints. X0 is taken at the start of the simulation.

        Parameters
        ----------
        model_reference: as required by constructor
        timepoints: list-float
            Time at which Jacobian is calculated
        suptitle: str
        is_plot: bool
        kwargs: dict
            Values used for simulation (start_time, end_time, num_point)
        """
        if isinstance(timepoints, float) or isinstance(timepoints, int):
            timepoints = [timepoints]
        ctlsb = cls(model_reference)
        rr_df = ctlsb.simulateRoadrunner(**kwargs)
        nrow = len(timepoints)
        ncol = len(rr_df.columns)
        fig, axes = plt.subplots(nrow, ncol, figsize=(15, 5))
        axes = np.reshape(axes, (nrow, ncol))
        for irow, timepoint in enumerate(timepoints):
            linear_df = ctlsb.simulateLinearSystem(timepoint=timepoint,
                  **kwargs)
            y_min = min(linear_df.min().min(), rr_df.min().min())
            if y_max is None:
                y_max = max(linear_df.max().max(), rr_df.max().max())
            for icol, column in enumerate(rr_df.columns):
                ax = axes[irow, icol]
                ax.plot(linear_df.index, linear_df[column], color="red")
                ax.plot(rr_df.index, rr_df[column], color="blue")
                ax.scatter(timepoint, y_min, s=40, marker="o", color="g")
                ax.set_ylim([y_min, y_max])
                if irow < nrow - 1:
                    ax.set_xticklabels([])
                if irow == 0:
                    ax.set_title(column, rotation=45)
                    if icol == 0:
                        ax.text(-3, 0.75*y_max, "Jacobian Time")
                        ax.legend(["linear", "nonlinear"])
                if icol > 0:
                    ax.set_yticklabels([])
                else:
                    ax.text(-2, y_max/2, "%2.1f" % timepoint)
        plt.suptitle(suptitle)
        if is_plot:
            plt.show()
        else:
            plt.close()
