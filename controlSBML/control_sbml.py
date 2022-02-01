"""LTI Control for SBML models"""

"""
TO DO:

1. Plot difference between time jacoabian at reference vs. Current.
2. Plot TOTAL residual SSQ vs. jacobian difference
"""

import control
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tellurium as te


ANT = "ant"
XML = "xml"


TYPE_MODEL = "type_model"  # libsbml model
TYPE_XML = "type_xml"  # XML string
TYPE_ANTIMONY = "type_xml"  # Antimony string
TYPE_FILE = "type_file" # File reference
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
        self.roadrunner = None  # Roadrunner object
        # Read the model
        if "RoadRunner" in str(type(model_reference)):
            self.roadrunner = model_reference
            model_reference = model_reference.getSBML()
        elif isinstance(model_reference, str):
            if model_reference[0:4] == "http":
                self.roadrunner = te.loadSBMLModel(model_reference)
            else:
                parts = model_reference.split(".")
                if len(parts) == 2:
                    if parts[1] == XML:
                        self.roadrunner = te.loadSBMLModel(model_reference)
                    elif parts[1] == ANT:
                        self.roadrunner = te.loadAntimonyModel(model_reference)
                    elif XML in model_reference.count:
                        self.roadrunner = te.loadSBMLModel(model_reference)
                    else:
                        # Assume string for antimony model
                        self.roadrunner = te.loada(model_reference)
                else:
                    self.roadrunner = te.loada(model_reference)
        else:
            raise ValueError("Invalid model reference")
        # Do the initializations
        self.boundary_species = self.roadrunner.getBoundarySpeciesIds()
        self._state_names = None
        self.antimony = self.roadrunner.getAntimony()
        self.roadrunner.reset()

    def _mkBoundarySpeciesFloating(self):
        for name in self.boundary_species:
            self.roadrunner.setBoundary(name, False)

    def _unmkBoundarySpeciesFloating(self):
        for name in self.boundary_species:
            self.roadrunner.setBoundary(name, True)

    @property
    def state_names(self):
        if self._state_names is None:
            self._mkBoundarySpeciesFloating()
            mat = self.roadrunner.getFullJacobian()
            self._state_names = list(mat.colnames)
            self._unmkBoundarySpeciesFloating()
        return self._state_names

    @property
    def jacobian(self):
        """
        Returns
        -------
        NamedArray with names for rows (rownames) and columns (colnames)
        Handle boundary species.
        """
        self._mkBoundarySpeciesFloating()
        mat = self.roadrunner.getFullJacobian()
        self._unmkBoundarySpeciesFloating()
        for idx, name in enumerate(mat.rownames):
            if name in self.boundary_species:
                mat[idx, :] = 0
        df = pd.DataFrame(mat, columns=mat.colnames, index=mat.rownames)
        return df

    @property
    def currentState(self):
        """
        Contructs vector of current state values (floating and boundary species)

        Returns
        -------
        Series
            index: str (state names)
        """
        values = list(self.get(self.state_names).values())
        return pd.Series(values, index=self.state_names)

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
              in zip(self.state_names, self.jacobian.columns)])
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

    def mkStateSpace(self, B=None, C=None, D=None):
        """
        Creates a control system object for the n X n jacobian.

        Parameters
        ----------
        B: np.array(n X p)
        C: np.array(q X n)
        D: np.array(q X p)

        Returns
        -------
        control.StateSpace
        """
        # Construct the matrices
        A = self.jacobian.values
        if B is None:
            B = np.repeat(0, A.shape[0])
            B = np.reshape(B, (A.shape[0], 1))
        if C is None:
            C = np.identity(A.shape[0])
        if D is None:
            D = B
        return control.StateSpace(A, B, C, D)

    def simulateLinearSystem(self, timepoint=0, start_time=0,
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
        sys = self.mkStateSpace()
        self.setTime(start_time)
        x0 = np.array(list(self.get(self.state_names).values()))
        x0 = np.reshape(x0, len(x0))
        self.setTime(cur_time)  # Restore the time
        # Run the linear simulation
        dt = (end_time - start_time)/num_point
        times = [n*dt for n in range(num_point)]
        times, y_vals = control.forced_response(sys, T=times, X0=x0)
        df = pd.DataFrame(y_vals.transpose(), index=times)
        df.columns = self.state_names
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

    @classmethod
    def evaluateAccuracy(cls, model_reference, timepoints, suptitle="",
         is_plot=True, **kwargs):
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
        
        Returns
        -------
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
            ymin = min(linear_df.min().min(), rr_df.min().min())
            ymax = max(linear_df.max().max(), rr_df.max().max())
            for icol, column in enumerate(rr_df.columns):
                ax = axes[irow, icol]
                ax.plot(linear_df.index, linear_df[column], color="red")
                ax.plot(rr_df.index, rr_df[column], color="blue")
                ax.scatter(timepoint, ymin, s=40, marker="o", color="g")
                ax.set_ylim([ymin, ymax])
                if irow < nrow - 1:
                    ax.set_xticklabels([])
                if irow == 0:
                    ax.set_title(column)
                    if icol == 0:
                        ax.text(-3, 0.75*ymax, "Jacobian Time")
                        ax.legend(["linear", "nonlinear"])
                if icol > 0:
                    ax.set_yticklabels([])
                else:
                    ax.text(-2, ymax/2, "%2.1f" % timepoint)
        plt.suptitle(suptitle)
        if is_plot:
            plt.show()
        else:
            plt.close()
