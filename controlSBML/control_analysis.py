"""LTI Control for SBML models. ControlAnalysis manipulates data and simulates."""

from controlSBML.control_base import ControlBase
import controlSBML.constants as cn
from controlSBML.options import Options

import control
from docstring_expander.expander import Expander
import pandas as pd


class ControlAnalysis(ControlBase):

    @staticmethod
    def _getSimulationParameters(**sim_opts):
        start_time = sim_opts[cn.O_START_TIME]
        end_time = sim_opts[cn.O_END_TIME]
        points_per_time = sim_opts[cn.O_POINTS_PER_TIME]
        num_point = int(points_per_time*(end_time - start_time)) + 1
        return start_time, end_time, num_point

    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def simulateLinearSystem(self, A_mat=None, B_mat=None, C_mat=None,
          timepoint=0, is_reduced=True, **kwargs):
        """
        Creates an approximation of the SBML model based on the Jacobian, and
        constructs predictions based on this Jacobian and the values of
        floating species at the start_time.

        Parameters
        ----------
        A_mat: np.array (n X n)
        B_mat: np.array (n X p)
        C_mat: np.array (q X n)
        timepoint: float
            Time at which Jacobian is taken
        #@expand

        Returns
        -------
        pd.dataframe
            columns: floating species
            index: time
        """
        options = Options(kwargs, [cn.SIM_DCT])
        sim_opts = options.parse()[0]
        start_time, end_time, num_point = self._getSimulationParameters(**sim_opts)
        cur_time = self.get(cn.TIME)
        self.setTime(timepoint)
        sys = self.makeStateSpace(A_mat=A_mat, B_mat=B_mat, C_mat=C_mat,
              is_reduced=is_reduced)
        self.setTime(start_time)
        x0 = self.getCurrentState(is_reduced=is_reduced)
        self.setTime(cur_time)  # Restore the time
        # Run the linear simulation
        dt = (end_time - start_time)/(num_point - 1)
        times = [start_time + n*dt for n in range(num_point)]
        times, y_vals = control.forced_response(sys, T=times, X0=x0)
        df = pd.DataFrame(y_vals.transpose(), index=times)
        df.columns = self.getSpeciesNames(is_reduced=is_reduced)
        return df

    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def simulateRoadrunner(self, **kwargs):
        """
        Runs a new roadrunner simulation.

        Parameters
        ----------
        #@expand

        Returns
        -------
        pd.dataframe
            columns: floating species
            index: time
        """
        options = Options(kwargs, [cn.SIM_DCT])
        sim_opts = options.parse()[0]
        start_time, end_time, num_point = self._getSimulationParameters(
              **sim_opts)
        #
        self.roadrunner.reset()
        data = self.roadrunner.simulate(start_time, end_time, num_point)
        columns = [c[1:-1] if c[0] =="[" else c for c in data.colnames]
        df = pd.DataFrame(data, columns=columns)
        df = df.set_index(cn.TIME)
        return df
