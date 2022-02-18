"""LTI Control for SBML models. ControlAnalysis manipulates data and simulates."""

from controlSBML.control_base import ControlBase
import controlSBML.constants as cn
from controlSBML.options import Options

import control
import numpy as np
import pandas as pd


class ControlAnalysis(ControlBase):

    @staticmethod
    def _getSimulationParameters(**sim_opts):
        start_time = sim_opts[cn.O_START_TIME]
        end_time = sim_opts[cn.O_END_TIME]
        points_per_time = sim_opts[cn.O_POINTS_PER_TIME]
        num_point = int(points_per_time*(end_time - start_time)) + 1
        return start_time, end_time, num_point

    def simulateLinearSystem(self, A_mat=None, timepoint=0, **kwargs):
        """
        Creates an approximation of the SBML model based on the Jacobian, and
        constructs predictions based on this Jacobian and the values of
        floating species at the start_time.

        Parameters
        ----------
        kwargs: dict
            cn.SIM_OPTS
        
        Returns
        -------
        pd.dataframe
            columns: floating species
            index: time
        """
        options = Options(dct=kwargs, default_dcts=[cn.SIM_DCT])
        sim_opts = options.parse()[0]
        start_time, end_time, num_point = self._getSimulationParameters(**sim_opts)
        cur_time = self.get(cn.TIME)
        self.setTime(timepoint)
        sys = self.makeStateSpace(A=A_mat)
        self.setTime(start_time)
        x0 = self.current_state
        self.setTime(cur_time)  # Restore the time
        # Run the linear simulation
        dt = (end_time - start_time)/(num_point - 1)
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
            cn.SIM_OPTS
        
        Returns
        -------
        pd.dataframe
            columns: floating species
            index: time
        """
        options = Options(dct=kwargs, default_dcts=[cn.SIM_DCT])
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
