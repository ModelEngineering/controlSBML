"""LTI Control for SBML models. ControlAnalysis manipulates data and simulates."""

from controlSBML.timeseries import Timeseries
from controlSBML.control_base import ControlBase
import controlSBML.constants as cn
from controlSBML.option_management.options import Options

import control
from docstring_expander.expander import Expander
import pandas as pd
import numpy as np


# Simulation parameters
class SimulationParameters(object):
    # Provides convenient access to simulation parameters in this module.

    def __init__(self, sim_opts):
        self.start_time = sim_opts[cn.O_START_TIME]
        self.end_time = sim_opts[cn.O_END_TIME]
        self.step_val = sim_opts[cn.O_STEP_VAL]
        self.A_df = sim_opts[cn.O_A_DF]
        self.B_df = sim_opts[cn.O_B_DF]
        self.C_df = sim_opts[cn.O_C_DF]
        self.points_per_time = sim_opts[cn.O_POINTS_PER_TIME]
        self.num_point = int(self.points_per_time*(
              self.end_time - self.start_time)) + 1


class ControlAnalysis(ControlBase):

    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def simulateLinearSystem(self, time=0, **kwargs):
        """
        Creates an approximation of the SBML model based on the Jacobian, and
        constructs predictions based on this Jacobian and the values of
        floating species at the start_time.

        Parameters
        ----------
        time: float
            Time at which Jacobian is taken
        #@expand

        Returns
        -------
        Timeseries
        """
        options = Options(kwargs, [cn.SIM_DCT])
        sim_opts = options.parse()[0]
        prms = SimulationParameters(sim_opts)
        cur_time = self.get(cn.TIME)
        self.setTime(time)
        if prms.B_df is None:
            prms.B_df = self._makeBDF()
        sys =self.makeStateSpace(A_mat=prms.A_df, B_mat=prms.B_df,
              C_mat=prms.C_df)
        X0 = self.state_ser
        self.setTime(prms.start_time)
        X0_vals = self.state_ser
        for name in X0.index:
            X0.loc[name] = X0_vals.loc[name]
        X0 = X0.values
        self.setTime(cur_time)  # Restore the time
        # Run the linear simulation
        dt = (prms.end_time - prms.start_time)/(prms.num_point - 1)
        times = [prms.start_time + n*dt for n in range(prms.num_point)]
        times, y_vals = control.forced_response(sys, T=times,
              X0=X0, U=prms.step_val)
        ts = Timeseries(np.transpose(y_vals),
              times=times, columns=self.output_names)
        return ts

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
        prms = SimulationParameters(sim_opts)
        #
        self.roadrunner.reset()
        data = self.roadrunner.simulate(prms.start_time, prms.end_time,
              prms.num_point)
        columns = [c[1:-1] if c[0] =="[" else c for c in data.colnames]
        return Timeseries(data, columns=columns)
