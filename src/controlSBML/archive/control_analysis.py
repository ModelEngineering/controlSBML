"""LTI Control for SBML models. ControlAnalysis manipulates data and simulates."""

from controlSBML.timeseries import Timeseries
from controlSBML.control_base import ControlBase
import controlSBML.constants as cn
from controlSBML import util
from controlSBML.option_management.options import Options

import control
from docstring_expander.expander import Expander
import pandas as pd
import numpy as np


class ControlAnalysis(ControlBase):

    class SimulationParameters(object):
    # Provides convenient access to simulation parameters in this module.

        def __init__(self, times=None, **sim_opts):
            self.start_time = sim_opts[cn.O_START_TIME]
            self.end_time = sim_opts[cn.O_END_TIME]
            self.step_val = sim_opts[cn.O_STEP_VAL]
            self.points_per_time = sim_opts[cn.O_POINTS_PER_TIME]
            self.num_point = int(self.points_per_time*(
                self.end_time - self.start_time)) + 1
            if times is None:
                self.times = util.makeSimulationTimes(self.start_time, self.end_time, self.num_point)
            else:
                self.times = list(times)

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
        prms = self.SimulationParameters(**sim_opts)
        #
        self.roadrunner.reset()
        data = self.roadrunner.simulate(prms.start_time, prms.end_time,
              prms.num_point)
        columns = [c[1:-1] if c[0] =="[" else c for c in data.colnames]
        return Timeseries(data, columns=columns)
    
    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def simulateLinearSystem(self, transfer_function, input_name, output_name, 
                             op_time=None, times=None, U=None, **kwargs):
        """
        Simulates the linear system described by the transfer function

        Parameters
        ----------
        transfer_function: control.TransferFunction
        input_name: str
        output_name: str
        op_time: float (time at which operating point values are obtained); if None, use steady state
        times: np.array (times for simulation)
        U: array (overrides step_val in sim_opts)
        #@expand

        Returns
        -------
        Timeseries
        """
        options = Options(kwargs, [cn.SIM_DCT])
        sim_opts = options.parse()[0]
        prms = self.SimulationParameters(times=times, **sim_opts)
        if U is None:
            U = float(prms.step_val)
        if op_time is None:
            if not self.setSteadyState():
                # If fail to set to steady state, then use the start time
                self.setTime(prms.start_time)
        else:
            self.setTime(op_time)
        initial_output = self.get(output_name)
        # Run the linear simulation
        self.setTime(prms.start_time)
        if isinstance(U, float) or isinstance(U, int):
            U = np.repeat(float(U), len(prms.times))
        _, y_vals = control.forced_response(transfer_function, T=prms.times, U=U)
        # Impulse response for initial condition of the output
        _, yi_vals = control.impulse_response(transfer_function, T=prms.times)
        ys = y_vals + yi_vals*initial_output
        df = pd.DataFrame({input_name: U, output_name: ys})
        ts = Timeseries(df, times=prms.times)
        return ts

    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def simulateMIMOLinearSystem(self, transfer_function_df,
                             op_time=None, times=None, U=None, **kwargs):
        """
        Simulates a MIMO systems specified as multiple SISO transfer functions.

        Parameters
        ----------
        transfer_function_df: pd.DataFrame / control.TransferFunction
            columns: output species
            index: input species
            values: control.TransferFunction
        op_time: float (time at which operating point values are obtained); if None, use steady state
        times: np.array (times for simulation)
        U: array (overrides step_val in sim_opts) or dictionary by input species
        #@expand

        Returns
        -------
        dict
            keys: input_name
            values: Timeseries with columnns for outputs
        """
        if isinstance(U, dict):
            U_dct = dict(U)
        else:
            U_dct = {n: U for n in transfer_function_df.index}
        # Do the simulations
        result_dct = {}
        for input_name in transfer_function_df.index:
            full_ts = None
            for output_name in transfer_function_df.columns:
                tf = transfer_function_df.loc[input_name, output_name]
                ts = self.simulateLinearSystem(tf, input_name=input_name,
                      output_name=output_name,
                      op_time=op_time, times=times, U=U_dct[input_name], **kwargs)
                if full_ts is None:
                    full_ts = ts
                else:
                    full_ts[output_name] = ts[output_name]
            result_dct[input_name] = full_ts
        return result_dct