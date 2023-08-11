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


# Simulation parameters
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
            self.times = times


class ControlAnalysis(ControlBase):

    @Expander(cn.KWARGS, cn.SIM_KWARGS)
    def simulateLinearSystem(self, transfer_function_df, input_name=None, output_name=None, 
                             op_time=None, times=None, U=None, **kwargs):
        """
        Creates an approximation of the SBML model using system identification to construt
        a transfer function.

        Parameters
        ----------
        transfer_function_df: pd.DataFrame / control.TransferFunction
            columns: output species
            index: input species
            values: control.TransferFunction
          if control.TransferFunction, then input_names and output_names must be specified
        input_name: str (name of input species if argument is a transfer function)
        output_name: str (name of output species if argument is a transfer function)
        op_time: float (time at which operating point values are obtained); if None, use steady state
        times: np.array (times for simulation)
        U: array (overrides step_val in sim_opts)
        #@expand

        Returns
        -------
        dict
            keys: input species
            values: Timeseries with columnns for outputs
        """
        if isinstance(transfer_function_df, control.TransferFunction):
            transfer_function_df = pd.DataFrame({output_name: [transfer_function_df]})
            transfer_function_df.index = [input_name]
            input_names = [input_name]
            output_names = [output_name]
        else:
            input_names = transfer_function_df.index
            output_names = transfer_function_df.columns
        options = Options(kwargs, [cn.SIM_DCT])
        sim_opts = options.parse()[0]
        prms = SimulationParameters(times=times, **sim_opts)
        if U is None:
            U = prms.step_val
        if op_time is None:
            self.setSteadyState()
        else:
            self.setTime(op_time)
        self.setTime(prms.start_time)
        # Run the linear simulation
        result_dct = {}
        for in_name in input_names:
            timeseries_dct = {}
            for out_name in output_names:
                tf = transfer_function_df.loc[in_name, out_name]
                _, y_vals = control.forced_response(tf, T=prms.times, U=U)
                # Impulse response for initial condition of the output
                _, yi_vals = control.impulse_response(tf, T=prms.times)
                y_vals = y_vals + yi_vals*self.get(out_name)
                timeseries_dct[out_name] = np.transpose(y_vals)
            result_dct[input_name] = Timeseries(pd.DataFrame(timeseries_dct), times=times)
        return result_dct

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
        prms = SimulationParameters(**sim_opts)
        #
        self.roadrunner.reset()
        data = self.roadrunner.simulate(prms.start_time, prms.end_time,
              prms.num_point)
        columns = [c[1:-1] if c[0] =="[" else c for c in data.colnames]
        return Timeseries(data, columns=columns)
