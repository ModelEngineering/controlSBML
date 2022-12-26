"""Simulate a system that may consist of subsystems."""

import controlSBML.constants as cn
import controlSBML as ctl
from controlSBML import util

import control
import numpy as np
import pandas as pd


def makeStateVector(sys, start_time=0):
    """
    Constructs the initial state vector recursively.

    Parameters
    ----------
    sys: inherits from control.InputOutputSystem
    start_time: float

    Returns
    -------
    list
    """
    x_lst = []
    if "InterconnectedSystem" in str(type(sys)):
        for sub_sys in sys.syslist:
            x_lst.extend(makeStateVector(sub_sys, start_time=start_time))
    elif isinstance(sys, ctl.NonlinearIOSystem):
        x_lst.extend(sys.makeStateSer().values)
    else:
        new_state = list(np.repeat(0, sys.nstates))
        x_lst.extend(new_state)
    result = [float(v) for v in x_lst]
    return result

def simulateSystem(sys, output_names=None, initial_x_vec=None, u_vec=None,
      start_time=cn.START_TIME, end_time=cn.END_TIME,
      points_per_time=cn.POINTS_PER_TIME):
    """
    Simulates the system. Provides default values for initial state by
    setting control systems to 0 and setting controlSBML.NonlinearIOSystem
    to their value time 0.

    Parameters
    ----------
    sys: inherits from control.InputOutputSystem
    output:names: list-str (names of the outputs)
    initial_x_vec: np.array or pd.Series
    u_vec: np.array or pd.Series
    start_time: float (start of simulation)
    end_time: float (end of simulation)
    points_per_time: int (number of simulation timepoints per second)

    Returns
    -------
    Timeseries
    """
    # Construc initial state vector if necesary
    if initial_x_vec is None:
        initial_x_vec = np.array(makeStateVector(sys, start_time=start_time))
    #
    times = util.makeSimulationTimes(start_time=start_time,
          end_time=end_time, points_per_time=points_per_time)
    if u_vec is not None:
        results = control.input_output_response(sys, times, X0=initial_x_vec,
            U=u_vec)
    else:
        results = control.input_output_response(sys, times, X0=initial_x_vec)
    output_mat = np.transpose(results.y)
    num_column = np.shape(output_mat)[1]
    is_int_columns = False
    if output_names is None:
        if isinstance(sys, ctl.NonlinearIOSystem):
            output_names = sys.output_names
        else:
            output_names = [str(n) for n in range(num_column)]
            is_int_columns = True
    ts = ctl.Timeseries(output_mat, columns=output_names, times=results.t)
    if is_int_columns:
        ts.columns = [int(c) for c in ts.columns]
    return ts
