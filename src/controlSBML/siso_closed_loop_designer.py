"""Designs a closed loop SISO System with a PID controller and Filter."""

import controlSBML.constants as cn
from controlSBML.option_management.option_manager import OptionManager
from controlSBML import util
from controlSBML.timeseries import Timeseries
from controlSBML.siso_closed_loop_system import SISOClosedLoopSystem

import control
import lmfit
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MAX_VALUE = 1e3  # Maximum value for a parameter
MIN_VALUE = 0  # Minimum value for a paramete
DEFAULT_INITIAL_VALUE = 1   # Default initial value for a
STEP_SIZE = 1
BELOW_MIN_MULTIPLIER = 1e-3
ABOVE_MAX_MULTIPLIER = 1e-3
PARAM_NAMES = ["kp", "ki", "kd", "kf"]
# Column names
COL_KP = "kp"
COL_KI = "ki"
COL_KD = "kd"
COL_KF = "kf"
COL_RMSE = "rmse"
COL_CLOSED_LOOP_SYSTEM = "closed_loop_system"
COL_CLOSED_LOOP_SYSTEM_TS = "closed_loop_system_ts"
COL_STEP_SIZE = "step_size"



def _calculateClosedLoopTf(sys_tf=None, kp=None, ki=None, kd=None, kf=None):
    # Construct the transfer functions
    is_none = True
    if sys_tf is None:
        raise ValueError("sys_tf must be defined")
    controller_tf = control.tf([0], [1])
    if kp is not None:
        is_none = False
        controller_tf += control.tf([kp], [1])
    if ki is not None:
        is_none = False
        controller_tf += control.tf([ki], [1, 0])
    if kd is not None:
        is_none = False
        controller_tf += control.tf([kd, 0], [1])
    # Filter
    if kf is not None:
        is_none = False
        filter_tf = control.TransferFunction([kf], [1, kf])
    else:
        filter_tf = 1
    if is_none:
        raise ValueError("At least one parameter must be defined")
    # Closed loop transfer function
    return control.feedback(controller_tf*sys_tf, filter_tf)


##################################################################
class SISOClosedLoopDesigner(object):

    def __init__(self, sys_tf, times=None, step_size=STEP_SIZE, is_history=True):
        """
        Args:
            sys_tf: control.TransferFunction (open loop system)
            is_history: bool (if True, then history is maintained)
        """
        self.sys_tf = sys_tf
        self.step_size = step_size
        if times is None:
            self.times = np.linspace(0, 5, 50)
        else:
            self.times = times
        # Internal state
        self.history = _History(self, is_history=is_history)
        # Outputs affected by self.set
        self.kp = None
        self.ki = None
        self.kd = None
        self.kf = None
        self.closed_loop_system = None
        self.closed_loop_system_ts = None
        # Outputs produced by self.design
        self.minimizer_result = None
        self.rmse = None # Root mean square of residuals
        # 
        self.history.add()

    @property
    def closed_loop_tf(self):
        return _calculateClosedLoopTf(sys_tf=self.sys_tf, kp=self.kp, ki=self.ki, kd=self.kd, kf=self.kf)
    
    def set(self, kp=None, ki=None, kd=None, kf=None):
        """
        Sets values of the design parameters

        Args:
            kp (float)
            ki (float)
            kd (float)
            kf (float)
        """
        for name in PARAM_NAMES:
            self.__setattr__(name, eval(name))
        self.closed_loop_system = None
        self.closed_loop_system_ts = None
        self.history.add()

    def get(self):
        """
        Provides design constants with non-None values

        Returns:
            dict: {name: value}
        """
        dct = {}
        for name in PARAM_NAMES:
            if self.__getattribute__(name) is not None:
                dct[name] = self.__getattribute__(name)
        return dct 

    def design(self, min_response=None, max_response=None, kp=False, ki=False, kd=False, kf=False):
        """
        Args:
            times (np.array): time points for the simulation
            min_response (float): minimum response value for the transfer function
            max_response (float): maximum response value for the transfer function
            kp, ki, kd, kf (bool, float): if True, the parameter is fitted. If float, then initial value.
        """
        def _calculateResiduals(params):
            """
            Args:
                params (lmfit.Parameters)
                sys_tf (control.TransferFunction)
            """
            # Calculate the closed loop transfer function
            kwargs = dict(sys_tf=self.sys_tf)
            kwargs.update(params.valuesdict())
            tf = _calculateClosedLoopTf(**kwargs)
            predictions = self.simulate(transfer_function=tf)
            if min_response is not None:
                predictions = np.array([v if v >= min_response else v*BELOW_MIN_MULTIPLIER for v in predictions])
            if max_response is not None:
                predictions = np.array([v if v <= max_response else v*ABOVE_MAX_MULTIPLIER for v in predictions])
            return STEP_SIZE - predictions
        # Construct the parameters
        params = lmfit.Parameters()
        for name in PARAM_NAMES:
            value = eval(name)
            if not isinstance(value, float):
                if not value:
                    continue
                value = DEFAULT_INITIAL_VALUE
            params.add(name, value=value, min=MIN_VALUE, max=MAX_VALUE)
        # Fit the parameters
        minimizer = lmfit.Minimizer(_calculateResiduals, params)
        self.minimizer_result = minimizer.leastsq()
        new_params = lmfit.Parameters()
        for name in PARAM_NAMES:
            if name in self.minimizer_result.params:
                self.__setattr__(name, self.minimizer_result.params[name].value)
                value = self.minimizer_result.params[name].value
                new_params.add(name, value=value, min=MIN_VALUE, max=MAX_VALUE)
        residuals = _calculateResiduals(new_params)
        self.rmse = np.sqrt(np.mean(residuals**2))
        self.history.add()

    def simulate(self, transfer_function=None, times=None):
        """
        Simulates the closed loop transfer function based on the parameters of the object.

        Args
            times (np.array): time points for the simulation
        Returns
            (np.array, np.array): predictions
        Raises
            ValueError: if there are no parameters defined for the closed loop transfer function
        """
        if transfer_function is None:
            transfer_function = self.closed_loop_tf
        if times is None:
            times = self.times
        _, predictions = control.forced_response(transfer_function, T=times, U=self.step_size)
        return predictions
    
    def plot(self, times=None, **kwargs):
        """
        Plots the step response if values are assigned to the closed loop parameters.

        Args:
            kwargs: arguments for OptionManager
        """
        predictions = self.simulate(times=times)
        df = pd.DataFrame({"time": times, "predictions": predictions})
        df["step_size"] = self.step_size
        ts = Timeseries(mat=df)
        if "is_plot" in kwargs:
            is_plot = kwargs["is_plot"]
        else:
            is_plot = True
        kwargs["is_plot"] = False
        plot_result = util.plotOneTS(ts, **kwargs)
        ax = plot_result.ax
        tf_text = util.simplifyTransferFunction(self.closed_loop_tf)
        ax.set_title(tf_text)
        if is_plot:
            plt.show()
        # Title lists values of the design parameters

    def makeNonlinearIOSystemClosedLoop(self, ctlsb, **kwargs):
        """
        Creates a SISOClosedLoopSystem using the parameters of the designer.

        Args:
            ctlsb: ControlSBML
            kwargs: arguments for SISOClosedLoopSystem
        Returns:
            control.Interconnect
        """
        siso = SISOClosedLoopSystem(ctlsb)
        new_kwargs = self.get()
        new_kwargs.update(kwargs)
        siso.makePIDClosedLoopSystem(**new_kwargs)
        self.closed_loop_system = siso.closed_loop_system
        self.closed_loop_system_ts = siso.makeStepResponse(start_time=0, end_time=100, step_size=500)
        self.history.add()
        util.plotOneTS(self.closed_loop_system_ts, markers=["", ""], figsize=(5, 5), xlabel="time", ax2=0)


class _History(object):
    # Maintains history of changes to design choices
    def __init__(self, designer, is_history=True):
        self.designer = designer
        self.is_history = is_history
        self._dct = None
        self.clear()

    def __len__(self):
        first = PARAM_NAMES[0]
        return len(self._dct[first])
    
    def clear(self):
        self._dct = {}
        for name in PARAM_NAMES:
            self._dct[name] = []
        self._dct[COL_CLOSED_LOOP_SYSTEM] = []
        self._dct[COL_CLOSED_LOOP_SYSTEM_TS] = []
        self._dct[COL_STEP_SIZE] = []
        self._dct[COL_RMSE] = []

    def add(self):
        if not self.is_history:
            return
        for name in PARAM_NAMES:
            self._dct[name].append(self.designer.__getattribute__(name))
        self._dct[COL_CLOSED_LOOP_SYSTEM].append(self.designer.closed_loop_system)
        self._dct[COL_CLOSED_LOOP_SYSTEM_TS].append(self.designer.closed_loop_system_ts)
        self._dct[COL_STEP_SIZE].append(self.designer.step_size)
        self._dct[COL_RMSE].append(self.designer.rmse)

    def undo(self):
        _ = self._dct.pop()

    def report(self):
        """
        Creates a dataframe of the history

        Returns:
            pd.DataFrame
        """
        df = pd.DataFrame(self._dct)
        return df

    def get(self, idx):
        """
        Returns the SISOClosedLoopDesigner at the specified index.

        Args:
            idx: int
        Returns:
            SISOClosedLoopDesigner
        """
        if idx > len(self) - 1:
            raise ValueError("idx must be less than %d" % len(self))
        # Construct entries for the desired history element
        dct = {}
        for name in self._dct.keys():
            dct[name] = self._dct[name][idx]
        designer = SISOClosedLoopDesigner(self.designer.sys_tf,
                                          times=self.designer.times,
                                          step_size=STEP_SIZE)
        for name in PARAM_NAMES:
            designer.__setattr__(name, dct[name])
        designer.closed_loop_system = dct[COL_CLOSED_LOOP_SYSTEM]
        designer.closed_loop_system_ts = dct[COL_CLOSED_LOOP_SYSTEM_TS]
        designer.rmse = dct[COL_RMSE]
        return designer