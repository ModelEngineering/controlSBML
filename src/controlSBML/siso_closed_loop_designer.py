"""
Designs a closed loop SISO System with a PID controller and Filter.

The design is done using transfer functions and is most appropriate if there is a fixed input species.
"""

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
LOWPASS_POLE = 1e4 # Pole for low pass filter
# Column names
COL_KP = "kp"
COL_KI = "ki"
COL_KD = "kd"
COL_KF = "kf"
COL_RESIDUAL_RMSE = "residual_rmse"
COL_CLOSED_LOOP_SYSTEM = "closed_loop_system"
COL_CLOSED_LOOP_SYSTEM_TS = "closed_loop_system_ts"
COL_STEP_SIZE = "step_size"


##################################################################
def _calculateClosedLoopTf(sys_tf=None, kp=None, ki=None, kd=None, kf=None, sign=-1):
    # Construct the transfer functions
    controller_tf = util.makePIDTransferFunction(kp=kp, ki=ki, kd=kd)
    # Filter
    if kf is not None:
        is_none = False
        filter_tf = control.TransferFunction([kf], [1, kf])
    else:
        filter_tf = 1
    # Closed loop transfer function
    forward_tf = sys_tf*controller_tf
    final_tf = control.feedback(forward_tf, filter_tf, sign=sign)
    # Ensure that the resulting transfer function is proper
    if len(final_tf.num[0][0]) == len(final_tf.den[0][0]):
        lowpass_tf = control.tf([LOWPASS_POLE], [1, LOWPASS_POLE])
        forward_tf = lowpass_tf*forward_tf
        final_tf = control.feedback(forward_tf, filter_tf, sign=sign)
    return final_tf


##################################################################
class SISOClosedLoopDesigner(object):

    def __init__(self, sys_tf, times=None, step_size=STEP_SIZE, is_history=True, sign=-1):
        """
        Args:
            sys_tf: control.TransferFunction (open loop system)
            is_history: bool (if True, then history is maintained)
            sign: int (if -1, then the residuals are multiplied by -1)
        """
        self.sys_tf = sys_tf
        self.step_size = step_size
        if times is None:
            self.times = np.linspace(0, 5, 50)
        else:
            self.times = times
        self.sign = sign
        # Internal state
        self.history = _History(self, is_history=is_history)
        # Outputs
        self.kp = None
        self.ki = None
        self.kd = None
        self.kf = None
        self.closed_loop_system = None
        self.closed_loop_system_ts = None
        self.residual_rmse = None
        self.minimizer_result = None
        #
        self._initializeDesigner()
        self.siso = None # SISOClosedLoopSystem
        # 
        self.history.add()

    @property
    def closed_loop_tf(self):
        return _calculateClosedLoopTf(sys_tf=self.sys_tf, kp=self.kp, ki=self.ki, kd=self.kd,
                                      kf=self.kf, sign=self.sign)
    
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
            value = eval(name)
            self.__setattr__(name, value)
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
                value = self.__getattribute__(name)
                dct[name] = value
        return dct 
    
    def _initializeDesigner(self):
        self.minimizer_result = None
        self.residual_rmse = None # Root mean square of residuals
        self.kp = None
        self.ki = None
        self.kd = None
        self.kf = None

    def design(self, min_response=None, max_response=None, kp=False, ki=False, kd=False, kf=False,
               residual_precision=5):
        """
        Args:
            times (np.array): time points for the simulation
            min_response (float): minimum response value for the transfer function
            max_response (float): maximum response value for the transfer function
            kp, ki, kd, kf (bool, float): if True, the parameter is fitted. If float, then initial value.
            residual_precision: int (number of decimal places for the residuals)
        """
        sample_period = np.mean(np.diff(self.times))
        self._initializeDesigner()
        def _calculateResiduals(params):
            """
            Args:
                params (lmfit.Parameters)
            """
            # Calculate the closed loop transfer function
            kwargs = dict(sys_tf=self.sys_tf)
            kwargs.update(params.valuesdict())
            kwargs.update({"sign": self.sign})
            tf = _calculateClosedLoopTf(**kwargs)
            dtf = control.c2d(tf, Ts=sample_period)
            _, predictions = self.simulate(transfer_function=dtf)
            if min_response is not None:
                predictions = np.array([v if v >= min_response else v*BELOW_MIN_MULTIPLIER for v in predictions])
            if max_response is not None:
                predictions = np.array([v if v <= max_response else v*ABOVE_MAX_MULTIPLIER for v in predictions])
            residuals = self.step_size - predictions
            residuals = np.round(residuals, residual_precision)
            return residuals
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
            if name in self.minimizer_result.params.valuesdict().keys():
                value = self.minimizer_result.params.valuesdict()[name]
                self.__setattr__(name, value)
                new_params.add(name, value=value, min=MIN_VALUE, max=MAX_VALUE)
        residuals = _calculateResiduals(new_params)
        self.residual_rmse = np.sqrt(np.mean(residuals**2))
        self.history.add()

    def simulate(self, transfer_function=None, period=None, times=None):
        """
        Simulates the closed loop transfer function based on the parameters of the object.

        Args
            transfer_function (control.TransferFunction): closed loop transfer function
            times (np.array): time points for the simulation
        Returns
            (np.array, np.array): times, predictions
        Raises
            ValueError: if there are no parameters defined for the closed loop transfer function
        """
        if transfer_function is None:
            transfer_function = self.closed_loop_tf
        if times is None:
            times = self.times
        if period is not None:
            U = np.sin(2*np.pi*times/period)
        else:
            U = np.repeat(1, len(times))
        U = U*self.step_size
        new_times, predictions = control.forced_response(transfer_function, T=times, U=U)
        return new_times, predictions
    
    def plot(self, times=None, period=None, **kwargs):
        """
        Plots the step response if values are assigned to the closed loop parameters.

        Args:
            kwargs: arguments for OptionManager
        """
        new_times, predictions = self.simulate(times=times, period=period)
        df = pd.DataFrame({"time": new_times, "predictions": predictions})
        df["step_size"] = self.step_size
        ts = Timeseries(mat=df)
        if "is_plot" in kwargs:
            is_plot = kwargs["is_plot"]
        else:
            is_plot = True
        kwargs["is_plot"] = False
        plot_result = util.plotOneTS(ts, **kwargs)
        ax = plot_result.ax
        param_dct = self.get()
        text = ["%s=%f " % (name, param_dct[name]) for name in param_dct.keys()]
        ax.set_title(" ".join(text))
        plt.rcParams['axes.titley'] = 0.9    # y is in axes-relative coordinates.
        if is_plot:
            plt.show()
        # Title lists values of the design parameters

    def evaluateNonlinearIOSystemClosedLoop(self, ctlsb, times=None, step_size=STEP_SIZE, 
                                            is_fixed_input_species=False,
                                            period=0, is_plot=True, **kwargs):
        """
        Creates a SISOClosedLoopSystem using the parameters of the designer.

        Args:
            ctlsb: ControlSBML
            times: np.array (times for simulation)
            step_size: float (step size for simulation)
            is_fixed_input_species: bool (if True, then the input species are fixed)
            kwargs: plot options
        Returns:
            control.Interconnect
        """
        if times is None:
            times = self.times
        start_time = times[0]
        end_time = times[-1]
        self.siso = SISOClosedLoopSystem(ctlsb, is_fixed_input_species=is_fixed_input_species)
        param_dct = {n: 0 for n in PARAM_NAMES}
        param_dct.update(self.get())
        self.siso.makePIDClosedLoopSystem(**param_dct)
        self.closed_loop_system = self.siso.closed_loop_system
        self.closed_loop_system_ts = self.siso.makeResponse(start_time=start_time, end_time=end_time, period=period,
                                                           step_size=step_size)
        self.history.add()
        plot_result = util.plotOneTS(self.closed_loop_system_ts, markers=["", ""],
                                     xlabel="time", ax2=0, is_plot=False, **kwargs)
        param_dct = self.get()
        text = ["%s=%f " % (name, param_dct[name]) for name in param_dct.keys()]
        plot_result.ax.set_title(" ".join(text))
        if is_plot:
            plt.show()


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
        self._dct[COL_RESIDUAL_RMSE] = []

    def add(self):
        if not self.is_history:
            return
        for name in PARAM_NAMES:
            self._dct[name].append(self.designer.__getattribute__(name))
        self._dct[COL_CLOSED_LOOP_SYSTEM].append(self.designer.closed_loop_system)
        self._dct[COL_CLOSED_LOOP_SYSTEM_TS].append(self.designer.closed_loop_system_ts)
        self._dct[COL_STEP_SIZE].append(self.designer.step_size)
        self._dct[COL_RESIDUAL_RMSE].append(self.designer.residual_rmse)

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
        designer.residual_rmse = dct[COL_RESIDUAL_RMSE]
        designer.history.add()
        return designer