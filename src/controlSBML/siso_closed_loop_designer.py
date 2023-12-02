"""
Designs a closed loop SISO System with a PID controller and Filter.

The design is done using transfer functions and is most appropriate if there is a fixed input species.
"""

import controlSBML.constants as cn
from controlSBML.option_management.option_manager import OptionManager
from controlSBML import util
from controlSBML import msgs 
from controlSBML.timeseries import Timeseries
from controlSBML.grid import Grid
from controlSBML.siso_design_evaluator import SISODesignEvaluator

import control
import lmfit
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MAX_VALUE = 1e3  # Maximum value for a parameter
MIN_VALUE = 0  # Minimum value for a paramete
DEFAULT_INITIAL_VALUE = 1   # Default initial value for a parameter
SETPOINT = 1
BELOW_MIN_MULTIPLIER = 1e-3
ABOVE_MAX_MULTIPLIER = 1e-3
PARAM_NAMES = ["kp", "ki", "kf"]
LOWPASS_POLE = 1e4 # Pole for low pass filter
# Column names
COL_RESIDUAL_MSE = "residual_rmse"
COL_CLOSED_LOOP_SYSTEM = "closed_loop_system"
COL_CLOSED_LOOP_SYSTEM_TS = "closed_loop_system_ts"
COL_SETPOINT = "setpoint"


##################################################################
def _calculateClosedLoopTransferFunction(open_loop_transfer_function=None, kp=None, ki=None, kd=None, kf=None, sign=-1):
    # Construct the transfer functions
    if open_loop_transfer_function is None:
        return None
    controller_tf = util.makePIDTransferFunction(kp=kp, ki=ki, kd=kd)
    # Filter
    if kf is not None:
        filter_tf = control.TransferFunction([kf], [1, kf])
    else:
        filter_tf = 1
    # Closed loop transfer function
    forward_tf = open_loop_transfer_function*controller_tf
    final_tf = control.feedback(forward_tf, filter_tf, sign=sign)
    # Ensure that the resulting transfer function is proper
    if len(final_tf.num[0][0]) == len(final_tf.den[0][0]):
        lowpass_tf = control.tf([LOWPASS_POLE], [1, LOWPASS_POLE])
        forward_tf = lowpass_tf*forward_tf
        final_tf = control.feedback(forward_tf, filter_tf, sign=sign)
    return final_tf


##################################################################
class SISOClosedLoopDesigner(object):

    def __init__(self, system, open_loop_transfer_function, times=None, setpoint=SETPOINT, is_steady_state=False,
                 is_history=True, sign=-1, input_name=None, output_name=None):
        """
        Constructs a SISOClosedLoopDesigner object. If the system has more than one input (output) and not input (output) is specified),
        then the first input (output) is used.
        Args:
            sbml_system: SBMLSystem
            sys_tf: control.TransferFunction (open loop system)
            is_steady_state: bool (if True, then the steady state is used)
            is_history: bool (if True, then history is maintained)
            sign: int (-1: negative feedback; +1: positive feedback)
            input_name: str (name of input species)
            output_name: str (name of output species)
        """
        self.system = system
        self.open_loop_transfer_function = open_loop_transfer_function
        self.setpoint = setpoint
        self.is_steady_state = is_steady_state
        if times is None:
            self.times = np.linspace(0, 5, 50)
        else:
            self.times = times
        if input_name is None:
            input_name = self.system.input_names[0]
        if output_name is None:
            output_name = self.system.output_names[0]
        self.input_name = input_name
        self.output_name = output_name
        # Calculated properties
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        self.num_point = len(self.times)
        self.sign = sign
        # Internal state
        self.history = _History(self, is_history=is_history)
        # Outputs
        self.kp = None
        self.ki = None
        self.kf = None
        self.closed_loop_system = None
        self.residual_mse = None
        self.minimizer_result = None
        #
        self._initializeDesigner()
        self.siso = None # SISOClosedLoopSystem
        # 
        self.history.add()

    @property
    def closed_loop_transfer_function(self):
        if self.open_loop_transfer_function is None:
            return None
        return _calculateClosedLoopTransferFunction(open_loop_transfer_function=self.open_loop_transfer_function, kp=self.kp, ki=self.ki,
                                      kf=self.kf, sign=self.sign)
    
    @property
    def closed_loop_timeseries(self):
        _, closed_loop_ts = self.simulateTransferFunction(transfer_function=self.closed_loop_transfer_function)
        return closed_loop_ts
    
    def set(self, kp=None, ki=None, kf=None):
        """
        Sets values of the design parameters

        Args:
            kp (float)
            ki (float)
            kd (float)
            kf (float)
        """
        value_dct = {cn.CP_KP: kp, "ki": ki, "kf": kf}
        for name, value in value_dct.items():
            if value is None:
                continue
            self.__setattr__(name, value)
        self.closed_loop_system = None
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
        self.residual_mse = None # Root mean square of residuals
        self.kp = None
        self.ki = None
        self.kd = None
        self.kf = None

    def _searchForFeasibleClosedLoopSystem(self, value_dct, fixeds, max_iteration=10):
        """
        Does a greedy search to find values of the parameters that are minimally stable.

        Args:
            value_dct: dict (name: value)
                key: name of parameter
                value: value of parameter or None
            fixeds: list-str (parameters whose values don't change)
            max_iteration: int (maximum number of iterations)
        Returns:
            dict: {name: value} or None (no stable result found)
                key: name of parameter
                value: value of parameter or None
        """
        evaluator = SISODesignEvaluator(self.system, input_name=self.input_name,
                                    output_name=self.output_name, setpoint=self.setpoint, times=self.times)
        MINIMAL_FACTOR = 0.01
        factor = 0.5
        def mult(dct, factor):
            new_dct = {}
            for name, value in dct.items():
                if name in fixeds:
                    new_dct[name] = value
                elif value is None:
                    new_dct[name] = None
                else:
                    new_dct[name] = factor*value
            return new_dct
        # Iterate to find values
        dct = dict(value_dct)
        last_stable_dct = None
        for _ in range(max_iteration):
            if evaluator.evaluate(**dct):
                if factor < MINIMAL_FACTOR:
                    break
                else:
                    # Try a bigger factor
                    factor = 1 + factor
                    last_stable_dct = dict(dct)
            else:
                if last_stable_dct is not None:
                    dct = dict(last_stable_dct)
                    break
                else:
                    factor = factor/2
            dct = mult(dct, factor)
        # Return the result
        if last_stable_dct is None:
            return None
        return dct

    def design(self, kp_spec=False, ki_spec=False, kf_spec=False, max_iteration=10,
               num_restart=5, min_value=MIN_VALUE, max_value=MAX_VALUE, is_grid_search=True,
               num_coordinate=3):
        """
        Design objective: Create a feasible system (stable, no negative inputs/outputs) that minimizes residuals.
        Args:
            kp_spec, ki_spec, kf_spec (bool, float): if True, the parameter is fitted. If float, then keeps at this value.
            num_restart: int (number of times to restart the minimizer)
            min_value: float/dict (parameter name: value)
            max_value: float/dict (parameter name: value)
            is_grid_search: restarts are done in different regions of the parameter space
            num_coordinate: int (number of coordinates for a parameter; minimum is 2)
        """
        def addAxis(grid, parameter_name, parameter_spec):
            """
            Adds the grid access based on the parameter specification.

            Args:
                grid: Grid
                parameter_name: str
                parameter_spec: None/bool/float
            """
            if parameter_spec is None:
                return
            if isinstance(parameter_spec, bool):
                if parameter_spec:
                    grid.addAxis(parameter_name, min_value=min_value, max_value=max_value, num_coordinate=num_coordinate)
                return
            if isinstance(parameter_spec, float):
                grid.addAxis(parameter_name, min_value=parameter_spec, max_value=parameter_spec, num_coordinate=1)
        #
        # Initial check
        if self.open_loop_transfer_function is not None:
            if (not util.isStablePoles(self.open_loop_transfer_function)) and (not util.isStableZeros(self.open_loop_transfer_function)):
                msg = "The open loop transfer function has unstable poles and zeros. Design may fail."
                msgs.warn(msg)
        # Initializations
        evaluator = SISODesignEvaluator(self.system, input_name=self.input_name,
                                    output_name=self.output_name, setpoint=self.setpoint, times=self.times)
        grid = Grid()
        addAxis(grid, cn.CP_KP, kp_spec)
        addAxis(grid, cn.CP_KI, ki_spec)
        addAxis(grid, cn.CP_KF, kf_spec)
        # Iterate to find values
        for point in grid.points:
            evaluator.evaluate(**point)
        # Record the result
        self.residual_mse = evaluator.residual_mse
        self.set(kp=evaluator.kp, ki=evaluator.ki, kf=evaluator.kf)
        if self.residual_mse is None:
            return msgs.warn("Could not find a feasible design.")
        else:
            self.history.add()

    def simulateTransferFunction(self, transfer_function=None, period=None):
        """
        Simulates the closed loop transfer function based on the parameters of the object.

        Args
            transfer_function (control.TransferFunction): closed loop transfer function
        Returns
            (np.array, np.array): times, predictions
        Raises
            ValueError: if there are no parameters defined for the closed loop transfer function
        """
        if transfer_function is None:
            if self.closed_loop_transfer_function is None:
                raise ValueError("No closed loop transfer function defined.")
            transfer_function = self.closed_loop_transfer_function
        if period is not None:
            U = np.sin(2*np.pi*self.times/period)
        else:
            U = np.repeat(1, len(self.times))
        U = U*self.setpoint
        new_times, predictions = control.forced_response(transfer_function, T=self.times, U=U)
        return new_times, predictions
    
    def plot(self, period=None, **kwargs):
        """
        Plots the step response if values are assigned to the closed loop parameters.

        Args:
            kwargs: arguments for OptionManager
        """
        new_times, predictions = self.simulateTransferFunction(period=period)
        df = pd.DataFrame({"time": new_times, "predictions": predictions})
        df["setpoint"] = self.setpoint
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

    def evaluate(self, is_plot=True, **kwargs):
        """
        Creates an SBMLSystem using the design parameters. Records the builder.

        Args:
            kwargs: plot options
        Returns:
            Timeseries (from simulation)
            AntimonyBuilder (from simulation)
        """
        param_dct = {n: None for n in PARAM_NAMES}
        param_dct.update(self.get())
        k_dct = {k: param_dct[k] for k in PARAM_NAMES}
        try:
            simulated_ts, antimony_builder = self.system.simulateSISOClosedLoop(setpoint=self.setpoint,
                                start_time=self.start_time, end_time=self.end_time, num_point=self.num_point,
                                is_steady_state=self.is_steady_state, inplace=False, **k_dct)
            success = True
        except Exception as exp:
            success = False
            msg = "Could not simulate the closed loop system."
            msgs.warn(msg)
        if success:
            if not "title" in kwargs:
                param_dct = self.get()
                text = ["%s=%f " % (name, param_dct[name]) for name in param_dct.keys()]
                title = " ".join(text)
            self.system.plotSISOClosedLoop(simulated_ts, self.setpoint, markers=["", ""], title=title, **kwargs)
            if is_plot:
                plt.show()
            self.history.add()
        return simulated_ts, antimony_builder

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
        self._dct[COL_SETPOINT] = []
        self._dct[COL_RESIDUAL_MSE] = []

    def add(self):
        if not self.is_history:
            return
        for name in PARAM_NAMES:
            self._dct[name].append(self.designer.__getattribute__(name))
        self._dct[COL_CLOSED_LOOP_SYSTEM].append(self.designer.closed_loop_system)
        self._dct[COL_SETPOINT].append(self.designer.setpoint)
        self._dct[COL_RESIDUAL_MSE].append(self.designer.residual_mse)

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
        designer = SISOClosedLoopDesigner(self.designer.system, self.designer.open_loop_transfer_function,
                                          times=self.designer.times,
                                          setpoint=SETPOINT, is_steady_state=self.designer.is_steady_state,
                                          is_history=self.designer.history.is_history, sign=self.designer.sign)
        for name in PARAM_NAMES:
            designer.__setattr__(name, dct[name])
        designer.closed_loop_system = dct[COL_CLOSED_LOOP_SYSTEM]
        designer.residual_mse = dct[COL_RESIDUAL_MSE]
        designer.history.add()
        return designer