"""Factory for making IOSystem objects"""

"""
Creates IOsystems to be used in constructing closed loop systems.
Inputs have the form in, in1, in2, ...
Outpus have the form out, out1, out2, ...

IOSystems created:
    makeStep: outputs a constant value
    makeSinusoid: creates an IOSystem that outputs a sinusoid
    makeAribtrarySignal: outputs a signal defined by a function
    makeFilter: creates an exponential filter
    makeStateFilter: creates an exponential filter for state variables
    makeAdder: creates an IOSystem that outputs the sum of the inputs
    makePassthur: creates an IOSystem that outputs the input
    makePIDController: creates a PID controller
    makeFullStateController: creates a PID controller

Factory can create a logger for the IOSystem using the is_log parameter. Logs are
in the factor property loggers.
"""

from controlSBML import constants as cn
from controlSBML import logger as lg
from controlSBML import util

import control
import numpy as np

IN = "in"
IN1 = "in1"
IN2 = "in2"
OUT = "out"
REF = "ref"
STATE = "state"


class IOSystemFactory(object):

    def __init__(self, name="factory", dt=1/cn.POINTS_PER_TIME, is_log=False):
        """
        Parameters
        ----------
        name: str
            Name of the InterconnectedSystem
        dt: float
            delta between simulation times
        is_log: bool
            Log activity of factory objects
        """
        self.name = name
        self.dt = dt
        self.registered_names = []  # List of names logged
        self.is_log = is_log
        self.loggers = []  # Loggers for created objects

    def _registerLogger(self, logger_name, item_names, logger=None):
        """Creates and records loggers in the factory."""
        is_duplicate = any([logger_name == l.logger_name for l in self.loggers])
        if is_duplicate:
            raise ValueError("Duplicate logger name: %s" % logger_name)
        if logger is None:
            logger = lg.Logger(logger_name, item_names)
        self.loggers.append(logger)
        return logger

    @staticmethod
    def _array2scalar(vec):
        """
        Converts an array with a single element to a scalar.

        Parameters
        ----------
        vec: float/array-like

        Returns
        -------
        float
        """
        try:
            val = vec[0]
            if len(vec) > 1:
                raise RuntimeError("Array has more than 1 element.")
        except:
            val = vec
        return val

    def initializeLoggers(self):
        """
        Initialize loggers.
        """
        for logger in self.loggers:
            logger.initialize()
    
    def add(self, logger, time, items):
        if self.is_log:
            logger.add(time, items)

    def report(self):
        """
        Creates of report from the logs of objects produced by this factory.

        Returns
        -------
        pd.DataFrame or None if not logging
        """
        if not self.is_log:
            return None
        if len(self.loggers) == 0:
            return None
        if len(self.loggers) == 1:
            return self.loggers[0].report()
        # Handle >1 elements
        logger = self.loggers[0]
        others = self.loggers[1:]
        merged_logger = logger.merge(self.name, others)
        return merged_logger.report()

    def makeNonlinearIOSystem(self, name, ctlsb):
        """
        Creates a "system" object.

        Parameters
        ----------

        Returns
        -------
        """
        system = ctlsb.makeNonlinearIOSystem(name)
        self._registerLogger(name, None, logger=system.updfcn_logger)
        return system

    def makePIDController(self, name, kp=2, ki=0, kd=0):
        """
        Creates a PID controller.
        NonlinearIOSystem with input IN and output OUT.
        States are:
            culmulative_err: cumulative control error
            last_err: last control error

        Parameters
        ----------
        name: str
            Name of the system
        kp: float
           proportional control constant
        ki: float
           integral control constant
        kd: float
           differential control constant.

        Returns
        -------
        control.NonlinearIOSystem
        """
        #
        logger_name = "%s_outfcn" % name
        logger = self._registerLogger(logger_name, [IN, "last_err",
              "acc_err", OUT])
        def updfcn(_, x_vec, u_vec, __):
            # u: float (error signal)
            # x_vec[0] - last_u
            # x_vec[1] - accumulated_u
            u_val = self._array2scalar(u_vec)
            daccumulated_u = u_val/self.dt
            dlast_u_val = (u_val - x_vec[0])/self.dt
            return dlast_u_val, daccumulated_u
        #
        def outfcn(time, x_vec, u_vec, __):
            # u: float (error signal)
            u_val = self._array2scalar(u_vec)
            last_u_val = x_vec[0]
            accumulated_u_val = x_vec[1]
            control_out =  kp*u_val  \
                          + ki*accumulated_u_val \
                          + kd*(u_val - last_u_val)/self.dt
            self.add(logger, time, [u_val, x_vec[0], x_vec[1], control_out])
            return control_out
        #
        return control.NonlinearIOSystem(
            updfcn, outfcn, inputs=[IN], outputs=[OUT],
            states=["last_err", "accumulated_err"],
            name=name)

    def makeFullStateController(self, controller_name, ss_sys, state_names, poles=-2, dcgain=1.0):
        """
        Creates a full state feedback controller for an SBML model
        where the system is linearized at the specified time.

        Parameters
        ----------
        controller_name: str
        ss_sys: State space linear system that approximates
        state_names: list-str
        dcgain: float
            factor for adjusting the reference input to get a closed loop
            transfer function of 1
        poles: list-float/float
            Desired poles; single float is the dominant pole
        time: float
            Time where system is lienarized

        Returns
        -------
        control.NonlinearIOSystem
           inputs:
               <state_variable> (excluding the input)
               ref: reference input
           outputs:
               out
        """
        num_state = len(state_names)
        controller_input_names = list(state_names)
        controller_input_names.insert(0, REF)  # first input
        is_distinct_poles = True
        if "len" in dir(poles):
            if len(set(poles)) < len(poles):
                is_distinct_poles = False
            poles = np.array(poles)
        else:
            # Ensure that poles are distinct
            poles = np.array([poles - 0.1*n for n in range(num_state)])
        if not is_distinct_poles:
            raise ValueError("Poles must be distinct. Not: %s" % str(poles))
        # Calculate the gain matrix
        kp_vec = control.place(ss_sys.A, ss_sys.B, poles)
        #
        def outfcn(_, __, u_vec, ___):
            # u_vec: list-float - reference, state variables
            ref = u_vec[0]/dcgain
            #arr = np.array(u_vec[1:])
            arr = u_vec[1:]
            output = ref - kp_vec.dot(arr)
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, inputs=controller_input_names, outputs=['out'],
            name=controller_name)
    
    def makeArbitrarySignal(self, name, signal_function=lambda t: 1, start_time=0, end_time=None):
        """
        Construct a NonlinearIOSystem that outputs the sum of
        The system has output OUT.

        Parameters
        ----------
        name: str
        signal_function: args: time; returns: float
        start_time: float (offset at which the signal begins)
        end_time: float (offset at which the signal ends)

        Returns
        -------
        NonlinearIOSystem
        """
        logger = self._registerLogger(name, [OUT])
        def outfcn(time, _, __, ___):
            """
            Creates a sine wave at the frequency specified.

            Parameters
            ----------
            time: float
            """
            output = 0.0
            if (time >= start_time) and (end_time is None or time <= end_time):
                output = signal_function(time)
            self.add(logger, time, [output])
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[OUT], inputs=[],
            name=name)
    
    def makeStep(self, name, step_size=1, start_time=0, end_time=None):
        """
        Outputs a constant value.

        Parameters
        ----------
        name: str
        constant: float
        start_time: float (offset at which the signal begins)
        end_time: float (offset at which the signal ends)

        Returns
        -------
        NonlinearIOSystem
        """
        return self.makeArbitrarySignal(name, lambda t: step_size, start_time=start_time, end_time=end_time)

    def makeSinusoid(self, name, amplitude=1, frequency=1, phase=0, dc_offset=0, start_time=0, end_time=None):
        """
        Construct a NonlinearIOSystem that outputs a sinusoid.
        The system has output OUT.

        Parameters
        ----------
        name: str
        amplitude: float
        frequency: float
        phase: float (phase in radians)
        dc_offset: float
        start_time: float (offset at which the signal begins)
        end_time: float (offset at which the signal ends)

        Returns
        -------
        NonlinearIOSystem
        """
        return self.makeArbitrarySignal(name, lambda t: amplitude*np.sin(frequency*t + phase) + dc_offset,
              start_time=start_time, end_time=end_time)

    def makeAdder(self, name, num_input=2, input_names=None):
        """
        Inputs two or more elements. Outputs their sum. Name is "sum".
        The inputs are IN1, IN2, ... or the names in input_names. The presence of
        a "-" in front of a name indicates that the input is subtracted.
        The output is OUT.

        Parameters
        ----------
        name: str
        num_input: int
        num_names: list-str

        Returns
        -------
        NonlinearIOSystem
        """
        if input_names is None:
            input_names = ["%s%d" % (IN, n) for n in range(1, num_input+1)]
        else:
            num_input = len(input_names)
        mult_arr = np.array([1 if not n.startswith("-") else -1 for n in input_names])
        input_names = [n.replace("-", "") for n in input_names]
        item_names = list(input_names)
        item_names.append(OUT)
        logger = self._registerLogger(name, item_names)
        def outfcn(time, _, u_vec, __):
            """
            Creates a sine wave at the frequency specified.

            Parameters
            ----------
            u: float, float
            """
            output = np.sum(u_vec*mult_arr)
            item_values = list(u_vec)
            item_values.append(output)
            self.add(logger, time, item_values)
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, inputs=input_names, outputs=[OUT], name=name)
   
    def makeFilter(self, name, constant=1):
        """
        Construct a NonlinearIOSystem for x' = -a*x + a*u, where u is the
        filter input and -a is the pole of the filer.
        The system has input 'in' and output 'out'.

        Parameters
        ----------
        filter_name: str
        constant: float (pole of the filter)

        Returns
        -------
        NonlinearIOSystem
        """
        logger = self._registerLogger(name, [IN, STATE, OUT])
        def updfcn(_, x_vec, u_vec, __):
            """
            Returns the derivative of the state.

            Parameters
            ----------
            time: float
            x_vec: float
            u_vec: float
            """
            x_val = self._array2scalar(x_vec)
            u_val = self._array2scalar(u_vec)
            dx_val = -constant*x_val + constant*u_val
            return dx_val
        #
        def outfcn(time, x_vec, u_vec, _):
            u_val = self._array2scalar(u_vec)
            x_val = self._array2scalar(x_vec)
            if np.isclose(constant, 0):
                output = u_val
            else:
                output = x_val
            self.add(logger, time, [u_val, x_val, output])
            return output
        #
        return control.NonlinearIOSystem(
            updfcn, outfcn, outputs=[OUT], inputs=[IN], states=[STATE],
            name=name)

    def makeStateFilter(self, filter_name, kf, input_system_name,
          output_system_name, state_names):
        """
        Construct a NonlinearIOSystem for a collection of exponential
        filters for named states. Provides the connections between the
        filter and its input and output. The constant is the exponent
        of an exponential. Note that if kf == 0, this is a passthru.

        Parameters
        ----------
        filter_name: str
        kf: float/list-float
            list must have same length as state_names
            if kf == 0, then this is a pass thru
        input_system_name: str
        output_system_name: str
        state_names: list-str

        Returns
        -------
        NonlinearIOSystem
            <filter_name>.<state_name>_in
            <filter_name>.<state_name>_out
        list-tuple-str (connections)
            input_system and output_system must have
              inputs with the same names as the states
        """
        def makePortName(state_name, sfx):
            return "%s_%s" % (state_name, sfx)
        def makeSignalName(state_name, sfx):
            return "%s.%s_%s" % (filter_name, state_name, sfx)
        def makeSystemSignalName(system_name, state_name):
            return "%s.%s" % (system_name, state_name)
        # Initializtions
        if util.isNumber(kf):
            kfs = np.repeat(kf, len(state_names))
        else:
            kfs = np.array(kf)
        input_dct = {}  # Name of the input to the filter
        output_dct = {}
        internal_dct = {}
        for state_name in state_names:
            input_dct[state_name] = makePortName(state_name, IN)
            output_dct[state_name] = makePortName(state_name, OUT)
            internal_dct[state_name] = makePortName(state_name, STATE)
        # Logger registration
        logging_names = list(input_dct.values())
        logging_names.extend(list(output_dct.values()))
        logging_names.extend(list(internal_dct.values()))
        logger = self._registerLogger(filter_name, logging_names)
        def updfcn(_, x_vec, u_vec, __):
            """
            Returns the derivative of the state.

            Parameters
            ----------
            time: float
            x_vec: float
                internal signal for all states
            u_vec: float
                input for all states

            Returns
            -------
            list-float
                change in state
            """
            dx_vals = [c*x + u for c, x, u in zip(kfs, x_vec, u_vec)]
            return dx_vals
        #
        def outfcn(time, x_vec, u_vec, _):
            """
            Parameters
            ----------
            time: float
            x_vec: float
                internal signal for all states
            u_vec: float
                input for all states

            Returns
            -------
            list-float
                outputs
            """
            outputs = [u if np.isclose(k, 0) else -k*x
                  for k, x, u in zip(kfs, x_vec, u_vec)]
            log_vals = list(u_vec)
            log_vals.extend(list(u_vec))
            log_vals.extend(outputs)
            log_vals.extend(list(x_vec))
            self.add(logger, time, log_vals)
            return outputs
        #
        inputs = list(input_dct.values())
        outputs = list(output_dct.values())
        internals = list(internal_dct.values())
        sys = control.NonlinearIOSystem(updfcn, outfcn, outputs=outputs,
            inputs=inputs, states=internals, name=filter_name)
        # Constrution connetions
        connections = []
        for state_name in state_names:
            # input system
            out_signal = makeSystemSignalName(input_system_name, state_name)
            in_signal = makeSignalName(state_name, IN)
            connections.append([in_signal, out_signal])
            # Output system
            in_signal = makeSystemSignalName(output_system_name, state_name)
            out_signal = makeSignalName(state_name, OUT)
            connections.append([in_signal, out_signal])
        #
        return sys, connections

    def makePassthru(self, name):
        """
        Makes a pass through system that outputs its input.

        Parameters
        ----------
        name: str
        input_name: str
        output_name: str

        Returns
        -------
        NonlinearIOSystem
        """
        logger = self._registerLogger(name, [IN, OUT])
        input_name = IN
        output_name = OUT
        def outfcn(time, _, u_vec, __):
            u_val = self._array2scalar(u_vec)
            output = u_val
            self.add(logger, time, [u_val, output])
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[output_name], inputs=[input_name], name=name)

    def makeMultiplier(self, name, factor=1):
        """
        Outputs a multiple of the input signal.

        Parameters
        ----------
        name: str
        factor: float

        Returns
        -------
        NonlinearIOSystem
        """
        logger = self._registerLogger(name, [IN, OUT])
        def outfcn(time, _, u_vec, __):
            u_val = self._array2scalar(u_vec)
            output = factor*u_val
            self.add(logger, time, [u_val, output])
            return output

        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[OUT], inputs=[IN], name=name)
