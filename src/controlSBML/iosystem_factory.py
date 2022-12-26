"""Factory for making IOSystem objects"""

"""
Creates IOsystems to be used in constructing closed loop systems.
Inputs have the form in, in1, in2, ...
Outpus have the form out, out1, out2, ...

IOSystems created:
    makeFilter: creates an exponential filter
    makeStateFilter: creates an exponential filter for state variables
    makeConstant: outputs a constant value
    makeAdder: creates an IOSystem that outputs the sum of the inputs
    makePassthur: creates an IOSystem that outputs the input
    makePIDController: creates a PID controller
    makeFullStateController: creates a PID controller
    makeSinusoid: creates an IOSystem that outputs a sinusoid
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

    def __init__(self, name="factory", dt=1/cn.POINTS_PER_TIME):
        """
        Parameters
        ----------
        name: str
            Name of the InterconnectedSystem
        dt: float
            delta between simulation times
        """
        self.name = name
        self.dt = dt
        self.registered_names = []  # List of names logged
        self.loggers = []  # Loggers for created objects

    def _registerLogger(self, logger_name, item_names, logger=None):
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

    def report(self):
        """
        Creates of report from the logs of objects produced by this factory.

        Returns
        -------
        pd.DataFrame
        """
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
        self._registerLogger(name, None, logger=system.logger)
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
        logger = self._registerLogger(name, [IN, "last_err",
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
            logger.add(time, [u_val, x_vec[0], x_vec[1], control_out])
            return control_out
        #
        return control.NonlinearIOSystem(
            updfcn, outfcn, inputs=[IN], outputs=[OUT],
            states=["last_err", "accumulated_err"],
            name=name)

    def makeFullStateController(self, name, ctlsb, poles=-2, time=0,
          dcgain=1.0):
        """
        Creates a full state feedback controller for an SBML model
        where the system is linearized at the specified time.

        Parameters
        ----------
        name: str
            Name of the system
        ctlsb: ControlSBML
            SISO system
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
        # Validity Checks
        if len(ctlsb.input_names) != 1:
            msg = "SBML model must have a single input. Has: %s"  \
                  % str(ctlsb.input_names)
            raise ValueError(msg)
        # Initializations
        state_space = ctlsb.makeStateSpace(time=time)
        controller_input_names = [n for n in ctlsb.state_names
             if not n in ctlsb.input_names]
        controller_input_names.insert(0, REF)  # first input
        num_state_input = len(controller_input_names) - 1
        is_distinct_poles = True
        try:
            _ = len(poles)
            if len(set(poles)) < len(poles):
                is_distinct_poles = False
            poles = np.array(poles)
        except:
            # Insert that poles are distinct
            #poles = np.array([poles + 0.1*n for n in range(num_state_input)])
            poles = np.array([poles - 0.1*n for n in 
                 range(ctlsb.num_state)])
        if not is_distinct_poles:
            raise ValueError("Poles must be distinct. Not: %s" % str(poles))
        # Calculate the gain matrix
        kp_vec = control.place(state_space.A, state_space.B, poles)
        def outfcn(_, __, u_vec, ___):
            # u_vec: list-float - reference, state variables
            ref = u_vec[0]/dcgain
            #arr = np.array(u_vec[1:])
            arr = u_vec
            output = ref - kp_vec.dot(arr)
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, inputs=controller_input_names, outputs=['out'],
            name=name)

    def makeSinusoid(self, name, amp, frequency):
        """
        Construct a NonlinearIOSystem that outputs a sinusoid.
        The system has output OUT.

        Parameters
        ----------
        name: str
        amp: float
        frequency: float

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
            output = amp*np.sin(time*frequency)
            logger.add(time, [output])
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[OUT], inputs=[],
            name=name)

    def makeAdder(self, name, num_input=2, input_names=None):
        """
        Inputs two or more elements. Outputs their sum. Name is "sum".
        The inputs are IN1, IN2, ... or the names in input_names.
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
            output = np.sum(u_vec)
            item_values = list(u_vec)
            item_values.append(output)
            logger.add(time, item_values)
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, inputs=input_names, outputs=[OUT], name=name)

    def makeFilter(self, filter_name, constant):
        """
        Construct a NonlinearIOSystem for x' = a*x + u, where u is the
        filter input.
        The system has input IN and output OUT.
        The output is normalized so that the DC gain of the filter is 1.
        A constant of 0 is a passthru.

        Parameters
        ----------
        filter_name: str
        constant: float
            e**expo_constant*time

        Returns
        -------
        NonlinearIOSystem
        """
        logger = self._registerLogger(filter_name, [IN, STATE, OUT])
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
            dx_val = constant*x_val + u_val
            return dx_val
        #
        def outfcn(time, x_vec, u_vec, _):
            u_val = self._array2scalar(u_vec)
            x_val = self._array2scalar(x_vec)
            if np.isclose(constant, 0):
                output = u_val
            else:
                output = -constant*x_val
            logger.add(time, [u_val, x_val, output])
            return output
        #
        return control.NonlinearIOSystem(
            updfcn, outfcn, outputs=[OUT], inputs=[IN], states=[STATE],
            name=filter_name)

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
            logger.add(time, log_vals)
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
            
            

    def makeConstant(self, name, constant):
        """
        Outputs a constant value.

        Parameters
        ----------
        name: str
        constant: float

        Returns
        -------
        NonlinearIOSystem
        """
        logger = self._registerLogger(name, [OUT])
        def outfcn(time, _,  __, ___):
            output = constant
            logger.add(time, [output])
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[OUT], inputs=[], name=name)

    def makePassthru(self, name, input_name=None, output_name=None):
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
        if input_name is None:
            input_name = IN
        if output_name is None:
            output_name = OUT
        def outfcn(time, _, u_vec, __):
            u_val = self._array2scalar(u_vec)
            output = u_val
            logger.add(time, [u_val, output])
            return output
        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[output_name], inputs=[input_name], name=name)

    def makeMultiplier(self, name, factor):
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
            logger.add(time, [u_val, output])
            return output

        #
        return control.NonlinearIOSystem(
            None, outfcn, outputs=[OUT], inputs=[IN], name=name)
