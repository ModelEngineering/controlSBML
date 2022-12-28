System Models Using SBML
========================

.. highlight:: python
   :linenothreshold: 5


``controlSBML`` provides ways to construct a system
model from an SBML model that can be used by the ``control`` package.
Two kinds of models can be created.
The first wraps the entire SBML model as a
``control.NonlinearIOSystem`` object.
This allows the ``control`` package to do the same kind of detailed
simulations that are provided by SBML simulators as as ``libroadrunner``
and ``Copasi``.
The second model is a linear approximation to the SBML model
that is encapsulated in a
``control.StateSpace`` object.

In the sequel, we use the following
Antiimony representation of an SBML model.

.. code-block:: python

    SIMPLE_MODEL = """
        species E1;
        J1: S1 -> S2; S1*E1
        S1 = 10; S2 = 0; E1 = 1;
    """

``ControlSBML`` Constructors
############################
We begin with constructing a ``ControlSBML`` object.
This requires having a representation of an SBML model.
This can be an XML file, an XML URL, an antimony file,
or a ``libroadrunner`` object.
We refer to any of these representations of the SBML model as
a *reference*.
The following constructs ``ControlSBML`` object for this model.

.. code-block:: python

    ctlsb = ctl.ControlSBML(SIMPLE_MODEL)

This object has a number of useful properties.

* ``ctlsb.antimony`` is the Antimony representation, even if the reference is an XML file or a ``libroadrunner`` object.
* ``ctlsb.roadrunner`` is the ``libroadrunner`` object.
* ``ctlsb.state_ser`` is a ``pandas`` ``Series`` that provides the values of the concentrations of floating species (the state variables).
* ``ctlsb.setTime(time)`` runs the simulation to a desired time. This updates ``ctlsb.state_ser``.

To create an object that is more useful for control analysis,
you should specify at least one *input* and at least one *output*.
An input is a chemical species that can be manipulated by
increasing or decreasing its value.
An output is a chemical species or reaction flux that can be
read.
An example of this constructor is:

.. code-block:: python

    ctlsb = ctl.ControlSBML(SIMPLE_MODEL,
        input_names=["E1"], output_names=["S2"])

Note that both ``input_names`` and ``output_names``
are lists of strings.
Inputs and outputs must be specified to make use of
subsequent capabilities of ``ControlSBML``.

Creating ``control.NonlinearIOSystem``
######################################

The following code fragment does a simulation using
the ``control`` package and plots
the output.

.. code-block:: python

    # Inputs and outputs
    input_names = ["E1"]
    output_names = ["S1", "S2"]
    # Create the ControlSBML and NonlinearIOSystem objects
    ctlsb = ctl.ControlSBML(SIMPLE_MODEL,
          input_names=input_names, output_names=output_names)
    sys = ctlsb.makeNonlinearIOSystem("simple")
    # Start simulation at time 0
    ctlsb.setTime(0)
    # Do the simulation using the controls package
    times = [0.1*n for n in range(51)]
    output = control.input_output_response(sys,
          times, X0=ctl.makeStateVector(sys))
    # Plot the result
    colors = ["orange", "green"]
    for idx in range(len(output_names)):
        plt.plot(times, output.outputs[idx,:], c=colors[idx])
    _ = plt.legend(output_names)
    plt.xlabel("time")
    plt.ylabel("concentration")

Below, we provide some details about this script.
In line 5,
``ctlsb`` is constructed to have the input ``E1`` and the outputs ``S1`` and ``S2``.
And, in line 7, ``sys`` is a ``control.NonlinearIOSystem`` object
that wraps the SBML model.
Line 9 sets the start time of the simulation to time 0
(which isn't necessary if ``sys`` is referenced only once).
The simulation using the ``control`` package is done
in line 12.
This requires specifying the times at which simulation results are
to be produced.
It also requires specifying the initial state of the ``NonlinearIOSystem``
object.
This state is obtained from the method ``ctl.makeStateVector``.
The output from the simulation is a two dimensional array since
there are two outputs.
``outputs[0, :]`` is ``S1``, and
``outputs[1, :]`` is ``S2``.
Below is the plot constructed by running this script.

.. image:: images/simple_model_plot.png
  :width: 400

Creating ``control.StateSpace``
###############################

We can construct a linear approximation of an SBML model by
selecting an operating point (the values of of state variables,
which are floating species) and estimating state values
using a linear approximation.
The operating point is defined by a simulation time at which
all state values are obtained.
The linear approximation is calculated using the
Jacobian evaluated at the selected time.

Once a ``ControlSBML`` object has been constructed, we use
the method ``makeStateSpace`` to create
a ``control.StateSpace`` object.
This is illustrated below.
