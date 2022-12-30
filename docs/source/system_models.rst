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
        J1: S1 -> S2; S1 + E1
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

``control.NonlinearIOSystem`` objects
#####################################

The method ``ctl.ControlSBML.makeNonlinearIOSystem`` constructs
a ``control.NonlinearIOSystem`` object for an SBML model.
The ``control.NonlinearIOSystem`` object provides a way to use
many of the capabilities of the control system library.
For example,


.. code-block:: python

    ctlsb = ctl.ControlSBML(SIMPLE_MODEL,
            input_names=["E1"], output_names=["S1", "S2"])
    sys = ctlsb.makeNonlinearIOSystem("simple")

``sys`` is a ``control.NonlinearIOSystem`` object,
and the following simulates this system for 5 s.

.. code-block:: python

    times = [0.1*n for n in range(51)]
    time_response = control.input_output_response(sys,
        times, X0=ctl.makeStateVector(sys))

``ctlsb`` is constructed to have the input ``E1`` and the outputs ``S1`` and ``S2``.
``sys`` is a ``control.NonlinearIOSystem`` object
that wraps the SBML model.
``time_response`` is a ``control.TimeResponseData`` object.
The simulation requires knowlege of initial values for all state variables,
which is provided by the method ``ctl.makeStateVector``.

Going a bit further, we introduce a plotting function for
the ``control.TimeResponseData`` object.

.. code-block:: python

    def plotTimeResponse(time_response):
        # Plots the results of running a simulation
        outputs = time_response.outputs
        times = time_response.time
        colors = ["orange", "green"]
        for idx in range(len(output_names)):
            if np.ndim(outputs) > 1:
                plt.plot(times, outputs[idx,:], c=colors[idx])
            else:
                plt.plot(times, outputs, c=colors[idx])
        _ = plt.legend(output_names)
        plt.xlabel("time")
        plt.ylabel("concentration")


We execute the statement below to plot the simulation results.

.. code-block:: python

    plotTimeResponse(time_response)

.. image:: images/simple_model_plot.png
  :width: 400

``control.StateSpace`` objects
##############################

A state space model is a linear system of differential equations.

.. math:: 
    
        \dot{\bf x}  &=  {\bf A} {\bf x} + {\bf B} {\bf u} \\
        {\bf y}      &=  {\bf C} {\bf x}

where:

.. math:: 

    {\bf x} \text{ has dimension }  n \times 1 \\
    {\bf u} \text{ has dimension }  p \times 1 \\
    {\bf y} \text{ has dimension }  q \times 1 \\
    {\bf A} \text{ has dimension }  n \times n \\
    {\bf B} \text{ has dimension }  n \times p \\
    {\bf C} \text{ has dimension }  q \times p \\


:math:`{\bf x}` is the state variable,
:math:`{\bf u}` is the forced input,
and :math:`{\bf y}` is the output.

``controlSBML`` uses the floating species
in the SBML model as state variables.
The :math:`{\bf u}` are the ``input_names``,
and the :math:`{\bf y}` are the ``output_names``.
A
linear approximation for an SBML model is constructed
using the Jacobian of the state variables at a specified operating point.
The operating point is a simulation time at which state variables are assigned their values
to calculate the Jacobian.

Once a ``ControlSBML`` object has been constructed,
the method ``makeStateSpace`` is used to create
a ``control.StateSpace`` object.
This is illustrated below to construct a ``control.StateSpace`` object using
time 0 as the operating point.

.. code-block:: python

    ctlsb = ctl.ControlSBML(SIMPLE_MODEL,
        input_names=["E1"], output_names=["S1", "S2"])
    state_space = ctlsb.makeStateSpace(time=0)

The resulting state space model is represented below:
:math:`{\bf A}` is in the upper left;
:math:`{\bf B}` is in the upper right;
and :math:`{\bf C}` is in the lower left.
Note that :math:`{\bf A}` is A
:math:`3 \times 3` matrix because
``E1``, ``S1``, and ``S2`` are floating species
and hence state variables.
:math:`{\bf B}` is a :math:`3 \times 1` matrix
because there is one forced input, ``E1``.
And, :math:`{\bf C}` is :math:`2 \times 3` because
there are two outputs, ``S1``, ``S2``.


.. image:: images/state_space_matrix.png
  :width: 200