Modeling Closed Loop
====================

.. highlight:: python
   :linenothreshold: 5

Here, we illustrate how to use ``controlSBML`` to model a closed loop system
for the SBML model.
The Antimony model is specified below.
Note that the use of a dollar sign (``$``) indicates
a fixed species whose concentration does not change.

.. code-block:: python

    FIXED_SPECIES_MODEL = """
        species E1;
        J1: $S1 -> S2; S1*E1
        J2: S2 -> ; k*S2
        $S1 = 10; S2 = 0; E1 = 1.0;
        k = 1
    """

If we simulate the above system, ``S2 = 10`` at steady state.
In this example, we construct a closed loop system
with the control objective that
``S2 = 5`` at steady state.
We will achieve this objective by properly
manipulating ``E1``.

Creating the NonlinearIOSystem object for the SBML model.

Setting the control objective.

Constructing the NonlinearIOSystem for the controller.

Specifying the closed loop system. Note that input names and output names
are available for interaconnection.

Evaluating the closed loop system.
