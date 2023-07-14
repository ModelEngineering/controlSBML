Concepts and Terms
==================

.. _concepts:

Control engineering is widely used in science and engineering
to improve performance and cost effectiveness.
Mechanical engineers use control engineering to build resilient
structures from inexpensive compoonents.
Electrical engineers use conrol engineering to create
circuits with lower power demands and larger ratios of signal to noise.
Computer engineers use control engineering to improve the
efficiency and quality of job scheduling in large computing complexes.

In biological systems, control engineering addresses applications such as:

* biostats (or chemostat) that maintain a constant population of microbes (chemical species);

* drug delivery systems that sense metabolite levels (e.g., glucose) and takes action to maintain these levels (e.g., administer insulin);

* prosthesises that are manipulated by body signals.

A first consideration in control systems is to establish
**control objectives**.
Common control objectives are:

* *stability*: Does reach a steady state value after a disturbance?

* *speed of response*: How quickly does the system settle after a disturbance?

* *oscillations*: Are there oscillations (which may be desired or undesired)?
  
* *bias*: Does the system deviate from the desired outputs?

Suppose we have a system that does not achieve one or more
of its control objectives.
One approach is to reengineer the system.
However, sometimes reengineering is either too expensive or it
is simply infeasible (i.e., changing a person' DNA).
Control engineering is largely about the second case.
We use the term **open loop system (OLS)** to refer to
the existing system.
Starting with the OLS,
control engineers add additional elements to create 
a **closed loop system (CLS)** that meets the control objectives.
These elements include:

* *Controllers* provide regulation (e.g., administration of insulin) in response to a biological signal (e.g., blood concentration of glucose).

* *Filters* reduce noise from the biological signal and the effector.

In order to proceed with control engineering,
the OLS must be properly enabled.
First, the OLS must expose **measured outputs**, the quantities
we want to control (e.g.,. glucose level).
Second, the OLS
must have one or more **control inputs** 
(e.g., administration of insulin) that effect changes in
measured outputs to achieve control objectives.

A central part of control engineering is predictive modeling.
A **system model**
predicts how a change in the input to a system affects its outputs.
System models may be numerical simulations or mathematical representations
(e.g., differential equations).
The system model of the OLS predicts
how control inputs affect measured outputs.
Models can also be constructed for controllers and filters.
With these models and knowledge of how elements are interconnected,
we
can construct a system model for the closed loop system
to evaluate its ability to meet control objectives.

Since modeling is time-consuming and knowledge-intenstive, we want
to reuse existing models.
Some biological systems are modeled using the Systems Biology Markup Language (SBML),
a community standard for biomedical models that provides a
way to describe and simulate behavior (e.g., how changes in enzyme levels
affect the concentrations of metabolites).
The
`CalTech control systems library <https://python-control.readthedocs.io/en/latest/intro.html>`_.
provides an extensive
collection of Python-based tools for control analsysis and design.
The ``controlSBML`` package provides ways for control engineers
to use SBML models in the control system library so that
existing biomedical models can used in control analysis and design.
