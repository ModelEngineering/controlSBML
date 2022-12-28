Concepts
========

Control engineering is a discipline that is applied to a wide variety
of fields, including mechanical, electrical, and computer engineering.
In biological systems, control engineering addresses applications such as:

* biostats (or chemostat) that maintain a constant population of microbes (chemical species);

* drug delivery systems that sense metabolite levels (e.g., glucose) and takes action to maintain these levels (e.g., administer insulin);

* prosthesises that are manipulated by body signals.

A first consideration in control systems is to establish
**control objectives**.
Common control objectives are:

* *stability*: does the system to a steady state value after a disturbance;

* *speed of response*: how quickly does the system settle after a disturbance (if they settle);

* *oscillations*: are their oscillations (which may be desired or undesired); and
  
* *bias*: does the system deviate from the desired outputs.

In control engineering, we start with an existing system that does
not achieve control objectives
(e.g., regulation of glucose concentrations).
The existing system is called the **open loop system (OLS)**.
Control engineering is of particular interest if the OLS
cannot easily be re-engineered to meet control objectives.
Control engineering creates a new system, the **closed loop system (CLS)**,
in which one or more new element is added.
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
The CalTech control systems library
(https://sourceforge.net/projects/python-control/) provides an extensive
collection of Python-based tools for control analsysis and design.
What is missing is an ability to combine SBML models with analysis
capabilities such as those in the control systems library.
