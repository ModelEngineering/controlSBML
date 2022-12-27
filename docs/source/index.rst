.. controlSBML documentation master file, created by
   sphinx-quickstart on Sun Dec 18 20:06:56 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``controlSBML``: Control engineering with SBML models
======================================================

Control engineering of biological systems
is of great interest to Bioengineers for applications such as:

* a biostat (or chemostat) that maintains a constant population of microbes (chemical species);

* a drug delivery system that senses metabolite levels (e.g., glucose) and takes action to maintain these levels (e.g., administer insulin);

* a prosthesis that is manipulated by body signals.

Very often biological systems are noisy and subject to disturbances that
affect their behavior.
Control engineering adds one or more elements
to regulate the orginal system through feedback:

* *Controllers* provide regulation (e.g., administration of insulin) in response to a biological signal (e.g., blood concentration of glucose).

* *Filters* reduce noise from the biological signal and the effector.

Control analysis and design proceeds by: (1) analyzing the behavior of the biological
system; (2) determining where to place new elements (and the characteristics
of these elements); and (3) evaluating the behavior of the resulting closed loop system.

Accomplishing the foregoing typically involves modeling the biological system
and new elements to be added and then using special purpose tools
for control analysis and design.
Many biological systems are modeled using the Systems Biology Markup Language (SBML)
a community standard for biomedical models that provides a
way to describe and simulate behavior (e.g., how changes in enzyme levels
affect the concentrations of metabolites).
Further,
the CalTech control systems library
(https://sourceforge.net/projects/python-control/) provides an extensive
collection of Python-based tools for control analsysis and design.

``controlSBML`` provides a bridge between the CalTech
control systems library.
At the lowest level, this bridge wraps an SBML model as an object known
to the controls package so that SBML models can be analyzed.
The approach is:

1. Construct a ``ControlSBML`` object for the SBML model.

2. From this object, obtain objects used in the CalTech package:

   a. ``NonlinearIOSystem`` object that wraps the entire SBML model; and
   
   b. ``StateSpace`` object that is a linear approximation to the SBML model

3. Do control analysis and design using the CalTech package.

Beyond this, ``controlSBML`` facilitates the construction of
models of controllers and filters to facilitate the evaluation of candidate
designs.
This is done through the use of *element factories* that create
``NonlinearIOSystem`` objects that can be used in the evaluation of
closed loop designs.

A third capability of ``controlSBML`` is to automate model construction
for an entire closed loop system.
At present, the package only supports models that have a single control input
and a single signal used in feedback (SISO systems).

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   basics
   closed_loop_system



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
