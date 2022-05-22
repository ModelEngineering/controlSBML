# controlSBML 
 The [Systems Biology Markup Language (SBML)](https://co.mbine.org/standards/sbml) is a community standard for representing
simulations of biological models, especially chemical reaction networks.
``controlSBML`` is a package that does control analysis and design of SBML models using
the [CalTech ``control`` package](http://python-control.sourceforge.net/manual/).

``controlSBML`` uses the following "systems characterization" of an SBML model:

* state variables are floating species;
* inputs are reaction fluxes;
* outputs are a subset of the state variables.

We note that reaction fluxes cannot be manipulated directly. Thus,implementation requires a mapping from a reaction flux to an *effector*, typically a chemical
species such as an enzyme that is specific to the reaction.
The keyword ``effector_dct`` refers to a python dictionary that does this mapping.
(See the ``makeNonlinearIOSystem`` method of ``ControlSBML``.)

``controlSBML`` interfaes to the ``control`` package by creating two kinds of ``control`` objects for an SBML model:

* ``NonlinearIOSystem``, which wraps the actual simulation of the SBML model; and
* ``StateSpace``, which is a linear MIMO model approximation of the SBML model at a specified time in the simulation.

These objects can be used in combination with other ``control`` objects to analyze system properties,
design controllers, and construct closed loop systems.

``controlSBML`` also provides analyses of SBML models as well as a means to manipulate
elements of the SBML model.

## Installation
``pip install controlSBML``

To find the current version:
```
import controlSBML as ctl
ctl.__version__
```

### Installing ``slycot`` (optional)
To get the full set of capabilities from the ``control`` package,
you need to install ``slycot``, a package that addresses requirements
for more sophisticated manipulation of MIMO models (e.g.,
solving algebraic Ricotti equations).

The best way is to install ``slycot`` is from binaries using ``anaconda``.
However,
at the present time, ``anaconda`` is difficult to set up on
Google Collaboratory.

Otherwise, you need to install from source. Below are the instructions
for Ubuntu.

* install cmake: ``pip install cmake --upgrade``
* install sikit-build: ``pip install scikit-build``
* install fortran: ``sudo apt-get install gfortran``
* provide path to fortran: ``export FC=`which gfortran``
* install BLAS: ``sudo apt-get install libatlas-base-dev``
* clone Sylcot: ``git clone --recurse-submodules https://github.com/python-control/Slycot.git``
* ``cd Slycot``
* ``python setup.py install`` 

## Constructing ControlSBML
```
ctlsb = ControlSBML(model_reference)
```
 where
``model_reference``can be any of the following:
* path to a local file with the extension ``.ant`` or ``.xml``
* a python string for an SBML or Antimony model representation
* a URL to an XML file

## Technical Summary

### ``ControlSBML``
#### Properties
* ``roadrunner`` is the roadrunner object for the SBML model.
* ``jacobian_df`` is the full Jacobian of the model at the current simulation time. Changes to ``roadrunner`` can change this property.
* ``antimony`` is the antimony representation of the model.
* ``A_df`` is the ${\bf A}$ matrix in the linear system represented as a pandas DataFrame.
* ``B_df`` is the ${\bf B}$ matrix in the linear system represented as a pandas DataFrame.
* ``C_df`` is the ${\bf C}$ matrix in the linear system represented as a pandas DataFrame.
* ``state_names`` is a list of floating species that constitute the state vector.
* ``input_names`` is a list of reaction names that specify the input vector.
* ``output_names`` is a list of floating species that constitute the output vector.

#### Methods
* ``setTime(new_time)`` resets ``roadrunner`` and runs a simulation from time zero to the time specified.
* ``makeStateSpace()`` creates a ``control`` state space object using the Jacobian of the model at the current simulation time.
* ``makeNonlinearIOSystem()`` creates a ``NonlinearIOSystem`` object that simulates the SBML model.
* ``plotTrueModel`` plots a roadrunner simulation.
* ``plotLinearApproximation`` plots the linear approximation provided by an ${\bf A}$ matrix. The default ${\bf A}$ matrix
is the Jacobian of the current ``RoadRunner`` instance. You can either specify a different ${\bf A}$ or set the ``RoadRunner`` instance to a different time.
* ``plotBode`` constructs bode plots for the SISO systems formed by all combinations of inputs and outputs.

### Other Functions
* ``makeTS`` creates a time series object from an array. A time series object is a DataFrame whose index in time (in integer milliseconds),
and the columns are variable names.
* ``simulateSystem`` provides a simplified way to simulate a system (which may be an ``interconnect``) by creating the times and initial,
returns a time series DataFrame.
* ``plotOneTS`` plots a single time series DataFrame.
* ``plotManyTS`` plots multiple time series DataFrames structured so that each column is a separate plot that compares the different times series.

## Example

[This notebook](https://github.com/ModelEngineering/controlSBML/blob/main/notebooks/Using-Control-SBML.ipynb) provides a simple example of using ``controlSBML``.

## Developer Notes
1. The package works with and without ``slycot``. So two virtual environments are needed for testing: ``ctl`` includes ``slycot`` and ``ctl_tst`` does not. Note
that continuous integration is done only *without* ``slycot``.

## Version History
* 0.2.16
* 0.2.15 5/21/2022
  * Fix bug in reverting the semantics of control input to be setting the species
    as a static.
* 0.2.14  5/11/2022
  * Provide callback for each manufactured IOsystemFactory
  * Reverted semantics of control input to a NonlinearIOSystem to be
setting the value rather than adding or subtracting a value.
* 0.2.13 5/9/2022
  * SISOClosedLoopSystem provides a way to construct a closed loop system
    for an SBML model. The system has a PID controller and a filter.
  * IOSysemFactory has a log
* 0.2.12 5/3/2022
  * IOSystemFactory creates IOSystem objects for Adder, Multiplier,
    Filter, PIDController, Sinusoid, Constant, Passthru
  * State inputs add to state, not replace the state value.
* 0.2.11 4/25/2022
  * Fix bug in calculating transfer function that incorrectly considered state
* 0.2.9 4/19/2022
  * Fluxes can be outputs
  * Construction of transfer function includes atol option for simplification
* 0.2.8 4/17/2022
  * Added options to plotTrueModel
  * Updated Using ControlSBML with an example of doing feedback 
* 0.2.7 4/11/2022
  * Species can be inputs
  * makeStateSpace, makeTransferFunction have timepoint argument
* 0.2.6 4/10/2020
  * Default for constructor: is_reduced=False
  * makeTransferFunction
* 0.2.5 4/8/2022
  * Improve performance by not recalculating Jacobian
  * Fix bugs related to implementation of is_reduced as applied on NonlinearIOSystem
* 0.2.4 3/31/2024 - Create reduced A matrix
  * mat2Df - fixed bug with printing column names
  * Create reduced A matrix for makeStateSpace so that A is non-singular
    Default output_names is all floating species
* 0.2.3, 3/22/2022 - Bug fix for mat2DF
* 0.2.2, 3/22/2022 - Bug fix
* 0.2.1, 3/22/2022 - Bug fix
* 0.2.0, 3/22/2022 - Version for class
  * ppMat - pretty print a matrix
  * plotMat - display a heatmap for a matrix
  * mat2TS - convert a matrix to a timeseries
  * mat2DF - convert a matrix to a dataframe
* 0.1.6, 3/16/2022
  * ``Using-Control-SBML.ipynb`` has an example of doing feedback control design
with ``controlSBML``.
  * ``control.NonlinearIOSystem`` wraps an SBML model. Can be used
    in the construction of systems using ``control.interconnect`` and in simulations using ``control.input_output_response``. One caveat is that this may work poorly for models implemented as SBML rate rules.
* 0.1.5, 3/5/2022.
  * More options for plotting and simulations
  * plotBode
  * Inputs are identified by reaction Ids
* 0.1.3, 2/13/2022. First release.
