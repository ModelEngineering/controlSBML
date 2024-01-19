# controlSBML
``controlSBML`` is a python packag that assists with control engineering of biological systems that are characterized by models in
 the [Systems Biology Markup Language (SBML)](https://co.mbine.org/standards/sbml), a community standard for representing
models of biological models, especially chemical reaction networks.
``controlSBML`` provides the following:
  * Read SBML models, viewing candidate inputs, and outputs, and running simulations.
  * System identification, including the creation of (SISO) transfer function objects from the python control systems library
  * Creating Antimony simulations of closed loop systems
  * Design of closed loop systems

Examples of usage are in this [directory](https://github.com/ModelEngineering/controlSBML/tree/main/examples)

## Installation
``pip install controlSBML``

To find the current version:
```
import controlSBML as ctl
ctl.__version__
```

## Version History
* 1.1.02 12/27/2023
    * Fix bug in dependencies. Update header documentation.
* 1.1.01 12/27/2023
    * Complete change in the architecture. Instead of using NonlinearIOSystems in the python control package, controlSBML generates Antimony code to implement staircase functions and closed loops.
    * Creation of a consistent API.
* 1.0.11 3/1/2023
    * Ensure that state variables are never negative.
* 1.0.10
    * Fixed bug with unequally spaced times
    * Fixed bug so can start a simulate at a time other than 0 and the correct initial state is obtained.
* 1.0.9 2/14/2023
    * Use temporary directory for plots created in tests
* 1.0.8 2/14/2023
    * Avoid error if Jacobian cannot be calculated.
    * Better handling of warnings
* 1.0.7 2/11/2023
    * Add missing dependency (lmfit)
* 1.0.6  2/11/2023
    * ControlSBML.makeSISOTransferFunctionBuilder creates a SISOTransferFunctionBuilder object. The plotStaircaseResponse method of SISOTransferFuntionBuilder indicates the controlability of an input for the output. fitTransferFunction fits a transfer function to the outputs produced by a staircase input.
    * plotStaircaseResponse shows effect of a staircase input on outputs
    * Remove cruft from effector_dct
    * Added plot option writefig which takes arguments, True, False, str (path)
* 1.0.5  1/22/2023
    * Fix bugs in NonlinearIOSystem relating to states and calculations in updfcn.
    * Changes to documentation
* 1.0.4
    * Fix bug so that makeStateSpace honors the time argument
    * Updated Sphinx documentation
    * Fix bug with __version__
* 1.0.3
    * Fix bug with file path in _version
* 1.0.1
    * src directory contains all pcakages
* 1.0.0
    * Redefined inputs as species adjustment (positive or negative)
    * ControlSBML.equals has an option to do a "quick check"
    * Deprecated the use of an effector dictionary
* 0.2.23 12/18/2022
    * Updates for using toml files
* 0.2.22
  * IOSystemFactor.makeStateFilter creates a vector of filters between 2 systems
  * SISOClosedLoopSystem.makeFullStateController creates multiple filters if kf != 0
* 0.2.21 5/30/2022
  * NonlinearIOSystem creates a logger.
  * SISOClosedLoopSystem.makeFullStateSystem has option for filters
  * Changed legend of step respoinse plot to use "reference" instead of "step"
* 0.2.20 5/27/2022
  * Fix bug in SISOClosedLoopSystem.evaluateControllability because scipy didn't
    handle nan values.
* 0.2.19 5/27/2022
  * Fix bug in SISOClosedLoopSystem.evaluateControllability because scipy didn't
    handle nan values.
* 0.2.18 5/26/2022
  * Fix small bug
* 0.2.17 5/26/2022
  * Deleted the callback_log implemented in 0.2.14.
* 0.2.16 5/24/2022
  * Refactored SISOClosedLoopSystem
  * Implemented SISOClosedLoopSystem.makeFullStateController
  * Fixed bug with makePIDController where ki, kd are ineffective.
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
