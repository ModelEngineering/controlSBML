# controlSBML
``controlSBML`` provides: (1) a bridge between SBML models and the python ``control`` package and (2) capabilities for the analysis and design of control systems for SBML models.

## Installation
``pip install controlSBML``

To find the current version:
```
import controlSBML as ctl
ctl.__version__
```

### Installing slycot (a ``control`` package dependency)
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
``ctlsb = ControlSBML(model_reference)`` where
``model_reference``can be any of the following:
* path to a local file with the extension ``.ant`` or ``.xml``
* a python string for an SBML or Antimony model representation
* a URL to an XML file

## Key properties
* ``roadrunner`` is the roadrunner object for the model
* ``jacobian`` is the full Jacobian of the model at the current simulation time. Changes to ``roadrunner`` can change this property.
* ``reduced_jacobian`` is the reduced Jacobian of the model at the current simulation time. It is reduced in the
sense that dependent species (those calculated from conservation laws) are removed as are all species
that are either independent of the others or do not affect the others.
* ``antimony`` is the antimony representation of the model

## Key methods
* ``getSpeciesNames`` returns a list of species names for either an unreduced or reduced model
* ``setTime(new_time)`` resets ``roadrunner`` and runs a simulation from time zero to the time specified.
* ``makeStateSpace()`` creates a ``control`` state space object using the Jacobian of the model at the current simulation time.
* ``plotTrueModel`` plots a roadrunner simulation.
* ``plotLinearApproximation`` plots the linear approximation provided by an ${\bf A}$ matrix. The default ${\bf A}$ matrix
is the Jacobian of the current ``RoadRunner`` instance. You can either specify a different ${\bf A}$ or set the ``RoadRunner`` instance to a different time.

## Example

[This notebook](https://github.com/ModelEngineering/controlSBML/blob/main/notebooks/Using-Control-SBML.ipynb) provides a simple example of using ``controlSBML``.

## Developer Notes
1. The package works with and without ``slycot``. So two virtual environments are needed for testing: ``ctl`` includes ``slycot`` and ``ctl_tst`` does not. Note
that continuous integration is done only *without* ``slycot``.

## Version History
* 0.1.4,
  * More options for plotting and simulations
* 0.1.3, 2/13/2022. First release.
