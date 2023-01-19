import numpy as np
import pandas as pd
import tellurium as te
import control
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")


LINEAR_MDL = """
S1 -> S2; k1*S1
S2 -> S3; k2*S2
S3 -> S4; k3*S3

k1 = 2
k2 = 1.5
k3 = 1
S1 = 10
S2 = 0
S3 = 0
S4 = 0
"""
LINEAR_RR = te.loada(LINEAR_MDL)
LINEAR_DATA = LINEAR_RR.simulate()
LINEAR_RR.plot(LINEAR_DATA)
LINEAR_STATE_NAMES = ["S1", "S2", "S3", "S4"]
LINEAR_PARAM_DCT = {"input_names": LINEAR_STATE_NAMES, "output_names": LINEAR_STATE_NAMES}

# We can get the derivative of a species concentration at the current simulation time
#timepoint = 0.5  # Time at which we want the derivative
#LINEAR_RR.reset()  # Start at time 0
#_ = LINEAR_RR.simulate(0, timepoint)  # Simulate to the timepoint
# Set values of state
LINEAR_RR["S1"] = 10
LINEAR_RR["S2"] = 0
LINEAR_RR["S3"] = 0
LINEAR_RR["S4"] = 0
# Get the derivatives
ds1 = LINEAR_RR["S1'"]
ds2 = LINEAR_RR["S2'"]
ds3 = LINEAR_RR["S3'"]
ds4 = LINEAR_RR["S4'"]

def updfcnSBML(timepoint, x_vec, __, param_dct):
    """
    Calculates the derivative of populations of hare (H) and lynx (L).

    Parameters
    ----------
    x_vec: array-float (S1, S2, S3, S4)
    param_dct: dict
        "roadrunner": ExtendedRoadrunner
        "state_names": list of states

    Returns
    -------
    list-float (same size as state)
    """
    state_names = param_dct["state_names"]
    roadrunner = param_dct["roadrunner"]
    roadrunner.reset()
    if timepoint < 2:
        # Too small of change
        dstates = np.repeat(0.0, len(state_names))
    else:
        # Run a simulation
        _ = roadrunner.simulate(0, timepoint, 2)
        #
        for idx, state_name in enumerate(state_names):
            roadrunner[state_name] = x_vec[idx]
        #
        dstates = []
        for state_name in state_names:
            name = "%s'" % state_name
            dvalue = roadrunner[name]
            dstates.append(dvalue)
    #
    return dstates

# Tests
param_dct = {"roadrunner": LINEAR_RR, "state_names": ["S1", "S2", "S3", "S4"]}
result = updfcnSBML(2, [10, 0, 0, 0], None, param_dct)
assert(len(result) == 4)
print("OK!")

def outfcnSBML(_, x_vec, __, param_dct):
    """
    Calculates the derivative of populations of hare (H) and lynx (L).

    Parameters
    ----------
    x_vec: array-float (S1, S2, S3, S4)
    param_dct: dict
        "state_names" list-str
        "output_names": list-str (subset of states)

    Returns
    -------
    list-float (same size as state)
    return x_vec
    """
    state_names = param_dct["state_names"]
    output_names = param_dct["output_names"]
    #
    output_vals = []
    for idx, state_name in enumerate(state_names):
        if state_name in output_names:
            output_vals.append(x_vec[idx])
    #
    return output_vals

# Tests
param_dct = {"roadrunner": LINEAR_RR, "state_names": ["S1", "S2", "S3", "S4"], "output_names": ["S2", "S3"]}
output_vals = outfcnSBML(None, [1, 1, 1, 1], None, param_dct)
assert(len(output_vals) == 2)
print("OK!")

output_names = ["S1", "S2", "S3", "S4"]
param_dct = {"roadrunner": LINEAR_RR, "state_names": ["S1", "S2", "S3", "S4"], "output_names": output_names}
linear_sbml_sys = control.NonlinearIOSystem(
    updfcnSBML, outfcnSBML, inputs=('u'), outputs=param_dct["output_names"],
    states=param_dct["state_names"], name='linear_sbml_sys')

# Simulate the system
X0 = [10, 0, 0, 0]               # Initial conditions for species concentrations
T = np.linspace(0, 5, 101)   # Simulation 70 years of time

# Simulate the system
t, y = control.input_output_response(linear_sbml_sys, T, 0, X0, params=param_dct)

output_names = param_dct["output_names"]
# Plot the response
plt.figure(1)
for idx in range(len(output_names)):
    plt.plot(t, y[idx])
plt.legend(output_names)
plt.show()
