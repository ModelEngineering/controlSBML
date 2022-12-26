# Model urls
from controlSBML._version import __version__
from controlSBML.control_sbml import ControlSBML
from controlSBML.constants import LegendSpec
from controlSBML.nonlinear_io_system import NonlinearIOSystem
from controlSBML.util import plotOneTS, plotManyTS, ppMat, mat2DF, plotMat,  \
      makeSimulationTimes, timeresponse2Timeseries
from controlSBML.sequential_model import SequentialModel
from controlSBML.simulate_system import simulateSystem, makeStateVector
from controlSBML.timeseries import Timeseries, TimeseriesSer
from controlSBML.iosystem_factory import IOSystemFactory
from controlSBML.siso_closed_loop_system import SISOClosedLoopSystem
mat2TS = Timeseries.mat2TS
BIOMODELS_DCT = {
      823: "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000823.2?filename=Varusai2018.xml",
      206: "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000206.2?filename=BIOMD0000000206_url.xml",
}

################FUNCTIONS#######################
def IOSystemFactory_CALLBACK_REPORT(call_name, time, x_vec, u_vec, dct, outputs):
    """
    Callback function used by IOSystemFactory that prints the inputs.

    Parameters
    ----------
    name: name of the component
    time: float
    x_vec: list-float
    u_vec: list-float
    dct; dict
        additional constant parameters
    outputs: list-float
    is_output: bool

    Returns
    -------
    str
    """
    line = call_name 
    line += "  %2.6f" % time
    line += "  %s" % str(x_vec)
    line += "  %s" % str(u_vec)
    line += "  %s" % str(outputs)
    line += "  %s" % str(dct)
    return line
