# Model urls
from controlSBML._version import __version__
from controlSBML.control_sbml import ControlSBML
from controlSBML.constants import LegendSpec
from controlSBML.nonlinear_io_system import NonlinearIOSystem
from controlSBML.util import plotOneTS, plotManyTS, ppMat, mat2DF, plotMat,  \
      makeSimulationTimes
from controlSBML.sequential_model import SequentialModel
from controlSBML.simulate_system import simulateSystem, makeStateVector
from controlSBML.timeseries import Timeseries, TimeseriesSer
from controlSBML.iosystem_factory import IOSystemFactory
mat2TS = Timeseries.mat2TS
BIOMODELS_DCT = {
      823: "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000823.2?filename=Varusai2018.xml",
      206: "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000206.2?filename=BIOMD0000000206_url.xml",
}
