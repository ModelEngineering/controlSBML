# Model urls
from controlSBML._version import __version__
from controlSBML.control_sbml import ControlSBML
from controlSBML.constants import LegendSpec
from controlSBML.nonlinear_io_system import NonlinearIOSystem
from controlSBML.util import plotOneTS, plotManyTS
from controlSBML.sequential_model import SequentialModel
from controlSBML.timeseries import Timeseries, TimeseriesSer
BIOMODELS_DCT = {
      823: "https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000823.2?filename=Varusai2018.xml"
}
