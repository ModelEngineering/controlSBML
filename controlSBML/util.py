import controlSBML.constants as cn
import controlSBML as ctl
from controlSBML.option_management.option_manager import OptionManager
from controlSBML.option_management.options import Options

from docstring_expander.expander import Expander
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import urllib.request


REPO_URL =  "https://github.com/ModelEngineering/controlSBML/raw/main"
MODULE_URLPAT =  REPO_URL + "/controlSBML/%s.py"  # Pattern for module URL
MODEL_URLPAT =  REPO_URL + "/models/%s"
MODEL_823_FILE = "biomodels_823.ant"
LOCAL_FILE = "local.txt"

def calculateMatrixDistance(mat1, mat2):
    """
    Calculates the distance between two matrices with the same shape.

    Parameters
    ----------
    mat1 - np.array
    mat2 - np.array
    
    Returns
    -------
    float
    """
    if np.shape(mat1) != np.shape(mat2):
        raise ValueError("Matrices must be the same shape.")
    return np.sum( (mat1 - mat2)**2)

def getSharedCodes(module_name="util"):
    """
    Obtains common codes from the github repository.

    Parameters
    ----------
    module_name: str
        name of the python module in the src directory
    """
    url = "https://github.com/ModelEngineering/controlSBML/raw/main/controlSBML/%s.py" % module_name
    local_python = "python.py"
    _, _ = urllib.request.urlretrieve(url=url, filename=local_python)
    with open(local_python, "r") as fd:
        code_str = "".join(fd.readlines())
    print(code_str)
    exec(code_str, globals())

def getModel(file_name=MODEL_823_FILE):
    """
    
    Parameters
    ----------
    file_name: str (filename with extension)
    
    Returns
    -------
    str
    """
    url = MODEL_URLPAT % file_name
    _, _ = urllib.request.urlretrieve(url=url, filename=LOCAL_FILE)
    with open(LOCAL_FILE, "r") as fd:
        model_str = "".join(fd.readlines())
    return model_str

@Expander(cn.KWARGS, cn.PLOT_KWARGS)
def plotOneTS(ts, **kwargs):
    """
    Plots a dataframe as multiple lines on a single plot. The
    columns are legends.

    Parameters
    ----------
    ts: TimeSeries
    #@expand
    """
    mgr = OptionManager(kwargs)
    mgr.plot_opts.set(cn.O_XLABEL, default="time")
    ax = mgr.plot_opts.get(cn.O_AX)
    if ax is None:
        _, ax = plt.subplots(1)
    ax.plot(ts.times, ts)
    legend_spec = cn.LegendSpec(ts.columns, crd=mgr.plot_opts[cn.O_LEGEND_CRD])
    mgr.plot_opts.set(cn.O_LEGEND_SPEC, default=legend_spec)
    mgr.doPlotOpts()
    mgr.doFigOpts()

@Expander(cn.KWARGS, cn.PLOT_KWARGS)
def plotManyTS(*tss, ncol=1, names=None, **kwargs):
    """
    Plots multiple Timeseries with the same columns so that each column
     is a different plot.

    Parameters
    ----------
    tss: sequence of Timeseries
    ncol: int (number of columns)
    names: list-str
        Names of the dataframes
    #@expand
    """
    mgr = OptionManager(kwargs)
    mgr.plot_opts.set(cn.O_XLABEL, default="time")
    nrow = int(np.ceil(len(tss[0].columns)/ncol))
    fig, axes = plt.subplots(nrow, ncol)
    axes = np.reshape(axes, (nrow, ncol))
    columns = list(tss[0].columns)
    for idx, col in enumerate(columns):
        new_mgr = mgr.copy()
        irow = int(np.floor(idx/ncol))
        icol = idx - irow*ncol
        ax = axes[irow, icol]
        new_mgr.plot_opts.set(cn.O_AX, default=ax)
        new_mgr.plot_opts.set(cn.O_TITLE, default=col)
        for ts_idx, ts in enumerate(tss):
            try:
                values = ts[col]
            except KeyError:
                raise ValueError("DataFrame %d lacks column %s"
                      % (ts_idx, col))
            ax.plot(ts.times, ts[col])
        if names is not None:
            legend_spec = cn.LegendSpec(names,
                   crd=new_mgr.plot_opts[cn.O_LEGEND_CRD])
            new_mgr.plot_opts.set(cn.O_LEGEND_SPEC, default=legend_spec)
        new_mgr.doPlotOpts()
    mgr.doFigOpts()

def makeSimulationTimes(start_time=cn.START_TIME, end_time=cn.END_TIME,
      points_per_time=cn.POINTS_PER_TIME):
    """
    Constructs the times for a simulation using the simulation options.

    Parameters
    ----------
    start_time: float
    end_time: float
    points_per_time: int
    
    Returns
    -------
    np.ndarray
    """
    num_point = int(points_per_time*(end_time - start_time))
    dt = (end_time - start_time)/num_point
    times = [start_time + dt*n for n in range(num_point)]
    times.append(end_time)  # Include the endpoint
    return np.array(times)

def mat2DF(mat, column_names=None, row_names=None):
    """
    Converts a numpy ndarray or array-like to a DataFrame.

    Parameters
    ----------
    mat: np.Array, NamedArray, DataFrame
    column_names: list-str
    row_names: list-str
    """
    if isinstance(mat, pd.DataFrame):
        df = mat
    else:
        if len(np.shape(mat)) == 1:
            mat = np.reshape(mat, (len(mat), 1))
        if column_names is None:
            if hasattr(mat, "colnames"):
                column_names = mat.colnames
        if column_names is not None:
            if len(column_names) == 0:
                column_names = None
        if row_names is None:
            if hasattr(mat, "rownames"):
                if len(mat.rownames) > 0:
                    row_names = mat.rownames
        if row_names is not None:
            if len(row_names) == 0:
                row_names = None
        df = pd.DataFrame(mat, columns=column_names, index=row_names)
    return df

def ppMat(mat, column_names=None, row_names=None, is_print=True):
    """
    Provides a pretty print for a matrix or DataFrame)

    Parameters
    ----------
    mat: np.Array, NamedArray, DataFrame
    column_names: list-str
    row_names: list-str
    """
    df = mat2DF(mat, column_names=column_names, row_names=row_names)
    if is_print:
        print(df)

def plotMat(mat, column_names=None, row_names=None, is_plot=True, **kwargs):
    """
    Creates a heatmap for the matrix.

    Parameters
    ----------
    mat: np.Array, NamedArray, DataFrame
    column_names: list-str
    row_names: list-str
    """
    df = mat2DF(mat, column_names=column_names, row_names=row_names)
    if is_plot:
        mgr = OptionManager(kwargs)
        ax = mgr.getAx()
        sns.heatmap(df, cmap="seismic", ax=ax)
        mgr.doPlotOpts()
        mgr.doFigOpts()

# TODO: Tests
def makeRoadrunnerSer(roadrunner, names):
    """
    Contructs a Series for the roadrunner names.

    Parameters
    ----------
    roadrunner: ExtendedRoadrunner
    names: list-str

    Returns
    -------
    pd.Series
    """
    values = list(getRoadrunnerValue(roadrunner, names).values())
    return pd.Series(values, index=names)

# TODO: Tests
def getRoadrunnerValue(roadrunner, names):
    """
    Provides the roadrunner values for a name. If no name,
    then all values are given.

    Parameters
    ----------
    roadrunner: ExtendedRoadrunner
    name: str/list-str

    Returns
    -------
    object/dict
    """
    if isinstance(names, str):
        return roadrunner[names]
    if names is None:
        names = roadrunner.keys()
    return {n: roadrunner[n] for n in names}

# TODO: Tests
def setRoadrunnerValue(roadrunner, name_dct):
    """
    Sets the values of names and values.

    Parameters
    ----------
    name_dct: dict
        key: str
        value: value
    """
    for name, value in name_dct.items():
        if isinstance(value, int):
            value = float(value)
        roadrunner[name] = value

def timeresponse2Timeseries(timeresponse, column_names=None):
    """
    Converts a control.Timeresponse object to a time series.

    Parameters
    ----------
    timeresponse: control.Timeresponse
    column_names: list-str
    
    Returns
    -------
    Timeseries
    """
    if len(np.shape(timeresponse.outputs)) == 1:
        df = pd.DataFrame({0: timeresponse.outputs})
    else:
        df = pd.DataFrame(timeresponse.outputs.transpose())
    if column_names is not None:
        df.columns = column_names
    df.index = timeresponse.t
    return ctl.Timeseries(df)
