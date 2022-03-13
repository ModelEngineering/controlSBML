import controlSBML.constants as cn
from controlSBML.option_management.option_manager import OptionManager

from docstring_expander.expander import Expander
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
def plotOneTS(df, **kwargs):
    """
    Plots a dataframe as multiple lines on a single plot. The
    columns are legends.

    Parameters
    ----------
    df: pd.DataFrame
    #@expand
    """
    mgr = OptionManager(kwargs)
    mgr.plot_opts.set(cn.O_XLABEL, default="time")
    plt.plot(df.index, df)
    legend_spec = cn.LegendSpec(df.columns, crd=mgr.plot_opts[cn.O_LEGEND_CRD])
    mgr.plot_opts.set(cn.O_LEGEND_SPEC, default=legend_spec)
    mgr.doPlotOpts()
    mgr.doFigOpts()

@Expander(cn.KWARGS, cn.PLOT_KWARGS)
def plotManyTS(*dfs, ncol=1, names=None, **kwargs):
    """
    Plots multiple dataframes with the same columns so that each column
    is a different plot.

    Parameters
    ----------
    dfs: sequence of pd.DataFrame
    ncol: int (number of columns)
    names: list-str
        Names of the dataframes
    #@expand
    """
    mgr = OptionManager(kwargs)
    mgr.plot_opts.set(cn.O_XLABEL, default="time")
    nrow = int(np.ceil(len(dfs)/ncol))
    fig, axes = plt.subplots(nrow, ncol)
    axes = np.reshape(axes, (nrow, ncol))
    columns = dfs[0].columns
    for idx, col in enumerate(columns):
        irow = int(np.floor(idx/ncol))
        icol = idx - irow*ncol
        ax = axes[irow, icol]
        mgr.plot_opts[cn.O_AX] = ax
        for df_idx, df in enumerate(dfs):
            try:
                values = df[col]
            except KeyError:
                raise ValueError("DataFrame %d lacks column %s"
                      % (df_idx, col))
            ax.plot(df.index, df[col])
        if names is not None:
            legend_spec = cn.LegendSpec(names,
                   crd=mgr.plot_opts[cn.O_LEGEND_CRD])
            mgr.plot_opts.set(cn.O_LEGEND_SPEC, default=legend_spec)
        mgr.doPlotOpts()
    mgr.doFigOpts()
