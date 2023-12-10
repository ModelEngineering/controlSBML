import controlSBML.constants as cn
import controlSBML as ctl
import control
from controlSBML.option_management.option_manager import OptionManager
from controlSBML import msgs

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
FIGURE_WRITER_IDX = 0  # Index of figures writter

############### FUNCTIONS ################

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

@Expander(cn.KWARGS, cn.ALL_KWARGS)
def plotOneTS(time_series, colors=None, markers=None, mgr=None, legend=True, **kwargs):
    """
    Plots a dataframe as multiple lines on a single plot. The
    columns are legends. If an OptionManger is provided,
    the caller must invoke doFig.

    Parameters
    ----------
    time_series: TimeSeries
    mgr: OptionsManager
    colors: list-str
    markers: list-str
    legend: bool (whether to show legend)
    #@expand

    Returns
    -------
    PlotResult
    """
    MARKERS = ["o", "s", "v", "^", "x", "+"]
    if colors is None:
        colors = ["green", "blue", "orange", "purple", "brown", "black"]
        [colors.extend(colors) for _ in range(5)]
    if markers is None:
        markers = MARKERS
        [markers.extend(markers) for _ in range(5)]
    elif markers == False:
        markers = np.repeat("", len(time_series.columns))
    if mgr is None:
        mgr = OptionManager(kwargs)
        is_fig = True
    else:
        is_fig = False
    mgr.plot_opts.set(cn.O_XLABEL, default="time")
    ax = mgr.plot_opts.get(cn.O_AX)
    if cn.O_FIGURE in kwargs.keys():
        fig = kwargs[cn.O_FIGURE]
    else:
        fig = None
    if ax is None:
        fig, ax = plt.subplots(1)
        mgr.plot_opts.set(cn.O_AX, ax)
    times = np.array(time_series.index)/cn.MS_IN_SEC
    sel_colors = list(colors)
    sel_markers = list(markers)
    for column in time_series.columns:
        ax.plot(times, time_series[column], color=sel_colors.pop(0), marker=sel_markers.pop(0))
    if legend:
        legend_spec = cn.LegendSpec(time_series.columns, crd=mgr.plot_opts[cn.O_LEGEND_CRD])
        mgr.plot_opts.set(cn.O_LEGEND_SPEC, default=legend_spec)
    mgr.doPlotOpts()
    if is_fig:
        mgr.doFigOpts()
    return PlotResult(ax=ax, ax2=None, fig=fig, time_series=time_series)

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

    Returns
    -------
    PlotResult
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
    return PlotResult(ax=axes, fig=fig)

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
    dt = 1.0/points_per_time
    dt_ms = int(cn.MS_IN_SEC*dt)
    start_ms = int(start_time*cn.MS_IN_SEC)
    end_ms = int(end_time*cn.MS_IN_SEC)
    times = np.arange(start_ms, end_ms + dt_ms, dt_ms)
    times = times/cn.MS_IN_SEC
    return np.array(times)

def makePIDTransferFunction(kp=None, ki=None, kd=None, name="controller", input_name=cn.IN, output_name=cn.OUT):
    """
    Constructs a PID transfer function.

    Args:
        kp: float
        ki: float
        kd: float

    Raises:
        ValueError: must provide at least one of kp, ki, kd

    Returns:
        control.TransferFunction
    """
    is_none = True
    controller_tf = control.tf([0], [1], name=name, inputs=input_name, outputs=output_name)
    if kp is not None:
        is_none = False
        controller_tf += control.tf([kp], [1])
    if ki is not None:
        is_none = False
        controller_tf += control.tf([ki], [1, 0])
    if kd is not None:
        is_none = False
        controller_tf += control.tf([kd, 0], [1])
    if is_none:
        raise ValueError("At least one of kp, ki, kd, kf must be defined")
    controller_tf.name = name
    controller_tf.input_labels = [input_name]
    controller_tf.output_labels = [output_name]
    return controller_tf

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
        try:
            roadrunner[name] = value
        except RuntimeError:
            msgs.warn("Could not set value for %s" % name)

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

def isNumber(item):
    return isinstance(item, float) or isinstance(item, int)

def simplifyTransferFunction(tf, tol=1e-3):
    """
    Simplifies a transfer function by removing terms that are close to 0
    and reducing the order of the numerator and polynomials.

    Parameters
    ----------
    tf: control.TransferFunction

    Returns
    -------
    control.TransferFunction
    """
    def prune(arr):
        arr = arr[0][0]
        arr_max = np.max(np.abs(arr))
        new_arr = np.array([n if np.abs(n) >= tol*arr_max else 0 for n in arr])
        return new_arr
    def isLastZero(arr):
        if np.isclose(arr[-1], 0):
            return True
        return False
    # Eliminate terms close to 0
    den = prune(tf.den)
    num = prune(tf.num)
    # Reduce the order of the transfer function
    den_len = len(den)
    num_len = len(num)
    min_len = min(den_len, num_len)
    for pos in range(min_len):
        if isLastZero(num) and isLastZero(den):
            num = num[:-1]
            den = den[:-1]
        else:
            break
    return control.TransferFunction(num, den)

def latexifyTransferFunction(tf, num_decimal=2):
    """
    Generates latex for the transfer function, rounding to the number of decimal digits.

    Args:
        tf (control.TransferFunction)
        num_decimal (int, optional): number of decimal digits
    """
    tf = simplifyTransferFunction(tf)
    numr = tf.num[0][0]
    denr = tf.den[0][0]
    # Construct the numerator string
    def latexifyPolynomial(polynomial, variable="s"):
        """Creates a latex representation of the polynomial in s
        Args:
            polynomial (list): high order coefficients first
        Returns:
            str
        """
        max_power = len(polynomial) - 1
        result_str = ""
        for idx, num in enumerate(polynomial):
            polynomial[idx] = round(num, num_decimal)
            power = max_power - idx
            if power == 0:
                poly_term = ""
            elif power == 1:
                poly_term = variable
            else:
                poly_term = "%s^{%s}" % (variable, power)
            if np.isclose(polynomial[idx], 1.0) and (len(poly_term) > 0):
                coefficient = ""
            else:
                coefficient = str(polynomial[idx])
            if not np.isclose(polynomial[idx], 0.0):
                if len(result_str) > 0:
                    if not "-" in coefficient:
                        result_str += " + "
                result_str += "%s %s" % (coefficient, poly_term)
        return result_str
    #
    numr_str = latexifyPolynomial(tf.num[0][0])
    denr_str = latexifyPolynomial(tf.den[0][0])
    if len(numr_str) > 0:
        if len(denr_str) > 0:
            latex = r'$\frac{%s}{%s}$' % (numr_str, denr_str)
        else:
            latex = r'$\infty$'
    else:
        latex = r'$0$'
    return latex

def setNoPlot(kwargs, default=False):
    """
    Sets the is_plot option to False.

    Parameters
    ----------
    is_plot: value of is_plot (default if not provided)
    kwargs: dict
    """
    new_kwargs = dict(kwargs)
    if cn.O_IS_PLOT in kwargs.keys():
        is_plot = kwargs[cn.O_IS_PLOT]
    else:
        is_plot = default
    new_kwargs[cn.O_IS_PLOT] = False
    return is_plot, new_kwargs

def cleanTimes(times):
    """
    Makes sure that int times are evenly spaced.

    Args:
        times: int
    Returns:
        np.array-int
    """
    expected_diff = np.mean(np.diff(times))
    expected_diff = int(np.round(expected_diff))
    new_times = np.repeat(None, len(times))
    new_times[0] = times[0]
    for idx in range(1, len(times)):
        new_times[idx] = new_times[idx - 1] + expected_diff
    arr = np.array(new_times)
    return arr.astype(float)

def calculateInitialValue(transfer_function):
    """
    Calculates the initial value of a transfer function.

    Args:
        transfer_function: control.TransferFunction
    Returns:
        float
    """
    numr = transfer_function.num[0][0]
    denr = transfer_function.den[0][0]
    if len(numr) < len(denr):
        return 0
    return numr[0]/denr[0]

def compareSingleArgumentFunctions(func1, func2, arg_min, arg_max, num_point=100, max_rmse=1e-3):
    """
    Compares the results of two functions that take a single argument.

    Args:
        func1: function
        func2: function
        arg_min: float
        arg_max: float
        num_point: int
        max_mse: float (maximum root mean square error if the same)
    """
    values = np.linspace(arg_min, arg_max, num_point)
    squared_diff_arr = np.array([(func1(v) - func2(v))**2 for v in values])
    mse = np.sum(squared_diff_arr)/len(values)
    if mse < max_rmse**2:
        return True
    else:
        return False
    
def allEqual(list1, list2):
    if len(list1) != len(list2):
        return False
    for l1, l2 in zip(list1, list2):
        if l1 != l2:
            return False
    return True
    
def makeRoadrunnerSymbolDct(roadrunner):
    """
    Creates a symbol dictory for the roadrunner model.

    Args:
        roadrunner: ExtendedRoadrunner
    Returns: dict
        key: str (symbol name)
        value: str (symbol type)
    """
    symbol_dct = {}
    for name in roadrunner.getFloatingSpeciesIds():
        symbol_dct[name] = cn.TYPE_FLOATING_SPECIES
    for name in roadrunner.model.getBoundarySpeciesIds():
        symbol_dct[name] = cn.TYPE_BOUNDARY_SPECIES
    try:
        stoichiometry_mat = roadrunner.getFullStoichiometryMatrix()
        for name in stoichiometry_mat.colnames:
            symbol_dct[name] = cn.TYPE_REACTION
    except Exception as exxp:
         msgs.warn("Could not get stoichiometry matrix. No reactions found.: %s" % exxp)
    for name in roadrunner.getGlobalParameterIds():
        symbol_dct[name] = cn.TYPE_PARAMETER
    for name in roadrunner.getAssignmentRuleIds():
        symbol_dct[name] = cn.TYPE_ASSIGNMENT
    #
    return symbol_dct

def isPositiveRealPart(complex_numbers):
    """
    Determines if any number in a vector has a positive real part.

    Args:
        complex_numbers: list-complex
    """
    for number in complex_numbers:
        if np.real(number) > 0:
            return True
    return False

def isStablePoles(transfer_function):
    """
    Checks if the transfer function is stable.

    Args:
        transfer_function: control.TransferFunction
    Returns:
        bool
    """
    return not isPositiveRealPart(transfer_function.poles())

def isStableZeros(transfer_function):
    """
    Checks if the transfer has zeros in the right half plane.

    Args:
        transfer_function: control.TransferFunction
    Returns:
        bool
    """
    return not isPositiveRealPart(transfer_function.zeros())

def roundToDigits(number:float, num_digits:int=2):
    """
    Rounds so that there are the specified number of non-zero digits.

    Args:
        number (float)
        num_digits (int, optional)
    """
    log_number = np.log10(number)
    if log_number > 0:
        return np.round(number, num_digits)
    required_decimal = -int(log_number) + num_digits
    return np.round(number, required_decimal)


############### CLASSES ################

class PlotResult(object):
    """
    Contains data returned from a plot method.

    Properties:
        time_series: Timeseries
            index: times in plot
            columns: variables plotted
        ax: Matplotlib.Axes (axis plotted)
    """

    def __init__(self, time_series=None, ax=None, ax2=None, fig=None):
        """
        Parameters
        ----------
        time_series: Timeseries
        ax: Matplotlib.Axes or array-Matplot.Axes
        ax2: Matplotlib.Axes or array-Matplot.Axes
        fig: Matplotlib.Figure
        """
        self.time_series = time_series
        self.ax = ax
        self.ax2 = ax2
        self.fig = fig

    def __repr__(self):
        """Avoid printer this object."""
        return("")
