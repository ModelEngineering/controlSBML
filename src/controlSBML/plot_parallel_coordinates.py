"""Does parallel coordinate plots for a dataframe."""

from colour import Color # type: ignore
from matplotlib import ticker # type: ignore
import numpy as np
import matplotlib.pyplot as plt # type: ignore
import pandas as pd # type: ignore


def plotParallelCoordinates(df, value_column=None, num_category=10, columns=None,
                            figsize=(15,5), round_digit=None, 
                            title=None, is_plot=True):
    """
    Does a parallel coordinate plots for with lines colored based on categorizing
    the value_column. One category is for nans.

    Args:
        df (pd.DataFrame): contains parameters to plot and the value_column
        value_column (object, optional): Column of dataframe that has the values to
            categorize. Defaults to first column.
        columns (str, optional): Columns in order
        num_category (int, optional): Number of categories to split the value_column,
            including the nan category
        title (str, optional): Title of the plot
        figsize (tuple, optional): Size of the plot
    """
    def fixupDF(df):
        """Cleans_ up the dataframe"""
        ser = df[value_column]
        for column in ["index", value_column, "Unnamed: 0"]:
            if column in df.columns:
                del df[column]
        df[C_DUMMY] = ser
        return df
    #
    C_DUMMY = "dummy"
    df = df.copy()
    if columns is None:
        columns = list(df.columns)
    columns = list(columns)
    if not value_column in columns:
        columns.append(value_column)
    df = df[columns]
    df[value_column] = df[value_column]
    scaled_df = df.copy()
    scaled_df = scaled_df.reset_index(drop=False)
    # Construct the categories for each row
    if value_column is None:
        value_column = columns[0]
    classify_vals = scaled_df[value_column]
    nan_idx = classify_vals.isnull()
    classify_arr = np.array(classify_vals)
    classify_arr[nan_idx] = np.mean(classify_arr[~nan_idx])
    value_column_str = str(value_column)

    # Create the categories and labels
    category_intervals = pd.cut(classify_arr, num_category)
    categories = pd.cut(classify_arr, num_category-1,
                        labels=range(num_category-1)).tolist()
    categories = np.array(categories)
    category_labels = list(range(num_category))

    # Add category for nans
    for idx, category in enumerate(categories):
        label = f"{category_intervals[idx].left} to {category_intervals[idx].right}"
        category_labels[category] = label
    category_labels[-1] = f"nan"
    categories[nan_idx] = num_category - 1

    # Create the colors for the categories
    red = Color("red")
    colors = list(red.range_to(Color("blue"),num_category-1))
    colors.reverse()
    colors.append("grey")

    # Final fixups
    scaled_df = fixupDF(scaled_df)
    df = fixupDF(df)
    columns = list(scaled_df.columns)
    x_vals = list(range(len(columns)))

    # Create  sublots along x axis
    _, axes = plt.subplots(1, len(columns), sharey=False, figsize=figsize)

    # Get min, max and range for each column
    # Normalize the data for each column
    min_max_range = {}
    for column in columns:
        min_max_range[column] = [df[column].min(), df[column].max(),
                                 np.ptp(df[column])]
        denom = np.ptp(df[column])
        if denom > 0:
            scaled_df[column] = np.true_divide(df[column] - df[column].min(), denom)
        else:
            scaled_df[column] = 0

    # Plot each row
    for i, ax in enumerate(axes):
        for idx in scaled_df.index:
            category = categories[idx]
            ax.plot(x_vals, scaled_df.loc[idx, columns], c=str(colors[category]))
            ax.set_xlim([x_vals[i], x_vals[i] + 1])
    
    # Set the tick positions and labels on y axis for each plot
    # Tick positions based on normalised data
    # Tick labels are based on original data
    def set_ticks_for_axis(dim, ax, ticks):
        min_val, max_val, val_range = min_max_range[columns[dim]]
        step = val_range / float(ticks-1)
        if round_digit is None:
            tick_labels = [min_val + step * i for i in range(ticks)]
        else:
            tick_labels = [np.round(min_val + step * i, round_digit) for i in range(ticks)]
        norm_min = scaled_df[columns[dim]].min()
        norm_range = np.ptp(scaled_df[columns[dim]])
        norm_step = norm_range / float(ticks-1)
        ticks = [round(norm_min + norm_step * i, 2) for i in range(ticks)]
        ax.yaxis.set_ticks(ticks)
        ax.set_yticklabels(tick_labels)

    for dim, ax in enumerate(axes):
        ax.xaxis.set_major_locator(ticker.FixedLocator([dim]))
        set_ticks_for_axis(dim, ax, ticks=6)
        ax.set_xticklabels([columns[dim]])
    
    # Move the final axis' ticks to the right-hand side
    axes[-1].remove()   # Remove the dummy
    axes = axes[:-1]
    x_vals = x_vals[:-1]
    columns = columns[:-1]
    ax = plt.twinx(axes[-2])
    dim = len(axes) - 1
    ax.xaxis.set_major_locator(ticker.FixedLocator([x_vals[-2], x_vals[-1]]))
    set_ticks_for_axis(dim, ax, ticks=6)
    ax.set_xticklabels([columns[-2], columns[-1]])

    # Remove space between subplots
    plt.subplots_adjust(wspace=0)
    axes[-1].remove()

    # Add legend to plot
    plt.legend(
        [plt.Line2D((0,1),(0,0), color=str(colors[v])) for v in range(num_category)],
             category_labels,
             bbox_to_anchor=(1.2, 1), loc=2, borderaxespad=0.)

    # Add title
    if title is None:
        title = f"Values by {value_column_str} category"
    else:
        title = f"{title}"
    plt.title(title)
    if is_plot:
        plt.show()