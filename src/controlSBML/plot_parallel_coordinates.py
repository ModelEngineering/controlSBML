"""Does parallel coordinate plots for a dataframe."""

from colour import Color # type: ignore
from matplotlib import ticker # type: ignore
import numpy as np
import matplotlib.pyplot as plt # type: ignore
import pandas as pd # type: ignore


def plotParallelCoordinates(df, value_column=None, num_category=10,
                            figsize=(15,5), is_plot=True):
    """
    Does a parallel coordinate plots for with lines colored based on categorizing
    the value_column.

    Args:
        df (pd.DataFrame): _description_
        value_column (object, optional): Column of dataframe that has the values to
            categorize. Defaults to first column.
        num_category (int, optional): Number of categories to split the value_column
    """
    columns = list(df.columns)
    df = df.copy()
    df = df.reset_index(drop=False)
    # Construct the categories for each row
    if value_column is None:
        value_column = columns[0]
    classify_vals = df[value_column]
    value_column_str = str(value_column)
    # Create the categories and labels
    category_intervals = pd.cut(classify_vals, num_category)
    categories = pd.cut(classify_vals, num_category, labels=range(num_category))
    category_labels = list(range(num_category))
    for idx, category in enumerate(categories):
        label = f"{category_intervals[idx].left} to {category_intervals[idx].right}"
        category_labels[category] = label
    # Create the colors for the categories
    red = Color("red")
    colors = list(red.range_to(Color("blue"),num_category))
    colors.reverse()
    # Clean up
    for column in ["index"]:
        del df[column]
    columns = list(df.columns)
    x_vals = list(range(len(columns)))
    #
    # Create  sublots along x axis
    #
    fig, axes = plt.subplots(1, len(columns), sharey=False, figsize=figsize)

    # Get min, max and range for each column
    # Normalize the data for each column
    min_max_range = {}
    for column in columns:
        min_max_range[column] = [df[column].min(), df[column].max(), np.ptp(df[column])]
        df[column] = np.true_divide(df[column] - df[column].min(), np.ptp(df[column]))

    # Plot each row
    for i, ax in enumerate(axes):
        for idx in df.index:
            category = categories[idx]
            ax.plot(x_vals, df.loc[idx, columns], c=str(colors[category]))
            if i < len(x_vals)-1:
                ax.set_xlim([x_vals[i], x_vals[i+1]])
    
    # Set the tick positions and labels on y axis for each plot
    # Tick positions based on normalised data
    # Tick labels are based on original data
    def set_ticks_for_axis(dim, ax, ticks):
        min_val, max_val, val_range = min_max_range[columns[dim]]
        step = val_range / float(ticks-1)
        tick_labels = [round(min_val + step * i, 2) for i in range(ticks)]
        norm_min = df[columns[dim]].min()
        norm_range = np.ptp(df[columns[dim]])
        norm_step = norm_range / float(ticks-1)
        ticks = [round(norm_min + norm_step * i, 2) for i in range(ticks)]
        ax.yaxis.set_ticks(ticks)
        ax.set_yticklabels(tick_labels)

    for dim, ax in enumerate(axes):
        ax.xaxis.set_major_locator(ticker.FixedLocator([dim]))
        set_ticks_for_axis(dim, ax, ticks=6)
        ax.set_xticklabels([columns[dim]])
    
    # Move the final axis' ticks to the right-hand side
    ax = plt.twinx(axes[-2])
    dim = len(axes) - 1
    ax.xaxis.set_major_locator(ticker.FixedLocator([x_vals[-2], x_vals[-1]]))
    set_ticks_for_axis(dim, ax, ticks=6)
    ax.set_xticklabels([columns[-2], columns[-1]])

    # Remove space between subplots
    axes[0].remove()
    plt.subplots_adjust(wspace=0)
    axes[-1].remove()

    # Add legend to plot
    plt.legend(
        [plt.Line2D((0,1),(0,0), color=str(colors[v])) for v in range(num_category)],
             category_labels,
             bbox_to_anchor=(1.2, 1), loc=2, borderaxespad=0.)

    plt.title(f"Values by {value_column_str} category")
    if is_plot:
        plt.show()