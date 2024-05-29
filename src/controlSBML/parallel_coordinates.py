"""Does parallel coordinate plots for a dataframe."""

from colour import Color # type: ignore
from matplotlib import ticker # type: ignore
import numpy as np
import matplotlib.pyplot as plt # type: ignore
import pandas as pd # type: ignore
from typing import Tuple, Dict, List

C_DUMMY = "dummy"
C_INDEX = "index"
C_UNNAMED = "Unnamed: 0"
LABEL_NAN = "nan"

class ParallelCoordinates(object):

    def __init__(self, df, value_column=None, num_category=10, columns=None,
                 figsize=(15,5), round_digit=None, title=None, 
                 noise_std=0.01, is_plot=True):
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
            noise_std (float, optional): Standard deviation of noise to add to the data 
        """
        # Save inputs
        self.unscaled_df = df.copy()
        if value_column is None:
            value_column = columns[0]
        self.value_column = value_column
        self.num_category = num_category
        self.noise_std = noise_std
        #
        self.figsize = figsize
        self.round_digit = round_digit
        self.title = title
        self.is_plot = is_plot
        self.values = self.unscaled_df[self.value_column]
        # Handle columns
        if columns is None:
            columns = list(self.unscaled_df.columns)
        self.columns = columns
        # Ensure that value_column is the last column
        if self.value_column in columns:
            self.columns.remove(value_column)
        self.columns.append(value_column)
        self.unscaled_df = self.unscaled_df[self.columns]
        # Manipulate data
        self.scaled_df = self._makeScaledDF()
        self.unscaled_df = self._fixupDF(self.unscaled_df)
        self.scaled_df = self._fixupDF(self.scaled_df)
        self.columns = list(self.scaled_df.columns)

    def _makeScaledDF(self)->pd.DataFrame:
        """
        Creates a dataframe that in which each column is scaled between 0 and 1.

        Returns:
            pd.DataFrame
        """
        scaled_df = self.unscaled_df.copy()
        scaled_df = self._fixupDF(scaled_df)
        # Adjust the value ranges
        for column in self.columns:
            ser = self.unscaled_df[column]
            min_val, range_val = ser.min(), np.ptp(ser)
            if range_val > 0:
                values = (ser - min_val)/range_val
            else:
                values = 0
            scaled_df[column] = values + np.random.normal(0, self.noise_std, len(ser))
        #
        return scaled_df

    def _fixupDF(self, df)->pd.DataFrame:
        """Cleans_ up the dataframe"""
        df = df.reset_index(drop=False)
        for column in [C_INDEX, C_UNNAMED]:
            if column in df.columns:
                del df[column]
        return df
    
    def _setAxisTicks(self, column:str, ax, num_tick:int=6):
        """
        Set the tick positions and labels on y axis for each plot
        Tick positions based on normalised data
        Tick labels are based on original data

        Args:
            column: str
            ax (matplotlib.Axis): _description_
            num_tick (int): Number of ticks
        """
        def min_step(ser):
            min_val = ser.min()
            range_val = np.ptp(ser)
            step = range_val/float(num_tick-1)
            return min_val, step
        #
        scaled_min, scaled_step = min_step(self.scaled_df[column])
        unscaled_min, unscaled_step = min_step(self.unscaled_df[column])
        if self.round_digit is None:
            tick_labels = [unscaled_min + unscaled_step * i for i in range(num_tick)]
            ticks = [scaled_min + scaled_step * i for i in range(num_tick)]
        else:
            tick_labels = [np.round(unscaled_min + unscaled_step * i, self.round_digit) for i in range(num_tick)]
            ticks = [np.round(scaled_min + scaled_step * i, self.round_digit) for i in range(num_tick)]
        is_ticks = [l != tick_labels[0] for l in tick_labels]
        is_ticks[0] = True
        ticks = [t for i, t in enumerate(ticks) if is_ticks[i]]
        tick_labels = [t for i, t in enumerate(tick_labels) if is_ticks[i]]
        ax.yaxis.set_ticks(ticks)
        ax.set_yticklabels(tick_labels)

    def _makeCategoriesAndLabels(self)->Tuple[np.array, np.array, np.array]:
        """
        Create the categories for each row.

        Returns:
            np.array - categories for rows
            np.array - untrimmed labels for each categories
            np.array - trimmed labels for each categories
        """
        NAN_CATEGORY = self.num_category - 1
        # Construct the categories for each row
        classify_vals = self.unscaled_df[self.value_column]
        nan_idx = np.array([np.isnan(v) or np.isinf(v) for v in classify_vals])
        classify_arr = np.array(classify_vals)
        classify_arr[nan_idx] = np.mean(classify_arr[~nan_idx])
        # Make the categories accounting for the nans
        categories = pd.cut(classify_arr, self.num_category-1,
                            labels=range(self.num_category-1)).tolist()
        categories = np.array(categories)
        categories[nan_idx] = NAN_CATEGORY
        # Create the categories and labels
        category_intervals = pd.cut(classify_arr, self.num_category)
        trimmed_category_labels = list(np.repeat("", self.num_category))
        untrimmed_category_labels = list(np.repeat("", self.num_category))
        for idx, category in enumerate(categories):
            left_value = category_intervals[idx].left
            right_value = category_intervals[idx].right
            label = f"{left_value} to {right_value}"
            untrimmed_category_labels[category] = label
            #
            if self.round_digit is not None:
                left_value = np.round(left_value, self.round_digit)
                right_value = np.round(right_value, self.round_digit)
                label = f"{left_value} to {right_value}"
                trimmed_category_labels[category] = label
            else:
                trimmed_category_labels[category] = untrimmed_category_labels[category]
        untrimmed_category_labels[NAN_CATEGORY] = LABEL_NAN
        trimmed_category_labels[NAN_CATEGORY] = LABEL_NAN
        #
        return categories, untrimmed_category_labels, trimmed_category_labels

    def plot(self):
        """
        Plot the parallel coordinates
        """
        categories, untrimmed_labels, trimmed_labels = self._makeCategoriesAndLabels()

        # Create the colors for the categories
        red = Color("red")
        colors = list(red.range_to(Color("blue"), self.num_category-1))
        colors.reverse()
        colors.append("grey")

        # Create  sublots along x axis
        _, axes = plt.subplots(1, len(self.columns), sharey=False, figsize=self.figsize)

        # Plot each row
        x_vals = list(range(len(self.columns)))
        for i, ax in enumerate(axes):
            for idx in self.scaled_df.index:
                category = categories[idx]
                ax.plot(x_vals, self.scaled_df.loc[idx, self.columns], c=str(colors[category]))
                ax.set_xlim([x_vals[i], x_vals[i] + 1])

        # Set the ticks for each axis
        for idx, ax in enumerate(axes):
            column = self.columns[idx]
            ax.xaxis.set_major_locator(ticker.FixedLocator([idx]))
            self._setAxisTicks(column, ax)
            ax.set_xticklabels([column])
        
        # Move the final axis' ticks to the right-hand side
        axes[-1].remove()   # Remove the dummy
        axes = axes[:-1]
        x_vals = x_vals[:-1]
        columns = self.columns[:-1]
        ax = plt.twinx(axes[-2])
        idx = len(axes) - 1
        ax.xaxis.set_major_locator(ticker.FixedLocator([x_vals[-2], x_vals[-1]]))
        self._setAxisTicks(columns[idx], ax)
        ax.set_xticklabels([columns[-2], columns[-1]])

        # Remove space between subplots
        plt.subplots_adjust(wspace=0)
        axes[-1].remove()

        # Add legend to plot
        valid_idxs = np.array([v for v in range(self.num_category) if v in categories])
        trimmed_category_arr = np.array(untrimmed_labels)[valid_idxs]
        plt.legend(
            [plt.Line2D((0,1),(0,0), color=str(colors[v])) for v in range(self.num_category)
                    if v in categories],
                trimmed_category_arr,
                bbox_to_anchor=(1.75, 1), loc=2, borderaxespad=0.)

        # Add title
        if self.title is None:
            value_column_str = str(self.value_column)
            title = f"Values by {value_column_str} category"
        else:
            title = f"{self.title}"
        plt.title(title)
        if self.is_plot:
            plt.show()

    @classmethod
    def plotParallelCoordinates(cls, df, value_column=None, num_category=10, columns=None,
                 figsize=(15,5), round_digit=None, title=None, is_plot=True):
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
        plotter = cls(df, value_column=value_column, num_category=num_category, columns=columns,
                    figsize=figsize, round_digit=round_digit, title=title, is_plot=is_plot)
        plotter.plot()