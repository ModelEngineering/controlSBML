"""Constructs bode plots for MIMO models."""

import control
import matplotlib.pyplot as plt
import numpy as np


def plotBode(sys, **kwargs):

freqs = np.array([np.pi*n/20 for n in range(1, 41)])
#freqs = np.array(range(1,100))
all_mag, all_phase, angle = sys.freqresp(freqs)
names = list(df.columns)
figure, axes = plt.subplots(len(df), len(df), figsize=(10, 10))
for idx1 in range(len(df)):
    for idx2 in range(len(df)):
        ax = axes[idx1, idx2]
        title = "%s->%s" % (names[idx1], names[idx2])
        arr = np.array([10e-5 if np.isclose(v, 0) else v for v in all_mag[idx1, idx2]])
        db = 20*np.log10(arr)
        db = np.array([v if v > -10 else -80 for v in db])
        y_limit = max(max(db), min(-db))
        log_angle = np.log(np.array(angle))
        ax.plot(log_angle, db)
        ax.plot([log_angle[0], log_angle[-1]], [0, 0], linestyle="dashed", color="grey")
        ax.set_ylim([-y_limit, y_limit])
        ax.set_ylim([-20, 20])
        ax.set_title(title)
        if idx1 < len(df) - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("log rads")
        if idx2 == 0:
            ax.set_ylabel("db")
        
