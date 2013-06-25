
import matplotlib.pyplot as plt
import numpy as np
import os
import snapshot
import collections

def setAxLinesBW(ax):
    """
    Take each Line2D in the axes, ax, and convert the line style to be
    suitable for black and white viewing.
    """
    MARKERSIZE = 3

    COLORMAP = {
        'b': {'marker': None, 'dash': (None,None)},
        'g': {'marker': None, 'dash': [5,5]},
        'r': {'marker': None, 'dash': [5,3,1,3]},
        'c': {'marker': None, 'dash': [1,3]},
        'm': {'marker': None, 'dash': [5,2,5,2,5,10]},
        'y': {'marker': None, 'dash': [5,3,1,2,1,10]},
        'k': {'marker': 'o', 'dash': (None,None)} #[1,2,1,10]}
        }

    for line in ax.get_lines():
        origColor = line.get_color()
        line.set_color('black')
        line.set_dashes(COLORMAP[origColor]['dash'])
        line.set_marker(COLORMAP[origColor]['marker'])
        line.set_markersize(MARKERSIZE)

def setFigLinesBW(fig):
    """
    Take each axes in the figure, and for each line in the axes, make the
    line viewable in black and white.
    """
    for ax in fig.get_axes():
        setAxLinesBW(ax)

def make_histogram(data, field, scalar=1.0):

    fig = plt.figure()
    ax  = fig.add_subplot(1, 1, 1)

    radii = data.get_bin_statistic("r")
    percentages = data.get_species_statistic("percentage")
    field_data = data.get_species_statistic(field)

    for species in field_data.iterkeys():
        if any(percentages[species] != 0.0):
            ax.plot(radii, field_data[species] * scalar, label=species)

    # update the view limits
    ax.set_xlim(0.0, radii[-1])
    ax.set_xlabel("Radius")
    ax.set_ylabel(field.capitalize())
    ax.legend(bbox_to_anchor=(0.95, 1.10), loc=2, borderaxespad=0.0)

    fig.savefig(field + str(data["snap_shot"]) + ".png", dpi=150)
    plt.close()

# Does the same thing as make_histogram, but takes a list of data sets and
# overlays the results on a single plot.
def make_histograms(datas, field, names=None, scalar=1.0):

    fig = plt.figure()
    ax  = fig.add_subplot(1, 1, 1)

    if names is None:
        names = ["" for i in datas]

    for data, name in zip(datas, names):
        radii = data.get_bin_statistic("r")
        percentages = data.get_species_statistic("percentage")
        field_data = data.get_species_statistic(field)

        for species in field_data.iterkeys():
            if any(percentages[species] > 0.05):
                error = 1.0 - (float(name) / 1000)
                ax.plot(radii, field_data[species] * scalar, label=species + "-" + str(error))

        # update the view limits
        ax.set_xlim(0.0, radii[-1])
        ax.set_xlabel("Radius (m)")
        ax.set_ylabel(field.capitalize() + ' (K)')

    setFigLinesBW(fig)
    ax.legend(bbox_to_anchor=(0.85, 1.10), loc=2, borderaxespad=0.0)

    fig.savefig(field + str(data["snap_shot"]) + ".pdf")
    plt.close()
