# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import numpy as np
from django.shortcuts import render
import pdb
# from django.http import HttpResponse


def prepare_planet_data(qs, precision=4):
    """converts numpy array containing all planet path info to list of plotly-readable strings

    Input:
            qs                  n_planets*1D-np.array
            precision           number of decimals kept during conversion (default=8)

    Output:
           x_list, y_list       lists of numpy_array in string form "[1.987, 2.234]"
    """
    x_list = []
    y_list = []
    N_planets = np.shape(qs)[0]
    for i in range(N_planets):
        q_planet = np.asarray(qs[i])
        x_list.append(numpy_to_plotlystring(q_planet[:, 0], precision=precision))
        y_list.append(numpy_to_plotlystring(q_planet[:, 1], precision=precision))

    return x_list, y_list


def numpy_to_plotlystring(numpy_array, precision=4):
    """converts 1Dnumpy array to string readable py plotly.js

    Input:
            numpy_array         1D-np.array
            precision           number of decimals kept during conversion (default=8)

    Output:
           string_array         data of numpy_array in string form "[1.987, 2.234]"
    """
    np.set_printoptions(threshold=np.nan)   # enables full printing
    n_data = np.shape(numpy_array)[0]
    # prevent array2string to cut string by ofsestting maximal width
    max_line_width=n_data*(8) # Output Format [#.####, ...], thus 7 chars per entry, added 100 for safekeeping
    string_array = np.array2string(numpy_array, separator=',', max_line_width=n_data*(8), precision=precision, formatter={'float_kind':lambda x: "%.4f" % x})
    np.set_printoptions(linewidth=75, threshold=1000)   # restores default
    return string_array

def cover(request):
    return render(request, 'pdf_print/base.html')


def pic_output(request):
    names = np.array(["Sun", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"])
    ms = np.array([1.00000597682,  0.00095486104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449, 7.692307692307693e-09])
    qs = np.array([[0,  0,  0],  [-3.5023653,  -3.8169847,  -1.5507963], [9.0755314,  -3.0458353,  -1.6483708], [8.3101420,  -16.2901086,  -7.2521278], [11.4707666,  -25.7294829,  -10.8169456], [-15.5387357,  -25.2225594,  -3.1902382]])
    vs = np.array([[0, 0, 0],  [0.00565429,  -0.00412490,  -0.00190589], [0.00168318,  0.00483525,  0.00192462], [0.00354178,  0.00137102,  0.00055029], [0.00288930,  0.00114527,  0.00039677], [0.00276725,  -0.00170702,  -0.00136504]])
    Tend = 20000
    nsteps = 2000

    from pdf_print.Solarclass import Solar

    # big bang
    Solar = Solar()

    # add planets
    Solar.add_planets(names,  ms,  qs,  vs)
    Solar.evolve(Tend,  nsteps)

    # retrieve planet paths
    qs = [planet.position for planet in Solar.planets]

    x,y = prepare_planet_data(qs)  # convert data in string format

    # prepare variable names uesd by plotly in HTML
    tracenames = []
    for i in range(Solar.n_planets):
        tracenames.append("trace%d" % i)

    context = {
        'x': x[3],
        'y': y[3],
        'ẗracenames' : tracenames
    }
    return render(request, 'pdf_print/pic.html', context)


def pdf_print(request):
    names = np.array(["Sun", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"])
    ms = np.array([1.00000597682,  0.00095486104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449, 7.692307692307693e-09])
    qs = np.array([[0,  0,  0],  [-3.5023653,  -3.8169847,  -1.5507963], [9.0755314,  -3.0458353,  -1.6483708], [8.3101420,  -16.2901086,  -7.2521278], [11.4707666,  -25.7294829,  -10.8169456], [-15.5387357,  -25.2225594,  -3.1902382]])
    vs = np.array([[0, 0, 0],  [0.00565429,  -0.00412490,  -0.00190589], [0.00168318,  0.00483525,  0.00192462], [0.00354178,  0.00137102,  0.00055029], [0.00288930,  0.00114527,  0.00039677], [0.00276725,  -0.00170702,  -0.00136504]])
    Tend = 20000
    nsteps = 2000
    cmap = "jet_r"
    return print_pdf(request, names, ms, qs, vs, Tend, nsteps, cmap)


def print_pdf(request, names, ms, qs, vs, Tend, nsteps, cmap):
    """performs evolution of gravitational problem with parameters

    Input:
            Planet Parameters:
            names           ... np.array of planet names
            ms              ... np.array of planet masses
            qs              ... np.array of planet initial positions
            vs              ... np.array of planet initial velocities
            tEnd            ... np.array of planet names
            nsteps          ... np.array of planet names

            Plot Parameters:
            cmap            ... matplolib colormap
    Output:
            matplotlib-plot showing evolution prited on browser

    """
    from pdf_print.Solarclass import Solar
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    import matplotlib.pyplot as plt
    import django


    # big bang
    Solar = Solar()

    # add planets
    Solar.add_planets(names,  ms,  qs,  vs)
    Solar.evolve(Tend,  nsteps)

    # retrieve planet paths
    qs = [planet.position for planet in Solar.planets]

    # plotting, inspired by http://scipy-cookbook.readthedocs.io/items/Matplotlib_Django.html
    fig = plt.Figure(figsize=(12, 8))
    ax = fig.add_subplot(111)
    cmap = plt.get_cmap(cmap)  # colormap for planet paths from input

    N = Solar.n_planets
    for i in range(N):
        q_planet = np.asarray(qs[i])
        color = cmap(float(i)/N)  # manually attribute colors to ensure same coloring of paths and endpoint
        ax.plot(q_planet[:, 0], q_planet[:, 1], "-", color=color)
        ax.plot(q_planet[-1, 0], q_planet[-1, 1], "o", color=color, label=Solar.planets[i].name)

    ax.legend()
    canvas = FigureCanvas(fig)
    response = django.http.HttpResponse(content_type='image/png')
    canvas.print_png(response)
    return response
