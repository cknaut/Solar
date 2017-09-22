# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import numpy as np
# import pdb
# from django.http import HttpResponse


def index(request):
    names = np.array(["Sun", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"])
    ms = np.array([1.00000597682,  0.00095486104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449, 7.692307692307693e-09])
    qs = np.array([[0,  0,  0],  [-3.5023653,  -3.8169847,  -1.5507963], [9.0755314,  -3.0458353,  -1.6483708], [8.3101420,  -16.2901086,  -7.2521278], [11.4707666,  -25.7294829,  -10.8169456], [-15.5387357,  -25.2225594,  -3.1902382]])
    vs = np.array([[0, 0, 0],  [0.00565429,  -0.00412490,  -0.00190589], [0.00168318,  0.00483525,  0.00192462], [0.00354178,  0.00137102,  0.00055029], [0.00288930,  0.00114527,  0.00039677], [0.00276725,  -0.00170702,  -0.00136504]])
    Tend = 40000
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
    matplotlib.use('Agg')
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
