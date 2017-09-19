# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import numpy as np
from django.shortcuts import render

from django.http import HttpResponse

from pdf_print.Solarclass import Solar


def index(request):
    
    # Initial Positions for Planets
    names = np.array(["Sun", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"])
    ms = np.array([1.00000597682,  0.00095486104043, 0.000285583733151, 0.0000437273164546, 0.0000517759138449, 7.692307692307693e-09])
    qs = np.array([[0,  0,  0],  [-3.5023653,  -3.8169847,  -1.5507963], [9.0755314,  -3.0458353,  -1.6483708], [8.3101420,  -16.2901086,  -7.2521278], [11.4707666,  -25.7294829,  -10.8169456], [-15.5387357,  -25.2225594,  -3.1902382]])
    vs = np.array([[0, 0, 0],  [0.00565429,  -0.00412490,  -0.00190589], [0.00168318,  0.00483525,  0.00192462], [0.00354178,  0.00137102,  0.00055029], [0.00288930,  0.00114527,  0.00039677], [0.00276725,  -0.00170702,  -0.00136504]])

    # big bang
    Solar_instance = Solar()

    # add planets
    Solar_instance.add_planets(names,  ms,  qs,  vs)
    Solar_instance.evolve(20000,  2000)
    Solar_instance.plot()
    return HttpResponse("Simulation completed.")

def pdf_view(request, path):
    with open(path, 'r') as pdf:
        response = HttpResponse(pdf.read(), mimetype='application/pdf')
        response['Content-Disposition'] = 'inline;filename=some_file.pdf'
        return response
    pdf.closed
