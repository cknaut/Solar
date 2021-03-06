# Generalized Hamiltonian N-Body Problem approximated by stoermer Verlet
import numpy as np
import numpy.linalg as linalg
import pdb
import sys
from time import time


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

# containg list of Planet objects as well as evolving functions and plotting functions
class Solar():

    def rhs(self, q, m):
        """Right hand side of ODE
        Input:  q   ...     np.array with positions
        Output: dp  ...     time-derivative of the velocities
        """
        # Number of bodies
        N = self.n_planets
        # Empty np.arrays for computed data
        dp = np.zeros_like(q)

        for k in range(N):
            dp_k = np.zeros(self.D)
            for i in range(N):
                if k != i:
                    qi_qk = q[i*self.D:(i+1)*self.D] - q[k*self.D:(k+1)*self.D]
                    qi_qk /= linalg.norm(qi_qk)**3
                    dp_k += m[i]*qi_qk
            dp_k *= self.G
            dp[k*self.D:(k+1)*self.D] = dp_k
        return dp

    def integrate_VV_one_step(self, h, m):
            """Integrate ODE with velocity verlet rule

            Input: y0     ... initial condition
                   h      ... timestep

            Output:
                    y ... solution of y after one propagation
            """

            qs = np.array([planet.position[-1] for planet in self.planets]).flatten()
            vs = np.array([planet.velocity[-1] for planet in self.planets]).flatten()

            # Velocity verlet rule
            q_new = qs + h*vs + 0.5*h**2*self.rhs(qs, m)
            v_new = vs + h / 2.0 * (self.rhs(qs, m) + self.rhs(q_new, m))

            # update position and velocities in planet class
            for i in range(self.n_planets):
                if self.planets[i].fix:
                    self.planets[i].changeposvel(self.planets[i].initpos, self.planets[i].initvel)
                else:
                    self.planets[i].changeposvel(q_new[i*self.D:(i * self.D) + self.D], v_new[i*self.D:(i * self.D) + self.D])

    def evolve(self, tEnd, nsteps):
        """performs evolution of gravitational problem by solving ODE System

        Input: tEnd          ... endtime
               nsteps        ... number of timesteps

        """
        tStart = self.timeline[-1]
        m = [planet.mass for planet in self.planets]
        h = (float(tEnd - tStart))/float(nsteps)

        for i in range(nsteps):
            self.integrate_VV_one_step(h, m)
            self.timeline.append(tStart + (i+1)*h)

    # initialize with empty list of planets
    # animate is used to invoc matplotlib animation
    def __init__(self, G=2.95912208286e-4, D=3):
        self.planets = []
        self.D = D  # number of spatial dimensions
        self.G = G  # gravitational constant
        self.timeline = [0]  # list of t-values, gets appended by evolving
        self.n_planets = 0

    # add planet objects to list, input np.array of matching size
    def add_planets(self, names, masses, initpos, initvels, fix):
        # pdb.set_trace()
        # sanity check: np.array as input
        try:
            n_add = np.shape(names)[0]
        except IndexError:
            sys.exit("ERROR: Input np.array to solar._init_")  # TODO: Djangofy

        # sanity check: number of arrays check out
        if not (n_add == np.shape(masses)[0] and n_add == np.shape(initpos)[0] and n_add == np.shape(initvels)[0]):
            sys.exit("ERROR: Number of planet mismatch in input data")  # TODO: Djangofy

        # sanity check: number of spatial dimensions
        if not (self.D == np.shape(initpos)[1] and self.D == np.shape(initvels)[1]):
            sys.exit("ERROR: Number of spatial dimensions mismatch, restrict to D=%d in input data" % self.D)  # TODO: Djangofy

        for i in range(n_add):
            self.planets.append(Planet(names[i], masses[i], initpos[i], initvels[i], fix[i]))
            self.n_planets += 1

'''
    def plot(self):
        """prints current path of planets
        """
        qs = [planet.position for planet in self.planets]
        fig = plt.figure(figsize=(12, 8))
        ax = fig.gca()
        for i in range(self.n_planets):
            q_planet = np.asarray(qs[i])
            plt.plot(q_planet[:, 0], q_planet[:, 1], "-", label=self.planets[i].name)

        plt.grid(True)
        plt.xlim(-50, 50)
        plt.ylim(-50, 50)
        ax.legend(loc="upper right")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        path = "solar_vv%d.pdf" % int(round(time()))
        plt.savefig(path)  # timestamp
        return path
'''


# contains information about planet
class Planet():

    def __init__(self, name, mass, initpos, initvel, fix):
        # Return a Planet object whose name is *name*, mass *mass*, initial position *initpos*, initial velocity *initvel*, current position *currentpos* and current velocity *currentvel*"""
        self.name = name
        self.fix = fix # 0, if position of planet fixed, 1 otherwise
        self.mass = mass
        self.initpos = initpos
        self.initvel = initvel
        self.position = [initpos]  # list of 1*3 dimensional numpy array containing position for ervery timestep
        self.velocity = [initvel]  # list of 1*3 dimensional numpy array containing velocity for ervery timestep

    def changeposvel(self, newpos, newvel):
        """updates position and velocity

        Input: newpos       ... np.array of new position
               newvel       ... np.array of new velocity
        """
        self.position.append(newpos)
        self.velocity.append(newvel)

    def x_as_plotly(self):
        """outputs x positions as strin "[1.222, 1.2226, .. ]"
        """
        xs = np.asarray(self.position)[:,0]
        return numpy_to_plotlystring(xs)

    def y_as_plotly(self):
        """outputs x positions as strin "[1.222, 1.2226, .. ]"
        """
        xs = np.asarray(self.position)[:,1]
        return numpy_to_plotlystring(xs)

    def z_as_plotly(self):
        """outputs x positions as strin "[1.222, 1.2226, .. ]"
        """
        xs = np.asarray(self.position)[:,2]
        return numpy_to_plotlystring(xs)
