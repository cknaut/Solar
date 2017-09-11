# Generalized Hamiltonian N-Body Problem approximated by stoermer Verlet
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
import pdb
from time import time

# containg Planet objects as well as evolving functions and plotting functions
class Universe():

    def rhs(self, q, D,m):
        """Right hand side

        Input: q ...  np.array with positions

        Output: dp ... time-derivative of the velocities
        """       
        # Number of bodies
        N = self.n_planets
        # Empty np.arrays for computed data
        dp = np.zeros_like(q)
        #######################################################
        #same as rhs, but without * m[k] 
        for k in range(N):
            dp_k = np.zeros(D)
            for i in range(N):
                if k != i:
                    qi_qk = q[i*D:(i+1)*D] - q[k*D:(k+1)*D]
                    qi_qk /= linalg.norm(qi_qk)**3
                    dp_k += m[i]*qi_qk
            dp_k *= self.G
            dp[k*D:(k+1)*D] = dp_k  
        #######################################################
        return dp


    def integrate_VV(self, y0, xStart, xEnd, steps, m, flag=False):
        r"""Integrate ODE with velocity verlet rule

        Input: y0     ... initial condition
               xStart ... start x
               xEnd   ... end   x
               steps  ... number of steps (h = (xEnd - xStart)/N)
               flag   ... flag == False return complete solution: (phi, phi', t)
                          flag == True  return solution at endtime only: phi(tEnd)

        Output: x ... variable
                y ... solution
        """
        x = np.zeros(steps+1)
        y = np.zeros((steps+1, np.size(y0)))
        #############################################################      
        D=3    
        q_0, v_0 = np.hsplit(y0, 2) 
        x[0]=xStart  
        N = self.n_planets
        h = (float(xEnd - xStart))/float(steps) 
        #setting initial positions and speed
        y[0, 0:(N*D)] = q_0
        y[0, (N*D):2*(N*D)] = v_0
        for i in range(steps):
            y[i+1,0:(N*D)] = y[i,0:(N*D)] + h*y[i, (N*D):2*(N*D)]+0.5*h**2*self.rhs(y[i,0:(N*D)],D,m) 
            y[i+1, (N*D):2*(N*D)] = y[i, (N*D):2*(N*D)] + h/2*(self.rhs(y[i,0:(N*D)],D,m)+self.rhs(y[i+1,0:(N*D)],D,m))             
            x[i+1]+=(i+1)*h 
        #############################################################
        if flag:
            return x[-1], y[-1][:]
        else:
            return x, y


    # initialize with empty list of planets
    def __init__(self, G=2.95912208286e-4):
        self.planets = []
        self.G = G # gravitational constant
        self.n_planets = 0

    # add planet object to list
    def add_planet(self, name, mass, initpos, initvel):        
        self.planets.append(Planet(name, mass, initpos, initvel))
        self.n_planets += 1 


    # performs evolution of gravitational *n_planets-body*, endtime *T* and *nrsteps* steps, returns arrays with results
    # TODO: Update Current Position of Planets 
    def evolve(self, T, nrsteps):

        D = 3 # number of spatial dimensions

        # initialize global arrays
        m = np.zeros(self.n_planets)
        q0 = np.zeros((self.n_planets,D))
        p0 = np.zeros((self.n_planets,D))
        y0 = np.zeros(self.n_planets*D)
        v0 = np.zeros((self.n_planets,D))

        # fill arrays
        for i in range(self.n_planets):
            planet = self.planets[i]
            m[i] = planet.mass
            q0[i] = planet.currentpos
            v0[i] = planet.currentvel
            p0[i] = v0[i]*m[i]

        q0 = q0.flatten() # flatten to 1*(D*n_planets array)
        p0 = p0.flatten() # flatten to 1*(D*n_planets array)
        v0 = v0.flatten() # flatten to 1*(D*n_planets array)

        y0 = np.hstack([q0, p0]) # stack to 2*(1*(D*n_planets array)) array
        y0vv = np.hstack([q0, v0]) # stack to 2*(1*(D*n_planets array)) array  
        t_vv, y_vv = self.integrate_VV(y0vv, 0, T, nrsteps, m,  False)
        return t_vv, y_vv


# contains information about planet
class Planet():
    
    def __init__(self, name, mass, initpos, initvel):
        # Return a Planet object whose name is *name*, mass *mass*, initial position *initpos*, initial velocity *initvel*, current position *currentpos* and current velocity *currentvel*"""
        self.name = name
        self.mass = mass
        self.initpos = initpos
        self.initvel = initvel
        self.currentpos = initpos
        self.currentvel = initvel
