# Generalized Hamiltonian N-Body Problem approximated by stoermer Verlet
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
import pdb
import sys
from time import time

# containg Planet objects as well as evolving functions and plotting functions
class Solar():

    def rhs(self, q, m):
        """Right hand side of ODE
        Input: q ...  np.array with positions
        Output: dp ... time-derivative of the velocities
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


    def integrate_VV(self, y0, tStart, tEnd, steps, m):
        r"""Integrate ODE with velocity verlet rule

        Input: y0     ... initial condition
               tStart ... start t
               tEnd   ... end   t
               steps  ... number of steps (h = (tEnd - tStart)/N)
               flag   ... flag == False return complete solution: (phi, phi', t)
                          flag == True  return solution at endtime only: phi(tEnd)

        Output: t ... variable
                y ... solution
        """
        t = np.zeros(steps+1)
        y = np.zeros((steps+1, np.size(y0)))
        #############################################################       
        q_0, v_0 = np.hsplit(y0, 2) 
        t[0]=tStart  
        N = self.n_planets
        h = (float(tEnd - tStart))/float(steps) 
        #setting initial positions and speed
        y[0, 0:(N*self.D)] = q_0
        y[0, (N*self.D):2*(N*self.D)] = v_0
        for i in range(steps):
            y[i+1,0:(N*self.D)] = y[i,0:(N*self.D)] + h*y[i, (N*self.D):2*(N*self.D)]+0.5*h**2*self.rhs(y[i,0:(N*self.D)],m) 
            y[i+1, (N*self.D):2*(N*self.D)] = y[i, (N*self.D):2*(N*self.D)] + h/2*(self.rhs(y[i,0:(N*self.D)],m)+self.rhs(y[i+1,0:(N*self.D)],m))             
            t[i+1]+=(i+1)*h 
        #############################################################
        return t, y

    def integrate_VV_one_step(self, y0, m):
            r"""Integrate ODE with velocity verlet rule

            Input: y0     ... initial condition
       

            Output:
                    y ... solution of y after one propagation
            """
            y = np.zeros_like(y0)
            y[0:(N*self.D)] = y0[0:(N*self.D)] + h*y0[(N*self.D):2*(N*self.D)]+0.5*h**2*self.rhs(y0[0:(N*self.D)],m) 
            y[(N*self.D):2*(N*self.D)] = y0[(N*self.D):2*(N*self.D)] + h/2*(self.rhs(y0[0:(N*self.D)],m)+self.rhs(y[0:(N*self.D)],m))            
            
            return t, y

    # initialize with empty list of planets
    def __init__(self, G=2.95912208286e-4, D=3):
        self.planets = []
        self.D = D # number of spatial dimensions
        self.G = G # gravitational constant
        self.n_planets = 0

    # add planet objects to list, input np.array of matching size
    def add_planets(self, names, masses, initpos, initvels):  
        #pdb.set_trace()
        
        # sanity check: np.array as input         
        try:
            n_add = np.shape(names)[0]
        except IndexError:
            sys.exit("ERROR: Input np.array to solar._init_") # TODO: Djangofy

        # sanity check: number of arrays check out
        if not (n_add == np.shape(masses)[0] and n_add == np.shape(initpos)[0] and n_add == np.shape(initvels)[0]):
             sys.exit("ERROR: Number of planet mismatch in input data") # TODO: Djangofy

        # sanity check: number of spatial dimensions
        if not (self.D == np.shape(initpos)[1] and self.D == np.shape(initvels)[1]):
             sys.exit("ERROR: Number of spatial dimensions mismatch, restrict to D=%d in input data" % self.D) # TODO: Djangofy
      
        for i in range(n_add):
            self.planets.append(Planet(names[i], masses[i], initpos[i], initvels[i]))
            self.n_planets += 1 
    

    # extracts mass array from current planet list
    def get_m(self):
        m = np.zeros(self.n_planets)
        for i in range(self.n_planets):
            m[i] = self.planets[i].mass
        return m


    # extracts yvv = [q,v]^T from current planet list
    def get_y(self):

        # initialize global arrays
        q = np.zeros((self.n_planets,self.D)) # positions
        v = np.zeros((self.n_planets,self.D)) # velocities
        p = np.zeros((self.n_planets,self.D)) # momentum


        # fill arrays
        for i in range(self.n_planets):
            planet = self.planets[i]
            q[i] = planet.currentpos
            v[i] = planet.currentvel

        return np.hstack([q.flatten(), v.flatten()]) # flatten and stack


        


    # performs evolution of gravitational *n_planets-body*, endtime *T* and *nrsteps* steps, returns arrays with results
    # TODO: Update Current Position of Planets 
    def evolve(self, T, nrsteps):
        y = self.get_y()
        m = self.get_m()
        t_vv, y_vv = self.integrate_VV(y, 0, T, nrsteps, m)
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
