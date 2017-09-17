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

    def integrate_VV_one_step(self, y0, h, m):
            """Integrate ODE with velocity verlet rule

            Input: y0     ... initial condition
                   h      ... timestep

            Output:
                    y ... solution of y after one propagation
            """
            y = np.zeros_like(y0)
            y[0:(self.n_planets*self.D)] = y0[0:(self.n_planets*self.D)] + h*y0[(self.n_planets*self.D):2*(self.n_planets*self.D)]+0.5*h**2*self.rhs(y0[0:(self.n_planets*self.D)],m) 
            y[(self.n_planets*self.D):2*(self.n_planets*self.D)] = y0[(self.n_planets*self.D):2*(self.n_planets*self.D)] + h/2*(self.rhs(y0[0:(self.n_planets*self.D)],m)+self.rhs(y[0:(self.n_planets*self.D)],m))            
            self.y_to_planets(y) # keep track of planet position on every iteration
            return y

    # TODO: Possible Speedup: Directly work with planet position as saved in planet classes, without y as intermediary
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
            y[i+1] = self.integrate_VV_one_step(y[i],h,m)                     
            t[i+1]+=(i+1)*h 
        #############################################################
        return t, y


    # initialize with empty list of planets
    # animate is used to invoc matplotlib animation
    def __init__(self, animate=False, G=2.95912208286e-4, D=3):
        self.animate = animate
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
            self.planets.append(Planet(names[i], masses[i], initpos[i], initvels[i], self.animate))
            self.n_planets += 1             
    

    # extracts mass array from current planet list
    def get_m(self):
        m = np.zeros(self.n_planets)
        for i in range(self.n_planets):
            m[i] = self.planets[i].mass
        return m


    # extracts yvv = [q,v]^T from current planet list
    def planet_to_y(self):

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

    # updates planets position according to y, inverse function of 
    def y_to_planets(self, y):
        q,v = np.hsplit(y,2)
        for i in range(self.n_planets):
            planet = self.planets[i]
            planet.changepos(q[i:i+3])
            planet.changevel(v[i:i+3])
           


    # performs evolution of gravitational *n_planets-body*, endtime *T* and *nrsteps* steps, returns arrays with results
    # TODO: Update Current Position of Planets 
    def evolve(self, T, nrsteps, retarray):
        """erforms evolution of gravitationa

        Input: T            ... endtime
               nrsteps      ... number of timesteps
                retarray    ... if TRUE: return array of planet path, if FALSE: Only update planet objects

        Output:
            t, y: time and planet path
        """
        y = self.planet_to_y()
        m = self.get_m()
        if retarray:
            t, y = self.integrate_VV(y, 0, T, nrsteps, m) # TODO: Possible Speedup: Inhibit creation of path arrays altogether, if retarray false
            return t, y
        else:
            self.integrate_VV(y, 0, T, nrsteps, m)
            return
            

# contains information about planet
class Planet():
    
    def __init__(self, name, mass, initpos, initvel, animate=False):
        # Return a Planet object whose name is *name*, mass *mass*, initial position *initpos*, initial velocity *initvel*, current position *currentpos* and current velocity *currentvel*"""
        self.name = name
        self.mass = mass
        self.initpos = initpos
        self.initvel = initvel
        self.position = [initpos] # list of 1*3 dimensional numpy array containing position for ervery timestep
        self.velocity = [initvel] # list of 1*3 dimensional numpy array containing velocity for ervery timestep
        self.animate = animate
        

    def changeposvel(self, newpos, newvel):
        """updates position and velocity 

        Input: newpos       ... np.array of new position
               newvel       ... np.array of new velocity
        """
        self.position.append(newpos)
        self.velocity.append(newpos)

 

    def changevel(self, newvel):
        self.currentvel = newvel
