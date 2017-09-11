# Generalized Hamiltonian N-Body Problem approximated by stoermer Verlet


#!/usr/bin/env python
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
from time import time


def rhs(q, D=3):
    """Right hand side

    Input: q ...  np.array with positions

    Output: dp ... time-derivative of the velocities
    """
    # Number of bodies
    N = q.size // D
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
        dp_k *= G
        dp[k*D:(k+1)*D] = dp_k  
    #######################################################
    return dp




def integrate_VV(y0, xStart, xEnd, steps, flag=False):
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
    N = np.size(q_0) // D
    h = (float(xEnd - xStart))/float(steps) 
    #setting initial positions and speed
    y[0, 0:(N*D)] = q_0
    y[0, (N*D):2*(N*D)] = v_0
    for i in range(steps):
        y[i+1,0:(N*D)] = y[i,0:(N*D)] + h*y[i, (N*D):2*(N*D)]+0.5*h**2*rhs(y[i,0:(N*D)]) 
        y[i+1, (N*D):2*(N*D)] = y[i, (N*D):2*(N*D)] + h/2*(rhs(y[i,0:(N*D)])+rhs(y[i+1,0:(N*D)]))             
        x[i+1]+=(i+1)*h 
    #############################################################
    if flag:
        return x[-1], y[-1][:]
    else:
        return x, y


# Unteraufgabe f)

G = 2.95912208286e-4

# Anfangswerte der Planeten
msun = 1.00000597682
qsun = np.array([0,0,0])
vsun = np.array([0,0,0])

mj = 0.00095486104043
qj = np.array([-3.5023653, -3.8169847, -1.5507963])

vj = np.array([0.00565429, -0.00412490, -0.00190589])

ms = 0.000285583733151
qs = np.array([9.0755314, -3.0458353, -1.6483708])
vs = np.array([0.00168318, 0.00483525, 0.00192462])

mu = 0.0000437273164546
qu = np.array([8.3101420, -16.2901086, -7.2521278])
vu = np.array([0.00354178, 0.00137102, 0.00055029])

mn = 0.0000517759138449
qn = np.array([11.4707666, -25.7294829, -10.8169456])
vn = np.array([0.00288930, 0.00114527, 0.00039677])

mp = 7.692307692307693e-09
qp = np.array([-15.5387357, -25.2225594, -3.1902382])
vp = np.array([0.00276725, -0.00170702, -0.00136504])

m = np.array([msun, mj, ms, mu, mn, mp])
q0 = np.hstack([qsun, qj, qs, qu, qn, qp])
p0 = np.hstack([msun*vsun, mj*vj, ms*vs, mu*vu, mn*vn, mp*vp])
y0 = np.hstack([q0, p0])

v0 = np.hstack([vsun, vj, vs, vu, vn, vp])
y0vv = np.hstack([q0, v0])

T = 20000
nrsteps = 2000


starttime = time()
t_vv, y_vv = integrate_VV(y0vv, 0, T, nrsteps, False)
endtime = time()
print('VV needed %f seconds for %f steps' % (endtime-starttime, nrsteps))

fig = plt.figure(figsize=(12,8))
ax = fig.gca()
plt.plot(y_vv[:,0], y_vv[:,1], "b-", label="Sonne")
plt.plot(y_vv[:,3], y_vv[:,4], "g-", label="Jupiter")
plt.plot(y_vv[:,6], y_vv[:,7], "r-", label="Saturn")
plt.plot(y_vv[:,9], y_vv[:,10], "c-", label="Uranus")
plt.plot(y_vv[:,12], y_vv[:,13], "m-", label="Neptun")
plt.plot(y_vv[:,15], y_vv[:,16], "k-", label="Pluto")

plt.plot(y_vv[-1,0], y_vv[-1,1], "bo")
plt.plot(y_vv[-1,3], y_vv[-1,4], "go")
plt.plot(y_vv[-1,6], y_vv[-1,7], "ro")
plt.plot(y_vv[-1,9], y_vv[-1,10], "co")
plt.plot(y_vv[-1,12], y_vv[-1,13], "mo")
plt.plot(y_vv[-1,15], y_vv[-1,16], "ko")
plt.grid(True)
plt.xlim(-50, 50)
plt.ylim(-50, 50)
ax.legend(loc="upper right")
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.savefig("solar_vv.pdf")

