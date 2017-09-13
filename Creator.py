# Tests Universe Class and Evolution Function

#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from Solar import Solar



names = np.array(["Sun","Jupiter","Uranus","Neptune","Pluto"])
ms = np.array([1.00000597682, 0.00095486104043,0.000285583733151,0.0000437273164546,0.0000517759138449,7.692307692307693e-09])
qs = np.array([[0,0,0],[-3.5023653, -3.8169847, -1.5507963],[9.0755314, -3.0458353, -1.6483708],[8.3101420, -16.2901086, -7.2521278],[11.4707666, -25.7294829, -10.8169456],[-15.5387357, -25.2225594, -3.1902382]])
ps = np.array([[0,0,0],[0.00565429, -0.00412490, -0.00190589],[0.00168318, 0.00483525, 0.00192462],[0.00354178, 0.00137102, 0.00055029],[0.00288930, 0.00114527, 0.00039677],[0.00276725, -0.00170702, -0.00136504]])
# Initial Positions for Planets
namesun = "Sun"
msun = 1.00000597682
qsun = np.array([0,0,0])
vsun = np.array([0,0,0])


namej = "Jupiter"
mj = 0.00095486104043
qj = np.array([-3.5023653, -3.8169847, -1.5507963])
vj = np.array([0.00565429, -0.00412490, -0.00190589])

names = "Saturn"
ms = 0.000285583733151
qs = np.array([9.0755314, -3.0458353, -1.6483708])
vs = np.array([0.00168318, 0.00483525, 0.00192462])

nameu = "Uranus"
mu = 0.0000437273164546
qu = np.array([8.3101420, -16.2901086, -7.2521278])
vu = np.array([0.00354178, 0.00137102, 0.00055029])

namen = "Neptune"
mn = 0.0000517759138449
qn = np.array([11.4707666, -25.7294829, -10.8169456])
vn = np.array([0.00288930, 0.00114527, 0.00039677])

namep = "Pluto"
mp = 7.692307692307693e-09
qp = np.array([-15.5387357, -25.2225594, -3.1902382])
vp = np.array([0.00276725, -0.00170702, -0.00136504])

# big bang
Solar = Solar()

Solar.add_planet(namesun, msun, qsun, vsun)
Solar.add_planet(namej, mj, qj, vj)
Solar.add_planet(names, ms, qs, vs)
Solar.add_planet(nameu, mu, qu, vu)
Solar.add_planet(namen, mn, qn, vn)
Solar.add_planet(namep, mp, qp, vp)

t_vv, y_vv = Solar.evolve(T = 20000, nrsteps = 2000)


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

