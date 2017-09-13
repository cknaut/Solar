# Tests Universe Class and Evolution Function

#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from Solar import Solar


# Initial Positions for Planets
names = np.array(["Sun","Jupiter","Saturn","Uranus","Neptune","Pluto"])
ms = np.array([1.00000597682, 0.00095486104043,0.000285583733151,0.0000437273164546,0.0000517759138449,7.692307692307693e-09])
qs = np.array([[0,0,0],[-3.5023653, -3.8169847, -1.5507963],[9.0755314, -3.0458353, -1.6483708],[8.3101420, -16.2901086, -7.2521278],[11.4707666, -25.7294829, -10.8169456],[-15.5387357, -25.2225594, -3.1902382]])
vs = np.array([[0,0,0],[0.00565429, -0.00412490, -0.00190589],[0.00168318, 0.00483525, 0.00192462],[0.00354178, 0.00137102, 0.00055029],[0.00288930, 0.00114527, 0.00039677],[0.00276725, -0.00170702, -0.00136504]])

# big bang
Solar = Solar(True)

# add planets
Solar.add_planets(names, ms, qs, vs)


t_vv, y_vv = Solar.evolve(T = 20000, nrsteps = 2000, retarray = True)


fig = plt.figure(figsize=(12,8))
ax = fig.gca()
plt.plot(y_vv[:,0], y_vv[:,1], "b-", label="Sonne")
plt.plot(y_vv[:,3], y_vv[:,4], "g-", label="Jupiter")
plt.plot(y_vv[:,6], y_vv[:,7], "r-", label="Saturn")
plt.plot(y_vv[:,9], y_vv[:,10], "c-", label="Uranus")
plt.plot(y_vv[:,12], y_vv[:,13], "m-", label="Neptun")
plt.plot(y_vv[:,15], y_vv[:,16], "k-", label="Pluto")

'''
plt.plot(y_vv[-1,0], y_vv[-1,1], "bo")
plt.plot(y_vv[-1,3], y_vv[-1,4], "go")
plt.plot(y_vv[-1,6], y_vv[-1,7], "ro")
plt.plot(y_vv[-1,9], y_vv[-1,10], "co")
plt.plot(y_vv[-1,12], y_vv[-1,13], "mo")
plt.plot(y_vv[-1,15], y_vv[-1,16], "ko")
'''

for i in range(6):
    plt.plot(Solar.planets[i].currentpos[0], Solar.planets[i].currentpos[1], "o") # TODO: Updating in Planet Class Falsch, ueberdenken, ob Planet class notwendig

plt.grid(True)
plt.xlim(-50, 50)
plt.ylim(-50, 50)
ax.legend(loc="upper right")
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.savefig("solar_vv.pdf")

