open("velocity_analysis_file.csv", "w").close()

import espressomd
import espressomd.io.writer.vtf
import numpy as np
import json 
import os
import matplotlib.pyplot as plt

def load_vectors(filename="velocity_analysis_file.csv"): 
    vectors = [] 
    
    with open(filename) as f: 
        for line in f: 
            line = line.strip() 
            
            if not line or line.startswith("#"): 
                continue 
        
            vec = np.array([float(x) for x in line.split(",")]) 
            vectors.append(vec) 

    return vectors

# tabulated data from [2014-wilson-jcp]
sWs = np.array([2.01, 2.10, 2.5, 3.0, 4.0, 6.0])
vWs = np.array([
0.65528,
0.73857,
0.87765,
0.93905,
0.97964,
0.99581,
])
wWs = np.array([
0.037336,
0.052035,
0.045466,
0.035022,
0.021634,
0.010159,
])
vWs2 = np.array([
0.63461,
0.59718,
0.49545,
0.41694,
0.31859,
0.21586,
])
vWs3 = np.array([
0.00498,
0.03517,
0.07393,
0.07824,
0.06925,
0.05078,
])
 
system = espressomd.System(box_l = [20.0, 20.0, 20.0])
system.time_step = 0.000001
system.cell_system.skin = 1
system.thermostat.set_brownian(kT = 0.0, gamma = 6.0 * np.pi, gamma_rotation = 8.0 * np.pi, seed = 1)
system.integrator.set_brownian_dynamics(2.0)

rs = [2.0, 2.01, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0]
 
for r in rs:

    #part1 = system.part.add(pos = [0.0, 0.0, 0.0], ext_force=[6*np.pi * (np.sqrt(3)/2), 6*np.pi * 0.5, 0.0])
    part1 = system.part.add(pos = [0.0, 0.0, 0.0])
    part2 = system.part.add(pos = [r, 0.0, 0.0])
    part3 = system.part.add(pos = [r/2, (np.sqrt(3)/2)*r, 0.0], ext_force=[0.0, -6*np.pi, 0.0])

    system.integrator.run(1)
    system.part.clear()

vectors = load_vectors() 
F = 6.0*np.pi
V0 = F/(6.0 * np.pi) 
W0 = V0 

VV0s = [abs(v[7]) / V0 for v in vectors]
WW0s = [np.linalg.norm(v[12:15] - v[9:12])/(2*W0) for v in vectors]
#parallel component of no force particles
VV02s = [abs(v[1] + v[4])/(2*V0) for v in vectors]
#orthogonal component of no force particles
VV03s = [(v[3] - v[0])/(2*V0) for v in vectors]


plt.figure()
plt.plot(sWs, vWs, "o:", label = "Wilson")
plt.plot(rs, VV0s, label = "RPY")
plt.title("Compare to fig. 2a in [2025-torres-arxiv]")
plt.xlabel("$s / a$")
plt.ylabel(r"$(V_x / V_0)$")
plt.legend()
plt.savefig("velocity_analysis1_one_force")
plt.figure()
plt.plot(sWs, wWs, ":o", label = "Wilson")
plt.plot(rs, WW0s, label = "RPY")
plt.title("Compare to fig. 2d in [2025-torres-arxiv]")
plt.xlabel("$s / a$")
plt.ylabel(r"$(\Omega_x/\Omega_0)$")
plt.legend()
plt.savefig("velocity_analysis2_one_force")
plt.figure()
plt.plot(sWs, vWs2, "o:", label = "Wilson")
plt.plot(rs, VV02s, label = "RPY")
plt.title("Compare to fig. 2b in [2025-torres-arxiv]")
plt.xlabel("$s / a$")
plt.ylabel(r"$(V_\parallel / V_0)$")
plt.legend()
plt.savefig("velocity_analysis3_one_force")
plt.figure()
plt.plot(sWs, vWs3, "o:", label = "Wilson")
plt.plot(rs, VV03s, label = "RPY")
plt.title("Compare to fig. 2c in [2025-torres-arxiv]")
plt.xlabel("$s / a$")
plt.ylabel(r"$(V_\perp / V_0)$")
plt.legend()
plt.savefig("velocity_analysis4_one_force")

print("Success!")