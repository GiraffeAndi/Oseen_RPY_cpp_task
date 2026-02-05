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
sWs = np.array([2.01, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0])
vWs = np.array([
1.79223228,
1.79070892,
1.77826951,
1.76311004,
1.74703222,
1.73072508,
1.71452988,
1.68314045,
1.65356957,
1.62599574,
1.6003915,
1.57664602,
1.55461831,
1.53416156
])
wWs = np.array([
0.15960749,
0.195478184,
0.209731929,
0.214756269,
0.215412007,
0.213577673,
0.210224103,
0.201041187,
0.190411103,
0.179481457,
0.16880121,
0.158630515,
0.149081378,
0.140187608
])
 
system = espressomd.System(box_l = [10.0, 10.0, 10.0])
system.time_step = 0.000001
system.cell_system.skin = 1
system.thermostat.set_brownian(kT = 0.0, gamma = 6.0 * np.pi, gamma_rotation = 8.0 * np.pi, seed = 1)
system.integrator.set_brownian_dynamics(2.0)

rs = [2.0 + 0.1 * i for i in range(11)]
 
for r in rs:

    part1 = system.part.add(pos = [0.0, 0.0, 0.0], ext_force = [0.0, 0.0, 1.0])
    part2 = system.part.add(pos = [r, 0.0, 0.0], ext_force = [0.0, 0.0, 1.0])
    part3 = system.part.add(pos = [r/2, (np.sqrt(3)/2) * r, 0], ext_force = [0.0, 0.0, 1.0])
    
    system.integrator.run(1)
    system.part.clear()

vectors = load_vectors() 

print(f"{vectors}")

V0 = 1.0/(6.0 * np.pi) 
W0 = V0 

VV0s = [(v[2] / V0) for v in vectors]
#Omega_x in source
WW0s = [(np.linalg.norm(v[-3]) / W0) for v in vectors]


plt.figure()
plt.plot(sWs, vWs, "o:", label = "Wilson")
plt.plot(rs, VV0s, label = "RPY")
plt.title("Compare to fig. 2 in [2014-wilson-jcp] or fig. 1a in [2025-torres-arxiv]")
plt.xlabel("$s / a$")
plt.ylabel(r"$M_{tt}\ (V_z / V_0)$")
plt.legend()
plt.savefig("velocity_analysis1")
plt.figure()
plt.plot(sWs, wWs, ":o", label = "Wilson")
plt.plot(rs, WW0s, label = "RPY")
plt.title("Compare to fig. 2 in [2014-wilson-jcp] or fig. 1b in [2025-torres-arxiv]")
plt.xlabel("$s / a$")
plt.ylabel(r"$M_{rt}\ (\Omega_x/\Omega_0)$")
plt.legend()
plt.savefig("velocity_analysis2")

print("Success!")