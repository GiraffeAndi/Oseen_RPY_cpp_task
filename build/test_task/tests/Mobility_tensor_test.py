open("tensor_analysis_file.csv", "w").close()

import espressomd
 
import espressomd.io.writer.vtf
 
import numpy as np
 
import json
 
import os

def load_diffusion_tensors(filename="tensor_analysis_file.csv"):
    tensors = []
    current = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith("#"):
                if current:
                    tensors.append(np.array(current, dtype=float))
                    current = []
            else:
                current.append([float(x) for x in line.split(",")])

        if current:
            tensors.append(np.array(current, dtype=float))

    return tensors
 
system = espressomd.System(box_l = [10.0, 10.0, 10.0])
 
system.time_step = 0.000001
 
"""Have to set the skin manually"""
 
system.cell_system.skin = 0.1
 
"""Setup for Brownian mechanics"""
 
system.thermostat.set_brownian(kT = 0.0, gamma = 6.0, gamma_rotation = 8.0, seed = 1)
 
system.integrator.set_brownian_dynamics(2.0)

rs = [2.0 + 0.1 * i for i in range(int((5.0 - 2.0) / 0.1) + 1)]
 
for r in rs:
    
    system.part.add(pos = [0.0, 0.0, 0.0])

    system.part.add(pos = [0.0, 0.0, r])
    
    system.integrator.run(1)

    system.part.clear()

tensors = load_diffusion_tensors()

mtt12zz = []

mtt12xx = []

mrr12zz = []

mrr12xx = []

mtr12xy = []

for M in tensors:
    
    mtt12xx.append(M[0, 3])

    mtt12zz.append(M[2, 5])

    mrr12zz.append(M[6, 9])

    mrr12xx.append(M[8, 11])

    mtr12xy.append(M[10, 0])

import matplotlib.pyplot as plt

plt.plot(rs, mtt12xx, label=r"$\perp$")

plt.plot(rs, mtt12zz, label=r"$\parallel$")

plt.legend()

plt.title("Compare to fig. 2.4b")

plt.ylim(0.02, 0.14)

plt.xlim(2., 5.)

ax = plt.gca()

ax.set_yticks([0.02, 0.05, 0.08, 0.11, 0.14])

ax2 = ax.twinx()

ax2.set_ylim(0.02, .14)

ax2.set_yticks([0.02, 0.05, 0.08, 0.11, 0.14])

plt.savefig("Tensor_analysis1")

plt.figure()

plt.plot(rs, mrr12xx, label=r"$\perp$")

plt.plot(rs, mrr12zz, label=r"$\parallel$")

plt.legend()

plt.title("Compare to fig. 2.4d")

plt.ylim(-0.02, 0.02)

plt.xlim(2., 5.)

ax = plt.gca()

ax.set_yticks([-0.02, -0.01, 0., 0.01, 0.02])

ax2 = ax.twinx()

ax2.set_ylim(-0.02, 0.02)

ax2.set_yticks([-0.02, -0.01, 0., 0.01, 0.02])

plt.savefig("Tensor_analysis2")

plt.figure()

plt.plot(rs, mtr12xy, label=r"$\perp$")

plt.legend()

plt.title("Compare to fig. 2.4f")

plt.ylim(0., 0.04)

plt.xlim(2., 5.)

ax = plt.gca()

ax.set_yticks([0., 0.01, 0.02, 0.03, 0.04])

ax2 = ax.twinx()

ax2.set_ylim(0., 0.04)

ax2.set_yticks([0., 0.01, 0.02, 0.03, 0.04])

plt.savefig("Tensor_analysis3")



