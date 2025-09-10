import espressomd

import espressomd.io.writer.vtf

import numpy as np

import json

system = espressomd.System(box_l = [10.0, 10.0, 10.0])

fp = open("trajectory_default_espresso.vtf", mode="w")

system.time_step = 0.1

"""Have to set the skin manually"""

system.cell_system.skin = 0.1

"""Triangle Setup"""

part1 = system.part.add(pos = [0.0, 0.0, 0.0], ext_force = [0.0, 0.0, 1.0])

part2 = system.part.add(pos = [2.5, 0.0, 0.0], ext_force = [0.0, 0.0, 1.0])

part3 = system.part.add(pos = [1.25, (np.sqrt(3)/2) * 2.5, 0], ext_force = [0.0, 0.0, 1.0])

"""Setup for Brownian mechanics"""

system.thermostat.set_brownian(kT = 0.0, gamma = 1.0, seed = 1)

system.integrator.set_brownian_dynamics()

"""VTF and json Output"""

trajectory_part1 = []

trajectory_part2 = []

trajectory_part3 = []

espressomd.io.writer.vtf.writevsf(system, fp)

espressomd.io.writer.vtf.writevcf(system, fp)

for i in range(50):

    system.integrator.run(1)

    trajectory_part1.append([float(x) for x in part1.pos])

    trajectory_part2.append([float(x) for x in part2.pos])

    trajectory_part3.append([float(x) for x in part3.pos])

    espressomd.io.writer.vtf.writevcf(system, fp)

fp.close()

trajectory = {
        "part1": trajectory_part1,

        "part2": trajectory_part2,

        "part3": trajectory_part3
        }

with open("trajectory_default_espresso.json", "w") as jf:

    json.dump(trajectory, jf, indent=2)

print("Finished. VTF and JSON files created.")






