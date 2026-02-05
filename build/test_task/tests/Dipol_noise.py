import espressomd
import espressomd.magnetostatics
import numpy as np
import json
import os
import vtf
 
system = espressomd.System(box_l = [20.0, 20.0, 20.0], periodicity=[False, False, False])

sep = 3.0
 
system.time_step = 0.01 
system.cell_system.skin = 0.1
actor = espressomd.magnetostatics.DipolarDirectSumCpu(prefactor=1.)
system.actors.add(actor)

"""Setup for Brownian mechanics"""
 
system.thermostat.set_brownian(kT = 1.0, gamma = 6.0 * np.pi, gamma_rotation = 8.0 * np.pi, seed = 1)
system.integrator.set_brownian_dynamics(2.0)
wca = system.non_bonded_inter[0, 0]
wca.wca.set_params(epsilon = 1.0, sigma = 2.0)
 
"""json Output"""
 
trajectory_part1 = []
trajectory_part2 = []
trajectory_part1_all = {}
trajectory_part2_all = {}
center_to_center_dist = []
center_to_center_dist_all = {}
scaled_average_part_dist = []

freq_list = [0.0 + 0.01 * i for i in range(301)]
H0 = 10.0

for freq in freq_list: 

    system.time = 0.0
    system.part.clear()

    center = np.array([10.0, 10.0, 10.0])
    part1 = system.part.add(pos = center + [-sep/2, 0.0, 0.0], type = 0, rotation=[True, True, True], dip=[-1.0, 0.0, 0.0])
    part2 = system.part.add(pos = center + [sep/2, 0.0, 0.0], type = 0, rotation=[True, True, True], dip=[2.0, 0.0, 0.0])

    center_to_center_dist = []
    trajectory_part1 = []
    trajectory_part2 = []


    for i in range(4000):
        
        #Rotating magnetic field
        t = system.time
        omega = 2.0 * np.pi * freq
        H_field = np.array([H0 * np.cos(omega * t), H0 * np.sin(omega * t), 0.0])
        #Torque update due to field
        
        for p in [part1, part2]:
            
            p.ext_torque = np.cross(p.dip, H_field)
            
        system.integrator.run(1)
        
        trajectory_part1.append([float(x) for x in part1.pos])
        trajectory_part2.append([float(x) for x in part2.pos])
        
        r1 = np.array(part1.pos)
        r2 = np.array(part2.pos)
        center_to_center_dist.append(np.linalg.norm(r2 - r1))

        trajectory_part1_all[freq] = trajectory_part1
        trajectory_part2_all[freq] = trajectory_part2
        center_to_center_dist_all[freq] = center_to_center_dist

 

import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 4, figsize=(16, 10), sharex=False)

fig.suptitle(
    r"Analysis of magnetic colloidal particles with external rotating magnetic field",
    fontsize=16
)

traj1 = np.array(trajectory_part1_all[0.05])
traj2 = np.array(trajectory_part2_all[0.05])
x1, y1 = traj1[:,0], traj1[:,1]
x2, y2 = traj2[:,0], traj2[:,1]

ax = axes[0, 0]
ax.plot(x1, y1)
ax.plot(x2, y2)
ax.text(0.90, 0.05, "(a)", transform=ax.transAxes, fontsize=12, va="bottom")
ax.set_xlabel(r"$X/a$", fontsize=13)
ax.set_ylabel(r"$Y/a$", fontsize=13)

traj1 = np.array(trajectory_part1_all[0.1])
traj2 = np.array(trajectory_part2_all[0.1])
x1, y1 = traj1[:,0], traj1[:,1]
x2, y2 = traj2[:,0], traj2[:,1]

ax = axes[0, 1]
ax.plot(x1, y1)
ax.plot(x2, y2)
ax.text(0.90, 0.05, "(b)", transform=ax.transAxes, fontsize=12, va="bottom")
ax.set_xlabel(r"$X/a$", fontsize=13)
ax.set_ylabel(r"$Y/a$", fontsize=13)

traj1 = np.array(trajectory_part1_all[0.2])
traj2 = np.array(trajectory_part2_all[0.2])
x1, y1 = traj1[:,0], traj1[:,1]
x2, y2 = traj2[:,0], traj2[:,1]

ax = axes[0, 2]
ax.plot(x1, y1)
ax.plot(x2, y2)
ax.text(0.90, 0.05, "(c)", transform=ax.transAxes, fontsize=12, va="bottom")
ax.set_xlabel(r"$X/a$", fontsize=13)
ax.set_ylabel(r"$Y/a$", fontsize=13)

traj1 = np.array(trajectory_part1_all[1.0])
traj2 = np.array(trajectory_part2_all[1.0])
x1, y1 = traj1[:,0], traj1[:,1]
x2, y2 = traj2[:,0], traj2[:,1]

ax = axes[0, 3]
ax.plot(x1, y1)
ax.plot(x2, y2)
ax.text(0.90, 0.05, "(d)", transform=ax.transAxes, fontsize=12, va="bottom")
ax.set_xlabel(r"$X/a$", fontsize=13)
ax.set_ylabel(r"$Y/a$", fontsize=13)

scaled = np.array(center_to_center_dist_all[0.05])
time = [0.0 + 0.01 + i for i in range(4000)]

ax = axes[1, 0]
ax.plot(time, scaled)
ax.text(0.90, 0.05, "(e)", transform=ax.transAxes, fontsize=12, va="bottom")
ax.set_xlabel("Runtime (s)", fontsize=13)
ax.set_ylabel(r"$r/a$", fontsize=13)

scaled = np.array(center_to_center_dist_all[0.1])

ax = axes[1, 1]
ax.plot(time, scaled)
ax.text(0.90, 0.05, "(f)", transform=ax.transAxes, fontsize=12, va="bottom")
ax.set_xlabel("Runtime (s)", fontsize=13)
ax.set_ylabel(r"$r/a$", fontsize=13)

scaled = np.array(center_to_center_dist_all[0.2])

ax = axes[1, 2]
ax.plot(time, scaled)
ax.text(0.90, 0.05, "(g)", transform=ax.transAxes, fontsize=12, va="bottom")
ax.set_xlabel("Runtime (s)", fontsize=13)
ax.set_ylabel(r"$r/a$", fontsize=13)

scaled = np.array(center_to_center_dist_all[1.0])

ax = axes[1, 3]
ax.plot(time, scaled)
ax.text(0.90, 0.05, "(h)", transform=ax.transAxes, fontsize=12, va="bottom")
ax.set_xlabel("Runtime (s)", fontsize=13)
ax.set_ylabel(r"$r/a$", fontsize=13)


plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig("magnetic_analysis.png")
    
