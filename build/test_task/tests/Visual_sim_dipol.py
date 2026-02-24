import espressomd
import espressomd.magnetostatics
import numpy as np
import json
import os
import vtf
 
system = espressomd.System(box_l = [100.0, 100.0, 100.0], periodicity=[False, False, False])
 
print("Please enter the desired separation in particle diameter(s) (sigma) between the particles.")
 
sep = float(input("> "))
 
res_path = os.path.join("results")
vtf_path = os.path.join(res_path, f"trajectory_dipol_noise_no_wca{sep}.vtf")
#json_path = os.path.join(res_path, f"trajectory_dipol_noise{sep}.json")
 
fp = open(vtf_path, mode="w")
 
system.time_step = 0.01 
system.cell_system.skin = 0.1
actor = espressomd.magnetostatics.DipolarDirectSumCpu(prefactor=1.)
system.actors.add(actor)
 
"""Triangle Setup"""

center = np.array([10.0, 10.0, 10.0])
part1 = system.part.add(pos = center + [-sep/2, 0.0, 0.0], type = 0, rotation=[True, True, True], dip=[-1.0, 0.0, 0.0])
part2 = system.part.add(pos = center + [sep/2, 0.0, 0.0], type = 0, rotation=[True, True, True], dip=[1.0, 0.0, 0.0])


"""Setup for Brownian mechanics"""
 
system.thermostat.set_brownian(kT = 1.0, gamma = 1, gamma_rotation = 1, seed = 1)
system.integrator.set_brownian_dynamics(2.0)
#wca = system.non_bonded_inter[0, 0]
#wca.wca.set_params(epsilon = 1.0, sigma = 2.0)
 
"""VTF and json Output"""
 
trajectory_part1 = []
trajectory_part2 = []

vtf._writevsf(system, fp, dipoles=True)
vtf._writevcf(system, fp, dipoles=True)

freq = 0.05
H0 = 10.0

for i in range(500):

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

    vtf._writevcf(system, fp, dipoles=True)
 
    
fp.close()
 
trajectory = {
    
    "part1": trajectory_part1,
    "part2": trajectory_part2
}
 
#with open(json_path, "w") as jf:
 
 #json.dump(trajectory, jf, indent = 2)
    
print("Finished. VTF file created.")

separ_path = os.path.join(res_path, "separations.json")
 
if os.path.exists(separ_path):
    
    with open(separ_path, "r") as tz:
        
        separations_data = json.load(tz)
 
 
else:
    
    separations_data = []
    
check_list = list(separations_data)
 
if sep not in separations_data:
    
    separations_data.append(sep)

if check_list != separations_data:
    
    with open(separ_path, "w") as t:
        
        json.dump(separations_data, t, indent = 2)
        
    print("Updated separations.json.")

else:
    
    print(f"{sep} already in separations.json, no changes made.")
