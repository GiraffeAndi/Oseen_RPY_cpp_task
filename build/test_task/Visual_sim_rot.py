import espressomd
import numpy as np
import json
import os
import vtf
 
system = espressomd.System(box_l = [20.0, 20.0, 20.0])
 
print("Please enter the desired separation in particle diameter(s) (sigma) between the particles.")
 
sep = float(input("> "))
 
res_path = os.path.join("results")
vtf_path = os.path.join(res_path, f"trajectory_rot_{sep}.vtf")
json_path = os.path.join(res_path, f"trajectory_rot_{sep}.json")
 
fp = open(vtf_path, mode="w")
 
system.time_step = 0.01 
system.cell_system.skin = 0.1
 
"""Triangle Setup"""

center = np.array([10.0, 10.0, 10.0])

part1 = system.part.add(pos = center + [-sep, 0.0, 0.0], type = 0, rotation=[True, True, True], ext_force=[0.5, 0.0, 0.0], ext_torque=[0.0, 0.0, 50.0], dip=[-1.0, 0.0, 0.0])
part2 = system.part.add(pos = center + [sep, 0.0, 0.0], type = 0, rotation=[True, True, True],  ext_force=[0.0, 0.5, 0.0], ext_torque=[0.0, 0.0, 30.0], dip=[1.0, 0.0, 0.0])
part3 = system.part.add(pos = center + [0, sep, 0], rotation=[True, True, True], type = 0, ext_force=[0.0, 0.0, 0.5], ext_torque=[20.0, 0.0, 0.0], dip=[0.0, 1.0, 0.0]) 

"""Setup for Brownian mechanics"""
 
system.thermostat.set_brownian(kT = 0.0, gamma = 1, gamma_rotation = 1, seed = 1)
system.integrator.set_brownian_dynamics(2.0)
wca = system.non_bonded_inter[0, 0]
wca.wca.set_params(epsilon = 1.0, sigma = 2.0)
 
"""VTF and json Output"""
 
trajectory_part1 = []
trajectory_part2 = []
trajectory_part3 = []

vtf._writevsf(system, fp, dipoles=True)
vtf._writevcf(system, fp, dipoles=True)

for i in range(500):
    
    system.integrator.run(1)

    trajectory_part1.append([float(x) for x in part1.pos])
    trajectory_part2.append([float(x) for x in part2.pos])
    trajectory_part3.append([float(x) for x in part3.pos])

    vtf._writevcf(system, fp, dipoles=True)
 
    
fp.close()
 
trajectory = {
    
    "part1": trajectory_part1,
    "part2": trajectory_part2,
    "part3": trajectory_part3
}
 
with open(json_path, "w") as jf:
 
 json.dump(trajectory, jf, indent = 2)
    
print("Finished. VTF and JSON files created.")

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
