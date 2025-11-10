import espressomd
 
import espressomd.io.writer.vtf
 
import numpy as np
 
import json
 
import os
 
system = espressomd.System(box_l = [10.0, 10.0, 10.0])
 
print("Please enter the desired separation in particle diameter(s) (sigma) between the particles.")
 
sep = float(input("> "))
 
res_path = os.path.join("results")
 
vtf_path = os.path.join(res_path, f"trajectory_{sep}.vtf")
 
json_path = os.path.join(res_path, f"trajectory_{sep}.json")
 
fp = open(vtf_path, mode="w")
 
system.time_step = 0.000001
 
"""Have to set the skin manually"""
 
system.cell_system.skin = 0.1
 
"""Triangle Setup"""
 
part1 = system.part.add(pos = [0.0, 0.0, 0.0], ext_force = [0.0, 0.0, 1.0])
 
part2 = system.part.add(pos = [sep, 0.0, 0.0], ext_force = [0.0, 0.0, 1.0])
 
part3 = system.part.add(pos = [sep/2, (np.sqrt(3)/2) * sep, 0], ext_force = [0.0, 0.0, 1.0])
 
"""Setup for Brownian mechanics"""
 
system.thermostat.set_brownian(kT = 0.0, gamma = 1.0, seed = 1)
 
system.integrator.set_brownian_dynamics(4.0)
 
"""VTF and json Output"""
 
trajectory_part1 = []
 
trajectory_part2 = []
 
trajectory_part3 = []
 
espressomd.io.writer.vtf.writevsf(system, fp)
 
espressomd.io.writer.vtf.writevcf(system, fp)
 
for i in range(150):
 
 trajectory_part1.append([float(x) for x in part1.pos])
 
 trajectory_part2.append([float(x) for x in part2.pos])
 
 trajectory_part3.append([float(x) for x in part3.pos])
 
 espressomd.io.writer.vtf.writevcf(system, fp)
 
 system.integrator.run(1)
    
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
