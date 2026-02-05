import numpy as np

# -------------------------------
# VTF/VSF writer helpers (UNCHANGED)
# -------------------------------
def _vtf_pid_map(system, types="all"):
    if not hasattr(types, "__iter__") or isinstance(types, str):
        types = [types]
    ids = [p.id for p in system.part if ("all" in types or p.type in types)]
    return dict(zip(ids, range(len(ids))))


def _writevsf(system, fp, *, types="all", radii=None, dipoles=False):
    radii = radii or {}
    vtf_index = _vtf_pid_map(system, types)
    lx, ly, lz = system.box_l
    fp.write(f"unitcell {lx} {ly} {lz}\n")
    n_written = 0
    dip_radii = []
    for pid, vtf_id in vtf_index.items():
        p = system.part.by_id(pid)
        rad = radii.get(p.type, 0.5)
        fp.write(f"atom {vtf_id} radius {rad} name {p.type} type {p.type}\n")
        n_written += 1
        if dipoles and getattr(p, "dipm", 0.0) > 0.0:
            dip_radii.append(rad)
    if dipoles:
        for i, rad in enumerate(dip_radii):
            fp.write(f"atom {n_written + i} radius {rad} name 999 type 999\n")


def _writevcf(system, fp, *, types="all", folded=False, dipoles=False):
    vtf_index = _vtf_pid_map(system, types)
    fp.write("\ntimestep indexed\n")
    dip_positions = []
    for pid, vtf_id in vtf_index.items():
        p = system.part.by_id(pid)
        pos = p.pos_folded if folded else p.pos
        fp.write(f"{vtf_id} {pos[0]} {pos[1]} {pos[2]}\n")
        if dipoles and getattr(p, "dipm", 0.0) > 0.0:
            tip = pos + np.array(p.dip) * 0.02
            dip_positions.append(tip)
    if dipoles:
        offset = len(vtf_index)
        for i, tip in enumerate(dip_positions):
            fp.write(f"{offset + i} {tip[0]} {tip[1]} {tip[2]}\n")
