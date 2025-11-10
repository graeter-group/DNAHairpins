import MDAnalysis as MDa
import numpy as np
import yaml
from pathlib import Path


def determine_end_to_end_indices(residue_indices_in, universe_in):
    hp_strand_5_prime_end = universe_in.residues[residue_indices_in[0]-1].atoms.select_atoms("name H5T")[0].index + 1
    nhp_strand_5_prime_end = universe_in.residues[residue_indices_in[1]-1].atoms.select_atoms("name H5T")[0].index + 1
    return hp_strand_5_prime_end, nhp_strand_5_prime_end

def determine_bonds(residue_indices_in_nhp, residue_indices_in_hp, universe_in, path_in):

    p_o5_prime_nhp = []
    o5prime_c5prime_nhp = []
    c5prime_c4prime_nhp = []
    c4prime_c3prime_nhp = []
    c3prime_o3prime_nhp = []
    o3prime_p_nhp = []

    for i in range(0, int(len(residue_indices_in_nhp)/2)):
        start = residue_indices_in_nhp[2*i]-1
        end = residue_indices_in_nhp[2*i+1]
        for residue_nhp_idx in range(start, end):
            cur_res = universe_in.residues[residue_nhp_idx]
            o5prime = cur_res.atoms.select_atoms("name O5'")[0].index + 1
            c5prime = cur_res.atoms.select_atoms("name C5'")[0].index + 1
            c4prime = cur_res.atoms.select_atoms("name C4'")[0].index + 1
            c3prime = cur_res.atoms.select_atoms("name C3'")[0].index + 1
            o3prime = cur_res.atoms.select_atoms("name O3'")[0].index + 1

            if residue_nhp_idx != end-1:
                next_p = universe_in.residues[residue_nhp_idx+1].atoms.select_atoms("name P")[0].index + 1
                o3prime_p_nhp.append((o3prime, next_p))

            if residue_nhp_idx != start:
                p = cur_res.atoms.select_atoms("name P")[0].index + 1
                p_o5_prime_nhp.append((p,o5prime))

            o5prime_c5prime_nhp.append((o5prime, c5prime))
            c5prime_c4prime_nhp.append((c5prime, c4prime))
            c4prime_c3prime_nhp.append((c4prime, c3prime))
            c3prime_o3prime_nhp.append((c3prime, o3prime))

    for bond_res, bond in enumerate(o5prime_c5prime_nhp):
        write_ndx_group(path_in, f"O5'-C5'.Res{bond_res+1}.NHP", bond)
    for bond_res, bond in enumerate(c5prime_c4prime_nhp):
        write_ndx_group(path_in, f"C5'-C4'.Res{bond_res+1}.NHP", bond)
    for bond_res, bond in enumerate(c4prime_c3prime_nhp):
        write_ndx_group(path_in, f"C4'-C3'.Res{bond_res+1}.NHP", bond)
    for bond_res, bond in enumerate(c3prime_o3prime_nhp):
        write_ndx_group(path_in, f"C3'-O3'.Res{bond_res+1}.NHP", bond)
    for bond_res, bond in enumerate(o3prime_p_nhp):
        write_ndx_group(path_in, f"O3'-P.Res{bond_res+1.5}.NHP", bond)
    for bond_res, bond in enumerate(p_o5_prime_nhp):
        write_ndx_group(path_in, f"P-O5'.Res{bond_res+2}.NHP", bond)


    p_o5_prime_hp = []
    o5prime_c5prime_hp = []
    c5prime_c4prime_hp = []
    c4prime_c3prime_hp = []
    c3prime_o3prime_hp = []
    o3prime_p_hp = []
    start = residue_indices_in_hp[0] -1
    end = residue_indices_in_hp[1]

    for residue_hp_idx in range(start, end):
        cur_res = universe_in.residues[residue_hp_idx]
        o5prime = cur_res.atoms.select_atoms("name O5'")[0].index + 1
        c5prime = cur_res.atoms.select_atoms("name C5'")[0].index + 1
        c4prime = cur_res.atoms.select_atoms("name C4'")[0].index + 1
        c3prime = cur_res.atoms.select_atoms("name C3'")[0].index + 1
        o3prime = cur_res.atoms.select_atoms("name O3'")[0].index + 1

        if residue_hp_idx != end - 1:
            next_p = universe_in.residues[residue_hp_idx + 1].atoms.select_atoms("name P")[0].index + 1
            o3prime_p_hp.append((o3prime, next_p))

        if residue_hp_idx != start:
            p = cur_res.atoms.select_atoms("name P")[0].index + 1
            p_o5_prime_hp.append((p, o5prime))

        o5prime_c5prime_hp.append((o5prime, c5prime))
        c5prime_c4prime_hp.append((c5prime, c4prime))
        c4prime_c3prime_hp.append((c4prime, c3prime))
        c3prime_o3prime_hp.append((c3prime, o3prime))


    for bond_res, bond in enumerate(o5prime_c5prime_hp):
        write_ndx_group(path_in, f"O5'-C5'.Res{bond_res + 1}.HP", bond)
    for bond_res, bond in enumerate(c5prime_c4prime_hp):
        write_ndx_group(path_in, f"C5'-C4'.Res{bond_res + 1}.HP", bond)
    for bond_res, bond in enumerate(c4prime_c3prime_hp):
        write_ndx_group(path_in, f"C4'-C3'.Res{bond_res + 1}.HP", bond)
    for bond_res, bond in enumerate(c3prime_o3prime_hp):
        write_ndx_group(path_in, f"C3'-O3'.Res{bond_res + 1}.HP", bond)
    for bond_res, bond in enumerate(o3prime_p_hp):
        write_ndx_group(path_in, f"O3'-P.Res{bond_res + 1.5}.HP", bond)
    for bond_res, bond in enumerate(p_o5_prime_hp):
        write_ndx_group(path_in, f"P-O5'.Res{bond_res + 2}.HP", bond)


def write_ndx_group(path, name, atom_indices_in):
    indices = np.asarray(atom_indices_in, dtype=int)
    with open(path, "a") as f:
        f.write(f"[ {name} ]\n")
        for i in range(0, len(indices), 15):
            f.write(" ".join(map(str, indices[i:i+15])) + "\n")
        f.write("\n")

with open("config.yaml") as f:
    cfg = yaml.safe_load(f)


data_dir = Path("../RawData/")
for gro_file in data_dir.glob("*.gro"):
    base_file_name = gro_file.stem
    print(base_file_name)
    system_name = base_file_name.split("_")[0]
    ndx_path_ete = f"../IndexFiles/{base_file_name}_EtE.ndx"
    ndx_path_bonds = f"../IndexFiles/{base_file_name}_BackboneBonds.ndx"

    u = MDa.Universe(gro_file)
    nhp_residue_indices = cfg["systems"][system_name]["nhpindices"][0]
    hp_residue_indices = cfg["systems"][system_name]["hpindices"][0]

    determine_bonds(nhp_residue_indices, hp_residue_indices, u, ndx_path_bonds)

    end_to_end_residue_indices = cfg["systems"][system_name]["endtoend"][0]
    end_to_end_indices = determine_end_to_end_indices(end_to_end_residue_indices, u)
    write_ndx_group(ndx_path_ete, "EndtoEndDistance", end_to_end_indices)