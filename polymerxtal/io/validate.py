"""
Functions to validate structure generated by Polymer Modeler.
"""
import numpy as np
import sys

try:
    from ovito.io import import_file

    use_ovito = True
except:  # noqa: E722
    from polymerxtal.data import atomic_radii

    use_ovito = False

from polymerxtal.polymod import readPDB


def readbond(bond_file):
    src = open(bond_file, "r")
    flag = 0
    bonds = []
    for line in src.readlines():
        ln = line.split()
        if ln:
            if ln[0] == "BONDS":
                flag = 1
                continue
            if flag:  # and eval(ln[1]) == 2:
                bonds.append([eval(ln[2]), eval(ln[3])])
    if not flag:
        raise Exception("%s file does not formatted correctly" % bond_file)
    src.close()
    return bonds


def validate_bonds(coords, path):
    """
    Validate all bonds are sucessfully connected

    Parameters
    ----------
    coords: dict or list of numpy-arrays
        atom coordinates
    path: str
        path to bonds.dat which contains bonds info

    Returns
    -------
    bool
        True if all the bonds are connected within 1.8 Angstrongs, False if otherwise
    """
    bonds = readbond(path)
    for bond in bonds:
        if np.linalg.norm(coords[bond[0] - 1] - coords[bond[1] - 1]) > 1.8:
            return False
    return True


def validate_coords(coords_path, bond_path):
    holder = readPDB(coords_path)
    bonds = readbond(bond_path)
    if use_ovito:
        pipeline = import_file(coords_path)
        types = pipeline.source.data.particles.particle_types
        for i in range(holder.num_atoms):
            for j in range(i + 1, holder.num_atoms):
                if ([i + 1, j + 1] not in bonds) and ([j + 1, i + 1] not in bonds):
                    if np.linalg.norm(holder.pos[i] - holder.pos[j]) <= (
                        types.type_by_name(holder.el_names[i]).radius
                        + types.type_by_name(holder.el_names[j]).radius
                    ):
                        return False
    else:
        for i in range(holder.num_atoms):
            for j in range(i + 1, holder.num_atoms):
                if ([i + 1, j + 1] not in bonds) and ([j + 1, i + 1] not in bonds):
                    if (holder.el_names[i] in atomic_radii) and (
                        holder.el_names[j] in atomic_radii
                    ):
                        if np.linalg.norm(holder.pos[i] - holder.pos[j]) <= (
                            atomic_radii[holder.el_names[i]]
                            + atomic_radii[holder.el_names[j]]
                        ):
                            return False
                    else:
                        if np.linalg.norm(holder.pos[i] - holder.pos[j]) < 1.5:
                            return False
    return True


def main(coords_path, bond_path):
    validate_coords(coords_path, bond_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
