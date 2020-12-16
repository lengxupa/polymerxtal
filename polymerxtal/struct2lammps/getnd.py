"""
This module generates LAMMPS data file using ovito
"""

import numpy as np

try:
    import ovito

    use_ovito = True
except:
    use_ovito = False


def getnd():
    if use_ovito:
        node = ovito.io.import_file("bonds_pdbfile.pdb")
        ovito.io.export_file(node, "bonds_old.lmpdat", "lammps_data", atom_style="full")
