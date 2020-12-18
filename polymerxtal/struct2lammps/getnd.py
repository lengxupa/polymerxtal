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
        node = ovito.io.import_file("./bonds/pdbfile.pdb")
        ovito.io.export_file(
            node, "./bonds/old.lmpdat", "lammps_data", atom_style="full"
        )
