"""
This module generates LAMMPS data file using ovito
"""

import numpy as np
import os

try:
    import ovito

    use_ovito = True
except:
    use_ovito = False

from polymerxtal.io import check_nanohub

def getnd():
    if use_ovito:
        node = ovito.io.import_file(".tmp/bonds/pdbfile.pdb")
        ovito.io.export_file(
            node, ".tmp/bonds/old.lmpdat", "lammps_data", atom_style="full"
        )
    
    elif check_nanohub():
    	os.system('ovitos getnd.py')
    
    else:
    	print("Cannot use function getnd - the package ovito is not installed or cannot be found.")

getnd()
