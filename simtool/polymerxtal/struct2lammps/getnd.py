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

use_nanohub = check_nanohub()

# Get the location of the current module
current_location = os.path.dirname(__file__)


def getnd():
    if use_ovito:
        node = ovito.io.import_file(".tmp/bonds/pdbfile.pdb")
        ovito.io.export_file(
            node, ".tmp/bonds/old.lmpdat", "lammps_data", atom_style="full"
        )

    elif use_nanohub:
        getnd_path = os.path.join(current_location, "getnd.py")
        os.system("ovitos %s" % getnd_path)

    else:
        print(
            "Cannot use function getnd - the package ovito is not installed or cannot be found."
        )


def main():
    getnd()


if __name__ == "__main__":
    main()
