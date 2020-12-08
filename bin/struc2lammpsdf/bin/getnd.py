import numpy as np, ovito

node = ovito.io.import_file('./bonds/pdbfile.pdb')
ovito.io.export_file(node, './bonds/old.lmpdat', 'lammps_data', atom_style='full')
