"""
Functions for calculating monomer properties.
"""
import os

from polymerxtal.io import run_polymod, read_latchout
from polymerxtal.polymod import readPDB

from .molecule import build_bond_list

# Get the location of the current module
current_location = os.path.dirname(__file__)


class Backbone:
    def __init__(self, name, path, head_index, tail_index, tail_torsion):
        self.name = name
        self.path = path
        self.head_index = head_index
        self.tail_index = tail_index
        self.tail_torsion = tail_torsion


def get_tail_torsion(mpath, head, tail):
    logfile = ".tmp/latch.log"
    infile = ".tmp/latch.in"
    des = open(infile, "w")
    des.write("bond_scale  1.2\n")
    des.write("monomer m1 %s %d %d\n" % (mpath, head, tail))
    des.write("torsion all fixed \n")
    des.write("chains 1\n")
    des.write("monomers 1\n")
    des.write("system_min 0 0 0\n")
    des.write("system_max 100 100 100\n")
    des.write("grid_size 5\n")
    des.write("stereo s pattern 1 m1\n")
    des.write("chain_stereo 1 s 1\n")
    des.write("sample none\n")
    des.write("energy none\n")
    des.write("write intermediate\n")
    des.write("log_file %s\n" % logfile)
    des.close()
    run_polymod(infile)
    latch_info = read_latchout(logfile)
    tail_torsion = eval(latch_info["monomer"][1]["info"][3][-1])
    return tail_torsion


class PolymerType:
    def __init__(self, name, path, backbone_atoms=[], side_atom=0):
        self.name = name
        self.path = path
        self.backbone_atoms = backbone_atoms
        self.side_atom = side_atom

    def __str__(self):
        return self.name

    def identify_backbone_path(self, find_inverse_path=False):
        h = readPDB(self.path)
        backbone_el = [h.el_names[i - 1] for i in self.backbone_atoms]

        if (
            find_inverse_path
            and (not self.side_atom)
            and backbone_el == backbone_el[::-1]
        ):
            raise Exception(
                "Polymer type %s does not have defects for head to head or tail to tail connections"
                % (self.name)
            )

        bonds = build_bond_list(h.pos, elements=h.el_names)

        center_head = self.backbone_atoms[0]
        center_tail = self.backbone_atoms[-1]
        if (
            self.side_atom
            and ((center_tail - 1, self.side_atom - 1) not in bonds)
            and ((self.side_atom - 1, center_tail - 1) not in bonds)
        ):
            raise ValueError("side atom must connect to tail backbone atom")

        head_atoms = []
        tail_atoms = []

        for key in bonds:
            atom1, atom2 = key
            if self.side_atom and (self.side_atom - 1 in [atom1, atom2]):
                continue
            elif center_head - 1 == atom1:
                if not center_tail - 1 == atom2:
                    head_atoms.append(atom2 + 1)
            elif center_head - 1 == atom2:
                if not center_tail - 1 == atom1:
                    head_atoms.append(atom1 + 1)
            elif center_tail - 1 == atom1:
                if not center_head - 1 == atom2:
                    tail_atoms.append(atom2 + 1)
            elif center_tail - 1 == atom2:
                if not center_head - 1 == atom1:
                    tail_atoms.append(atom1 + 1)

        backbones = {}
        for tail in tail_atoms:
            for head in head_atoms:
                tail_torsion = get_tail_torsion(self.path, head, tail)
                while tail_torsion < 0:
                    tail_torsion += 360
                while tail_torsion >= 360:
                    tail_torsion -= 360
                tail_torsion = (
                    60
                    if 0 <= tail_torsion < 120
                    else 180
                    if 120 <= tail_torsion < 240
                    else 300
                )
                m_name = "h_m" + str(tail) + "T" + str(tail_torsion)
                backbones[m_name] = Backbone(
                    m_name, self.path, head, tail, tail_torsion
                )
                if find_inverse_path:
                    tail_torsion = get_tail_torsion(self.path, tail, head)
                    while tail_torsion < 0:
                        tail_torsion += 360
                    while tail_torsion >= 360:
                        tail_torsion -= 360
                    tail_torsion = (
                        60
                        if 0 <= tail_torsion < 120
                        else 180
                        if 120 <= tail_torsion < 240
                        else 300
                    )
                    m_name = "t_m" + str(tail) + "T" + str(tail_torsion)
                    backbones[m_name] = Backbone(
                        m_name, self.path, tail, head, tail_torsion
                    )

        return backbones


monomer_path = os.path.join(current_location, "..", "data", "pdb", "monomer")

polymer_types = {
    "PE": PolymerType(
        "PE", os.path.join(monomer_path, "PE.pdb"), backbone_atoms=[2, 3]
    ),
    "PP": PolymerType(
        "PP", os.path.join(monomer_path, "PP.pdb"), backbone_atoms=[2, 3], side_atom=8
    ),
    "PS": PolymerType(
        "PS", os.path.join(monomer_path, "PS.pdb"), backbone_atoms=[2, 3], side_atom=8
    ),
    "POM": PolymerType(
        "POM", os.path.join(monomer_path, "POM.pdb"), backbone_atoms=[2, 3]
    ),
    "PTFE": PolymerType(
        "PTFE", os.path.join(monomer_path, "PTFE.pdb"), backbone_atoms=[2, 3]
    ),
    "PVC": PolymerType(
        "PVC", os.path.join(monomer_path, "PVC.pdb"), backbone_atoms=[2, 3], side_atom=8
    ),
    "PAN": PolymerType(
        "PAN", os.path.join(monomer_path, "PAN.pdb"), backbone_atoms=[3, 2], side_atom=5
    ),
}
