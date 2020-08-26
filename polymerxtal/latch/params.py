import random
from datetime import datetime

from .scan import Scanner

class Params:
    def __init__(self):
        v0 = [0.0, 0.0, 0.0]
        self.system_min = v0
        self.system_max = v0
        self.system_size = v0
        self.domain_size = v0
        self.bond_scale = 0.0
        self.temperature = 0.1 # not zero
        self.density = 0.0
        self.grid_size = 0.0
        self.energy_cutoff = 0.0
        self.self_avoid_cutoff = 0.0
        self.total_volume = 0.0
        self.excluded_volume = 0.0
        self.torsion_step = 0.0
        self.backbone_bond_length = 0.0
        self.num_monomers_stddev = 0.0
        self.chain_length_histo_bin = 0.0
        self.rng_seed = random.seed(datetime.now())
        self.num_chains = 0
        self.num_stereo = 0
        self.num_monomers = 0
        self.num_domains_x = 1
        self.num_domains_y = 1
        self.num_domains_z = 1
        self.total_domains = 0
        self.num_configs = 30
        self.bond_cutoff = 0
        self.num_delta_steps = 0
        self.max_monomer_atoms = 0
        self.recalculate_positions = 0
        self.recalculate_neighbors = 0
        self.sample_monte_carlo = 0
        self.explicit_volume = 0
        self.inverted_volume = 0
        self.isolate_chains = 0
        self.write_wrapped_pdb = 0
        self.write_unwrapped_pdb = 0
        self.write_wrapped_xyz = 0
        self.write_unwrapped_xyz = 0
        self.write_chain_length_histo = 0
        self.write_chain_length = 0
        self.write_torsion_histo = 0
        self.write_intermediate = 0

    def readParams(self, path):
        s = createScanner(path)
        done = 0
        while not done:
            pass


