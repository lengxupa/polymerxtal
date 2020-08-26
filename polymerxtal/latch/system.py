class PolymerSystem:
    def __init__(self, chains=0, domains=0, rngs=0, pending_atoms=0):
        self.chains = chains
        self.domains = domains
        self.rngs = rngs
        self.pending_atoms = pending_atoms

    def __str__(self):
        return f'chains: {self.chains}\ndomains: {self.domains}\nrngs: {self.rngs}\npending_atoms: {self.pending_atoms}'