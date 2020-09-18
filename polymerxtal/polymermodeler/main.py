# ============================================================================
# main.py -- Linear Amorphous Thermoplastic CHain builder top level functions
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import numpy as np

from .config import MAX_BONDS
from .params import Params
from .system import PolymerSystem
from .timer import Timer
from .utils import CA2CC, AMU2GRAM

# File scope
params = Params()
polysys = PolymerSystem()
log_file = {'name': ''}
status_file = {'name': ''}
total_monomers = 0
update_count = 0
timer = Timer()


# ============================================================================
# logMessage()
# ----------------------------------------------------------------------------
# Result: report a message to the user
# ============================================================================
def logMessage(msg):
    if log_file['name'] and msg:
        log_file['file'].write(msg + "\n")


# ============================================================================
# cleanup()
# ----------------------------------------------------------------------------
# Result: free all top level structures
# ============================================================================
def cleanup():
    polysys.cleanupSystem(params)  # before freeParams()
    params.freeParams()
    if log_file['name']:
        log_file['file'].close()
    if status_file['name']:
        status_file['file'].close()


# ============================================================================
# showVersion()
# ----------------------------------------------------------------------------
# Result: write version and build info to the FILE f
# ============================================================================
def showVersion(f):
    f.write("OpenMP threading disabled\n")
    f.write("MAX_BONDS: %d\n" % MAX_BONDS)


# ============================================================================
# main()
# ----------------------------------------------------------------------------
# Result: read input, run simulation, write output, kitchen sink...;
# return EXIT_SUCCESS or exit() with EXIT_FAILURE via choke()
# ============================================================================
def main(args):
    done = 0
    infile = ''
    usage = [
        "Usage: %s [flags] infile\n", "Flags:\n", "   -v   Show version and build info, then exit\n",
        "   -h   Show help message, then exit\n"
    ]
    mini = np.zeros(3)
    maxi = np.zeros(3)
    m = []
    s = []
    f = []

    timer.startTimer()

    for p in args:
        if p:
            if '-' == p[0]:
                if len(p) > 2 and p[2]:
                    # multi-char flag
                    raise SyntaxError("\nUnknown flag: \"%s\"\n" % p)
                if len(p) == 2:
                    if p[1] == 'h':
                        print(''.join(usage))
                        done += 1
                    elif p[1] == 'v':
                        print(''.join(usage))
                        done += 1
                    else:
                        raise SyntaxError("\nUnknown flag: \"%s\"\n" % p[1])
                        break
            else:
                if infile:
                    raise SyntaxError("\nInput file specified twice: %s, %s\n" % (infile, p))
                infile = p

    if done:
        cleanup()
        return 0
    if not infile:
        raise SyntaxError("\nMissing input file\n")

    # Read input
    params.readParams(infile)

    # Maximum monomer mass
    monomer_mass = 0.0
    monomer_count = 0
    m = params.known_monomers
    while m and hasattr(m, 'next'):
        if m.central_mass > monomer_mass:
            monomer_mass = m.central_mass
        monomer_count += 1
        m = m.next
    if 0 == monomer_count:
        raise ValueError("No monomers specified")

    # System size and chain distribution
    # XXX params.excluded_volume does not consider volume overlaps
    if params.explicit_volume:
        params.total_volume = 1.0
        if params.system_max[0] < params.system_min[0]:
            raise ValueError("System minimum x (%f) > system maximum x (%f)" %
                             (params.system_min[0], params.system_max[0]))
        params.system_size[0] = params.system_max[0] - params.system_min[0]
        params.total_volume *= params.system_size[0]
        if params.system_max[1] < params.system_min[1]:
            raise ValueError("System minimum y (%f) > system maximum y (%f)" %
                             (params.system_min[1], params.system_max[1]))
        params.system_size[1] = params.system_max[1] - params.system_min[1]
        params.total_volume *= params.system_size[1]
        if params.system_max[2] < params.system_min[2]:
            raise ValueError("System minimum z (%f) > system maximum z (%f)" %
                             (params.system_min[2], params.system_max[2]))
        params.system_size[2] = params.system_max[2] - params.system_min[2]
        params.total_volume *= params.system_size[2]  # total box volume
        if params.inverted_volume:
            params.total_volume = params.excluded_volume
        else:
            params.total_volume -= params.excluded_volume
        # Now params.total_volume holds the true polymer volume
        if 0 == params.num_monomers:
            params.num_monomers = int(
                (params.density * CA2CC(params.total_volume)) / float(params.num_chains * AMU2GRAM(monomer_mass)))
    else:
        params.total_volume = CC2CA(
            float(params.num_chains * params.num_monomers) * AMU2GRAM(monomer_mass) / params.density)  # polymer volume
        params.system_min[0] = 0.0
        params.system_min[1] = 0.0
        params.system_min[2] = 0.0
        params.system_max[0] = np.cbrt(params.total_volume + params.excluded_volume)
        params.system_max[1] = params.system_max[0]
        params.system_max[2] = params.system_max[0]
        params.system_size[0] = params.system_max[0]
        params.system_size[1] = params.system_max[1]
        params.system_size[2] = params.system_max[2]
    params.half_system_size[0] = 0.5 * params.system_size[0]
    params.half_system_size[1] = 0.5 * params.system_size[1]
    params.half_system_size[2] = 0.5 * params.system_size[2]

    # Spatial domain parameters
    params.total_domains = params.num_domains_x * params.num_domains_y * params.num_domains_z
    params.domain_size[0] = params.system_size[0] / float(params.num_domains_x)
    params.domain_size[1] = params.system_size[1] / float(params.num_domains_y)
    params.domain_size[2] = params.system_size[2] / float(params.num_domains_z)

    # Log file
    if params.log_file:
        log_file['file'] = open(params.log_file, "w")

    # Status file
    if params.status_file:
        status_file['file'] = open(params.status_file, "w")

    # Report parameters
    log_file['file'].write("\n========================================")
    log_file['file'].write("========================================\n")
    showVersion(log_file['file'])
    log_file['file'].write("========================================")
    log_file['file'].write("========================================\n\n")
    params.reportParams(log_file, log_file['file'])


#int
#main(int argc, char **argv)
#{
#   int done = 0;
#   int monomer_count;
#   int max_monomers;
#   int total_atoms;
#   int i, j, k, n;
#   char *p;
#   char *infile = NULL;
#   const char usage[] =
#      "Usage: %s [flags] infile\n"
#      "Flags:\n"
#      "   -v   Show version and build info, then exit\n"
#      "   -h   Show help message, then exit\n";
#   Real monomer_mass;
#   Real monomer_mean;
#   Real monomer_stddev;
#   Real total_mass;
#   Real g;
#   Vector min, max;
#   Monomer *m;
#   Stereo *s;
#   FILE *f;
