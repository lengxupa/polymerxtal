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

from .params import Params
from .system import PolymerSystem
from .timer import Timer
from .stdio import FILE

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
    #polysys.buildsystem(params)

    #Unfinished


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
