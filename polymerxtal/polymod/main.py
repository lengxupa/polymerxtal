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
import random

from .build import *  # noqa: F403
from .chain import Chain, createChain
from .config import MAX_BONDS
from .domain import createDomain  # Bin, , Domain

# from .monomer import Monomer, Bond
from .params import Params
from .stdio import *  # noqa: F403

# from .stereo import Stereo
from .system import PolymerSystem
from .timer import Timer
from .types import *  # noqa: F403
from .utils import *  # noqa: F403

# from .zmatrix import ZMatrix

# File scope
polysys = PolymerSystem()
log_file = FILE()  # noqa: F405
status_file = FILE()  # noqa: F405
timer = Timer()
params = Params()


# ============================================================================
# logMessage()
# ----------------------------------------------------------------------------
# Result: report a message to the user
# ============================================================================
def logMessage(msg):
    if log_file.path and msg:
        log_file.printf(msg)
        log_file.printf("\n")


# ============================================================================
# cleanup()
# ----------------------------------------------------------------------------
# Result: free all top level structures
# ============================================================================
def cleanup(params):
    polysys.cleanupSystem(params)  # before freeParams()
    del params
    if log_file.path and stdout != log_file:  # noqa: F405
        log_file.fclose()
    if status_file.path and stdout != status_file:  # noqa: F405
        status_file.fclose()


# ============================================================================
# showVersion()
# ----------------------------------------------------------------------------
# Result: write version and build info to the FILE f
# ============================================================================
def showVersion(f):
    # fprintf(f, "%s version %s\n", prog_name, prog_version);
    # fprintf(f, "Compiled with \"%s\"\n", prog_cc_line);
    # fprintf(f, "Linked with \"%s\"\n", prog_ld_line);
    # ifdef HAVE_OPENMP
    # fprintf(f, "OpenMP threading enabled\n");
    # else
    f.printf("OpenMP threading disabled\n")
    # endif
    # fprintf(f, "Real type: %s\n", REAL_STR);
    f.printf("MAX_BONDS: %d\n" % MAX_BONDS)


# ============================================================================
# main()
# ----------------------------------------------------------------------------
# Result: read input, run simulation, write output, kitchen sink...;
# return EXIT_SUCCESS or exit() with EXIT_FAILURE via choke()
# ============================================================================
def main(args):
    done = 0
    infile = ""
    # usage = (
    #    "Usage: %s [flags] infile\n"
    #    + "Flags:\n"
    #    + "   -v   Show version and build info, then exit\n"
    #    + "   -h   Show help message, then exit\n"
    # )
    mini = np.zeros(3)
    maxi = np.zeros(3)
    total_monomers = 0

    def SHOW_USAGE(f):
        pass

    timer.startTimer()
    polysys = PolymerSystem()
    params = Params()

    for p in args:
        if p:
            if "-" == p[0]:
                if len(p) > 2 and p[2]:
                    # multi-char flag
                    SHOW_USAGE(stderr)  # noqa: F405
                    raise SyntaxError('\nUnknown flag: "%s"\n' % p)
                if len(p) == 2:
                    if p[1] == "h":
                        SHOW_USAGE(stdout)  # noqa: F405
                        done += 1
                    elif p[1] == "v":
                        showVersion(stdout)  # noqa: F405
                        done += 1
                    else:
                        SHOW_USAGE(stderr)  # noqa: F405
                        raise SyntaxError('\nUnknown flag: "%s"\n' % p[1])
                        break
            else:
                if infile:
                    SHOW_USAGE(stderr)  # noqa: F405
                    raise SyntaxError(
                        "\nInput file specified twice: %s, %s\n" % (infile, p)
                    )
                infile = p

    if done:
        cleanup()
        return 0
    if not infile:
        SHOW_USAGE(stderr)  # noqa: F405
        raise SyntaxError("\nMissing input file\n")

    # Read input
    params.readParams(infile)

    # Maximum monomer mass
    monomer_mass = 0.0
    monomer_count = 0
    m = params.known_monomers
    while m and hasattr(m, "next"):
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
            raise ValueError(
                "System minimum x (%f) > system maximum x (%f)"
                % (params.system_min[0], params.system_max[0])
            )
        params.system_size[0] = params.system_max[0] - params.system_min[0]
        params.total_volume *= params.system_size[0]
        if params.system_max[1] < params.system_min[1]:
            raise ValueError(
                "System minimum y (%f) > system maximum y (%f)"
                % (params.system_min[1], params.system_max[1])
            )
        params.system_size[1] = params.system_max[1] - params.system_min[1]
        params.total_volume *= params.system_size[1]
        if params.system_max[2] < params.system_min[2]:
            raise ValueError(
                "System minimum z (%f) > system maximum z (%f)"
                % (params.system_min[2], params.system_max[2])
            )
        params.system_size[2] = params.system_max[2] - params.system_min[2]
        params.total_volume *= params.system_size[2]  # total box volume
        if params.inverted_volume:
            params.total_volume = params.excluded_volume
        else:
            params.total_volume -= params.excluded_volume
        # Now params.total_volume holds the true polymer volume
        if 0 == params.num_monomers:
            params.num_monomers = int(
                (params.density * CA2CC(params.total_volume))  # noqa: F405
                / float(params.num_chains * AMU2GRAM(monomer_mass))  # noqa: F405
            )
    else:
        params.total_volume = CC2CA(  # noqa: F405
            float(params.num_chains * params.num_monomers)
            * AMU2GRAM(monomer_mass)  # noqa: F405
            / params.density
        )  # polymer volume
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
    params.total_domains = (
        params.num_domains_x * params.num_domains_y * params.num_domains_z
    )
    params.domain_size[0] = params.system_size[0] / float(params.num_domains_x)
    params.domain_size[1] = params.system_size[1] / float(params.num_domains_y)
    params.domain_size[2] = params.system_size[2] / float(params.num_domains_z)

    # Log file
    log_file = stdout  # noqa: F405
    if params.log_file:
        log_file = openFile(params.log_file, "w")  # noqa: F405

    # Status file
    status_file = stdout  # noqa: F405
    if params.status_file:
        status_file = openFile(params.status_file, "w")  # noqa: F405 F841  # noqa:

    # Report parameters
    log_file.printf("\n========================================")
    log_file.printf("========================================\n")
    showVersion(log_file)
    log_file.printf("========================================")
    log_file.printf("========================================\n\n")
    params.reportParams(log_file)

    # Random number generators: one per Domain
    polysys.rngs = {}
    for i in range(params.total_domains):
        polysys.rngs[i] = params.rng_seed + i

    # Domains
    polysys.domains = {}
    n = 0
    for i in range(params.num_domains_z):
        for j in range(params.num_domains_y):
            for k in range(params.num_domains_x):
                mini[0] = params.system_min[0] + float(k) * params.domain_size[0]
                maxi[0] = mini[0] + params.domain_size[0]
                mini[1] = params.system_min[1] + float(j) * params.domain_size[1]
                maxi[1] = mini[1] + params.domain_size[1]
                mini[2] = params.system_min[2] + float(i) * params.domain_size[2]
                maxi[2] = mini[2] + params.domain_size[2]
                polysys.domains[n] = createDomain(n, mini, maxi, params)
                n += 1

    # Chains
    max_monomers = 0
    total_atoms = 0
    total_mass = 0.0
    for i in range(params.num_chains):
        polysys.chains[i] = Chain()
        # Choose a length in Monomers
        if params.num_monomers_stddev > 0.0:
            n = int(
                random.gauss(float(params.num_monomers), params.num_monomers_stddev)
            )
        else:
            n = params.num_monomers
        # Choose a Stereo option for this Chain
        j = selectWeight(  # noqa: F405
            params.chain_stereo_weights, params.num_stereo, polysys.rngs[0]
        )
        s = params.chain_stereo[j]
        s.selection_count += 1
        k = n * (params.max_monomer_atoms - 2) + 2
        total_mass += n * monomer_mass
        if s.term and hasattr(s.term, "next"):
            n += 1
            k += s.term.num_atoms - 2
        if n > max_monomers:
            max_monomers = n
        # global total_monomers
        total_monomers += n
        polysys.chains[i] = createChain(s, n, k, not params.recalculate_positions)
        total_atoms += k

        # Choose a Domain
        polysys.chains[i].init_domain = int(random.random() * params.total_domains)
        polysys.chains[i].domain = polysys.chains[i].init_domain
    log_file.printf("\nMaximum atoms: %d\n" % total_atoms)
    log_file.printf(
        "Rough estimate of density: %f g/cm^3\n\n"
        % (AMU2GRAM(total_mass) / CA2CC(params.total_volume))  # noqa: F405
    )

    # Build Chains
    # ifdef HAVE_OPENMP
    # pragma omp parallel default(shared) num_threads(params.total_domains)
    #   {  /* Each thread builds in a specific Domain */
    #      buildChains(&sys, &params, omp_get_thread_num());
    #   }
    # else
    buildChains(polysys, params, 0)  # noqa: F405
    # endif

    # Final log outputs
    total_mass = 0.0
    total_atoms = 0
    monomer_mean = 0.0
    for i in range(params.num_chains):
        if not polysys.chains[i].dead:
            total_mass += polysys.chains[i].mass
            total_atoms += polysys.chains[i].curr_atom
            monomer_mean += float(polysys.chains[i].curr_monomer)
    monomer_mean /= float(params.num_chains)
    monomer_stddev = 0.0
    for i in range(params.num_chains):
        if not polysys.chains[i].dead:
            g = float(polysys.chains[i].curr_monomer) - monomer_mean
            monomer_stddev += g * g
    if params.num_chains > 1:
        monomer_stddev /= float(params.num_chains - 1)  # bias corrected
    log_file.printf("Monomer selection:\n")
    m = params.known_monomers
    while m and hasattr(m, "next"):
        log_file.printf(
            "   %s: %.2f %%\n"
            % (m.name, 100.0 * float(m.selection_count) / float(total_monomers))
        )
        m = m.next
    log_file.printf("\nStereochemistry selection:\n")
    s = params.known_stereo
    while s and hasattr(s, "next"):
        log_file.printf(
            "   %s: %.2f %%\n"
            % (s.name, 100.0 * float(s.selection_count) / float(params.num_chains))
        )
        s = s.next
    log_file.printf("Mean polymerization: %f\n" % monomer_mean)
    log_file.printf("   Standard deviation: %f\n" % np.sqrt(monomer_stddev))
    log_file.printf("Total atoms: %d\n" % total_atoms)
    log_file.printf(
        "Final density: %f g/cm^3\n"
        % (AMU2GRAM(total_mass) / CA2CC(params.total_volume))  # noqa: F405
    )

    # PDB output
    if params.write_wrapped_pdb:
        f = openFile("chains_wrapped.pdb", "w")  # noqa: F405
        n = 1
        for i in range(params.num_chains):
            if not polysys.chains[i].dead:
                polysys.chains[i].writeChainPDB(n, 1, f, params)
        f.fclose()
    if params.write_unwrapped_pdb:
        f = openFile(".tmp/chains_unwrapped.pdb", "w")  # noqa: F405
        n = 1
        for i in range(params.num_chains):
            if not polysys.chains[i].dead:
                polysys.chains[i].writeChainPDB(n, 0, f, params)
        f.fclose()

    # XYZ output
    if params.write_wrapped_xyz:
        f = openFile("chains_wrapped.xyz", "w")  # noqa: F405
        f.printf("%d\nChains\n" % total_atoms)
        for i in range(params.num_chains):
            if not polysys.chains[i].dead:
                polysys.chains[i].writeChainXYZ(1, f, params)
        f.fclose()
    if params.write_unwrapped_xyz:
        f = openFile("chains_unwrapped.xyz", "w")  # noqa: F405
        f.printf("%d\nChains\n" % total_atoms)
        for i in range(params.num_chains):
            if not polysys.chains[i].dead:
                polysys.chains[i].writeChainXYZ(0, f, params)
        f.fclose()

    # Chain length output
    if params.write_chain_length_histo:
        max_len = 0.0

        for i in range(params.num_chains):
            nm = polysys.chains[i].num_monomers - 1
            if (not polysys.chains[i].dead) and polysys.chains[i].length[nm] > max_len:
                max_len = polysys.chains[i].length[nm]
        num_bins = int(np.ceil(max_len / params.chain_length_histo_bin))
        bins = {}
        for i in range(num_bins):
            bins[i] = 0
        for i in range(params.num_chains):
            nm = polysys.chains[i].num_monomers - 1
            if not polysys.chains[i].dead:
                bins[
                    int(polysys.chains[i].length[nm] / params.chain_length_histo_bin)
                ] += 1
        f = openFile("chain_length_histo.dat", "w")  # noqa: F405
        f.printf("# length (A)    count\n")
        for i in range(num_bins):
            f.printf("%f  %d\n" % ((i + 1) * params.chain_length_histo_bin, bins[i]))
        f.fclose()
        del bins
    if params.write_chain_length:
        length = {}
        weight = {}

        for i in range(max_monomers):
            length[i] = 0.0
            weight[i] = 0.0
        for i in range(params.num_chains):
            if polysys.chains[i].dead:
                continue
            for j in range(polysys.chains[i].num_monomers):
                length[j] += polysys.chains[i].length[j] * polysys.chains[i].weight[j]
                weight[j] += polysys.chains[i].weight[j]
        f = openFile("mean_chain_length.dat", "w")  # noqa: F405
        f.printf("# monomers  mean length (A)\n")
        for i in range(max_monomers):
            f.printf(
                "%d  %f\n" % (i + 1, length[i] / weight[i] if weight[i] > 0.0 else 0.0)
            )
        f.fclose()
        del length
        del weight

    # Torsion output
    if params.write_torsion_histo:
        count = {}

        for i in range(360):
            count[i] = 0
        for i in range(params.num_chains):
            if polysys.chains[i].dead:
                continue
            for j in range(360):
                count[j] += polysys.chains[i].torsion_count[j]
        f = openFile("torsion_histo.dat", "w")  # noqa: F405
        f.printf("# angle (deg)   count\n")
        for i in range(360):
            f.printf("%d  %d\n" % (i, count[i]))
        f.fclose()
        del count

    # Intermediate files for LAMMPS inputs */
    # N.B. This block of code is largely duplicated in pdb2im.c; it should be
    # abstracted to accept arrays of positions, bonds ...
    if params.write_intermediate:
        pos = np.zeros(3)

        f = openFile(".tmp/atom_type.dat", "w")  # noqa: F405
        writeAtomTypesDreiding(f)  # noqa: F405
        f.fclose()
        f = openFile(".tmp/bond_type.dat", "w")  # noqa: F405
        writeBondTypesDreiding(f)  # noqa: F405
        f.fclose()
        na = 0
        for i in range(params.num_chains):
            na += polysys.chains[i].curr_atom
        f = openFile(".tmp/atoms.dat", "w")  # noqa: F405
        f.printf("%d atoms\n\nATOMS\n\n" % na)
        na = 1
        for i in range(params.num_chains):
            if polysys.chains[i].dead:
                continue
            for j in range(polysys.chains[i].curr_atom):
                pos = polysys.chains[i].zm.getPosition(j, pos)
                foldPosition(  # noqa: F405
                    pos, params.system_min, params.system_max, params.system_size
                )
                f.printf(
                    "%d  %d  %d  %f  %.6f  %.6f  %.6f\n"
                    % (
                        na,
                        i + 1,
                        getAtomTypeIndex(  # noqa: F405
                            polysys.chains[i].zm.entries[j].type
                        ),  # noqa: F405
                        0.0,  # XXX charge
                        pos[0],
                        pos[1],
                        pos[2],
                    )
                )
                na += 1
        f.fclose()
        f = openFile(".tmp/bonds.dat", "w")  # noqa: F405
        f.printf("BONDS\n\n")
        na = 0
        nb = 1
        offset = 0
        for i in range(params.num_chains):
            if polysys.chains[i].dead:
                continue
            for j in range(polysys.chains[i].curr_atom):
                if polysys.chains[i].zm.entries[j].bond_index > -1:
                    type1 = polysys.chains[i].zm.entries[j].type
                    k = polysys.chains[i].zm.entries[j].bond_index
                    type2 = polysys.chains[i].zm.entries[k].type
                    f.printf(
                        "%d  %d  %d  %d\n"
                        % (
                            nb,
                            getBondTypeIndex(type1, type2),  # noqa: F405
                            na + 1,
                            offset + polysys.chains[i].zm.entries[j].bond_index + 1,
                        )
                    )
                    nb += 1
                na += 1
            b = polysys.chains[i].extra_bonds
            while b and hasattr(b, "next"):
                type1 = polysys.chains[i].zm.entries[b.index1].type
                type2 = polysys.chains[i].zm.entries[b.index2].type
                f.printf(
                    "%d  %d  %d  %d\n"
                    % (
                        nb,
                        getBondTypeIndex(type1, type2),  # noqa: F405
                        offset + b.index1 + 1,
                        offset + b.index2 + 1,
                    )
                )
                nb += 1
                b = b.next
            offset = na
        f.fclose()

    # Finish
    log_file.printf("Total build and analysis time: %g s\n" % timer.getElapsedTime())
    cleanup(params)
    del params
    return 0
