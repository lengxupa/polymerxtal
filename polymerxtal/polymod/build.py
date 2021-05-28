# ============================================================================
# build.py -- buildChains(), the main simulation function
# ----------------------------------------------------------------------------
# Author: Benjamin P. Haley, Tongtong Shen, Purdue University
# Copyright (c) 2012 Purdue University
# ----------------------------------------------------------------------------
# See the LICENSE file for information on usage and redistribution of this
# file and for a DISCLAIMER OF ALL WARRANTIES.
# ============================================================================

import numpy as np

# from .chain import Chain
from .config import REAL_MAX
from .domain import OccAtom  # Domain,
from .energy import Energy

# from .exclude import ExclCylinder, ExclSlab, ExclSphere
from .monomer import Monomer, Torsion
from .params import Params

# from .stereo import Stereo
# from .system import PolymerSystem
from .utils import *  # noqa: F403
from .vector import getNearestSqDist

# from .zmatrix import ZMatrix

params = Params()
update_count = 0


# ============================================================================
# updateStatus()
# ----------------------------------------------------------------------------
# Result: update the addition of a monomer to a Chain; write status update
# after all Chains have updated
# ============================================================================
def updateStatus():
    curr_monomers = 0
    # pragma omp critical (status)

    # printf("thread %d: status update %d of %d\n", omp_get_thread_num(),
    #        update_count+1, params.total_domains);
    # fflush(stdout);

    global update_count
    update_count += 1
    if update_count == params.total_domains:
        for i in range(params.num_chains):
            curr_monomers += polysys.chains[i].curr_monomer  # noqa: F405
            pct = int(
                100.0 * float(curr_monomers) / float(total_monomers)  # noqa: F405
            )  # noqa: F405
            # ifdef RAPPTURE_STATUS
            # fprintf(status_file, "=RAPPTURE-PROGRESS=>%d\n", pct);
            # else
            BARLEN = 30
            fputc("|", status_file)  # noqa: F405
            maxi = int(float(pct * BARLEN) * 0.01)
            for i in range(BARLEN):
                fputc("=" if i < maxi else " ", status_file)  # noqa: F405
            status_file.printf(  # noqa: F405
                "| %3d%% (%f s)   " % (pct, timer.getElapsedTime())  # noqa: F405
            )  # noqa: F405
            if stdout == status_file:  # noqa: F405
                fputc("\r" if pct < 100 else "\n", status_file)  # noqa: F405
            else:
                status_file.rewind()  # noqa: F405
            # endif
            update_count = 0


# ============================================================================
# interpolateTorsion()
# ----------------------------------------------------------------------------
# Result: set the energy, *E, and the probability, *p, associated with the
# specified torsion angle, using a linear interpolation in
# m->torsion_angles[ti]
# ============================================================================
def interpolateTorsion(angle, m, ti, E, p):
    i = 0
    found = 0

    while (not found) and i < m.num_torsions[ti] - 1:
        if m.torsion_angles[ti][i] <= angle and m.torsion_angles[ti][i + 1] > angle:
            found += 1
        else:
            i += 1
    if found:
        dang = m.torsion_angles[ti][i + 1] - m.torsion_angles[ti][i]
        dE = m.torsion_energies[ti][i + 1] - m.torsion_energies[ti][i]
        dp = m.torsion_probs[ti][i + 1] - m.torsion_probs[ti][i]
    else:
        dang = m.torsion_angles[ti][0] - m.torsion_angles[ti][i] + 360.0
        dE = m.torsion_energies[ti][0] - m.torsion_energies[ti][i]
        dp = m.torsion_probs[ti][0] - m.torsion_probs[ti][i]
    if dang > 0.0:
        E = m.torsion_energies[ti][i] + dE * (angle - m.torsion_angles[ti][i]) / dang
        p = m.torsion_probs[ti][i] + dp * (angle - m.torsion_angles[ti][i]) / dang
    else:
        E = m.torsion_energies[ti][i]  # noqa: F841
        p = m.torsion_probs[ti][i]  # noqa: F841


# ============================================================================
# setTorsions()
# ----------------------------------------------------------------------------
# Result: set the last num_torsions backbone torsion angles of the tail
# monomer, of type m, in the Chain c; the tail monomer backbone indices
# start at zm_offset in c->zm, and the corresponding angles, energies,
# and probabilities are found in
# m->torsion_{probs,angles,energies}[torsion_offset+i]; updates the sum
# of bonded interactions, Eb, and the weight of bonded interactions, wb.
# ============================================================================
def setTorsions(
    c, m, num_torsions, zm_offset, torsion_offset, rng, Eb, wb, torsion_step
):

    for i in range(num_torsions):
        if torsion_step > 0.0:  # modify existing choice
            if m.torsions[i + torsion_offset] == Torsion.TORSION_FIXED:
                # no change
                if Eb:
                    Eb += 0.0
                if wb:
                    wb *= 1.0
            elif m.torsions[i + torsion_offset] == Torsion.TORSION_FREE:
                sign = -1.0 if random.random() < 0.5 else 1.0  # noqa: F405
                c.zm.entries[zm_offset + i].torsion_angle += sign * torsion_step
                if c.zm.entries[zm_offset + i].torsion_angle >= 360.0:
                    c.zm.entries[zm_offset + i].torsion_angle -= 360.0
                if c.zm.entries[zm_offset + i].torsion_angle < 0.0:
                    c.zm.entries[zm_offset + i].torsion_angle += 360.0
                if Eb:
                    Eb += 0.0
                if wb:
                    wb *= 1.0
            elif m.torsions[i + torsion_offset] == Torsion.TORSION_ENERGY:
                sign = -1.0 if random.random() < 0.5 else 1.0  # noqa: F405
                c.zm.entries[zm_offset + i].torsion_angle += sign * torsion_step
                if c.zm.entries[zm_offset + i].torsion_angle >= 360.0:
                    c.zm.entries[zm_offset + i].torsion_angle -= 360.0
                if c.zm.entries[zm_offset + i].torsion_angle < 0.0:
                    c.zm.entries[zm_offset + i].torsion_angle += 360.0
                interpolateTorsion(
                    c.zm.entries[zm_offset + i].torsion_angle,
                    m,
                    i + torsion_offset,
                    E,  # noqa: F405
                    p,  # noqa: F405
                )
                if Eb:
                    Eb += E  # noqa: F405
                if wb:
                    wb *= m.torsion_prob_min[i + torsion_offset] / p  # noqa: F405

        else:  # initial choice
            if m.torsions[i + torsion_offset] == Torsion.TORSION_FIXED:
                if m.torsion_angles[i + torsion_offset][0] > -REAL_MAX:
                    c.zm.entries[zm_offset + i].torsion_angle = m.torsion_angles[
                        i + torsion_offset
                    ][0]
            elif m.torsions[i + torsion_offset] == Torsion.TORSION_FREE:
                c.zm.entries[zm_offset + i].torsion_angle = (
                    random.random() * 360.0  # noqa: F405
                )  # noqa: F405
            elif m.torsions[i + torsion_offset] == Torsion.TORSION_ENERGY:
                n = selectWeight(  # noqa: F405
                    m.torsion_probs[i + torsion_offset],
                    m.num_torsions[i + torsion_offset],
                    rng,
                )
                c.zm.entries[zm_offset + i].torsion_angle = m.torsion_angles[
                    i + torsion_offset
                ][n]
                if Eb:
                    Eb += m.torsion_energies[i + torsion_offset][n]
                if wb:
                    wb *= (
                        m.torsion_prob_min[i + torsion_offset]
                        / m.torsion_probs[i + torsion_offset][n]
                    )


# ============================================================================
# rejectConfig()
# ----------------------------------------------------------------------------
# Result: return 1 if any atom in the tail monomer of c intrudes on an
# excluded region, else return 0
# ============================================================================
def rejectConfig(c, indices, num_indices, store_positions, start_unset, p):
    pos = np.zeros(3)

    if (
        (not p.excluded_cylinders)
        and (not p.excluded_slabs)
        and (not p.excluded_spheres)
    ):
        return 0
    for i in range(num_indices):
        pos = c.zm.getPosition(indices[i], pos)
        # Cache positions to avoid size-of-monomer recursion for getPosition();
        # note that this is NOT the final position...
        if store_positions:
            c.zm.setPosition(indices[i], pos)
        foldPosition(pos, p.system_min, p.system_max, p.system_size)  # noqa: F405
        ec = p.excluded_cylinders
        while ec and hasattr(ec, "next"):
            if ec.invert != ec.insideExclCylinder(pos):
                return 1
            ec = ec.next
        es = p.excluded_slabs
        while es and hasattr(es, "next"):
            if es.invert != es.insideExclSlab(pos):
                return 1
            es = es.next
        esph = p.excluded_spheres
        while esph and hasattr(esph, "next"):
            if esph.invert != esph.insideExclSphere(pos):
                return 1
            esph = esph.next
    if store_positions:
        # Unset temporarily cached positions
        for i in range(start_unset, num_indices):
            c.zm.clearPosition(indices[i])
    return 0


# ============================================================================
# fillNeighbors()
# ----------------------------------------------------------------------------
# Result: fill oa_nbrs[] with pointers to the top of OccAtom lists in
# unique neighbor bins; oa_nbrs[i] will be NULL if the ith neighbor is a
# duplicate of a bin already in oa_nbrs[]
# ============================================================================
def fillNeighbors(s, d, nbrs, domain_edge, oa_nbrs):
    nd = [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]

    # nd[] represents the neighbor domain indices, the indices of
    # d->nbr_domains[], representing the Domains that hold the neighboring
    # Bins; d->nbr_domains[0] is d (and nbrs[0] is index of the Bin that
    # holds the current atom under investigation)

    if domain_edge[MIN_X_EDGE] and (not domain_edge[MIN_Y_EDGE]):  # noqa: F405
        nd[1] = 4
    elif (not domain_edge[MIN_X_EDGE]) and domain_edge[MIN_Y_EDGE]:  # noqa: F405
        nd[1] = 2
    elif domain_edge[MIN_X_EDGE] and domain_edge[MIN_Y_EDGE]:  # noqa: F405
        nd[1] = 1

    if domain_edge[MIN_Y_EDGE]:  # noqa: F405
        nd[2] = 2

    if domain_edge[MAX_X_EDGE] and (not domain_edge[MIN_Y_EDGE]):  # noqa: F405
        nd[3] = 5
    elif (not domain_edge[MAX_X_EDGE]) and domain_edge[MIN_Y_EDGE]:  # noqa: F405
        nd[3] = 2
    elif domain_edge[MAX_X_EDGE] and domain_edge[MIN_Y_EDGE]:  # noqa: F405
        nd[3] = 3

    if domain_edge[MIN_X_EDGE]:  # noqa: F405
        nd[4] = 4

    if domain_edge[MAX_X_EDGE]:  # noqa: F405
        nd[5] = 5

    if domain_edge[MAX_Y_EDGE] and (not domain_edge[MIN_X_EDGE]):  # noqa: F405
        nd[6] = 7
    elif (not domain_edge[MAX_Y_EDGE]) and domain_edge[MIN_X_EDGE]:  # noqa: F405
        nd[6] = 4
    elif domain_edge[MAX_Y_EDGE] and domain_edge[MIN_X_EDGE]:  # noqa: F405
        nd[6] = 6

    if domain_edge[MAX_Y_EDGE]:  # noqa: F405
        nd[7] = 7

    if domain_edge[MAX_Y_EDGE] and (not domain_edge[MAX_X_EDGE]):  # noqa: F405
        nd[8] = 7
    elif (not domain_edge[MAX_Y_EDGE]) and domain_edge[MAX_X_EDGE]:  # noqa: F405
        nd[8] = 5
    elif domain_edge[MAX_Y_EDGE] and domain_edge[MAX_X_EDGE]:  # noqa: F405
        nd[8] = 8

    if domain_edge[MIN_Z_EDGE]:  # noqa: F405
        for i in range(9, 18):
            nd[i] = nd[i - 9] + 9

    if domain_edge[MAX_Z_EDGE]:  # noqa: F405
        for i in range(18, 27):
            nd[i] = nd[i - 18] + 18

    for i in range(27):
        oa_nbrs[i] = s.domains[d.nbr_domains[nd[i]]].bins[nbrs[i]].oa_list
        for j in range(i):
            if nbrs[i] == nbrs[j] and nd[i] == nd[j]:
                # Bin list already in oa_nbrs[]
                oa_nbrs[i] = OccAtom()
                break  # out of j loop


# ============================================================================
# sumInteractions()
# ----------------------------------------------------------------------------
# Result: return the sum of non-bonded interactions of all the atoms in the
# tail monomer of s->chains[chain_index] with all surrounding atoms
# ============================================================================
def sumInteractions(
    s,
    d,
    chain_index,
    indices,
    num_indices,
    start_unset,
    store_positions,
    p,
    energy_func,
):
    E = 0.0
    ecut2 = p.energy_cutoff * p.energy_cutoff
    c = s.chains[chain_index]
    prev_bin = -1
    prev_domain = -1
    ipos = np.zeros(3)
    jpos = np.zeros(3)
    oa_nbrs = {}
    for i in range(27):
        oa_nbrs[i] = OccAtom()
        oa_nbrs[i].create()
    # Timer t;

    # startTimer(&t);
    # Loop over each tail atom
    for i in range(num_indices):
        itype = c.zm.entries[indices[i]].type.element_index

        # printf("Domain %d: getting position %d of chain %d in sumInteractions\n",
        #        d->index, indices[i]+1, chain_index+1);
        # fflush(stdout);

        ipos = c.zm.getPosition(indices[i], ipos)
        # Cache positions to avoid size-of-monomer recursion for getPosition();
        # note that this is NOT the final position...
        if store_positions:
            c.zm.setPosition(indices[i], ipos)
        foldPosition(ipos, p.system_min, p.system_max, p.system_size)  # noqa: F405
        if (
            ipos[0] > d.max[0]
            or ipos[0] < d.min[0]
            or ipos[1] > d.max[1]
            or ipos[1] < d.min[1]
            or ipos[2] > d.max[2]
            or ipos[2] < d.min[2]
        ):
            # ipos is in another Domain ...
            domain = hashBin(  # noqa: F405
                ipos, p.system_min, p.domain_size, p.num_domains_x, p.num_domains_y
            )
            if domain != prev_domain:
                nd = s.domains[domain]
                bin = hashBin(  # noqa: F405
                    ipos, nd.min, nd.bin_size, nd.num_bins_x, nd.num_bins_y
                )  # noqa: F405
                if p.recalculate_neighbors:
                    getNeighborIndices(  # noqa: F405
                        bin,
                        d.nbr_bins,
                        d.domain_edge,
                        nd.num_bins_x,
                        nd.num_bins_y,
                        nd.num_bins_z,
                    )
                    fillNeighbors(s, nd, d.nbr_bins, d.domain_edge, oa_nbrs)
                else:
                    fillNeighbors(
                        s, nd, nd.bins[bin].nbr_bins, nd.bins[bin].domain_edge, oa_nbrs
                    )
                prev_domain = domain
        else:
            # ipos is in d
            bin = hashBin(  # noqa: F405
                ipos, d.min, d.bin_size, d.num_bins_x, d.num_bins_y
            )  # noqa: F405
            if bin != prev_bin:
                if p.recalculate_neighbors:
                    getNeighborIndices(  # noqa: F405
                        bin,
                        d.nbr_bins,
                        d.domain_edge,
                        d.num_bins_x,
                        d.num_bins_y,
                        d.num_bins_z,
                    )
                    fillNeighbors(s, d, d.nbr_bins, d.domain_edge, oa_nbrs)
                else:
                    fillNeighbors(
                        s, d, d.bins[bin].nbr_bins, d.bins[bin].domain_edge, oa_nbrs
                    )
                prev_bin = bin

        # Interactions with other atoms in this monomer, separated by more than
        # p->bond_cutoff bonds, not yet added to spatial bins
        for j in range(i):
            if not c.zm.isBonded(indices[i], indices[j], p.bond_cutoff):
                jtype = c.zm.entries[indices[j]].type.element_index
                jpos = c.zm.getPosition(indices[j], jpos)
                foldPosition(  # noqa: F405
                    jpos, p.system_min, p.system_max, p.system_size
                )  # noqa: F405
                r2 = getNearestSqDist(jpos, ipos, p.system_size, p.half_system_size)
                if r2 < ecut2:
                    E += energy_func(itype, jtype, r2)

        # Loop over spatial bins: oa_nbrs[0] is the center bin that holds the
        # tail atom; atoms in these bins could be in other chains or from
        # previous monomers in the same chain
        for j in range(27):
            oa = oa_nbrs[j]
            while oa and hasattr(oa, "next"):
                contrib = 1
                if p.isolate_chains and oa.chain != chain_index:
                    contrib = 0
                if oa.chain == chain_index and c.zm.isBonded(
                    indices[i], oa.atom, p.bond_cutoff
                ):
                    contrib = 0
                if contrib:
                    cn = s.chains[oa.chain]
                    jtype = cn.zm.entries[oa.atom].type.element_index

                    # printf("Domain %d: getting position %d of chain %d in sumInteractions\n",
                    # d->index, oa->atom+1, oa->chain+1);
                    # fflush(stdout);

                    jpos = cn.zm.getPosition(oa.atom, jpos)
                    foldPosition(  # noqa: F405
                        jpos, p.system_min, p.system_max, p.system_size
                    )  # noqa: F405
                    r2 = getNearestSqDist(jpos, ipos, p.system_size, p.half_system_size)
                    if r2 < ecut2:
                        E += energy_func(itype, jtype, r2)
                oa = oa.next
    if store_positions:
        # Unset temporarily cached positions
        for i in range(start_unset, num_indices):
            c.zm.clearPosition(indices[i])
    return E


# ============================================================================
# buildChains()
# ----------------------------------------------------------------------------
# Result: build all Chains in s->domains[index] to completion
# ============================================================================
def buildChains(s, p, domain_index):
    system_done = 0
    atom_indices = {}
    prev_angles = {}
    for i in range(p.max_monomer_atoms):
        atom_indices[i] = 0
        prev_angles[i] = 0.0
    kT_inv = 1.0 / (kB * p.temperature)  # noqa: F405
    m = Monomer()
    m.create()
    rng = s.rngs[domain_index]
    d = s.domains[domain_index]
    pos = np.zeros(3)
    # Vector p1, p3;

    while not system_done:
        for i in range(p.num_chains):
            # Add any pending atoms to d
            # pragma omp critical (pending_atom)
            s.getPendingAtoms(d, p)

            # Choose a Chain randomly
            chain_index = int(random.random() * p.num_chains)  # noqa: F405
            if domain_index != s.chains[chain_index].domain:
                continue  # on to the next Chain
            if (
                s.chains[chain_index].curr_monomer == s.chains[chain_index].num_monomers
                or s.chains[chain_index].dead
            ):
                break  # test system_done and restart dead Chains

            # The tail atom of s->chains[chain_index] is in
            # s->domains[domain_index]; only this instance of buildChains()
            # will access s->chains[chain_index]
            c = s.chains[chain_index]
            store_positions = (not p.recalculate_positions) and c.zm.num_positions > 2

            # Add a Monomer to the Chain ZMatrix
            m_prev = m
            nm = c.curr_monomer
            if (
                nm == c.num_monomers - 1
                and c.stereo.term
                and hasattr(c.stereo.term, "next")
            ):
                m = c.stereo.term
            else:
                m = c.stereo.getNextMonomer(rng)
            # pragma omp atomic
            m.selection_count += 1
            c.addMonomer(m, p.backbone_bond_length, 1)

            # printf("Domain %d: added monomer %d to chain %d\n", domain_index,
            #        c->curr_monomer, chain_index);
            # fflush(stdout);
            n0 = 0
            n1 = 0
            if 1 == c.curr_monomer:
                # Choose a position in Domain d for the first monomer
                for j in range(m.num_atoms):
                    atom_indices[j] = j
                num_test_indices = m.num_atoms
                i_mon = 3
                torsion_offset = 3
                # Positions for atoms 0 and 1 are set, below, so, if we cache
                # temp positions while summing hard core interactions, unset
                # all positions starting at index 2.
                start_unset = 2
                k = 0
                # Force single chain backbone to lie along x
                pos[0] = d.min[0] + random.random() * (  # noqa: F405
                    d.max[0] - d.min[0]
                )  # noqa: F405
                pos[1] = d.min[1] + random.random() * (  # noqa: F405
                    d.max[1] - d.min[1]
                )  # noqa: F405
                pos[2] = d.min[2] + random.random() * (  # noqa: F405
                    d.max[2] - d.min[2]
                )  # noqa: F405
                c.zm.setPosition(0, pos)
                # Select a random orientation of the atom 0 - atom 1 bond;
                # by convention, getPosition() orients this bond along x
                # if the atom 1 position is not set.
                th = 2.0 * np.pi * random.random()  # noqa: F405
                ph = np.arccos(2.0 * random.random() - 1.0)  # noqa: F405
                r = c.zm.entries[1].bond_length
                g = r * np.sin(ph)
                pos[0] += g * np.cos(th)
                pos[1] += g * np.sin(th)
                pos[2] += r * np.cos(ph)
                c.zm.setPosition(1, pos)
                Ei = sumInteractions(
                    s,
                    d,
                    chain_index,
                    atom_indices,
                    num_test_indices,
                    start_unset,
                    store_positions,
                    p,
                    Energy.energySelfAvoid,
                )
                k += 1
                while k < p.num_configs and (
                    Ei > 0.0
                    or rejectConfig(
                        c,
                        atom_indices,
                        num_test_indices,
                        store_positions,
                        start_unset,
                        p,
                    )
                ):
                    pos[0] = d.min[0] + random.random() * (  # noqa: F405
                        d.max[0] - d.min[0]
                    )  # noqa: F405
                    pos[1] = d.min[1] + random.random() * (  # noqa: F405
                        d.max[1] - d.min[1]
                    )  # noqa: F405
                    pos[2] = d.min[2] + random.random() * (  # noqa: F405
                        d.max[2] - d.min[2]
                    )  # noqa: F405
                    c.zm.setPosition(0, pos)
                    # Select a random orientation of the atom 0 - atom 1 bond;
                    # by convention, getPosition() orients this bond along x
                    # if the atom 1 position is not set.
                    th = 2.0 * np.pi * random.random()  # noqa: F405
                    ph = np.arccos(2.0 * random.random() - 1.0)  # noqa: F405
                    r = c.zm.entries[1].bond_length
                    g = r * np.sin(ph)
                    pos[0] += g * np.cos(th)
                    pos[1] += g * np.sin(th)
                    pos[2] += r * np.cos(ph)
                    c.zm.setPosition(1, pos)
                    Ei = sumInteractions(
                        s,
                        d,
                        chain_index,
                        atom_indices,
                        num_test_indices,
                        start_unset,
                        store_positions,
                        p,
                        Energy.energySelfAvoid,
                    )
                    k += 1
                if k == p.num_configs:
                    # Unable to find a good configuration

                    #   logMessage("Rejecting chain %d: ", chain_index+1);
                    #   logMessage(" unable to position monomer 1");

                    bad_chain += 1  # noqa: F405
                    break  # out of inner build loop
            else:
                # Not the first Monomer added to c
                i_mon = c.i_monomer[nm]
                torsion_offset = 2
                atom_indices[0] = c.zm.entries[i_mon].bond_index
                for j in range(m.num_atoms - 2):
                    atom_indices[j + 1] = i_mon + j
                num_test_indices = m.num_atoms - 1
                # Position of former tail atom (atom_indices[0]) was updated
                # in addMonomer(), so, if we cache temp positions while
                # summing interactions, unset all positions starting at
                # index 1.
                start_unset = 1
            # Now, i_mon is the index of the first torsion to adjust and
            # torsion_offset is the index of the corresponding m->torsions[]
            # entry; atom_indices[] holds the indices of the atoms whose
            # positions should be tested.
            num_torsions = m.num_bb - torsion_offset

            # Set initial torsions for new monomer
            Eb = 0.0
            wt = 1.0
            setTorsions(c, m, num_torsions, i_mon, torsion_offset, rng, Eb, wt, 0.0)

            # Sample configurations
            if p.sample_monte_carlo:
                Ei = Eb + sumInteractions(
                    s,
                    d,
                    chain_index,
                    atom_indices,
                    num_test_indices,
                    start_unset,
                    store_positions,
                    p,
                    p.energy_func,
                )
                Emax = Echosen = Ei
                for j in range(p.num_configs):
                    # Loop over configurations; Ei holds energy of the inital
                    # configuration
                    for k in range(num_torsions):
                        prev_angles[k] = c.zm.entries[i_mon + k].torsion_angle
                    Eb = 0.0
                    # Previously this added +/- torsion_step to the initial
                    # choice; now we choose from the available distribution
                    # of angles each time...

                    # &
                    # setTorsions(c, m, num_torsions, i_mon, torsion_offset, rng,
                    #             &Eb, NULL, p->torsion_step);

                    setTorsions(
                        c, m, num_torsions, i_mon, torsion_offset, rng, Eb, 0, 0.0
                    )
                    Ef = Eb + sumInteractions(
                        s,
                        d,
                        chain_index,
                        atom_indices,
                        num_test_indices,
                        start_unset,
                        store_positions,
                        p,
                        p.energy_func,
                    )
                    if Ef > Emax:
                        Emax = Ef
                    dE = Ef - Ei
                    if dE > 0 and np.exp(-dE * kT_inv) < random.random():  # noqa: F405
                        # Reject configuration
                        for k in range(num_torsions):
                            c.zm.entries[i_mon + k].torsion_angle = prev_angles[k]
                    else:
                        # Keep new configuration
                        Ei = Ef
                        Echosen = Ef
                wt = np.exp((Echosen - Emax) * kT_inv)
            bad_chain = rejectConfig(
                c, atom_indices, num_test_indices, store_positions, start_unset, p
            )
            n1 += 1
            while n1 < p.num_configs and bad_chain:
                if 1 == c.curr_monomer:
                    # Choose a position in Domain d for the first monomer
                    for j in range(m.num_atoms):
                        atom_indices[j] = j
                    num_test_indices = m.num_atoms
                    i_mon = 3
                    torsion_offset = 3
                    # Positions for atoms 0 and 1 are set, below, so, if we cache
                    # temp positions while summing hard core interactions, unset
                    # all positions starting at index 2.
                    start_unset = 2
                    k = 0
                    # Force single chain backbone to lie along x
                    pos[0] = d.min[0] + random.random() * (  # noqa: F405
                        d.max[0] - d.min[0]
                    )  # noqa: F405
                    pos[1] = d.min[1] + random.random() * (  # noqa: F405
                        d.max[1] - d.min[1]
                    )  # noqa: F405
                    pos[2] = d.min[2] + random.random() * (  # noqa: F405
                        d.max[2] - d.min[2]
                    )  # noqa: F405
                    c.zm.setPosition(0, pos)
                    # Select a random orientation of the atom 0 - atom 1 bond;
                    # by convention, getPosition() orients this bond along x
                    # if the atom 1 position is not set.
                    th = 2.0 * np.pi * random.random()  # noqa: F405
                    ph = np.arccos(2.0 * random.random() - 1.0)  # noqa: F405
                    r = c.zm.entries[1].bond_length
                    g = r * np.sin(ph)
                    pos[0] += g * np.cos(th)
                    pos[1] += g * np.sin(th)
                    pos[2] += r * np.cos(ph)
                    c.zm.setPosition(1, pos)
                    Ei = sumInteractions(
                        s,
                        d,
                        chain_index,
                        atom_indices,
                        num_test_indices,
                        start_unset,
                        store_positions,
                        p,
                        Energy.energySelfAvoid,
                    )
                    while k < p.num_configs and (
                        Ei > 0.0
                        or rejectConfig(
                            c,
                            atom_indices,
                            num_test_indices,
                            store_positions,
                            start_unset,
                            p,
                        )
                    ):
                        k += 1
                        pos[0] = d.min[0] + random.random() * (  # noqa: F405
                            d.max[0] - d.min[0]
                        )  # noqa: F405
                        pos[1] = d.min[1] + random.random() * (  # noqa: F405
                            d.max[1] - d.min[1]
                        )  # noqa: F405
                        pos[2] = d.min[2] + random.random() * (  # noqa: F405
                            d.max[2] - d.min[2]
                        )  # noqa: F405
                        c.zm.setPosition(0, pos)
                        # Select a random orientation of the atom 0 - atom 1 bond;
                        # by convention, getPosition() orients this bond along x
                        # if the atom 1 position is not set.
                        th = 2.0 * np.pi * random.random()  # noqa: F405
                        ph = np.arccos(2.0 * random.random() - 1.0)  # noqa: F405
                        r = c.zm.entries[1].bond_length
                        g = r * np.sin(ph)
                        pos[0] += g * np.cos(th)
                        pos[1] += g * np.sin(th)
                        pos[2] += r * np.cos(ph)
                        c.zm.setPosition(1, pos)
                        Ei = sumInteractions(
                            s,
                            d,
                            chain_index,
                            atom_indices,
                            num_test_indices,
                            start_unset,
                            store_positions,
                            p,
                            Energy.energySelfAvoid,
                        )
                    k += 1
                    if k == p.num_configs:
                        # Unable to find a good configuration

                        #   logMessage("Rejecting chain %d: ", chain_index+1);
                        #   logMessage(" unable to position monomer 1");

                        bad_chain += 1
                        break  # out of inner build loop
                else:
                    # Not the first Monomer added to c
                    i_mon = c.i_monomer[nm]
                    torsion_offset = 2
                    atom_indices[0] = c.zm.entries[i_mon].bond_index
                    for j in range(m.num_atoms - 2):
                        atom_indices[j + 1] = i_mon + j
                    num_test_indices = m.num_atoms - 1
                    # Position of former tail atom (atom_indices[0]) was updated
                    # in addMonomer(), so, if we cache temp positions while
                    # summing interactions, unset all positions starting at
                    # index 1.
                    start_unset = 1
                # Now, i_mon is the index of the first torsion to adjust and
                # torsion_offset is the index of the corresponding m->torsions[]
                # entry; atom_indices[] holds the indices of the atoms whose
                # positions should be tested.
                num_torsions = m.num_bb - torsion_offset

                # Set initial torsions for new monomer
                Eb = 0.0
                wt = 1.0
                setTorsions(c, m, num_torsions, i_mon, torsion_offset, rng, Eb, wt, 0.0)

                # Sample configurations
                if p.sample_monte_carlo:
                    Ei = Eb + sumInteractions(
                        s,
                        d,
                        chain_index,
                        atom_indices,
                        num_test_indices,
                        start_unset,
                        store_positions,
                        p,
                        p.energy_func,
                    )
                    Emax = Echosen = Ei
                    for j in range(p.num_configs):
                        # Loop over configurations; Ei holds energy of the inital
                        # configuration
                        for k in range(num_torsions):
                            prev_angles[k] = c.zm.entries[i_mon + k].torsion_angle
                        Eb = 0.0
                        # Previously this added +/- torsion_step to the initial
                        # choice; now we choose from the available distribution
                        # of angles each time...

                        # setTorsions(c, m, num_torsions, i_mon, torsion_offset, rng,
                        #             &Eb, NULL, p->torsion_step);

                        setTorsions(
                            c,
                            m,
                            num_torsions,
                            i_mon,
                            torsion_offset,
                            rng,
                            Eb,
                            NULL,  # noqa: F405
                            0.0,
                        )
                        Ef = Eb + sumInteractions(
                            s,
                            d,
                            chain_index,
                            atom_indices,
                            num_test_indices,
                            start_unset,
                            store_positions,
                            p,
                            p.energy_func,
                        )
                        if Ef > Emax:
                            Emax = Ef
                        dE = Ef - Ei
                        if (
                            dE > 0
                            and np.exp(-dE * kT_inv) < random.random()  # noqa: F405
                        ):  # noqa: F405
                            # Reject configuration
                            for k in range(num_torsions):
                                c.zm.entries[i_mon + k].torsion_angle = prev_angles[k]
                        else:
                            # Keep new configuration
                            Ei = Ef
                            Echosen = Ef
                    wt = np.exp((Echosen - Emax) * kT_inv)
                bad_chain = rejectConfig(
                    c, atom_indices, num_test_indices, store_positions, start_unset, p
                )
                n1 += 1
            if n1 == p.num_configs and bad_chain and c.curr_monomer > 2:
                # Adjust previous monomer
                j = c.i_monomer[c.curr_monomer - 2]
                d.removeDeadChains(i, j)
                while j < c.i_monomer[c.curr_monomer - 1]:
                    c.zm.clearPosition(j)
                    j += 1
                i_mon = c.i_monomer[c.curr_monomer - 2]
                atom_indices[0] = c.zm.entries[i_mon].bond_index
                k = 1
                for j in range(i_mon, c.i_monomer[c.curr_monomer - 1]):
                    atom_indices[k] = j
                    k += 1
                num_test_indices = k
                torsion_offset = 2
                num_torsions = m_prev.num_bb - torsion_offset
                start_unset = 1
                n1 = 0
                # Set initial torsions for previous monomer
                Eb = 0.0
                wt = 1.0
                setTorsions(
                    c, m_prev, num_torsions, i_mon, torsion_offset, rng, Eb, wt, 0.0
                )
                # Sample configurations
                if p.sample_monte_carlo:
                    Ei = Eb + sumInteractions(
                        s,
                        d,
                        chain_index,
                        atom_indices,
                        num_test_indices,
                        start_unset,
                        store_positions,
                        p,
                        p.energy_func,
                    )
                    Emax = Echosen = Ei
                    for j in range(p.num_configs):
                        # Loop over configurations; Ei holds energy of the
                        # inital configuration
                        for k in range(num_torsions):
                            prev_angles[k] = c.zm.entries[i_mon + k].torsion_angle
                        Eb = 0.0
                        # Previously this added +/- torsion_step to the initial
                        # choice; now we choose from the available distribution
                        # of angles each time...

                        # setTorsions(c, m_prev, num_torsions, i_mon,
                        #            torsion_offset, rng, &Eb, NULL,
                        #            p->torsion_step);

                        setTorsions(
                            c,
                            m_prev,
                            num_torsions,
                            i_mon,
                            torsion_offset,
                            rng,
                            Eb,
                            0.0,
                            0.0,
                        )
                        Ef = Eb + sumInteractions(
                            s,
                            d,
                            chain_index,
                            atom_indices,
                            num_test_indices,
                            start_unset,
                            store_positions,
                            p,
                            p.energy_func,
                        )
                        if Ef > Emax:
                            Emax = Ef
                        dE = Ef - Ei
                        if (
                            dE > 0
                            and np.exp(-dE * kT_inv) < random.random()  # noqa: F405
                        ):  # noqa: F405
                            # Reject configuration
                            for k in range(num_torsions):
                                c.zm.entries[i_mon + k].torsion_angle = prev_angles[k]
                        else:
                            # Keep new configuration
                            Ei = Ef
                            Echosen = Ef
                    wt = np.exp((Echosen - Emax) * kT_inv)
                bad_chain = rejectConfig(
                    c, atom_indices, num_test_indices, store_positions, start_unset, p
                )
                n1 += 1
                while n1 < p.num_configs and bad_chain:
                    # Set initial torsions for previous monomer
                    Eb = 0.0
                    wt = 1.0
                    setTorsions(
                        c, m_prev, num_torsions, i_mon, torsion_offset, rng, Eb, wt, 0.0
                    )
                    # Sample configurations
                    if p.sample_monte_carlo:
                        Ei = Eb + sumInteractions(
                            s,
                            d,
                            chain_index,
                            atom_indices,
                            num_test_indices,
                            start_unset,
                            store_positions,
                            p,
                            p.energy_func,
                        )
                        Emax = Echosen = Ei
                        for j in range(p.num_configs):
                            # Loop over configurations; Ei holds energy of the
                            # inital configuration
                            for k in range(num_torsions):
                                prev_angles[k] = c.zm.entries[i_mon + k].torsion_angle
                            Eb = 0.0
                            # Previously this added +/- torsion_step to the initial
                            # choice; now we choose from the available distribution
                            # of angles each time...

                            # setTorsions(c, m_prev, num_torsions, i_mon,
                            #            torsion_offset, rng, &Eb, NULL,
                            #            p->torsion_step);

                            setTorsions(
                                c,
                                m_prev,
                                num_torsions,
                                i_mon,
                                torsion_offset,
                                rng,
                                Eb,
                                0.0,
                                0.0,
                            )
                            Ef = Eb + sumInteractions(
                                s,
                                d,
                                chain_index,
                                atom_indices,
                                num_test_indices,
                                start_unset,
                                store_positions,
                                p,
                                p.energy_func,
                            )
                            if Ef > Emax:
                                Emax = Ef
                            dE = Ef - Ei
                            if (
                                dE > 0
                                and np.exp(-dE * kT_inv) < random.random()  # noqa: F405
                            ):  # noqa: F405
                                # Reject configuration
                                for k in range(num_torsions):
                                    c.zm.entries[i_mon + k].torsion_angle = prev_angles[
                                        k
                                    ]
                            else:
                                # Keep new configuration
                                Ei = Ef
                                Echosen = Ef
                        wt = np.exp((Echosen - Emax) * kT_inv)
                    bad_chain = rejectConfig(
                        c,
                        atom_indices,
                        num_test_indices,
                        store_positions,
                        start_unset,
                        p,
                    )
                    n1 += 1
                if bad_chain:
                    # Not able to reconfigure the previous monomer; set n0 to
                    # the value that will break out of the outer build loop
                    n0 = p.num_configs - 1
                else:
                    # Previous monomer is now ok, but keep bad_chain set so
                    # that the outer loop will reconfigure the tail, given the
                    # new state of the previous monomer
                    bad_chain += 1
                    # Found good, new configuration for previous monomer
                    if (not p.recalculate_positions) and c.zm.num_positions > 2:
                        for j in range(
                            c.i_monomer[c.curr_monomer - 2],
                            c.i_monomer[c.curr_monomer - 1],
                        ):
                            # Cache final unwrapped positions
                            pos = c.zm.getPosition(j, pos)
                            c.zm.setPosition(j, pos)
                            # Add folded positions to grid
                            foldPosition(  # noqa: F405
                                pos, p.system_min, p.system_max, p.system_size
                            )  # noqa: F405
                            if (
                                pos[0] > d.max[0]
                                or pos[0] < d.min[0]
                                or pos[1] > d.max[1]
                                or pos[1] < d.min[1]
                                or pos[2] > d.max[2]
                                or pos[2] < d.min[2]
                            ):
                                # pos is in another Domain
                                # pragma omp critical (pending_atom)
                                s.addPendingAtom(chain_index, j)
                            else:  # pos is in d
                                d.addAtom(pos, chain_index, j)
            # end adjust previous monomer
            # else:
            # TODO torsion delta ...
            n0 += 1
            while n0 < p.num_configs and bad_chain:
                n1 = 0
                if 1 == c.curr_monomer:
                    # Choose a position in Domain d for the first monomer
                    for j in range(m.num_atoms):
                        atom_indices[j] = j
                    num_test_indices = m.num_atoms
                    i_mon = 3
                    torsion_offset = 3
                    # Positions for atoms 0 and 1 are set, below, so, if we cache
                    # temp positions while summing hard core interactions, unset
                    # all positions starting at index 2.
                    start_unset = 2
                    k = 0
                    # Force single chain backbone to lie along x
                    pos[0] = d.min[0] + random.random() * (  # noqa: F405
                        d.max[0] - d.min[0]
                    )  # noqa: F405
                    pos[1] = d.min[1] + random.random() * (  # noqa: F405
                        d.max[1] - d.min[1]
                    )  # noqa: F405
                    pos[2] = d.min[2] + random.random() * (  # noqa: F405
                        d.max[2] - d.min[2]
                    )  # noqa: F405
                    c.zm.setPosition(0, pos)
                    # Select a random orientation of the atom 0 - atom 1 bond;
                    # by convention, getPosition() orients this bond along x
                    # if the atom 1 position is not set.
                    th = 2.0 * np.pi * random.random()  # noqa: F405
                    ph = np.arccos(2.0 * random.random() - 1.0)  # noqa: F405
                    r = c.zm.entries[1].bond_length
                    g = r * np.sin(ph)
                    pos[0] += g * np.cos(th)
                    pos[1] += g * np.sin(th)
                    pos[2] += r * np.cos(ph)
                    c.zm.setPosition(1, pos)
                    Ei = sumInteractions(
                        s,
                        d,
                        chain_index,
                        atom_indices,
                        num_test_indices,
                        start_unset,
                        store_positions,
                        p,
                        Energy.energySelfAvoid,
                    )
                    k += 1
                    while k < p.num_configs and (
                        Ei > 0.0
                        or rejectConfig(
                            c,
                            atom_indices,
                            num_test_indices,
                            store_positions,
                            start_unset,
                            p,
                        )
                    ):
                        pos[0] = d.min[0] + random.random() * (  # noqa: F405
                            d.max[0] - d.min[0]
                        )  # noqa: F405
                        pos[1] = d.min[1] + random.random() * (  # noqa: F405
                            d.max[1] - d.min[1]
                        )  # noqa: F405
                        pos[2] = d.min[2] + random.random() * (  # noqa: F405
                            d.max[2] - d.min[2]
                        )  # noqa: F405
                        c.zm.setPosition(0, pos)
                        # Select a random orientation of the atom 0 - atom 1 bond;
                        # by convention, getPosition() orients this bond along x
                        # if the atom 1 position is not set.
                        th = 2.0 * np.pi * random.random()  # noqa: F405
                        ph = np.arccos(2.0 * random.random() - 1.0)  # noqa: F405
                        r = c.zm.entries[1].bond_length
                        g = r * np.sin(ph)
                        pos[0] += g * np.cos(th)
                        pos[1] += g * np.sin(th)
                        pos[2] += r * np.cos(ph)
                        c.zm.setPosition(1, pos)
                        Ei = sumInteractions(
                            s,
                            d,
                            chain_index,
                            atom_indices,
                            num_test_indices,
                            start_unset,
                            store_positions,
                            p,
                            Energy.energySelfAvoid,
                        )
                        k += 1
                    if k == p.num_configs:
                        # Unable to find a good configuration

                        #   logMessage("Rejecting chain %d: ", chain_index+1);
                        #   logMessage(" unable to position monomer 1");

                        bad_chain += 1
                        break  # out of inner build loop
                else:
                    # Not the first Monomer added to c
                    i_mon = c.i_monomer[nm]
                    torsion_offset = 2
                    atom_indices[0] = c.zm.entries[i_mon].bond_index
                    for j in range(m.num_atoms - 2):
                        atom_indices[j + 1] = i_mon + j
                    num_test_indices = m.num_atoms - 1
                    # Position of former tail atom (atom_indices[0]) was updated
                    # in addMonomer(), so, if we cache temp positions while
                    # summing interactions, unset all positions starting at
                    # index 1.
                    start_unset = 1
                # Now, i_mon is the index of the first torsion to adjust and
                # torsion_offset is the index of the corresponding m->torsions[]
                # entry; atom_indices[] holds the indices of the atoms whose
                # positions should be tested.
                num_torsions = m.num_bb - torsion_offset

                # Set initial torsions for new monomer
                Eb = 0.0
                wt = 1.0
                setTorsions(c, m, num_torsions, i_mon, torsion_offset, rng, Eb, wt, 0.0)

                # Sample configurations
                if p.sample_monte_carlo:
                    Ei = Eb + sumInteractions(
                        s,
                        d,
                        chain_index,
                        atom_indices,
                        num_test_indices,
                        start_unset,
                        store_positions,
                        p,
                        p.energy_func,
                    )
                    Emax = Echosen = Ei
                    for j in range(p.num_configs):
                        # Loop over configurations; Ei holds energy of the inital
                        # configuration
                        for k in range(num_torsions):
                            prev_angles[k] = c.zm.entries[i_mon + k].torsion_angle
                        Eb = 0.0
                        # Previously this added +/- torsion_step to the initial
                        # choice; now we choose from the available distribution
                        # of angles each time...

                        # setTorsions(c, m, num_torsions, i_mon, torsion_offset, rng,
                        #             &Eb, NULL, p->torsion_step);

                        setTorsions(
                            c,
                            m,
                            num_torsions,
                            i_mon,
                            torsion_offset,
                            rng,
                            Eb,
                            NULL,  # noqa: F405
                            0.0,
                        )
                        Ef = Eb + sumInteractions(
                            s,
                            d,
                            chain_index,
                            atom_indices,
                            num_test_indices,
                            start_unset,
                            store_positions,
                            p,
                            p.energy_func,
                        )
                        if Ef > Emax:
                            Emax = Ef
                        dE = Ef - Ei
                        if (
                            dE > 0
                            and np.exp(-dE * kT_inv) < random.random()  # noqa: F405
                        ):  # noqa: F405
                            # Reject configuration
                            for k in range(num_torsions):
                                c.zm.entries[i_mon + k].torsion_angle = prev_angles[k]
                        else:
                            # Keep new configuration
                            Ei = Ef
                            Echosen = Ef
                    wt = np.exp((Echosen - Emax) * kT_inv)
                bad_chain = rejectConfig(
                    c, atom_indices, num_test_indices, store_positions, start_unset, p
                )
                n1 += 1
                while n1 < p.num_configs and bad_chain:
                    if 1 == c.curr_monomer:
                        # Choose a position in Domain d for the first monomer
                        for j in range(m.num_atoms):
                            atom_indices[j] = j
                        num_test_indices = m.num_atoms
                        i_mon = 3
                        torsion_offset = 3
                        # Positions for atoms 0 and 1 are set, below, so, if we cache
                        # temp positions while summing hard core interactions, unset
                        # all positions starting at index 2.
                        start_unset = 2
                        k = 0
                        # Force single chain backbone to lie along x
                        pos[0] = d.min[0] + random.random() * (  # noqa: F405
                            d.max[0] - d.min[0]
                        )  # noqa: F405
                        pos[1] = d.min[1] + random.random() * (  # noqa: F405
                            d.max[1] - d.min[1]
                        )  # noqa: F405
                        pos[2] = d.min[2] + random.random() * (  # noqa: F405
                            d.max[2] - d.min[2]
                        )  # noqa: F405
                        c.zm.setPosition(0, pos)
                        # Select a random orientation of the atom 0 - atom 1 bond;
                        # by convention, getPosition() orients this bond along x
                        # if the atom 1 position is not set.
                        th = 2.0 * np.pi * random.random()  # noqa: F405
                        ph = np.arccos(2.0 * random.random() - 1.0)  # noqa: F405
                        r = c.zm.entries[1].bond_length
                        g = r * np.sin(ph)
                        pos[0] += g * np.cos(th)
                        pos[1] += g * np.sin(th)
                        pos[2] += r * np.cos(ph)
                        c.zm.setPosition(1, pos)
                        Ei = sumInteractions(
                            s,
                            d,
                            chain_index,
                            atom_indices,
                            num_test_indices,
                            start_unset,
                            store_positions,
                            p,
                            energySelfAvoid,  # noqa: F405
                        )
                        while k < p.num_configs and (
                            Ei > 0.0
                            or rejectConfig(
                                c,
                                atom_indices,
                                num_test_indices,
                                store_positions,
                                start_unset,
                                p,
                            )
                        ):
                            k += 1
                            pos[0] = d.min[0] + random.random() * (  # noqa: F405
                                d.max[0] - d.min[0]
                            )  # noqa: F405
                            pos[1] = d.min[1] + random.random() * (  # noqa: F405
                                d.max[1] - d.min[1]
                            )  # noqa: F405
                            pos[2] = d.min[2] + random.random() * (  # noqa: F405
                                d.max[2] - d.min[2]
                            )  # noqa: F405
                            c.zm.setPosition(0, pos)
                            # Select a random orientation of the atom 0 - atom 1 bond;
                            # by convention, getPosition() orients this bond along x
                            # if the atom 1 position is not set.
                            th = 2.0 * np.pi * random.random()  # noqa: F405
                            ph = np.arccos(2.0 * random.random() - 1.0)  # noqa: F405
                            r = c.zm.entries[1].bond_length
                            g = r * np.sin(ph)
                            pos[0] += g * np.cos(th)
                            pos[1] += g * np.sin(th)
                            pos[2] += r * np.cos(ph)
                            c.zm.setPosition(1, pos)
                            Ei = sumInteractions(
                                s,
                                d,
                                chain_index,
                                atom_indices,
                                num_test_indices,
                                start_unset,
                                store_positions,
                                p,
                                Energy.energySelfAvoid,
                            )
                        k += 1
                        if k == p.num_configs:
                            # Unable to find a good configuration

                            #   logMessage("Rejecting chain %d: ", chain_index+1);
                            #   logMessage(" unable to position monomer 1");

                            bad_chain += 1
                            break  # out of inner build loop
                    else:
                        # Not the first Monomer added to c
                        i_mon = c.i_monomer[nm]
                        torsion_offset = 2
                        atom_indices[0] = c.zm.entries[i_mon].bond_index
                        for j in range(m.num_atoms - 2):
                            atom_indices[j + 1] = i_mon + j
                        num_test_indices = m.num_atoms - 1
                        # Position of former tail atom (atom_indices[0]) was updated
                        # in addMonomer(), so, if we cache temp positions while
                        # summing interactions, unset all positions starting at
                        # index 1.
                        start_unset = 1
                    # Now, i_mon is the index of the first torsion to adjust and
                    # torsion_offset is the index of the corresponding m->torsions[]
                    # entry; atom_indices[] holds the indices of the atoms whose
                    # positions should be tested.
                    num_torsions = m.num_bb - torsion_offset

                    # Set initial torsions for new monomer
                    Eb = 0.0
                    wt = 1.0
                    setTorsions(
                        c, m, num_torsions, i_mon, torsion_offset, rng, Eb, wt, 0.0
                    )

                    # Sample configurations
                    if p.sample_monte_carlo:
                        Ei = Eb + sumInteractions(
                            s,
                            d,
                            chain_index,
                            atom_indices,
                            num_test_indices,
                            start_unset,
                            store_positions,
                            p,
                            p.energy_func,
                        )
                        Emax = Echosen = Ei
                        for j in range(p.num_configs):
                            # Loop over configurations; Ei holds energy of the inital
                            # configuration
                            for k in range(num_torsions):
                                prev_angles[k] = c.zm.entries[i_mon + k].torsion_angle
                            Eb = 0.0
                            # Previously this added +/- torsion_step to the initial
                            # choice; now we choose from the available distribution
                            # of angles each time...

                            # setTorsions(c, m, num_torsions, i_mon, torsion_offset, rng,
                            #             &Eb, NULL, p->torsion_step);

                            setTorsions(
                                c,
                                m,
                                num_torsions,
                                i_mon,
                                torsion_offset,
                                rng,
                                Eb,
                                NULL,  # noqa: F405
                                0.0,
                            )
                            Ef = Eb + sumInteractions(
                                s,
                                d,
                                chain_index,
                                atom_indices,
                                num_test_indices,
                                start_unset,
                                store_positions,
                                p,
                                p.energy_func,
                            )
                            if Ef > Emax:
                                Emax = Ef
                            dE = Ef - Ei
                            if (
                                dE > 0
                                and np.exp(-dE * kT_inv) < random.random()  # noqa: F405
                            ):  # noqa: F405
                                # Reject configuration
                                for k in range(num_torsions):
                                    c.zm.entries[i_mon + k].torsion_angle = prev_angles[
                                        k
                                    ]
                            else:
                                # Keep new configuration
                                Ei = Ef
                                Echosen = Ef
                        wt = np.exp((Echosen - Emax) * kT_inv)
                    bad_chain = rejectConfig(
                        c,
                        atom_indices,
                        num_test_indices,
                        store_positions,
                        start_unset,
                        p,
                    )
                    n1 += 1
                if n1 == p.num_configs and bad_chain and c.curr_monomer > 2:
                    # Adjust previous monomer
                    j = c.i_monomer[c.curr_monomer - 2]
                    d.removeDeadChains(i, j)
                    while j < c.i_monomer[c.curr_monomer - 1]:
                        c.zm.clearPosition(j)
                        j += 1
                    i_mon = c.i_monomer[c.curr_monomer - 2]
                    atom_indices[0] = c.zm.entries[i_mon].bond_index
                    k = 1
                    for j in range(i_mon, c.i_monomer[c.curr_monomer - 1]):
                        atom_indices[k] = j
                        k += 1
                    num_test_indices = k
                    torsion_offset = 2
                    num_torsions = m_prev.num_bb - torsion_offset
                    start_unset = 1
                    n1 = 0
                    # Set initial torsions for previous monomer
                    Eb = 0.0
                    wt = 1.0
                    setTorsions(
                        c, m_prev, num_torsions, i_mon, torsion_offset, rng, Eb, wt, 0.0
                    )
                    # Sample configurations
                    if p.sample_monte_carlo:
                        Ei = Eb + sumInteractions(
                            s,
                            d,
                            chain_index,
                            atom_indices,
                            num_test_indices,
                            start_unset,
                            store_positions,
                            p,
                            p.energy_func,
                        )
                        Emax = Echosen = Ei
                        for j in range(p.num_configs):
                            # Loop over configurations; Ei holds energy of the
                            # inital configuration
                            for k in range(num_torsions):
                                prev_angles[k] = c.zm.entries[i_mon + k].torsion_angle
                            Eb = 0.0
                            # Previously this added +/- torsion_step to the initial
                            # choice; now we choose from the available distribution
                            # of angles each time...

                            # setTorsions(c, m_prev, num_torsions, i_mon,
                            #            torsion_offset, rng, &Eb, NULL,
                            #            p->torsion_step);

                            setTorsions(
                                c,
                                m_prev,
                                num_torsions,
                                i_mon,
                                torsion_offset,
                                rng,
                                Eb,
                                0.0,
                                0.0,
                            )
                            Ef = Eb + sumInteractions(
                                s,
                                d,
                                chain_index,
                                atom_indices,
                                num_test_indices,
                                start_unset,
                                store_positions,
                                p,
                                p.energy_func,
                            )
                            if Ef > Emax:
                                Emax = Ef
                            dE = Ef - Ei
                            if (
                                dE > 0
                                and np.exp(-dE * kT_inv) < random.random()  # noqa: F405
                            ):  # noqa: F405
                                # Reject configuration
                                for k in range(num_torsions):
                                    c.zm.entries[i_mon + k].torsion_angle = prev_angles[
                                        k
                                    ]
                            else:
                                # Keep new configuration
                                Ei = Ef
                                Echosen = Ef
                        wt = np.exp((Echosen - Emax) * kT_inv)
                    bad_chain = rejectConfig(
                        c,
                        atom_indices,
                        num_test_indices,
                        store_positions,
                        start_unset,
                        p,
                    )
                    n1 += 1
                    while n1 < p.num_configs and bad_chain:
                        # Set initial torsions for previous monomer
                        Eb = 0.0
                        wt = 1.0
                        setTorsions(
                            c,
                            m_prev,
                            num_torsions,
                            i_mon,
                            torsion_offset,
                            rng,
                            Eb,
                            wt,
                            0.0,
                        )
                        # Sample configurations
                        if p.sample_monte_carlo:
                            Ei = Eb + sumInteractions(
                                s,
                                d,
                                chain_index,
                                atom_indices,
                                num_test_indices,
                                start_unset,
                                store_positions,
                                p,
                                p.energy_func,
                            )
                            Emax = Echosen = Ei
                            for j in range(p.num_configs):
                                # Loop over configurations; Ei holds energy of the
                                # inital configuration
                                for k in range(num_torsions):
                                    prev_angles[k] = c.zm.entries[
                                        i_mon + k
                                    ].torsion_angle
                                Eb = 0.0
                                # Previously this added +/- torsion_step to the initial
                                # choice; now we choose from the available distribution
                                # of angles each time...

                                # setTorsions(c, m_prev, num_torsions, i_mon,
                                #            torsion_offset, rng, &Eb, NULL,
                                #            p->torsion_step);

                                setTorsions(
                                    c,
                                    m_prev,
                                    num_torsions,
                                    i_mon,
                                    torsion_offset,
                                    rng,
                                    Eb,
                                    0.0,
                                    0.0,
                                )
                                Ef = Eb + sumInteractions(
                                    s,
                                    d,
                                    chain_index,
                                    atom_indices,
                                    num_test_indices,
                                    start_unset,
                                    store_positions,
                                    p,
                                    p.energy_func,
                                )
                                if Ef > Emax:
                                    Emax = Ef
                                dE = Ef - Ei
                                if (
                                    dE > 0
                                    and np.exp(-dE * kT_inv)
                                    < random.random()  # noqa: F405
                                ):  # noqa: F405
                                    # Reject configuration
                                    for k in range(num_torsions):
                                        c.zm.entries[
                                            i_mon + k
                                        ].torsion_angle = prev_angles[k]
                                else:
                                    # Keep new configuration
                                    Ei = Ef
                                    Echosen = Ef
                            wt = np.exp((Echosen - Emax) * kT_inv)
                        bad_chain = rejectConfig(
                            c,
                            atom_indices,
                            num_test_indices,
                            store_positions,
                            start_unset,
                            p,
                        )
                        n1 += 1
                    if bad_chain:
                        # Not able to reconfigure the previous monomer; set n0 to
                        # the value that will break out of the outer build loop
                        n0 = p.num_configs - 1
                    else:
                        # Previous monomer is now ok, but keep bad_chain set so
                        # that the outer loop will reconfigure the tail, given the
                        # new state of the previous monomer
                        bad_chain += 1
                        # Found good, new configuration for previous monomer
                        if (not p.recalculate_positions) and c.zm.num_positions > 2:
                            for j in range(
                                c.i_monomer[c.curr_monomer - 2],
                                c.i_monomer[c.curr_monomer - 1],
                            ):
                                # Cache final unwrapped positions
                                pos = c.zm.getPosition(j, pos)
                                c.zm.setPosition(j, pos)
                                # Add folded positions to grid
                                foldPosition(  # noqa: F405
                                    pos, p.system_min, p.system_max, p.system_size
                                )
                                if (
                                    pos[0] > d.max[0]
                                    or pos[0] < d.min[0]
                                    or pos[1] > d.max[1]
                                    or pos[1] < d.min[1]
                                    or pos[2] > d.max[2]
                                    or pos[2] < d.min[2]
                                ):
                                    # pos is in another Domain
                                    s.addPendingAtom(chain_index, j)
                                else:  # pos is in d
                                    d.addAtom(pos, chain_index, j)
                # end adjust previous monomer
                # else:
                # TODO torsion delta ...
                n0 += 1
            if n0 == p.num_configs and bad_chain:
                c.dead += 1
                continue  # on to the next Chain

            # Found good configuration for the current chain
            if (not p.recalculate_positions) and c.zm.num_positions > 2:
                # Cache final unwrapped positions to ease recursive burden in
                # subsequent calls to getPosition()
                for j in range(c.i_monomer[nm], c.curr_atom):

                    # printf("Domain %d: getting position %d of chain %d in final lock\n",
                    #        d->index, j+1, chain_index+1);
                    # fflush(stdout);

                    pos = c.zm.getPosition(j, pos)
                    c.zm.setPosition(j, pos)

            # Record torsion selections
            for j in range(num_torsions):
                c.torsion_count[int(c.zm.entries[i_mon + j].torsion_angle)] += 1

            # Update length and weight
            c.length[nm] = c.getChainLength()
            c.weight[nm] = wt  # FIXME

            # Add folded tail atom positions to grid
            for j in range(c.i_monomer[nm], c.curr_atom):

                # printf("Domain %d: getting position %d of chain %d in grid add\n",
                #        d->index, j+1, chain_index+1);
                # fflush(stdout);

                pos = c.zm.getPosition(j, pos)
                foldPosition(  # noqa: F405
                    pos, p.system_min, p.system_max, p.system_size
                )  # noqa: F405
                if (
                    pos[0] > d.max[0]
                    or pos[0] < d.min[0]
                    or pos[1] > d.max[1]
                    or pos[1] < d.min[1]
                    or pos[2] > d.max[2]
                    or pos[2] < d.min[2]
                ):
                    # pos is in another Domain */

                    # printf("Domain %d: atom %d of chain %d changed Domain\n",
                    #        d->index, j+1, chain_index+1);
                    # fflush(stdout);

                    # pragma omp critical (pending_atom)
                    s.addPendingAtom(chain_index, j)
                else:  # pos is in d
                    d.addAtom(pos, chain_index, j)

            # Report progress
            updateStatus()

            # Update c->domain; if c->domain != domain_index, another thread will
            # continue building c

            # printf("Domain %d: getting position %d of chain %d in domain update\n",
            #        d->index, c->tail_index+1, chain_index+1);
            # fflush(stdout);

            pos = c.zm.getPosition(c.tail_index, pos)
            foldPosition(pos, p.system_min, p.system_max, p.system_size)  # noqa: F405
            c.domain = hashBin(  # noqa: F405
                pos, p.system_min, p.domain_size, p.num_domains_x, p.num_domains_y
            )
        # End loop over Chains

        # Remove dead Chains
        for i in range(p.num_chains):
            if s.chains[i].dead:
                d.removeDeadChains(i, 0)
                if s.chains[i].domain == domain_index:
                    s.chains[i].resetChain()

        # Check exit condition: all Chains have reached num_monomers
        chains_done = 0
        for i in range(p.num_chains):
            # Previously tested s->chains[i]->dead here, but dead Chains should
            # be restarted, not counted against completed Chains...
            if s.chains[i].curr_monomer == s.chains[i].num_monomers:
                chains_done += 1
        if chains_done == p.num_chains:
            system_done += 1
    # end while !system_done
    # bruteForceCheck(s, d, p, 1.3);
    atom_indices = {}
    prev_angles = {}
