# This script simulates data with msprime, as described in the Skov et al., 2018 paper (see https://doi.org/10.1371/journal.pgen.1007641 and https://github.com/LauritsSkov/Introgression-detection/tree/master/Simulating%20data). WARNING: I (EKM) have not carefully checked this script -- please verify that it works before using in a study!!!

import argparse
import msprime as msp
import numpy as np
import sys


def combine_segs(segs, get_segs = False):
    """
    Taken from:
    https://github.com/LauritsSkov/Introgression-detection/tree/master/Simulating%20data
    """
    merged = np.empty([0, 2])
    if len(segs) == 0:
        if get_segs:
            return([])
        else:
            return(0)
    sorted_segs = segs[np.argsort(segs[:, 0]), :]
    for higher in sorted_segs:
        if len(merged) == 0:
            merged = np.vstack([merged, higher])            
        else:
            lower = merged[-1, :]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1, :] = (lower[0], upper_bound) 
            else:
                merged = np.vstack([merged, higher])
    if get_segs:
        return(merged)
    else:
        return(np.sum(merged[:, 1] - merged[:, 0])/seq_len)


def simulate_under_skov2018(nceu, nyri, nsite, seed):
    """
    Taken from 
    https://github.com/LauritsSkov/Introgression-detection/tree/master/Simulating%20data
    """
    # Generation time, mutation rate and recomination rate
    gen_time = 29.0 
    rec_rate = 1.2e-8 
    mu = 1.2e-8

    # Population sizes
    N_DE = 5000
    N_AF = 27122
    N_ancestral_eurasians = 5000
    N_ancestral_humans = 7000
    N_neanderthal_human_ancestral = 18296
    N_Europe = 3899
    N_Asia = 5054

    N_bottleneck_eurasians = 1305
    N_bottleneck_nonAfricans = 250

    # Split Times
    T_human_archaic = 656908 / gen_time
    T_bottleneck_nonAfricans = 62041 / gen_time
    T_eurasians = 57100 / gen_time
    T_europa_asia = 41997 / gen_time

    # Bottleneck duration
    T_bottlenecks = 100

    # Admixture times and proportions from archaic humans
    admixtureproportion = 5 #(in percent)
    T_GF_DE = 54000 / gen_time


    # Number of samples
    n_outgroup = nyri # 100
    n_ingroup = nceu  # 2
    samples = [msp.Sample(0, 0)] * n_outgroup + [msp.Sample(1, 0)]*n_ingroup

    population_configurations = [
        msp.PopulationConfiguration(initial_size = N_AF),
        msp.PopulationConfiguration(initial_size = N_Europe),
        msp.PopulationConfiguration(initial_size = N_Asia),
        msp.PopulationConfiguration(initial_size = N_DE)   
    ]

    demographic_events = [
        # European and Asian populations merge
        msp.MassMigration(
            time = T_europa_asia, source = 2, destination = 1, proportion = 1.0),
        msp.PopulationParametersChange(
            time = T_europa_asia, initial_size = N_ancestral_eurasians, 
            growth_rate = 0, population_id = 1),

        # Archaic admixture into Europeans
        msp.MassMigration(
            time = T_GF_DE, source = 1, destination = 3,
            proportion = admixtureproportion/100.0),

        # Eurasians experience bottleneck
        msp.PopulationParametersChange(
            time = T_eurasians - T_bottlenecks, 
            initial_size = N_bottleneck_eurasians, growth_rate = 0, population_id = 1),
        msp.PopulationParametersChange(
            time = T_eurasians, 
            initial_size = N_ancestral_eurasians, growth_rate = 0, population_id = 1),

        # Out of africa bottleneck in Eurasians and merge in Africans
        msp.PopulationParametersChange(
            time = T_bottleneck_nonAfricans - T_bottlenecks, 
            initial_size = N_bottleneck_nonAfricans, growth_rate = 0, population_id = 1),
        msp.MassMigration(
            time = T_bottleneck_nonAfricans, source = 1, destination = 0, proportion = 1.0),
        msp.PopulationParametersChange(
            time = T_bottleneck_nonAfricans, 
            initial_size = N_ancestral_humans, growth_rate = 0, population_id = 0),

        # Archiac and modern humans merge 650,000 years ago
        msp.MassMigration(
            time = T_human_archaic, source = 3, destination = 0, proportion = 1.0),
        msp.PopulationParametersChange(
            time = T_human_archaic, 
            initial_size = N_neanderthal_human_ancestral, 
            growth_rate = 0, population_id = 0), 
    ]

    ts = msp.simulate(
        samples = samples,
        population_configurations = population_configurations,
        demographic_events = demographic_events,
        length = nsite,
        recombination_rate = rec_rate,
        random_seed = seed,
        mutation_rate = mu,
        record_migrations = True
    )

    return ts


def write_eigenstrat(ts, nsites):
    chrm = '1'
    Testpopulation = ts.get_samples(1)
    AF_ids = ts.get_samples(0)

    # Write SNP positions
    with open("out.snp", 'w') as fo:
        for variant in ts.variants():
            pos = int(variant.site.position)
            dist = str(pos / nsites)
            pos = str(pos)
            row = '\t'.join([chrm + ':' + pos, chrm, dist, pos, 'A', 'G'])
            fo.write(row + '\n')

    # Write genotype output
    AF_ids = ts.get_samples(0)
    with open('out.1.geno', 'w') as fo:
        for variant in ts.variants():
            row = ''.join([str(x) for x in variant.genotypes[AF_ids]]) + '\n'
            fo.write(row)

    Testpopulation = ts.get_samples(1)
    with open('out.ADMIXED.geno', 'w') as fo:
        for variant in ts.variants():
            row = ''.join([str(x) for x in variant.genotypes[Testpopulation]]) + '\n'
            fo.write(row)

    with open('out.ADMIXED.ind', 'w') as fo:
        for i, x in enumerate(Testpopulation):
            row = "0:" + str(i) + "\tU\tADMIXED\n"
            fo.write(row)


def get_true_introgressed_segments(ts):
    """
    Taken from:
    https://github.com/LauritsSkov/Introgression-detection/tree/master/Simulating%20data
    """
    # Keep track of which segments are actually introgressed (in this case from pop 3 into pop 1)
    Testpopulation = ts.get_samples(1)
    AF_ids = ts.get_samples(0)

    de_seg = {i: [] for i in Testpopulation}

    for mr in ts.migrations():
        if mr.source == 1 and mr.dest == 3:
            for tree in ts.trees(leaf_lists=True):
                if mr.left > tree.get_interval()[0]:
                    continue
                if mr.right <= tree.get_interval()[0]:
                    break
                for l in tree.leaves(mr.node):
                    if l in Testpopulation:
                        #print l, mr
                        de_seg[l].append(tree.get_interval())

    true_de_segs = [combine_segs(np.array(de_seg[i]), True) for i in sorted(de_seg.keys())]
    return true_de_segs


def write_true_introgressed_segs(ts, true_de_segs):
    # This produces a .bed file with four columns
    # chr1
    # Start site of archaic segment
    # End site of archaic segment
    # Haplotype ID (0-100 for this sim)
    #
    # Need to turn it into a format where 
    # Matrix 0/1
    # Column is haplotype
    # Row is SNP -- requires us to extract SNPs...

    #chrom = 1
    #with open('archaic_segments_chr.bed', 'w') as out:
    #    for haplotype, archaic_segments in enumerate(true_de_segs):
    #        for archaic_segment in archaic_segments:
    #            out.write('chr{}\t{}\t{}\t{}\n'.format(chrom, int(archaic_segment[0]), int(archaic_segment[1]), haplotype))

    with open('out.ADMIXED.anc', 'w') as fo:
        for variant in ts.variants():
            row = ''
            pos = int(variant.site.position)
            for haplotype, archaic_segments in enumerate(true_de_segs):
                keep = '0'
                for archaic_segment in archaic_segments:
                    start = int(archaic_segment[0])
                    end = int(archaic_segment[1])
                    if (start <= pos) and (pos <= end):
                        keep = '1'
                row = row + keep
            fo.write(row + '\n')


def main(args):
    """
    """
    ts = simulate_under_skov2018(args.nceu, args.nyri, args.nsites, args.seed)
    #ts = simulate_under_skov2018(100, 100, 1000000, 12345)

    write_eigenstrat(ts, args.nsites)
    #write_eigenstrat(ts, 1000000)

    true_de_segs = get_true_introgressed_segments(ts)
    write_true_introgressed_segs(ts, true_de_segs)

    sys.stdout.write("# done")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-nc", "--nceu",
                        type=int, default=100,
                        help="No. haplotypes in CEU (target) population")
    parser.add_argument("-ny", "--nyri",
                        type=int, default=100,
                        help="No. haplotypes in YRI (reference) population")
    parser.add_argument("-ns", "--nsites",
                        type=int, default=5000000,
                        help="No. sites per locus (window size in ArchIE)")
    parser.add_argument("-x", "--seed",
                        type=int,
                        help="Random seed")
    main(parser.parse_args())
