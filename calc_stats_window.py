import argparse
import sys
import math
import csv
import copy
import scipy

import numpy as np
import scipy.stats as sp
import sklearn.cluster
import sklearn.metrics.pairwise

from bisect import bisect_left

import warnings
warnings.filterwarnings("ignore")

def min_dist_ref(genotypes, ref_genotypes, focal_idx):
    f = np.array(genotypes)[:,focal_idx]
    t_r = np.transpose(np.array(ref_genotypes))
    d = []
    for r in t_r:
        r_np = np.array(r)
        arr = np.array([f, r_np])
        d.append(np.max(sklearn.metrics.pairwise.pairwise_distances(arr)))
    return(np.min(d))

def distance_vector(genotypes, focal_idx):
    dist = sklearn.metrics.pairwise.pairwise_distances(np.transpose(genotypes))
    focal_dist = dist[focal_idx]
    focal_dist.sort()
    m = np.mean(focal_dist)
    var = np.var(focal_dist)
    skew = sp.skew(focal_dist)
    kurtosis = sp.kurtosis(focal_dist)
    return([m] + [var] + [skew] + [kurtosis])

def num_private(clean_geno, focal_idx):
    """"input is a genotype matrix with (rows=samples), (columns=snps) and SNPs in reference removed"""
    focal_hap = clean_geno[focal_idx]
    return(np.sum(focal_hap))


def parse_snp(snp):
    """Returns a list of SNPs"""
    pos = []
    with open(snp, 'r') as f:
        for line in f:
            pos.append(int(line.split("\t")[3]))
    return(pos)

def parse_geno(geno_f):
    """Returns a list of lists of genotypes (2D array) (sample x snp)"""
    geno = []
    with open(geno_f, 'r') as f:
            for line in f:
                ll = list(line)[:-1]
                li = [int(i) for i in ll]
                geno.append(li)
    return(geno)

def parse_anc(anc_f, n_samples):
    """Returns a list n_samples long of the number of archaic bases per sample"""
    anc = []
    with open(anc_f, 'r') as f:
        for line in f:
            anc.append(list(line)[:-1])
    return(anc)

def parse_ind(ind_f):
    inds = []
    with open(ind_f, 'r') as f:
        for line in f:
            inds.append(line.split("\t")[0])
            inds.append(line.split("\t")[0]) # do it twice so that each individual shows up twice (diploid)
    return(inds)

def label(bases, snp_list, n_sites, arch_thresh, not_arch_thresh):
    archaic_anc = 0
    prev_pos = 0
    previous_anc = 0
    for idx,base in enumerate(bases):
        s = snp_list[idx]
        if base == "1" and previous_anc == "1":
            archaic_anc = archaic_anc + (s-previous_pos)
        previous_pos = s
        previous_anc = base
    prop = archaic_anc/n_sites
    if prop > arch_thresh:
        return([1,0,0, prop])
    elif prop < not_arch_thresh:
        return([0,1,0, prop])
    else:
        return([0,0,1, prop])

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script to calculate stats for ArchIE. Prints stats to stdout.")
    parser.add_argument("-s", "--snp", action="store", required=True, help=".snp file with positions (eigenstrat format)")
    parser.add_argument("-i", "--ind", action="store", required=True, help=".ind file of IDS (eigenstrat format)")
    parser.add_argument("-a", "--admix", action="store", required=True, help=".geno file for admixed/test population (eigenstrat format)")
    parser.add_argument("-r", "--reference", action="store", required=False, help=".geno file for reference population (eigenstrat format).")
    parser.add_argument("-c", "--chrom", action="store", required=True, help="Chromosome")
    parser.add_argument("-b", "--begin", action="store", required=True, help="Beginning position")
    parser.add_argument("-e", "--end", action="store", required=True, help="End position")
    parser.add_argument("-w", "--window", action="store", required=True, help="Window length (eg. 50000)")
    parser.add_argument("-z", "--step", action="store", required=True, help="Step size for the window (eg. 10000)")
    parser.add_argument("-l", "--label", action="store_true", required=False, help="Enable flag if you are using simulated data and would like to label it as archaic/not")
    parser.add_argument("--anc", action="store", required=False, help=".anc file to output a label (needed with -l)")
    args = parser.parse_args()

    try:
        mutation_positions = parse_snp(args.snp)
    except FileNotFoundError:
        print("Error reading SNP file. Check file path.", file=sys.stderr)
        sys.exit(1)

    try:
        inds = parse_ind(args.ind)
    except FileNotFoundError:
        print("Error reading IND file. Check file path.", file=sys.stderr)
        sys.exit(1)

    try:
        genotypes = parse_geno(args.admix)
    except FileNotFoundError:
        print("Error reading admix .geno file. Check file path.", file=sys.stderr)
        sys.exit(1)

    try:
        ref_geno = parse_geno(args.reference)
    except FileNotFoundError:
        print("Error reading reference .geno file. Check file path.", file=sys.stderr)
        sys.exit(1)

    step = int(args.step)
    window = int(args.window)
    START = int(args.begin)
    END = int(args.end)
    n_sites = END - START
    n_samples = len(genotypes[0])

    if args.label:
        try:
            arch = parse_anc(args.anc, n_samples)
        except FileNotFoundError:
            print("Error reading anc file. Check file path.", file=sys.stderr)
            sys.exit(1)

    header=["Dist.Mean", "Dist.Var", "Dist.Skew", "Dist.Kurtosis", "Min.Dist", "Num.Private", "Haplotype.num", "Chrom", "Start", "End"]
    if args.label:
        header = header + ["Archaic", "Not.archaic", "In.between", "Proportion.archaic"]
    print(*header, sep="\t")

    starts = list(range(START, END, step))
    for start in starts:
        end = start+window
        start_idx = bisect_left(mutation_positions, start)
        end_idx = bisect_left(mutation_positions, end)
        curr_genotypes = genotypes[start_idx:end_idx+1]
        curr_ref_geno = ref_geno[start_idx:end_idx+1]
        curr_mutation_positions = mutation_positions[start_idx:end_idx+1]

        t_ref = list(map(list, zip(*curr_ref_geno)))
        t_geno =  list(map(list, zip(*curr_genotypes)))
        pos_to_remove = set() # contains indexes to remove
        haps = []
        for idx, hap in enumerate(t_ref):
            for jdx, site in enumerate(hap):
                if site == 1:
                    pos_to_remove.add(jdx)

        for idx, hap in enumerate(t_geno):
            haps.append([v for i, v in enumerate(hap) if i not in pos_to_remove])

        print(str(start)+"\t"+str(len(curr_genotypes[0]))+"\t"+str(len(curr_mutation_positions)), file=sys.stderr)

        #individual level

        for focal_idx in range(0, n_samples):
            dist = distance_vector(curr_genotypes, focal_idx)
            min_d = [min_dist_ref(curr_genotypes, curr_ref_geno, focal_idx)]
            n_priv = [num_private(np.array(haps), focal_idx)]

            output = dist + min_d + n_priv + [inds[focal_idx]] + [args.chrom] + [start] + [end] # need to put scalar elements into array to concatenate
            if args.label:
                focal_arch = [row[focal_idx] for row in arch ]
                lab = label(focal_arch, mutation_positions, n_sites, 0.7, 0.3)
                output = output + lab

            print(*output, sep="\t") # print stats to standard out
