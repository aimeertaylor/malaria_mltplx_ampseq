import allel
import os
import numpy as np
import pandas as pd
import sys
import csv

# file output prefix:
prefix = sys.argv[1]
# path to vcf:
vcf_path = sys.argv[2]
# path to tabix for vcf:
tabix_path = sys.argv[3]
# amplicon coordinates file
coords_file = sys.argv[4]


def removeMissingStats(region, region_length):
    """ Calculates haplotype frequencies, haplotypic diversity, and nucleotide
    diversity for a given region of a vcf. Subjects with missing data for
    1+ variant will be removed before calculation.

    Args:
      region: a haplotypeArray covering coordinates of interest.

    Returns:
      List containing:
        - Number of individuals with no missing variant data in this region;
          only these individuals' data is used in further calculations.
        - Haplotypic diversity for region. Returns 0 if only 1 subject present.
        - Nucleotide diversity.
        - List of haplotype frequencies.
    """
    # remove missing in region
    keep_subject = np.ones(region.shape[1], dtype=int)
    for i in range(1, region.shape[1]):
        if -1 in region[:, i]:
            keep_subject[i] = 0
    region_complete = region.compress(condition=keep_subject, axis=1)
    # calculate haplotype frequencies
    freqs = region_complete.distinct_frequencies()
    nh = region_complete.n_haplotypes
    if nh == 1:
        hap_div = 0
    else:
        hap_div = allel.haplotype_diversity(region_complete)
    # calculate nucleotide diversity specifically on nonmissing region
    ac = region_complete.count_alleles()
    diffs = allel.mean_pairwise_difference(ac, fill=0)
    pi = np.sum(diffs)/region_length
    return [nh, hap_div, pi, freqs]


amplicons_stats = []
amplicons_freqs = []

with open(coords_file) as f:
    f.readline()
    for line in f:
        (chrom, start, stop, name) = line.strip().split("\t")
        coords = chrom + ":" + start + "-" + stop
        print(name + ": " + coords)
        amp_length = int(stop) - int(start)
        amp_vcf = allel.read_vcf(input=vcf_path,
                                 tabix=tabix_path,
                                 region=coords,
                                 fields="calldata/GT")
        if amp_vcf is None:
            print("no variants")
            nh = 0
            hap_div = 0
            pi = 0
            freqs = {1}
        else:
            gt = allel.GenotypeArray(amp_vcf["calldata/GT"])
            # remove any het calls, convert to haploid
            gt.mask = gt.is_het()
            gt_hom_only = gt.fill_masked(value=-1)
            gt_hap_array = gt_hom_only.haploidify_samples()
            # get stats on each amplicon and add to lists
            nh, hap_div, pi, freqs = removeMissingStats(gt_hap_array,
                                                        amp_length)
        amplicons_stats.append((name, nh, hap_div, pi))
        amplicons_freqs.append(freqs)
    # write stats to file
    file_stats = os.getcwd() + "/" + prefix + "_hap_div_stats.csv"
    with open(file_stats, 'w') as out:
        csv_out = csv.writer(out)
        csv_out.writerow(['name', 'nh', 'hd', 'pi'])
        csv_out.writerows(amplicons_stats)
    # write freqs to separate file
    file_freqs = os.getcwd() + "/" + prefix + "_hap_div_freqs.csv"
    with open(file_freqs, 'w') as out:
        csv_out = csv.writer(out)
        csv_out.writerows(amplicons_freqs)
