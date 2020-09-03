import allel
import os
import numpy as np
import pandas as pd
import sys

# file output prefix:
prefix = sys.argv[1]
# path to vcf:
vcf_path = sys.argv[2]


def removeMissingStats(region):
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
    pi = np.sum(diffs)/200
    return [nh, hap_div, pi, freqs]


chromosomes = allel.read_vcf(vcf_path, fields=['CHROM'])
chromosomes_list = np.unique(chromosomes['variants/CHROM'])
for chrom in chromosomes_list:
    print(chrom)
    # read in that chromosome data only
    callset = allel.read_vcf(vcf_path, region=chrom, fields='*')
    gt = allel.GenotypeArray(callset["calldata/GT"])
    # remove any het calls, convert to haploid
    gt.mask = gt.is_het()
    gt_hom_only = gt.fill_masked(value=-1)
    gt_hap_array = gt_hom_only.haploidify_samples()
    # remove individuals with missing data and calculate stats
    n_list, w, n = allel.windowed_statistic(pos=callset["variants/POS"],
                                            values=gt_hap_array,
                                            statistic=removeMissingStats,
                                            size=200,
                                            step=50,
                                            start=1)
    df = pd.DataFrame(list(zip(n_list, w, n)),
                      columns=["n_list", "windows_n", "n_var_n"])
    file_name = os.getcwd() + "/" + prefix + chrom + "_hap_div.csv"
    df.to_csv(file_name, header=True)
