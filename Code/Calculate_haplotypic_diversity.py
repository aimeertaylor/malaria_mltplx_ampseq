import allel
import numpy as np
import pandas as pd
import os
import time

# Input path to file directory:
os.chdir("/Volumes/seq_plasmodium/emilylav/Mali-R01/")
# Input the file names for each chromosome's vcf for the region of interest.
# vcf's should have no indels and have telomeric regions filtered out.
vcf_list = ["subtelomeric_no_indels_pf3k_Mali_Pf3D7_01_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_02_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_03_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_04_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_05_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_06_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_07_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_08_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_09_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_10_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_11_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_12_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_13_v3.recode.vcf",
            "subtelomeric_no_indels_pf3k_Mali_Pf3D7_14_v3.recode.vcf"]


def calculateHapDiv(region):
    """ Calculates haplotypic diversity for a given region of a vcf.
    Subjects with missing data for 1+ variant will be removed
    before calculation.

    Args:
      region: a haplotypeArray covering coordinates of interest.

    Returns:
      Haplotypic diversity based on subjects with complete data for every
      variant in the region. Returns 0 if only 1 subject present.
    """
    # remove missing in region
    keep_subject = np.ones(region.shape[1], dtype=int)
    for i in range(1, region.shape[1]):
        if -1 in region[:, i]:
            keep_subject[i] = 0
    region_complete = region.compress(condition=keep_subject, axis=1)
    # calculate haplotype frequencies
    freqs = region_complete.distinct_frequencies()
    n = region_complete.n_haplotypes
    homozygosity = sum(freqs**2)
    if n == 1:
        hap_div = 0
    else:
        hap_div = (1 - homozygosity) * (n/(n-1))
    return hap_div


for chrom in vcf_list:
    callset = allel.read_vcf(chrom)
    # create genotype array
    gt = allel.GenotypeArray(callset["calldata/GT"])
    # remove any het calls
    gt.mask = gt.is_het()
    gt_hom_only = gt.fill_masked(value=-1)
    gt_hap_array = gt_hom_only.haploidify_samples()
    H, windows, n_var = allel.windowed_statistic(pos=callset["variants/POS"],
                                                 values=gt_hap_array,
                                                 statistic=calculateHapDiv,
                                                 size=200, step=50,
                                                 start=1)
    df = pd.DataFrame(list(zip(H, windows, n_var)),
                      columns=["H", "window", "n_variants"])
    file_name = os.getcwd() + "/" + chrom[:-4] + "_hap_div.csv"
    df.to_csv(file_name, header=True)
