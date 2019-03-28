#!/bin/bash

source /broad/software/scripts/useuse
reuse R-3.5

cd /seq/plasmodium/emilylav/subtelomeric/new

Rscript /seq/plasmodium/emilylav/subtelomeric/Calculate_haplotype_diversity.R
