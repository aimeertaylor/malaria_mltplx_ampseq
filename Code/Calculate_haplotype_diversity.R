#!/usr/bin/env Rscript

suppressPackageStartupMessages({
library(pegas)
library(plyr)
library(dplyr)
})

# set up df with lengths and subtelomeric coordinates of chromosomes
chromosomes.df <- data.frame("chromosome" = c("Pf3D7_01_v3", "Pf3D7_02_v3", 
                                              "Pf3D7_03_v3", "Pf3D7_04_v3", 
                                              "Pf3D7_05_v3", "Pf3D7_06_v3", 
                                              "Pf3D7_07_v3", "Pf3D7_08_v3", 
                                              "Pf3D7_09_v3", "Pf3D7_10_v3", 
                                              "Pf3D7_11_v3", "Pf3D7_12_v3", 
                                              "Pf3D7_13_v3", "Pf3D7_14_v3"),
                             "start" = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                                         1, 1, 1),
                             "stop" = c(640851, 947102, 1067971, 1200490, 
                                        1343557, 1418242, 1445207, 1472805,
                                        1541735, 1687656, 2038340, 2271494, 
                                        2925236, 3291936),
                             "included.start" = c(126317, 120523, 105912, 
                                                  174828, 86819, 82608, 111472,
                                                  106592, 127415, 117241, 
                                                  140047, 98977, 129716, 71367),
                             "included.stop" = c(481383, 782133, 1001583, 
                                                 1066629, 1296200, 1296332, 
                                                 1319083, 1298329, 1373723, 
                                                 1514417, 1933151, 2126572, 
                                                 2807160, 3128120))

# set up dfs to add results of CalculateHapDivSlidingWindow to 
subtelomeric.Senegal.overall.df <- data.frame("chromosome" = character(), 
                                              "hap.div" = numeric(),
                                              "start" = numeric(),
                                              "stop" = numeric())

subtelomeric.FG.overall.df <- data.frame("chromosome" = character(), 
                                         "hap.div" = numeric(),
                                         "start" = numeric(),
                                         "stop" = numeric())

#####
# function definitions

CalculateHapDivOnRegion <- function(df, start, stop) {
  # Calculates haplotype diversity of samples in a data frame between positions 
  # start and stop, for samples with no missing data in region.
  #
  # Args:
  #  df: data frame containing all haplotypes and a POS column (position).
  #  start: beginning of region for calculation (coordinate, not row number).
  #  stop: end coordinate of region for calculation.
  #
  # Returns:
  #  List of the haplotype diversity of given haplotypes in df with positions 
  # between start and stop, as well as the number of non-missing samples. 
  temp <- as.data.frame(t(df %>% 
                            filter(POS >= start & POS < stop) %>% 
                            select(-CHROM, -POS, -ID, -REF, -ALT, -QUAL, -INFO, -FILTER, -FORMAT)))
  if (dim(temp)[1] == 0) {
    return(0)
  } else { 
    na.temp <- na.omit(temp)  # make sure that removing these samples is the way to go
    if (dim(temp)[1] == 0) {
      return(0)
    } else {
      test <- apply(na.temp, 1, paste, collapse = "")
      counts <- as.data.frame(table(test)) 
      n <- sum(counts$Freq)
      counts$proportion <- counts$Freq/n
      homozygosity <- sum(counts$proportion^2)
      hap_div <- (1 - homozygosity) * (n/(n-1))
      
      result <- list(hap_div, n)
      return(result) 
    }
  }
}

CalculateHapDivSlidingWindow <- function(df_overall, win.size, jump, start.pos,
                                         stop.pos) {
  # Calls CalculateHapDivOnRegion and applies on sliding window over a given region.
  # Args: 
  #  df_overall: data frame containing haplotypes and coordinates.
  #  win.size: size (in coordinate unit) of window.
  #  jump: how far to move between windows (in coordinate unit).
  #  start.pos: beginning of region to apply all sliding windows.
  #  stop.pos: end of region to apply all sliding windows.
  # Returns: 
  #  data frame with columns "start", "stop", "hap.div", and "n" per window
  df.to.return <- data.frame("start" = seq(start.pos, stop.pos, by = jump),
                             "stop" = 0,
                             "hap.div" = 0,
                             "n" = 0)
  df.to.return$stop <- df.to.return$start + win.size
  for (i in 1:dim(df.to.return)[1]) {
    temp.start <- df.to.return$start[i]
    temp.stop <- df.to.return$stop[i]
    result <- CalculateHapDivOnRegion(df_overall, temp.start, temp.stop)
    df.to.return$hap.div[i] <- result[[1]]
    df.to.return$n[i] <- result[[2]]
  }
  return(df.to.return)
}


GenerateHapDivDFPerRegion <- function(region, window.length, jump.length,
                                      df.to.append, 
                                      save.incrementally = FALSE) {
  # Wrapper function to generate haplotype diversity data per chromosome
  # per region, using CalculateHapDivSlidingWindow. 
  #
  # Args:
  #  region: character string with name of region (in file name).
  #  window.length: number of base pairs within sliding window
  #  jump.length: number of base pairs to move between sliding windows.
  #  df.to.append: data frame to add generated data to, with columns chromosome,
  #  hap.div, start, and stop. 
  #  save.incrementally: option to save data per chromosome/region individually,
  #  as an .RData file, instead of just at the end of all loops. 
  #
  # Returns:
  #  Data frame that was given as arg, with all generated data. 
  for (n in 1:14) {
    # set up internal variables
    file.name <- paste0("/seq/plasmodium/emilylav/subtelomeric/", 
                        region, "-", chromosomes.df$chromosome[n], ".vcf.gz")
    sw.name <- paste0(region, ".sw.", chromosomes.df$chromosome[n])
    
    # load data from VCF and filter to only SNPs 
    current.loci <- VCFloci(file.name)
    current.snps <- is.snp(current.loci)
    current.variants <- read.vcf(file.name, which.loci = which(current.snps))
    colnames(current.variants) <- 1:length(current.variants)
    current.loci.snps <- current.loci[current.snps,]
    
    # replace things that aren't ref/ref or alt/alt with NA
    # this assumes biallelic SNPs only! 
    for (i in 1:dim(current.loci.snps)[1]) {
      current.variants[[i]] <- as.character(current.variants[[i]])
      ref <- paste0(current.loci.snps$REF[i],"/", current.loci.snps$REF[i])
      alt <- paste0(current.loci.snps$ALT[i],"/", current.loci.snps$ALT[i])
      for (j in 1:dim(current.variants)[1]) {
        # go through each sample per variant in current.variants
        if (current.variants[[i]][j] == ref) {
          current.variants[[i]][j] <- "0"
        } else if (current.variants[[i]][j] == alt) {
          current.variants[[i]][j] <- "1"
        } else {
          current.variants[[i]][j] <- NA
        }
      }
    }
    
    # combine rewritten allele data with metadata in current.loci.snps
    current.all <- cbind(current.loci.snps, as.data.frame(t(current.variants)))
    
    # call helper functions to calculate haplotype diversity
    df <- CalculateHapDivSlidingWindow(df_overall = current.all, 
                                       win.size = window.length, 
                                       jump = jump.length, 
                                       start = chromosomes.df$start[n], 
                                       stop = chromosomes.df$stop[n])
    df$chromosome <- chromosomes.df$chromosome[n]
    
    # add new data to df with all chromosome data
    df.to.append <- rbind(df.to.append, df)
    
    if (save.incrementally) {
      data.name <- paste(region, chromosomes.df$chromosome[n], window.length, 
                        jump.length, ".RData", sep = "_")
      save(df, file = data.name)
    }
    
  }
  return(df.to.append)
}

#####
# call function per region

subtelomeric.Senegal.df <- GenerateHapDivDFPerRegion(region = "Senegal", 
                                                     window.length = 200,
                                                     jump.length = 50,
                                                     df.to.append = subtelomeric.Senegal.overall.df,
                                                     save.incrementally = TRUE)

save(subtelomeric.Senegal.df, file = "subtelomeric.Senegal.df.RData")

subtelomeric.FG.overall.df <- GenerateHapDivDFPerRegion(region = "FG", 
                                                        window.length = 200,
                                                        jump.length = 50,
                                                        df.to.append = subtelomeric.FG.overall.df,
                                                        save.incrementally = TRUE)

save(subtelomeric.FG.overall.df, file = "subtelomeric.FG.overall.df.RData")