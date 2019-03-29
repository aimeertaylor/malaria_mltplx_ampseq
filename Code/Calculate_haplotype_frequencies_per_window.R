suppressPackageStartupMessages({
  library(pegas)
  library(plyr)
  library(dplyr)
})

setwd("/Volumes/seq_plasmodium/emilylav/subtelomeric")
load("overall.subtelomeric.RData")

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

# set up dfs to add results of CalculateHaplotypeFrequencies to 
hap.frequencies.Senegal <- data.frame("start" = as.numeric(),
                                      "stop" = as.numeric(), 
                                      "S.hap.div" = as.numeric(),
                                      "S.n" = as.numeric(),
                                      "chromosome" = as.character(),
                                      "FG.hap.div" = as.numeric(),
                                      "FG.n" = as.numeric(),
                                      "haplotypes" = as.character(),
                                      "hap.frequencies" = as.numeric())

hap.frequencies.FG <- data.frame("start" = as.numeric(),
                                 "stop" = as.numeric(), 
                                 "S.hap.div" = as.numeric(),
                                 "S.n" = as.numeric(),
                                 "chromosome" = as.character(),
                                 "FG.hap.div" = as.numeric(),
                                 "FG.n" = as.numeric(),
                                 "haplotypes" = as.character(),
                                 "hap.frequencies" = as.numeric())

CalculateHaplotypeFrequencies <- function(region, df.to.append, df.of.windows, 
                                      save.incrementally = FALSE) {
  # Wrapper function to calculate haplotype frequency data per chromosome
  # per region, for all samples with complete data over that window.  
  #
  # Args:
  #  region: character string with name of region (in file name).
  #  df.to.append: data frame to add generated data to, with columns chromosome,
  #  S.hap.div, FG.hap.div, S.n, FG.n, haplotypes, hap.frequencies, start, stop.
  #  df.of.windows: data frame with coordinates of windows over which to
  #  calculate the haplotypes and freqs. Columns start and stop needed.
  #  save.incrementally: option to save data per chromosome/region individually,
  #  as an .RData file, instead of just at the end of all loops. 
  #
  # Returns:
  #  Data frame that was given as arg, with all generated data. 
  for (n in 1:14) {
    # set up internal variables
    file.name <- paste0(#"/seq/plasmodium/emilylav/subtelomeric/", # add this if running on server!
                        region, "-", chromosomes.df$chromosome[n], ".vcf.gz")
    #sw.name <- paste0(region, ".sw.", chromosomes.df$chromosome[n])
    
    # load data and filter to only SNPs 
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
    
    # pull current chromosome window data from overall dataset
    temp.windows <- df.of.windows %>% 
      filter(chromosome == chromosomes.df$chromosome[n]) 
    
    for (k in 1:nrow(temp.windows)) {
      start <- temp.windows$start[k]
      stop <- temp.windows$stop[k]
      temp.per.window <- as.data.frame(t(current.all %>% 
                                           filter(POS >= start & POS < stop) %>% 
                                           select(-CHROM, -POS, -ID, -REF, -ALT, -QUAL, -INFO, -FILTER, -FORMAT)))
      # include only samples with complete data for that window
      na.temp.per.window <- na.omit(temp.per.window) 
      # create haplotypes as character strings of alleles
      haps <- apply(na.temp.per.window, 1, paste, collapse = "")
      # calculate frequencies of haplotypes
      counts <- as.data.frame(table(haps)) 
      n <- sum(counts$Freq)
      counts$proportion <- counts$Freq/n
      
      # add lists of haplotypes and their respective frequencies to df
      temp.windows$haplotypes[k] <- list(as.character(counts$haps))
      temp.windows$hap.freqs[k] <- list(counts$proportion)
    }
    df.to.append <- rbind(df.to.append, temp.windows) 
    
    if (save.incrementally) {
      data.name <- paste("hap_freqs_", region, chromosomes.df$chromosome[n], 
                         window.length, jump.length, ".RData", sep = "_")
      save(temp.windows, file = data.name)
    }
  }
}

#####
# call function per region
hap.frequencies.Senegal <- CalculateHaplotypeFrequencies(region = "Senegal", 
                                                         df.to.append = hap.frequencies.Senegal,
                                                         df.of.windows = subtelomeric.overall.df,
                                                         save.incrementally = TRUE)

save(hap.frequencies.Senegal, file = "hap.frequencies.Senegal.RData")

hap.frequencies.FG <- CalculateHaplotypeFrequencies(region = "FG", 
                                                    df.to.append = hap.frequencies.FG,
                                                    df.of.windows = subtelomeric.overall.df,
                                                    save.incrementally = TRUE)

save(hap.frequencies.FG, file = "hap.frequencies.FG.RData")
