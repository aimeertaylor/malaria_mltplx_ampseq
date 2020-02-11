suppressPackageStartupMessages({
  library(ggplot2)
  library(plyr)
  library(dplyr)
  library(RColorBrewer)
})

rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

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

#####
setwd("~/Desktop/malaria_mltplx_ampseq/RData/haplotype_frequencies")
load("hap.frequencies.FG.RData")
load("hap.frequencies.Senegal.RData")

# 2019-09-26
# limit windows to only those with >2 haplotypes to plot without 0.5 pileup
for (i in 1:nrow(hap.frequencies.FG)) {
  hap.frequencies.FG$FG.n.haps[[i]] <- length(hap.frequencies.FG$haplotypes[[i]])
}

for (j in 1:nrow(hap.frequencies.Senegal)) {
  hap.frequencies.Senegal$S.n.haps[[j]] <- length(hap.frequencies.Senegal$haplotypes[[j]])
}

over_2_haps_FG <- hap.frequencies.FG %>%
  filter(FG.n.haps > 2, FG.n >= 10) %>%
  select(-haplotypes, -hap.freqs)

over_2_haps_S <- hap.frequencies.Senegal %>%
  filter(S.n.haps > 2, S.n >= 10) %>%
  select(-haplotypes, -hap.freqs)

over_2_haps_all <- plyr::join(over_2_haps_FG, over_2_haps_S, type = "inner")

# plot jitter comparison of subtelomeric data for FG vs Senegal
g1 <- ggplot(over_2_haps_all, aes(x = FG.hap.div, y = S.hap.div, color = chromosome)) +
  geom_jitter(alpha = 0.5, size = 1) +
  facet_wrap(~ chromosome)

ggsave(g1, filename = "~/Desktop/hap.div.jitter.all.chromosomes.png", width = 10, height = 10)

# plot 2d bin comparison of subtelomeric data for FG vs Senegal
g2 <- ggplot(over_2_haps_all, aes(x = FG.hap.div, y = S.hap.div)) +
  stat_bin2d(bins = 25) +
  scale_fill_gradientn(colors = r, trans = "log") +
  facet_wrap( ~ chromosome)

ggsave(g2, filename = "~/Desktop/hap.div.2dbin.all.chromosomes.png", width = 10, height = 10)

# plot hap.div scores along chromosome coordinates - hard to see in this combined file
g3 <- ggplot(no.hets.comparison.tidy, aes(x = start, y = hap.div, color = region)) +
  geom_point(alpha = 0.3) + 
  facet_wrap( ~ chromosome)

ggsave(g3, filename = "hap.div.by.coordinate.all.chromosomes.png", width = 10, height = 10)

# generate a single plot per chromosome
for (i in 1:14) {
  temp.tidy <- subtelomeric.overall.df.tidy %>% 
    filter(chromosome == chromosomes.df$chromosome[i], metric == "S.hap.div" | metric == "FG.hap.div")
  
  g <- ggplot(temp.tidy, aes(x = start, y = value, color = metric)) +
    geom_point(alpha = 0.2) +
    xlim(chromosomes.df$start[i], chromosomes.df$stop[i])
  
  temp.filename <- paste0("hap.div.by.coordinate.", chromosomes.df$chromosome[i], ".png")
  ggsave(g, filename = temp.filename, width = 8, height = 5)
}

