# 2019-09 
# comparing antigens of interest with haplotypically diverse windows

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

setwd("~/Desktop/malaria_mltplx_ampseq/RData")
load("overall.subtelomeric.RData")

hap_div_S_FG <- subtelomeric.overall.df %>% 
  #select(-haplotypes, -hap.freqs) %>%
  unite("window", c("chromosome", "start", "stop"), sep = "_", remove = TRUE)

# format into BED file to intersect with gene IDs gff:
hap_div_windows_BED <- data.frame("chrom" = subtelomeric.overall.df$chromosome,
                           "chromStart" = subtelomeric.overall.df$start,
                           "chromEnd" = subtelomeric.overall.df$stop)
write_tsv(x = hap_div_windows_BED, path = "~/Desktop/hap_div_windows.BED")

# then run on command line: 
# bedtools intersect -a hap_div_windows.BED -b PlasmoDB-25_Pfalciparum3D7.GENES.gff -wo > hap_div_windows_in_genes.tsv

# read back in the intersection of windows in coding regions
windows_intersection <- read_delim("~/Desktop/hap_div_windows_in_genes.tsv", 
                                   "\t", escape_double = FALSE, col_names = FALSE, 
                                   trim_ws = TRUE)

# tidy data to get a df with all the info we want per window, including gene info
windows_intersection <- windows_intersection %>% 
  select(-X4, -X5, -X6, -X9, -X11, -X12, -X14, -X16, -X17, -X18, -X19) %>%
  separate(X13, c(NA, "gene_ID"), sep = "=") %>%
  separate(X15, c(NA, "size"), sep = "=") %>%
  unite("window", X1:X3, sep = "_", remove = TRUE)

colnames(windows_intersection) <- c("window", "gene_start", "gene_end", "gene_strand", 
                                    "gene_ID", "gene_size", "overlap_size")

windows_intersection <- join(windows_intersection, hap_div_S_FG)

# next, make a df with the top hap.div (HD) scoring window per gene in Senegal
max_HD_Senegal_per_gene <- windows_intersection %>% 
  group_by(gene_ID) %>% 
  slice(which.max(S.hap.div)) %>%
  select(-FG.hap.div, -FG.n)

colnames(max_HD_Senegal_per_gene) <- c("S_window", "gene_start", 
                                       "gene_end", "gene_strand", "gene_ID", 
                                       "gene_size", "S_overlap_size", "S_hap_div", 
                                       "S_n")
# and in FG data
max_HD_FG_per_gene <- windows_intersection %>% 
  group_by(gene_ID) %>% 
  slice(which.max(FG.hap.div)) %>%
  select(-S.hap.div, -S.n)

colnames(max_HD_FG_per_gene) <- c("FG_window", "gene_start", "gene_end", 
                                  "gene_strand", "gene_ID", "gene_size", 
                                  "FG_overlap_size", "FG_hap_div", "FG_n")

# now can input a list of gene IDs and pull the top windows per region per gene
# this is just a working list of gene IDs, can always change later
genes <- c("PF3D7_0423800", "PF3D7_0731500", "PF3D7_1201300", "PF3D7_0102700", 
           "PF3D7_0424100", "PF3D7_0902800", "PF3D7_1216600", "PF3D7_1133400", 
           "PF3D7_1335900", "PF3D7_0304600", "PF3D7_1335100", "PF3D7_0930300", 
           "PF3D7_0902900", "PF3D7_0508000", "PF3D7_0207000", "PF3D7_0207700", 
           "PF3D7_1016900", "PF3D7_1449000", "PF3D7_0214900", "PF3D7_0207900", 
           "PF3D7_0620400", "PF3D7_1222600", "PF3D7_0511400", "PF3D7_1035400", 
           "PF3D7_1136200", "PF3D7_0721700", "PF3D7_0828800", "PF3D7_0404900", 
           "PF3D7_1012200", "PF3D7_0905400", "PF3D7_1321900", "PF3D7_0105400.1", 
           "PF3D7_1334400", "PF3D7_0206800", "PF3D7_0721100", "PF3D7_0707300", 
           "PF3D7_0616500", "PF3D7_1342500", "PF3D7_0624400", "PF3D7_0422100", 
           "PF3D7_0207600", "PF3D7_0206900.1", "PF3D7_1242000", "PF3D7_1463900", 
           "PF3D7_0212600", "PF3D7_0108700", "PF3D7_1028700", "PF3D7_0623000", 
           "PF3D7_0612700", "PF3D7_1033200", "PF3D7_0317100", "PF3D7_1352500", 
           "PF3D7_1116700", "PF3D7_0620000", "PF3D7_1036000", "PF3D7_0207800", 
           "PF3D7_0713700", "PF3D7_0526900", "PF3D7_0513700", "PF3D7_0104000", 
           "PF3D7_1351800.2", "PF3D7_1310500", "PF3D7_1017100", "PF3D7_0405900", 
           "PF3D7_0702900", "PF3D7_0406200", "PF3D7_1420700", "PF3D7_0606800", 
           "PF3D7_1334600", "PF3D7_1334300", "PF3D7_0103900", "PF3D7_1404900", 
           "PF3D7_0419700", "PF3D7_1035800", "PF3D7_1222300", "PF3D7_0511600", 
           "PF3D7_1115600", "PF3D7_0729200", "PF3D7_0831800", "PF3D7_1372200", 
           "PF3D7_1343700", "PF3D7_0709000", "PF3D7_0417200", "PF3D7_0810800",
           "PF3D7_0523000")

genes_df <- data.frame("gene_ID" = genes)

# pull the top window per region for each gene and write csv
genes_df <- join(genes_df, max_HD_Senegal_per_gene) %>%
  join(max_HD_FG_per_gene) %>%
  select(gene_ID, gene_start, gene_end, gene_strand, gene_size, S_window, S_hap_div, 
         S_overlap_size, S_n, FG_window, FG_hap_div, FG_overlap_size, FG_n)
  # re-orders columns to group by gene info, S info, FG info

write.csv(genes_df, file = "~/Desktop/genes_hap_div.csv")

#####
# pulling top windows in both regions to see what genes they're in
load("top.windows.subtelomeric.RData")

# format into BED file to intersect with gene IDs gff:
top_windows_BED <- data.frame("chrom" = top.subtelomeric.overall.df$chromosome,
                                  "chromStart" = top.subtelomeric.overall.df$start,
                                  "chromEnd" = top.subtelomeric.overall.df$stop)
write_tsv(x = top_windows_BED, path = "~/Desktop/top_windows.BED")

# then run on command line: 
# bedtools intersect -a top_windows.BED -b PlasmoDB-25_Pfalciparum3D7.GENES.gff -wo > hap_div_windows_in_genes.tsv

# read back in the intersection of windows in coding regions
top_windows_intersection <- read_delim("~/Desktop/top_windows_in_genes.tsv", 
                                   "\t", escape_double = FALSE, col_names = FALSE, 
                                   trim_ws = TRUE)

# tidy data to get a df with all the info we want per window, including gene info
top_windows_intersection <- top_windows_intersection %>% 
  select(-X4, -X5, -X6, -X9, -X11, -X12, -X14, -X16, -X17, -X18, -X19) %>%
  separate(X13, c(NA, "gene_ID"), sep = "=") %>%
  separate(X15, c(NA, "size"), sep = "=") %>%
  unite("window", X1:X3, sep = "_", remove = TRUE)

colnames(top_windows_intersection) <- c("window", "gene_start", "gene_end", "gene_strand", 
                                    "gene_ID", "gene_size", "overlap_size")

# format intial top.windows input to have "window" field
top_windows_S_FG <- top.subtelomeric.overall.df %>% 
  #select(-haplotypes, -hap.freqs) %>%
  unite("window", c("chromosome", "start", "stop"), sep = "_", remove = TRUE)

top_windows_intersection <- join(top_windows_intersection, top_windows_S_FG)

write.csv(top_windows_intersection, file = "~/Desktop/top_windows_intersection.csv")

#####
# want a list of all genes so that I can filter to ones with any a highest window 
# score in FG and S that is not necessarily the same exact window
genes_all <- data.frame("gene_ID" = max_HD_FG_per_gene$gene_ID)
max_HD_per_region_per_gene <- join(genes_all, max_HD_Senegal_per_gene,) %>%
  join(max_HD_FG_per_gene)
                                  
top_genes_max_HD <- max_HD_per_region_per_gene %>%
  filter(FG_hap_div > 0.5, S_hap_div > 0.5)
  #separate(FG_window, c("chromosome", "FG_start", "FG_stop"), sep = "_", remove = FALSE) %>%
  #select(-FG_start, -FG_stop)
top_genes_max_HD$chromosome <- substr(top_genes_max_HD$S_window, 1, 11)

# plot jitter comparison of subtelomeric data for FG vs Senegal
ggplot(top_genes_max_HD, aes(x = FG_hap_div, y = S_hap_div, color = chromosome)) +
  geom_jitter(alpha = 0.5, size = 1) +
  facet_wrap(~ chromosome) +
  xlim(0, 1) + 
  ylim(0, 1)

#ggsave(g1, filename = "~/Desktop/hap.div.jitter.all.chromosomes.png", width = 10, height = 10)

# plot 2d bin comparison of subtelomeric data for FG vs Senegal
ggplot(top_genes_max_HD, aes(x = FG_hap_div, y = S_hap_div)) +
  stat_bin2d(bins = 25) +
  scale_fill_gradientn(colors = r, trans = "log") +
  facet_wrap( ~ chromosome)

#ggsave(g2, filename = "~/Desktop/hap.div.2dbin.all.chromosomes.png", width = 10, height = 10)

ggplot(top_genes_max_HD, aes(x = chromosome)) +
  geom_histogram(stat = "count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "number of genes")

### pulling top window per gene if window > 0.7
top_Senegal_windows <- max_HD_Senegal_per_gene %>%
  filter(S_hap_div > 0.7)

write.csv(top_Senegal_windows, file = "~/Desktop/top_Senegal_windows.csv")

# pulling top window per gene from Helb et al., 2015
Helb_et_al_DataS1_600antigens <- read_delim("~/Desktop/Helb_et_al_DataS1_600antigens.bed", 
                                            "\t", escape_double = FALSE, trim_ws = TRUE)


genes_df <- data.frame("gene_ID" = Helb_et_al_DataS1_600antigens$gene)

# pull the top window per region for each gene and write csv
genes_df <- join(genes_df, max_HD_Senegal_per_gene)
  select(gene_ID, gene_start, gene_end, gene_strand, gene_size, S_window, S_hap_div, 
         S_overlap_size, S_n, FG_window, FG_hap_div, FG_overlap_size, FG_n)
# re-orders columns to group by gene info, S info, FG info

write.csv(genes_df, file = "~/Desktop/genes_Helb_hap_div.csv")
