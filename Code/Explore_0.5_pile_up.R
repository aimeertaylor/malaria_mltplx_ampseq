#################################################################
# In this script I'm trying to understand the reason for the pile
# up of diversity = 0.5 in the FG heterozygone

# To-do
# Expore the effect of any MAF cut offs
# Consider filter or CIs around diversity to account for n-artifacts
# Consider effect of overlapping windows
# Re generate plots of diversity against K and SNPno for all windows
# Simulate repulsive point process
#################################################################

#================================================================
# Exploring diversity vs n (no. of samples with data per window)
# Focus on FG
#================================================================
rm(list = ls())
load('~/Documents/malaria_mltplx_ampseq/RData/haplotype_diversity/subtelomeric.FG.overall.df.RData')

# Chrom numeric
chrom = as.numeric(gsub('Pf3D7_', '', gsub('_v3', '', subtelomeric.FG.overall.df$chromosome)))

# Plot of n against diversity
par(mfrow = c(4,4), mar = c(2,2,1,1))
for(c in 1:14){
  ind = which(chrom == c)
  plot(x = subtelomeric.FG.overall.df$n[ind], 
       y = subtelomeric.FG.overall.df$hap.div[ind],
       ylim = c(0,1), xlim = range(subtelomeric.FG.overall.df$n), 
       col = adjustcolor('black', alpha.f = 0.25), pch = 20, 
       xaxt = 'n', yaxt = 'n', main = sprintf("Chromosome %s", c))
  title(ylab = 'Diversity (0 to 1)', xlab = 'n (0 to 108)', line = 0)
}
n_count = table(subtelomeric.FG.overall.df$n)

# Visualise pile up at 0.5 when n > 100 and not
# Does not appear to be linked to low n
par(mfrow = c(4,4), mar = c(1,1,1,1))
for(c in 1:14){
  ind = which(chrom == c) 
  counts = table(round(subtelomeric.FG.overall.df$hap.div[ind], 1))
  maxY = sort(counts, decreasing = T)[2]
  plot(x = as.numeric(names(counts)), y = counts, type = 'l', 
       ylim = c(0,maxY), xlim = c(0,1), xaxt = 'n', yaxt = 'n')
  title(xlab = 'Diversity (0 to 1)', ylab = 'Counts (relative)', line = 0)
  
  ind = which((chrom == c) & (subtelomeric.FG.overall.df$n > 100))
  counts = table(round(subtelomeric.FG.overall.df$hap.div[ind], 1))
  lines(x = as.numeric(names(counts)), y = counts, type = 'l', col = 'blue')
  
  abline(v = 0.5, lty = 'dotted')
  mtext(text = sprintf("Chromosome %s", c), side = 4, line = -1, cex = 0.75)
}

# Zoom in on 0.5s: 
ind = which(subtelomeric.FG.overall.df$hap.div > 0.49 & subtelomeric.FG.overall.df$hap.div < 0.51)
sum(ind, na.rm = T)
subtelomeric.FG.overall.df[ind,]


#================================================================
# Exploring diversity as a function of SNP count and cardinality
# For both data sets cardinality << 2^SNPno.
#================================================================
rm(list = ls())
load('../RData/FromEmily/top.windows.hap.frequencies.FG.RData', )
load('../RData/FromEmily/top.windows.hap.frequencies.Senegal.RData')
X = list('French Guiana' = top.windows.hap.frequencies.FG, 'Senegal' = top.windows.hap.frequencies.Senegal)
Z = lapply(X, function(z){
  z$K = apply(z, 1, function(x)length(x$haplotypes))
  z$SNPno. = apply(z, 1, function(x){length(strsplit(x$haplotypes[[1]][1],split='')[[1]])})
  return(z)
})
 
par(mfrow = c(2,3))
MaxK = max(sapply(Z, function(x)max(x$K)))
MaxSNPno =  max(sapply(Z, function(x)max(x$SNPno.)))
for(i in 1:2){
  
  ind = ifelse(i == 1, 5, 3) # Extract diversity column
    
  # Plot diversity against SNP count
  plot(x = Z[[i]]$SNPno., y = Z[[i]][,ind], pch = 20, col = adjustcolor('black', alpha.f = 0.2),
       panel.first = grid(), ylab = "Haplotype diversity", xlab = 'SNP count per hapotype', 
       main = names(X)[i], las = 1, xlim = c(1,MaxSNPno))
  
  # Plot diversity against cardinality
  plot(x = Z[[i]]$K, y = Z[[i]][,ind], pch = 20, col = adjustcolor('black', alpha.f = 0.2),
       panel.first = grid(), ylab = "Haplotype diversity", xlab = 'Cardinality', 
       main = names(X)[i], las = 1, xlim = c(1,MaxK))
  
  # Plot cardinalty against SNP count
  plot(x = Z[[i]]$SNPno., y = Z[[i]]$K, pch = 20, col = adjustcolor('black', alpha.f = 0.2),
       panel.first = grid(), ylab = "Cardinality", xlab = 'SNP count per hapotype', 
       main = names(X)[i], las = 1, ylim = c(1,MaxK), xlim = c(1,MaxSNPno))
}

# Many of the K = 2 have high MAF (understandable given > 0.5 diversity only)  
Kcounts = sapply(Z, function(x){z = table(x$K); z/sum(z)}) # Twice
par(mfrow = c(1,2))
sapply(Kcounts, barplot) # plots of window cardinaloty

