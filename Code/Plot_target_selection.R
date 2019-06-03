###################################################
# Previous crude experiments among markers all with cardinality > 2 
# suggest that approaches ranked in order of increasing error around estimates of IBD were:
#
# Maximal diversity 
# Sample at random
# Maximise distance
# Maximise log(distance) * effective cardinality
#
# The latter approach was based on the observation that log(distance) is roughly
# the same order of magnitude as effective cardinality and so both are given 
# roughly equal weighting. It maximises distances beteen ghost points, 
# not seen points
#
# Morever, it may be preferable to preference distance over cardinality? 
# If we used the log trick, we'd need to first check log(distances)
# are always on off the same magnitudes as effective cardinality
# 
# NEED TO FIND A WAY TO EVALUATE WHETHER ONE METHOD BETTER ANOTHEr
# COULD JUST SIMULATE. VARIATION IN APPROACH
###################################################
load('../RData/data_FG_S.RData') # Load data
source('./target_selection.R') # functions to select targets

data_$Med_Keff = apply(data_[c('Keff_FG', 'Keff_S')], 1, median) # Median eff. card.
data_$Med_div = 1-(1/data_$Med_Keff) # Med. diversity
data_$Med_div_01 = (data_$Med_div-min(data_$Med_div))/(max(data_$Med_div)-min(data_$Med_div)) # Normalise for colour res.

# Function to map number in zero one to colour
col_mapper <- function(num){
  x = colorRamp(rainbow(100, end = 0.8))(num)
  rgb(x[1], x[2], x[3], maxColorValue=255)
}

# Copy and pasted from plasmodb (there are better ways to do this)
chr_lengths <- c(640851, 947102, 1067971, 1200490, 1343557, 
                 1418242, 1445207, 1472805, 1541735, 1687656, 
                 2038340, 2271494, 2925236, 3291936)

# Static distance target strategies
panel_strategies = c('Sample at random', 
                     'Maximise effective cardinality', 
                     'Maximise log(distance) * effective cardinality', 
                     'Maximise distance')
m = 100 # Number of markers we'd like


#=============================================================
# Plot all possible targets
#=============================================================
par(mfrow = c(3,2))
finite_dt = data_$dt[is.finite(data_$dt)]
plot(y = data_$chrom, x = data_$pos, pch = 4, las = 1, bty = 'n', ylim = c(1,14), 
     yaxt = 'n', ylab = 'Chromosome', xlab = 'Position', main = 'All possible targets', 
     col = unlist(sapply(data_$Med_div_01, col_mapper)))
axis(side = 2, at = 1:14, las = 1, tick = F)
for(chr in 1:14){segments(y0 = chr, y1 = chr, x0 = 1, x1 = chr_lengths[chr])}
title(xlab = bquote("Stdev. finite"~italic(d[t]) == .(round(sd(finite_dt)))), 
      line = -2, cex = 0.75, adj = 1)
title(xlab = bquote("Med. finite"~italic(d[t]) == .(round(median(finite_dt)))), 
      line = -3, cex = 0.75, adj = 1)
title(xlab = bquote("Mean eff"~italic(K[t]) == .(round(mean(c(data_$Keff_FG, data_$Keff_S))))), 
      line = -4, cex = 0.75, adj = 1)

#=============================================================
# Plot positions for previous target selection functions 
# Variance doesn't capture spatial clustering
#=============================================================
for(panel_strategy in panel_strategies){
  
  Targets = target_selection_static_dt(data_, m, panel_strategy)
  
  # Recalculate distance
  Targets$dt <- c(diff(Targets$pos), Inf)
  pos_change_chrom <- 1 + which(diff(Targets$chrom) != 0) # find places where chromosome changes
  Targets$dt[pos_change_chrom-1] <- Inf
  finite_dt = Targets$dt[is.finite(Targets$dt)]
  
  print(finite_dt)
  
  plot(y = Targets$chrom, x = Targets$pos, pch = 20, las = 1, bty = 'n', ylim = c(1,14), 
       yaxt = 'n', ylab = 'Chromosome', xlab = 'Position', 
       col = unlist(sapply(Targets$Med_div_01, col_mapper)), 
       main = sprintf("Select %s: %s", m, panel_strategy))
  axis(side = 2, at = 1:14, las = 1, tick = F)
  for(chr in 1:14){
    segments(y0 = chr, y1 = chr, x0 = 1, x1 = chr_lengths[chr])
  }
  title(xlab = bquote("Stdev. finite"~italic(d[t]) == .(round(sd(finite_dt)))), 
        line = -2, cex = 0.75, adj = 1)
  title(xlab = bquote("Med. finite"~italic(d[t]) == .(round(median(finite_dt)))), 
        line = -3, cex = 0.75, adj = 1)
  title(xlab = bquote("Mean eff"~italic(K[t]) == .(round(mean(c(Targets$Keff_FG, Targets$Keff_S))))), 
        line = -4, cex = 0.75, adj = 1)
}



#=============================================================
# Plot positions dynamic dt function 
# Variance doesn't capture spatial clustering
#=============================================================
Targets = target_selection_dynmic_dt(data_, m)

# Recalculate distance
Targets$dt <- c(diff(Targets$pos), Inf)
pos_change_chrom <- 1 + which(diff(Targets$chrom) != 0) # find places where chromosome changes
Targets$dt[pos_change_chrom-1] <- Inf
finite_dt = Targets$dt[is.finite(Targets$dt)]

plot(y = Targets$chrom, x = Targets$pos, pch = 20, las = 1, bty = 'n', ylim = c(1,14), 
     yaxt = 'n', ylab = 'Chromosome', xlab = 'Position', 
     col = unlist(sapply(Targets$Med_div_01, col_mapper)), 
     main = "Matérn Type III + Strauss esque")

axis(side = 2, at = 1:14, las = 1, tick = F)
for(chr in 1:14){
  segments(y0 = chr, y1 = chr, x0 = 1, x1 = chr_lengths[chr])
}

title(xlab = bquote("Stdev. finite"~italic(d[t]) == .(round(sd(finite_dt)))), 
      line = -2, cex = 0.75, adj = 1)
title(xlab = bquote("Med. finite"~italic(d[t]) == .(round(median(finite_dt)))), 
      line = -3, cex = 0.75, adj = 1)
title(xlab = bquote("Mean eff"~italic(K[t]) == .(round(mean(c(Targets$Keff_FG, Targets$Keff_S))))), 
      line = -4, cex = 0.75, adj = 1)


