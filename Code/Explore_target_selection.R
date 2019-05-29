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
load('~/Documents/BroadLaptop/TM_border/QuantLocalPfConnIBD/')
load('../RData/data_FG_S.RData') # Load data

data_$Med_Keff = apply(data_[c('Keff_FG', 'Keff_S')], 1, median) # Median eff. card.
data_$Med_div = 1-(1/data_$Med_Keff) # Med. diversity
data_$Med_div_01 = (data_$Med_div-min(data_$Med_div))/(max(data_$Med_div)-min(data_$Med_div)) # Normalise for colour res.

# Function to map number in zero one to colour
col_mapper <- function(num){
  x = colorRamp(rainbow(100, end = 0.8))(num)
  rgb(x[1], x[2], x[3], maxColorValue=255)
}


#=============================================================
# Plot all 
#=============================================================
# Copy and pasted from plasmodb (there are better ways to do this)
chr_lengths <- c(640851, 947102, 1067971, 1200490, 1343557, 
                 1418242, 1445207, 1472805, 1541735, 1687656, 
                 2038340, 2271494, 2925236, 3291936)

par(mfrow = c(1,1))
plot(y = data_$chrom, x = data_$pos, pch = 4, las = 1, bty = 'n', ylim = c(1,14), 
     yaxt = 'n', ylab = 'Chromosome', xlab = 'Position', main = 'All data', 
     col = unlist(sapply(data_$Med_div_01, col_mapper)))
axis(side = 2, at = 1:14, las = 1, tick = F)
for(chr in 1:14){segments(y0 = chr, y1 = chr, x0 = 1, x1 = chr_lengths[chr])}



#=============================================================
# Previous target selection functions 
#=============================================================
panel_strategies = c('Sample at random', 
                     'Maximise effective cardinality', 
                     'Maximise log(distance) * effective cardinality', 
                     'Maximise distance')

# The function for target selection, assuming data has effective cardinality 
target_selection = function(data_, m, panel_strategy){
  
  # Create distances
  data_$dt <- c(diff(data_$pos), Inf)
  pos_change_chrom <- 1 + which(diff(data_$chrom) != 0) # find places where chromosome changes
  data_$dt[pos_change_chrom-1] <- Inf
  
  # Check distance and effective cardinality are the same magnitude
  Med_Keff*log(data_$dt)
  
  # Order by distance, eff. cardinality, log(distances) * eff. cardinality
  Order_dt = sort.int(log(data_$dt), decreasing = T, index.return = T)$ix
  Order_Keff = sort.int(Med_Keff, decreasing = T, index.return = T)$ix
  Order_logdt_Keff = sort.int(Med_Keff*log(data_$dt), decreasing = T, index.return = T)$ix
  
  if(panel_strategy == "Maximise distance"){
    inds = Order_dt[1:m]} 
  if(panel_strategy == "Maximise effective cardinality"){
    inds = Order_Keff[1:m]}
  if(panel_strategy == "Maximise log(distance) * effective cardinality"){
    inds = Order_logdt_Keff[1:m]}
  if(panel_strategy == 'Sample at random'){
    inds = sample(Mmax, m, replace = F)}
  
  return(data_[sort(inds),])
}

#=============================================================
# Plot positions for previous target selection functions 
# Variance doesn't capture spatial clustering
#=============================================================
par(mfrow = c(3,2))
m = 100 # Number of markers we'd like
for(panel_strategy in panel_strategies){
  Targets = target_selection(data_, m, panel_strategy)
  
  # Recalculate distance
  Targets$dt <- c(diff(Targets$pos), Inf)
  pos_change_chrom <- 1 + which(diff(Targets$chrom) != 0) # find places where chromosome changes
  Targets$dt[pos_change_chrom-1] <- Inf
  finite_dt = Targets$dt[is.finite(Targets$dt)]
  
  print(finite_dt)
  
  plot(y = Targets$chrom, x = Targets$pos, pch = 20, las = 1, bty = 'n', ylim = c(1,14), 
       yaxt = 'n', ylab = 'Chromosome', xlab = 'Position', 
       col = unlist(sapply(Targets$Med_div_01, col_mapper)), 
       main = panel_strategy)
  axis(side = 2, at = 1:14, las = 1, tick = F)
  for(chr in 1:14){
    segments(y0 = chr, y1 = chr, x0 = 1, x1 = chr_lengths[chr])
  }
  title(xlab = bquote("Stdev. finite"~italic(d[t]) == .(round(sd(finite_dt)))), 
        line = -2, cex = 0.75, adj = 1)
  title(xlab = bquote("Med. finite"~italic(d[t]) == .(round(median(finite_dt)))), 
        line = -3, cex = 0.75, adj = 1)
  
}




#=============================================================
# New target selection function: can we do better? 
# http://www1.cmc.edu/pages/faculty/MHuber/Research/talks/huber_talk_2011h.pdf
#=============================================================
target_selection_new = function(data_, m, R = 1000){
  
  #========================================================
  # First step inspired by Matérn Type III process
  #========================================================
  # Sample with probability proportional to diversity
  # Each naturally has a birthday in [0, ∞) - its position
  ind_propose = sort(sample(Mmax, m, prob = Med_Keff, replace = F)) 
  
  # Run time forward, only allowing birth if not within R of older 
  ind_acceptd = ind_propose[1] 
  for(t in 2:m){
    ind_b4 = tail(ind_acceptd,1)
    data_t = data_[c(ind_b4,ind_propose[t]),] # subset data
    dt <- diff(data_t$pos) # calculate dt
    dt[diff(data_t$chrom) != 0] <- Inf # Make inf if chrom changed
    if(dt < R){
      next()
    } else {
      ind_acceptd = c(ind_acceptd, ind_propose[t])
    }
  }
  
  #========================================================
  # Second stepp inspired by Strauss process
  #========================================================
  m_remain = m - length(ind_acceptd)
  
  while(m_remain > 0){
    
    ind_remain = (1:nrow(data_))[-ind_acceptd]
    ind_propose = sort(sample(ind_remain, m_remain, prob = Med_Keff[ind_remain], replace = F)) 
    
    for(t in 1:m_remain){
      
      diffs = ind_acceptd - ind_propose[t]
      diff_sign_diffs = diff(sign(diffs))
      
      if(all(diff_sign_diffs == 0)){
        # If the all accepted are larger
        if(all(diffs > 0)){
          # subset data
          data_t = data_[c(ind_propose[t],ind_acceptd[1]),] 
          dt <- diff(data_t$pos) # Two distances 
          dt[diff(data_t$chrom) != 0] <- Inf 
          if(dt < R){ # If close to either of its neighbours, exclude
            next()} else {
              ind_acceptd = c(ind_propose[t], ind_acceptd)}
        } 
        
        # If the all accepted are smaller
        if(all(diffs < 0)){
          # subset data
          data_t = data_[c(ind_propose[t],ind_acceptd[1]),] 
          dt <- diff(data_t$pos) # Two distances 
          dt[diff(data_t$chrom) != 0] <- Inf 
          if(dt < R){ # If close to either of its neighbours, exclude
            next()} else {
              ind_acceptd = c(ind_acceptd, ind_propose[t])}
        } 
      } else {
        # Find closest neighbours in ind_accepted
        n1_ind = which(diff_sign_diffs != 0)
        n2_ind = n1_ind + 1
        
        # subset data
        data_t = data_[c(ind_acceptd[n1_ind],ind_propose[t],ind_acceptd[n2_ind]),] 
        dts <- diff(data_t$pos) # Two distances 
        dts[diff(data_t$chrom) != 0] <- Inf 
        if(any(dts < R)){ # If close to either of its neighbours, exclude
          next()} else {
            ind_acceptd = c(ind_acceptd, ind_propose[t])}
      }
      ind_acceptd = sort(ind_acceptd)
    }
    m_remain = m - length(ind_acceptd)
    print(m_remain)
  }
  return(data_[sort(ind_acceptd),])
}


#=============================================================
# Plot positions for previous target selection functions 
# Variance doesn't capture spatial clustering
#=============================================================
Targets = target_selection_new(data_, m)

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

