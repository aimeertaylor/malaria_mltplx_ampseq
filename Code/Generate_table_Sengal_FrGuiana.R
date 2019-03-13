###################################################################
# This script simulates data using diversity estimates from 
# French Guiana and Senegal
#
# To-do: 
# - calculate example CIs
###################################################################

## Set up
rm(list = ls())
set.seed(1)
library(ggplot2)
library(dplyr)
library(Rcpp)
library(doParallel)
library(doRNG)
library(MCMCpack) # For dirichlet
source("~/Dropbox/IBD_IBS/PlasmodiumRelatedness/Code/simulate_data.R")
sourceCpp("~/Dropbox/IBD_IBS/PlasmodiumRelatedness/Code/hmmloglikelihood.cpp")
registerDoParallel(cores = detectCores()-2)
epsilon <- 0.001 # Fix epsilon throughout
rs <- seq(from = 0.01, to = 0.99, length.out = 5)
ks <- c(1,5,10,50)
samplesizes = c(96*c(0.25, 1:3)) # add maximum sample size when working with real data
set.seed(1) # for reproducibility
nrepeats <- 100 # number of repeats 
kfixed <- 12 
rfixed <- 0.5 
panel_strategies = c('Sample at random', 
                     'Maximise effective cardinality', 
                     'Maximise log(distance) * effective cardinality', 
                     'Maximise distance')
countries = c('French Guiana', 'Senegal')
RUN = T

## Mechanism to generate Ys given fs, distances, k, r, epsilon
simulate_Ys_hmm <- function(frequencies, distances, k, r, epsilon){
  Ys <- simulate_data(frequencies, distances, k = k, r = r, epsilon, rho = 7.4 * 10^(-7))
  return(Ys)
}

## Mechanism to compute MLE given fs, distances, Ys, epsilon
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon){
  ndata <- nrow(frequencies)
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optim(par = c(kfixed, 0.5), fn = function(x) - ll(x[1], x[2]))
  rhat <- optimization$par
  return(rhat)
}



#====================================================
# Data 
#====================================================
# Format diversity data from French Guiana and Senegal datasets
data_ = read.csv(file = '~/Dropbox/IBD_IBS/Original_data/high.hap.div.csv', as.is = T) # This is old data without frequencies
load('~/Documents/BroadLaptop/AmpSeq/RData/FromEmily/top.windows.hap.frequencies.FG.RData', )
load('~/Documents/BroadLaptop/AmpSeq/RData/FromEmily/top.windows.hap.frequencies.Senegal.RData')

# Check identity of first 5 columns before creating data_
if(identical(top.windows.hap.frequencies.FG[,1:5], top.windows.hap.frequencies.Senegal[,1:5])){
  data_ = top.windows.hap.frequencies.FG[,1:5] 
  i <- sapply(data_, is.factor) 
  data_[i] <- lapply(data_[i], as.character) # Convert factors to characters
} else {
  stop('Common data are not identical as expected')
}
Mmax <- nrow(data_) # Total number of markers

# Extract chromosome
data_$chrom = as.numeric(do.call(rbind, strsplit(data_$chromosome, '_'))[,2]) 
all(data_$chrom == cummax(data_$chrom)) # Check order

# Extract position as median
data_$pos = apply(data_[,c('start','stop')] , 1, function(x){median(as.numeric(x))}) 
reordered_data_list = list()
for(chr in unique(data_$chrom)){
  x = data_$pos[chr == data_$chrom] # Extract positions per chromosome 
  print(all(x == cummax(x))) # If not all True, not all are monotonically increasing
  inds = sort.int(x, index.return = T)$ix # Sort positions
  reordered_data_list[[chr]] = data_[chr == data_$chrom, ][inds, ]
}

# Check now all in order: yes
data_ = do.call(rbind, reordered_data_list)
all(sapply(unique(data_$chrom), function(chr){
  x = data_$pos[chr == data_$chrom] 
  all(x == cummax(x)) 
}))


# Create distances
data_$dt <- c(diff(data_$pos), Inf)
pos_change_chrom <- 1 + which(diff(data_$chrom) != 0) # find places where chromosome changes
data_$dt[pos_change_chrom-1] <- Inf

# Extract frequencies 
MaxHapCountFG = max(sapply(top.windows.hap.frequencies.FG$hap.freqs, length))
MaxHapCountS = max(sapply(top.windows.hap.frequencies.Senegal$hap.freqs, length))
frequencies_S <- array(0, dim = c(Mmax, MaxHapCountS))
frequencies_FG <- array(0, dim = c(Mmax, MaxHapCountFG))
for(i in 1:Mmax){
  FreqFG_i = top.windows.hap.frequencies.FG$hap.freqs[[i]]
  FreqS_i = top.windows.hap.frequencies.Senegal$hap.freqs[[i]]
  frequencies_FG[i,1:length(FreqFG_i)] <- FreqFG_i
  frequencies_S[i,1:length(FreqS_i)] <- FreqS_i
}

#++++++++++++++++++++++++++++
# Diversities appear to be off = raise in meeting 
# Due to the sample-estimate correction - think about
myDiv_S = 1-rowSums(frequencies_S^2)
myDiv_FG = 1-rowSums(frequencies_FG^2)
par(mfrow = 2:1, pty = 'm')
plot(myDiv_S,data_$S.hap.div, pch = 20, col = adjustcolor(1,alpha.f = 0.5), 
     xlab = 'Calculated Aimee', ylab = 'Emily Provided')
abline(a = 0, b = 1, col = 'blue')
plot(myDiv_FG,data_$FG.hap.div, pch = 20, col = adjustcolor(1,alpha.f = 0.5),
     xlab = 'Calculated Aimee', ylab = 'Emily Provided')
abline(a = 0, b = 1, col = 'blue')
#++++++++++++++++++++++++++++

# Calculate effective cardinalities frequencies 
data_$Keff_FG = 1/(1-myDiv_FG)
hist(data_$Keff_FG); range(data_$Keff_FG)
data_$Keff_S = 1/(1-myDiv_S)
hist(data_$Keff_S); range(data_$Keff_S)

# Create an ordering of makers for panel strategy
# Note that if we select top m diverse regions, there is a lot overlap
# due to correlation between high Keff and distance. 
Med_Keff = apply(data_[c('Keff_FG', 'Keff_S')], 1, median)
par(mfrow = c(1,2))
plot(Med_Keff, data_$dt)
plot(Med_Keff, log(data_$dt))
# log(distances is of the same order of magnitude as Keff)

# There are many ways to make this more nuanced
Order_dt = sort.int(log(data_$dt), decreasing = T, index.return = T)$ix
Order_Keff = sort.int(Med_Keff, decreasing = T, index.return = T)$ix
Order_dt_Keff = sort.int(Med_Keff*log(data_$dt), decreasing = T, index.return = T)$ix

freqs_pos_inds = function(m, panel_strategy){
  if(panel_strategy == "Maximise distance"){
    inds = Order_dt[1:m]} 
  if(panel_strategy == "Maximise effective cardinality"){
    inds = Order_Keff[1:m]}
  if(panel_strategy == "Maximise log(distance) * effective cardinality"){
    inds = Order_dt_Keff[1:m]}
  if(panel_strategy == 'Sample at random'){
    inds = sample(Mmax, m, replace = F)
  }
  return(sort(inds))
}

# What do the average Keffs look like? 
sapply(samplesizes, function(m){
  inds = freqs_pos_inds(m,  panel_strategies[1])
  colMeans(data_[inds,c('Keff_FG', 'Keff_S', 'dt')])
})

sapply(samplesizes, function(m){
  inds = freqs_pos_inds(m,  panel_strategies[2])
  colMeans(data_[inds,c('Keff_FG', 'Keff_S', 'dt')])
})

sapply(samplesizes, function(m){
  inds = freqs_pos_inds(m,  panel_strategies[3])
  colMeans(data_[inds,c('Keff_FG', 'Keff_S', 'dt')])
})

# What do the distances look like? 
par(mfrow  = c(2,2))
for(m in samplesizes){
  inds_unordered = freqs_pos_inds(m, panel_strategies[1])
  inds_ordered = freqs_pos_inds(m, panel_strategies[2])
  plot(data_$dt[inds_unordered], 
       data_$dt[inds_ordered], 
       log = 'xy')
}

#=====================================
# Calculate RMSE 
#=====================================
tables_many_repeats_m <- array(dim = c(length(samplesizes), length(rs), length(panel_strategies), length(countries)),
                               dimnames = list(samplesizes, rs, panel_strategies, countries))
inds_store <- list() # Store inds  

if(RUN){
  system.time(
    for(panel_strategy in panel_strategies){
      for (m in samplesizes){
        
        inds = freqs_pos_inds(m, panel_strategy)
        inds_store[[panel_strategy]][[as.character(m)]] = inds
        freq_FG <- frequencies_FG[inds, ]
        freq_S <- frequencies_S[inds, ]
        
        # Recalculate distances
        distances <- c(diff(data_$pos[inds]), Inf)
        pos_change_chrom <- 1 + which(diff(data_$chrom[inds]) != 0) # find places where chromosome changes
        distances[pos_change_chrom-1] <- Inf
        
        for (r in rs){
          rhats_FG <- foreach (irepeat = 1:nrepeats, .combine = rbind) %dorng% {
            Ys <- simulate_Ys_hmm(freq_FG, distances, kfixed, r, epsilon)
            c(compute_rhat_hmm(freq_FG, distances, Ys, epsilon))
          }
          
          rhats_S <- foreach (irepeat = 1:nrepeats, .combine = rbind) %dorng% {
            Ys <- simulate_Ys_hmm(freq_S, distances, kfixed, r, epsilon)
            c(compute_rhat_hmm(freq_S, distances, Ys, epsilon))
          }
          
          tables_many_repeats_m[as.character(m), as.character(r), panel_strategy, 'French Guiana'] <- sqrt(mean((rhats_FG[,2] - r)^2))
          tables_many_repeats_m[as.character(m), as.character(r), panel_strategy, 'Senegal'] <- sqrt(mean((rhats_S[,2] - r)^2))
        }
      }
    }
  )
  save(tables_many_repeats_m, inds_store, file = '../RData/Tables_rs_FG_S.RData')
}











