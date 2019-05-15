###########################################################
# This script explores the impact of clustering markers on 
# a single chromosome on the error in estimates of r and k 
# 
# Estimates of r appear to get progressively worse
# Simularly for estimates of k, except under complete 
# independence (i.e. if it were possible to have say 100
# markers on 100 chromosomes) or complete clustering (i.e.
# all markers on a single chromosome)
###########################################################

## Set up
rm(list = ls())
set.seed(1)
library(Rcpp)
library(doParallel)
library(doRNG)
require(kableExtra) # to print as table
source("./simulate_data_and_IBD.R")
sourceCpp("~/Dropbox/IBD_IBS/PlasmodiumRelatedness/Code/hmmloglikelihood.cpp")
registerDoParallel(cores = detectCores()-2)
epsilon <- 0.001 # Fix epsilon throughout
set.seed(1) # for reproducibility
nrepeats <- 1000 # number of repeats 
kfixed <- 2
rfixed <- 0.5
load("~/Dropbox/IBD_IBS/PlasmodiumRelatedness/RData/hmmInput_freqs.RData") 
Data_ <- hmmInput_freqs$TM_WGS[,c(1,2)] # Using Thai positions as a template (Why no chromosome 14???)
m <- 100
K <- 2
frequencies <- matrix(0.5, nrow = m, ncol = 2) # all 0.5
pcts = seq(0,100,10)
chrom_clust = 1 # Chromosome onto which to cluster markers
chrom_clust_ind = which(Data_$chrom == chrom_clust)
RUN = F
PDF = T


## Mechanism to generate Ys given fs, distances, k, r, epsilon
simulate_Ys_hmm <- function(frequencies, distances, k, r, epsilon){
  Ys_IBDs <- simulate_data(frequencies, distances, k = k, r = r, epsilon, rho = 7.4 * 10^(-7))
  return(Ys_IBDs)
}

## Mechanism to compute MLE given fs, distances, Ys, epsilon
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon){
  ndata <- nrow(frequencies)
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho = 7.4 * 10^(-7))
  optimization <- optim(par = c(kfixed, 0.5), fn = function(x) - ll(x[1], x[2]))
  rhat <- optimization$par
  return(rhat) 
}

if(RUN){
  #===================================================
  # Generate distances for different spacing strategies 
  #===================================================
  # first sample different disatances
  Distances = lapply(pcts, function(pct){
    indices_clust = sample(chrom_clust_ind, m*pct/100)
    indices_other = sample(nrow(Data_), m - length(indices_clust))
    indices = sort(c(indices_clust, indices_other))
    data_ <- Data_[indices,]
    data_$dt <- c(diff(data_$pos), Inf)
    pos_change_chrom <- 1 + which(diff(data_$chrom) != 0) # find places where chromosome changes
    data_$dt[pos_change_chrom-1] <- Inf
    data_$dt
  })
  
  #===================================================
  # Generate results for different spacing strategies 
  #===================================================
  Results_dts = lapply(Distances, function(distances){
    pars <- foreach (irepeat = 1:nrepeats, .combine = rbind) %dorng% {
      Ys_IBDs = simulate_Ys_hmm(frequencies, distances, kfixed, rfixed, epsilon)
      compute_rhat_hmm(frequencies, distances, Ys_IBDs$Ys, epsilon)}
  })
  names(Results_dts) = pcts
  
  #===================================================
  # Generate results under independence and append
  #===================================================
  distances = rep(Inf, length = m)
  Results_inf <- foreach (irepeat = 1:nrepeats, .combine = rbind) %dorng% {
    Ys_IBDs = simulate_Ys_hmm(frequencies, distances, kfixed, rfixed, epsilon)
    compute_rhat_hmm(frequencies, distances, Ys_IBDs$Ys, epsilon)}
  Results_dts[['-10']] = Results_inf 
  save(Distances, Results_dts, file = '../RData/Results_dts.RData')
}


#===================================================
# Plot error associated with different strategies
#===================================================
if(!RUN){load('../RData/Results_dts.RData')}
if(PDF){pdf('../Plots/Plots_chromosomally_clustered.pdf')}
require(RColorBrewer)

par(mfrow = c(2,1), mar = c(4,4,1,1))
rmse_r = sapply(Results_dts, function(x){sqrt(median((x[,2] - rfixed)^2))})
rmse_k = sapply(Results_dts, function(x){sqrt(median((x[,1] - kfixed)^2))})
plot(x = c(pcts, -10), y = rmse_r, pch = 16, 
     xlab = sprintf('Percent of %s markers targeted to chromosome %s', m, chrom_clust), 
     ylab = 'Root median square error')  
axis(side = 1, at = -10, labels = 'Independence', las = 2, cex.axis = 0.5)

plot(x = c(pcts, -10), y = rmse_k, pch = 16, 
     xlab = sprintf('Percent of %s markers targeted to chromosome %s', m, chrom_clust), 
     ylab = 'Root median square error')
axis(side = 1, at = -10, labels = 'Independence', las = 2, cex.axis = 0.5)

#===================================================
# Plot estimates associated with different strategies
#===================================================
require(RColorBrewer)
cols <- brewer.pal(length(Distances), 'Spectral')

#---------------------------------------------------
# Histograms of r
#---------------------------------------------------
par(mfrow = c(3,4), mar = c(4,4,1,1))
# First plot result assumint independence
hist(Results_dts[[length(Results_dts)]][,2], col = 'gray', yaxt = 'n', freq = F, 
     xlim = c(0,1), xlab = '', ylab = '', main = 'Independence')
abline(v = 0.5, lwd = 2)

# Second plot results for different percent of markers targeted to chrom_clust
for(i in 1:(length(Results_dts)-1)){
  hist(Results_dts[[i]][,2], col = cols[i], yaxt = 'n', freq = F, 
       xlim = c(0,1), main = sprintf('%s percent',names(Results_dts)[i]), xlab = '', ylab = '')
  abline(v = 0.5, lwd = 2)}

#---------------------------------------------------
# Plots of k
#---------------------------------------------------
par(mfrow = c(3,4), mar = c(4,4,1,1))
# First plot result assumint independence
plot(sort(Results_dts[[length(Results_dts)]][,1]), col = 'gray', main = 'Independence', 
     pch = 20, ylim = c(0,20), ylab = 'Estimate of k')

# Second plot results for different percent of markers targeted to chrom_clust
for(i in 1:(length(Results_dts)-1)){
  plot(sort(Results_dts[[i]][,1]), col = cols[i], main = sprintf('%s percent',names(Results_dts)[i]), 
       pch = 20, ylim = c(0,20), ylab = 'Estimate of k')}


#---------------------------------------------------
# Changes in distance with clustering
#---------------------------------------------------
# Proportion NA doesn't change much with clustering: 
par(mfrow = c(2,1))
plot(sapply(Distances, function(x){mean(is.infinite(x))}), ylab = 'Proportion Inf')
plot(sapply(Distances, function(x){mean(x[!is.infinite(x)])}), ylab = 'Mean finite distance')

if(PDF){dev.off()}
