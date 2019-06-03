###################################################################
# This script simulates data using targets selected from 
# French Guiana and Senegal
# Is it better to pick targets dynamically or statically? 
###################################################################

## Set up
rm(list = ls())
set.seed(1)
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
ks <- c(1,10,100)
samplesizes = c(96*c(0.25, 1:3)) # add maximum sample size when working with real data
set.seed(1) # for reproducibility
nrepeats <- 1000 # number of repeats 
kfixed <- 8
rfixed <- 0.5 
panel_strategies = c('static', 'dynamic')
countries = c('French Guiana', 'Senegal')
RUN = T
load(file = '../RData/freq_FG_S.RData')
load('../RData/data_FG_S.RData')
source('./target_selection.R')
Mmax = nrow(data_)
data_$Med_Keff = apply(data_[,c('Keff_FG', 'Keff_S')], 1, median)
data_$wID = 1:Mmax

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

tables_many_repeats_m <- array(dim = c(length(samplesizes), length(rs), length(panel_strategies), length(countries)),
                               dimnames = list(samplesizes, rs, panel_strategies, countries))
inds_store <- list() # Store inds 

#=====================================
# Calculate RMSE 
#=====================================
if(RUN){
  system.time(
    for(panel_strategy in panel_strategies){
      for (m in samplesizes){
        
        # Choose strategy 
        target_selection = ifelse(panel_strategy == 'dynamic', 
                                  target_selection_dynmic_dt, 
                                  target_selection_static_dt)
        
        # Extract window ID
        inds = target_selection(data_, m)$wID
        
        # Extract frequencies
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
        
        # Store selected window IDs
        inds_store[[panel_strategy]][[as.character(m)]] = inds
      }
    }
  )
  save(tables_many_repeats_m, inds_store, file = '../RData/Tables_rs_FG_S_dynamic_vs_static.RData')
}


