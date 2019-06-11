###################################################################
# This script generates RMSE based on parasite pairs whose data are
# simulated according to the 
###################################################################
## Set up
rm(list = ls())
set.seed(1) # for reproducibility
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
nrepeats <- 100 # number of repeats 
kfixed <- kinit <- 8 
rfixed <- rinit <- 0.5 
rho <- 7.4 * 10^(-7)
load('../RData/sanger_amp_data_for_IBDsim.RData') # Data
Countries = names(frqs_per_country)
RUN = T

## Mechanism to generate Ys given fs, distances, k, r, epsilon
simulate_Ys_hmm <- function(frequencies, distances, k, r, epsilon){
  Ys <- simulate_data(frequencies, distances, k = k, r = r, epsilon, rho)
  return(Ys)
}

## Mechanism to compute MLE given fs, distances, Ys, epsilon
compute_rhat_hmm <- function(frequencies, distances, Ys, epsilon){
  ndata <- nrow(frequencies)
  ll <- function(k, r) loglikelihood_cpp(k, r, Ys, frequencies, distances, epsilon, rho)
  optimization <- optim(par = c(kinit, rinit), fn = function(x) - ll(x[1], x[2]))
  rhat <- optimization$par
  return(rhat)
}

#=============================================================
# Simulate data
#=============================================================
distances <- amp_data$dt

if(RUN){
  system.time(
    for(Country in Countries){
      frequencies = frqs_per_country[[Country]]
      for(r in rs){
        for(k in ks){
          estimands <- foreach(irepeat = 1:nrepeats, .combine = list) %dorng% {
            Ys <- simulate_Ys_hmm(frequencies, distances, k, r, epsilon)
            c(compute_rhat_hmm(frequencies, distances, Ys, epsilon))
          }
          RMSE_r = sqrt(mean((estimands[,2] - r)^2))
          RMSE_k = sqrt(mean((estimands[,2] - r)^2))
          # Add confidence intervals
        }
      }
    }
  )
}


