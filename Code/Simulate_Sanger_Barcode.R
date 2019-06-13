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
nboot <- 100 # number of parametric bootstrap iterations for
kfixed <- kinit <- 8 
rfixed <- rinit <- 0.5 
rho <- 7.4 * 10^(-7)
load('../RData/sanger_amp_data_for_IBDsim.RData') # Data
Countries = names(frqs_per_country)
RUN = F

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
# (1.6 hr on MacBook Air nrepeats and nboot = 100) 
#=============================================================
distances <- amp_data$dt
PPair_results <- list()
RMSEr_results <- array(NA, dim = c(length(ks),length(rs)), dimnames = list(ks, rs))
RMSEk_results <- RMSEr_results

if(RUN){
  system.time(
    for(Country in Countries){
      frequencies = frqs_per_country[[Country]]
      for(r in rs){
        for(k in ks){
          
          PPair_results[[sprintf('k=%s r=%s',k,r)]] = vector('list', length = nrepeats)
          
          for(i in 1:nrepeats){
            
            # Simulate data
            Ys <- simulate_Ys_hmm(frequencies, distances, k, r, epsilon)
            
            # Estimate r and k for simulated data 
            rkhat <- compute_rhat_hmm(frequencies, distances, Ys, epsilon)
            
            # Generate parametric bootstrap draws
            rkhats_boot = foreach(iboot = 1:nboot, .combine = rbind) %dopar% {
              Ys_boot <- simulate_Ys_hmm(frequencies, distances, k = rkhat[1], r = rkhat[2], epsilon)
              rkhat_boot <- compute_rhat_hmm(frequencies, distances, Ys_boot, epsilon)
              rkhat_boot
            }
            
            # Calculate 95% CIs
            rkhats_CIs <- apply(rkhats_boot, 2, quantile, probs = c(0.025, 0.975))
            
            # Concatenate the results and name 3 x 2 matrix
            to_return_i = rbind(rkhat, rkhats_CIs)
            colnames(to_return_i) = c('k','r')
            PPair_results[[sprintf('k=%s r=%s',k,r)]][[i]] = to_return_i
          }
          
          # Extract point estimates 
          rkhats = sapply(PPair_results[[sprintf('k=%s r=%s',k,r)]], function(x)x['rkhat',])
          
          # Calculate RMSE
          RMSE_k = sqrt(mean((rkhats['k',] - k)^2)) # Single number
          RMSE_r = sqrt(mean((rkhats['r',] - r)^2)) # Single number
        
          # Store RMSE
          RMSEk_results[as.character(k), as.character(r)] = RMSE_k
          RMSEr_results[as.character(k), as.character(r)] = RMSE_r
          
          writeLines(sprintf('\nFinished: r = %s, k = %s',r,k))
        }
      }
    }
  )
}

#save(PPair_results, RMSEr_results, RMSEk_results, file = '../RData/Sanger_Amplicon_SimResults.Rdata')
