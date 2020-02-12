###################################################################
# This script generates RMSE based on parasite pairs whose data are
# simulated according to the GTseq targets listed on google drive 
# with and without CSP/TRAP/SERA2/AMA1
#
# (I need to get the frequencies for these targets from specific 
# populations)
###################################################################
## Set up
rm(list = ls())


# I haven't got contact with googledrive to work yet but could be 
# an later option to prevent version control issues 
# library("googledrive") 

# Notes from sanger barcode version
# Folling inports  "amp_data" and "frqs_per_country"
load('../RData/sanger_amp_data_for_IBDsim.RData') # Data

# ===============================================
# Ideally this should be done in a separate file e.g. Process_Sanger_Barcode
# Import targets: 
Targets = read.csv('../RawData/GTseq_all_targets_08Nov19_Primers_UDI-AS - target coordinates.csv',
         row.names = "amplicon_name", stringsAsFactors = F) 

amp_data = data.frame(Amplicon_name = rownames(Targets), 
                      Chr = Targets$chromosome, 
                      Start = Targets$start, 
                      Stop = Targets$end, 
                      pos = floor(Targets$start + 0.5*(Targets$end - Targets$start)), # position as floor of median 
                      chrom = as.numeric(sapply(strsplit(Targets$chromosome, split = "_"), function(x)x[2])))

# Re-order then add dt
amp_data$dt = apply(data_[,c('start','stop')] , 1, function(x){median(as.numeric(x))}) 

# Extract position as median
reordered_data_list = list()
for(chr in sort(unique(amp_data$chrom))){
  x = amp_data$pos[chr == amp_data$chrom] # Extract positions per chromosome 
  print(all(x == cummax(x))) # If not all True, not all are monotonically increasing
  inds = sort.int(x, index.return = T)$ix # Sort positions
  reordered_data_list[[chr]] = amp_data[chr == amp_data$chrom, ][inds, ]
}

# Check now all in order: yes
amp_data = do.call(rbind, reordered_data_list)
all(sapply(unique(amp_data$chrom), function(chr){
  x = amp_data$pos[chr == amp_data$chrom] 
  all(x == cummax(x)) 
}))

# Create distances
amp_data$dt <- c(diff(amp_data$pos), Inf)
pos_change_chrom <- 1 + which(diff(amp_data$chrom) != 0) # find places where chromosome changes
amp_data$dt[pos_change_chrom-1] <- Inf

# Make a fake frqs_per_country for now
frqs_per_country 
# ===============================================

# Code for simulating under HMM
set.seed(1) # for reproducibility
library(dplyr)
library(Rcpp)
library(doParallel)
library(doRNG)
library(MCMCpack) # For dirichlet
source("~/Dropbox/IBD_IBS/PlasmodiumRelatedness/Code/simulate_data.R")
sourceCpp("~/Dropbox/IBD_IBS/PlasmodiumRelatedness/Code/hmmloglikelihood.cpp")
registerDoParallel(cores = detectCores()-2)

# Specify options
epsilon <- 0.001 # Fix epsilon throughout
rs <- seq(from = 0.01, to = 0.99, length.out = 5)
ks <- c(1,10,50)
nrepeats <- 10 # number of repeats 
nboot <- 10 # number of parametric bootstrap iterations for
kfixed <- kinit <- 8 
rfixed <- rinit <- 0.5 
rho <- 7.4 * 10^(-7)

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
# (1.6 hr on MacBook Air nrepeats and nboot = 100) 
#=============================================================
distances <- amp_data$dt


if(RUN){
  system.time(
    for(Country in Countries){
      
      frequencies = frqs_per_country[[Country]]
      PPair_results <- list()
      RMSEr_results <- array(NA, dim = c(length(ks),length(rs)), dimnames = list(ks, rs))
      RMSEk_results <- RMSEr_results
      
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
      save(PPair_results, RMSEr_results, RMSEk_results, 
           sprintf(file = '../RData/Sanger_Amplicon_SimResults_%s.Rdata', Country))
    }
  )
}

