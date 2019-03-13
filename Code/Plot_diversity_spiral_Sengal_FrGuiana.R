#=================================================
# Spiral plot of diversity
#=================================================
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

pdf('../Plots/DiversitySpiral.pdf', height = 15, width = 10)

# Generate a plot of diversity over the genome
library(RColorBrewer) # For colours
library(plotrix) # For plotting data radially
par(mfrow = c(2,1))

col_mapper <- function(num){
  y <- colorRamp(rev(brewer.pal(9, "Spectral")))(num)
  z <- rgb(y[,1], y[,2], y[,3], maxColorValue = 255)
  return(z)
}

plt.lns <- seq(1, nrow(data_), 1)
angles <- seq(0, 5*360, length=nrow(data_))%%360

# Normalise diversity
MIN = min(data_$FG.hap.div, data_$S.hap.div)
MAX = max(data_$FG.hap.div, data_$S.hap.div)
div_FG <- (data_$FG.hap.div-MIN)/(MAX-MIN)
div_S <- (data_$S.hap.div-MIN)/(MAX-MIN)

polar.plot(plt.lns, polar.pos=angles+20, main = 'Unpermutated data', 
           labels ='', show.grid.labels = FALSE, show.grid = F, line.col = NA, 
           rp.type = "symbols", point.symbols=20, 
           point.col = col_mapper(div_FG))
polar.plot(plt.lns, polar.pos=angles-20, add = T, 
           labels ='', show.grid.labels = FALSE, show.grid = F, line.col = NA, 
           rp.type = "symbols", point.symbols=20, 
           point.col = col_mapper(div_S))

legend('topright', bty = 'n', pch = 20, col = col_mapper(seq(0,1,0.2)), title = 'Diversity', 
       legend = format((seq(0,1,0.2)*(MAX-MIN))+MIN, digits = 2, drop0trailing = F))



# Generate a plot of diversity over the genome after perturbing windows
ptdata = data_[base::sample(1:nrow(data_), size = nrow(data_), replace = F), ]
plt.lns <- seq(1, nrow(ptdata), 1)
angles <- seq(0, 5*360, length=nrow(ptdata))%%360
MIN = min(ptdata$FG.hap.div, ptdata$S.hap.div)
MAX = max(ptdata$FG.hap.div, ptdata$S.hap.div)
div_FG <- (ptdata$FG.hap.div-MIN)/(MAX-MIN)
div_S <- (ptdata$S.hap.div-MIN)/(MAX-MIN)

polar.plot(plt.lns, polar.pos=angles+20, main = 'Randomly permuted', 
           labels ='', show.grid.labels = FALSE, show.grid = F, line.col = NA, 
           rp.type = "symbols", point.symbols=20, 
           point.col = col_mapper(div_FG))
polar.plot(plt.lns, polar.pos=angles-20, add = T, 
           labels ='', show.grid.labels = FALSE, show.grid = F, line.col = NA, 
           rp.type = "symbols", point.symbols=20, 
           point.col = col_mapper(div_S))

legend('topright', bty = 'n', pch = 20, col = col_mapper(seq(0,1,0.2)), title = 'Diversity', 
       legend = format((seq(0,1,0.2)*(MAX-MIN))+MIN, digits = 2, drop0trailing = F))
dev.off()
