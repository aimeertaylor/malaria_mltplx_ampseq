###################################################################
# This script plots RMSE based on microhaplotype data from 
# French Guiana and Senegal
###################################################################
load('../RData/Tables_rs_FG_S.RData')
PDF = T
if(PDF){pdf('../Plots/Table_Senegal_FrGuiana.pdf', height = 15, width = 10)}

# Plot
require(RColorBrewer)
par(mfrow = c(4,2))

for(panel_strategy in panel_strategies){
  for(country in countries){
    
    Table <- tables_many_repeats_m[,,panel_strategy,country]
    ms <- as.numeric(rownames(Table))
    rs <- as.numeric(colnames(Table))
    cols <- brewer.pal(length(rs), 'Spectral')
    max_m <- max(ms)
    
    # Calculate average effective Ks
    if(country == 'French Guiana'){
      av_Keff = sapply(inds_store[[panel_strategy]],function(inds){
        av_Keff = mean(data_$Keff_FG[inds])}) 
    } else {
      av_Keff = sapply(inds_store[[panel_strategy]],function(inds){
        av_Keff = mean(data_$Keff_S[inds])}) 
    }
    
    # Hack
    av_Keff = av_Keff[!is.na(av_Keff)]
    
    plot(NULL, 
         ylim = c(0, pmax(max(Table), 0.35)), 
         xlim = range(ms), bty = 'n', xaxt = 'n', las = 2, 
         ylab = 'Root mean squared error', 
         xlab = 'Number of markers (effective cardinality)', 
         main = paste0(panel_strategy, " (", country, ")"), cex.main = 0.75)
    axis(side = 1, at = ms, labels = paste0(ms, " (", round(av_Keff,1), ")"))
    for(j in 1:ncol(Table)){
      lines(y = Table[,j], x = ms, pch = 20, panel.first = grid(), 
            type = 'b', col= cols[j])
    }
    legend('topright', bty = 'n', pch = 20, col = cols, legend = format(round(rs,2), drop0trailing = F), 
           title = expression('Relatedness,'~italic(r)))
  }
}
if(PDF){dev.off()}