###################################################################
# This script plots RMSE based on microhaplotype data from 
# French Guiana and Senegal
###################################################################
load('../RData/Tables_rs_FG_S_dynamic_vs_static.RData')
PDF = T
countries = dimnames(tables_many_repeats_m)[[4]]
panel_strategies = dimnames(tables_many_repeats_m)[[3]]
load(file = '../RData/freq_FG_S.RData')
load('../RData/data_FG_S.RData')

if(PDF){pdf('../Plots/Table_Senegal_FrGuiana_dynamic_vs_static.pdf', height = 7, width = 7)}

# Plot
require(RColorBrewer)
par(mfrow = c(2,2))

  for(country in countries){
    
    for(panel_strategy in panel_strategies){
    
    Table <- tables_many_repeats_m[,,panel_strategy,country]
    ms <- as.numeric(rownames(Table))
    rs <- as.numeric(colnames(Table))
    cols <- rainbow(length(rs),end = 0.8)
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
         main = paste0(panel_strategy, " (", country, ")"), cex.main = 1.2)
    
    
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