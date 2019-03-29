# tidying data and finding top windows
suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
})

setwd("/Volumes/seq_plasmodium/emilylav/subtelomeric")
load(subtelomeric.Senegal.overall.df.RData)
load(subtelomeric.FG.overall.df.RData)

# combine data from both regions
colnames(subtelomeric.Senegal.overall.df) <- c("start", "stop", "S.hap.div", "S.n", "chromosome")
colnames(subtelomeric.FG.overall.df) <- c("start", "stop", "FG.hap.div", "FG.n", "chromosome")
subtelomeric.overall.df <- join(subtelomeric.Senegal.overall.df, subtelomeric.FG.overall.df)

# filter to top windows: >0.5 scores in both regions
top.subtelomeric.overall.df <- subtelomeric.overall.df %>% filter(S.hap.div > 0.5 & FG.hap.div > 0.5)

# convert both to tidy data format
subtelomeric.overall.df.tidy <- gather(subtelomeric.overall.df, metric, value, -start, -stop, -chromosome)
top.subtelomeric.overall.tidy <- gather(top.subtelomeric.overall.df, metric, value, -start, -stop, -chromosome)

save(subtelomeric.overall.df, subtelomeric.overall.df.tidy, file = "overall.subtelomeric.RData")
save(top.subtelomeric.overall.df, top.subtelomeric.overall.tidy, file = "top.windows.subtelomeric.RData")
