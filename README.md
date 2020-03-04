# malaria_mltplx_ampseq
All code-related files pertaining to the Neafsey lab malaria mltplx ampseq project

## To-do:
- make code into R package so others (Emily) can take lead
- fix spelling mistake colombia
- simulate CTSA on its own (to compare)
- simulate 24 SNP barcode (to compare)
- simulate Sanger mali (to compare and so code runs for all countries and panels combined)
- simulate paragon panel (to compare)
- Bronwyn asked: what is the minimally informative step (e.g. at what point does the GTseq panel become useless)? Evaluate in light of results on low parasitaemia samples
- summarise number of markers per panel as well as diversity / effective cardinality
- add effective to cardinality
- finish visualisations for different countries and for RMSE
- Angela suggested evaluating panels on a more applied versus statistical level (confidence interval length and RMSE) by enumerating what proportion of samples pairs have r estimates that are statistically distinguishable from some threshold (e.g. rhat > 0.25) that might be used in an application. This type of analysis ought to be done by simulation of the whole genome, followed by marker selection per panel of interest. Otherwise, panels with more markers are liable to have less mendelian variation and thus fewer distinguishable sample pairs simply because there are fewer simulated samples pairs with realised IBD > 0.25 (we want to compare the informativeness of the panels, not mendelian variation in realised IBD).

## Viz ideas for paper: 
- add existing plots to Rmd and combine panels
- add chrom lengths and CTSA annotations to position plot (static)
- add dynamic annotations to effective cardinality and diversity plots (Rshiny) 
