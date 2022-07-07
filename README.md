# malaria_mltplx_ampseq
Initial code-related files pertaining to the Neafsey lab malaria mltplx ampseq project besides the R package https://github.com/artaylor85/paneljudge (created April 2020).

*Important*: all commits prior to and inc. 11th April 2020 pertain to work done without use of the package paneljudge. 

# To-do (also see google doc. of ms)
- simulate paragon panel 
- simulate CTSA on its own 
- simulate 24 SNP barcode 
- simulate Sanger Mali (to compare and so code runs for all countries and panels combined)
- Bronwyn asked: what is the minimally informative panel (e.g. at what point does a panel become useless)? Evaluate in light of results on low parasitaemia samples, when realistic drop out will be apparent.  
- Angela suggested evaluating panels on a more applied versus statistical level (confidence interval length and RMSE) by enumerating what proportion of samples pairs have r estimates that are statistically distinguishable from some threshold (e.g. rhat > 0.25) that might be used in an application. This type of analysis ought to be done by simulation of the whole genome, followed by marker selection per panel of interest. Otherwise, panels with more markers are liable to have less Mendelian variation and thus fewer distinguishable sample pairs simply because there are fewer simulated samples pairs with realised IBD > 0.25 (we want to compare the informativeness of the panels, not Mendelian variation in realised IBD).




