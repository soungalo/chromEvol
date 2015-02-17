# chromEvol
chromEvol is a free, open-source program for analyzing changes in chromosome-number along a phylogeny. Given a phylogeny and chromosome number data, the program infers the location and type of chromosome number transitions along the tree.  
Features		
- Test hypotheses about chromosome numbe evolution
- Determine ploidy levels for taxa
- Find the chromosome base-number of a phylogeny

Chromosome number is a remarkably dynamic feature of eukaryotic evolution. Chromosome numbers can change by a duplication of the whole genome (a process termed polyploidy), or by single chromosome changes (ascending dysploidy via, e.g., chromosome fission or descending dysploidy via, e.g., chromosome fusion).
Of the various mechanisms of chromosome number change, polyploidy has received significant attention because of the impact such an event may have on the organism.
ChromEvol implements a series of likelihood models for the evolution of chromosome numbers. By comparing the fit of the different models to biological data, it may be possible to gain insight regarding the pathways by which the evolution of chromosome number proceeds. For each model, the program estimates the rates for the possible transitions assumed by the model, infers the set of ancestral chromosome numbers, and estimates the location along the tree for which polyploidy events (and other chromosome number changes) occurred. For further methodological details, see the publications and manual on the Downloads page.
ChromEvol was first introduced in Mayrose et al in 2010. Further options and new evolutionary models were added in Glick & Mayrose in 2014. Version 2.0 is now available for download, along with an external computational pipeline which fascilitates easy ploidy inference.
When citing the chromEvol program, please refer to: Mayrose et al, 2010 and Glick & Mayrose, 2014.

Software website: http://www.tau.ac.il/~itaymay/cp/chromEvol/index.html
