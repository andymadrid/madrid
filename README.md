# madrid: Methylation Age Determination via Read-Informed Data (MADRID) epigenetic clock

## Welcome
Hey, thanks for stopping on by! This is baby's first R package, so bear with me here.

The package is pretty straightforward; all you really need to start is a bsseq object.

At its core, all it really does is filter your bsseq object for the select few CpGs used to predict age from peripheral blood. There's a function to guess the sex of your samples, but that is very much in the beta testing stage of things. Feel free to test it out, just know that your results may vary. To increase its accuracy, more samples are needed to enhance its predictive power.

## Vignette

Check out the tutorial of the general functionality of the package at:

http://htmlpreview.github.io/?https://raw.githubusercontent.com/andymadrid/madrid/main/vignettes/introduction.html

## Sequencing-based data to array-based clocks

Additional functionality has been added to the predictAge() function to now allow it to predict age from several array-based clocks, including Horvath's (original), Hannum's, Levine's (PhenoAge), Horvath's (Skin-Blood), Lin's, and DunedinPACE. Essentially, it takes a matrix of methylation values, formats them a bit, filters to only those CpGs on the arrays (27k, 450k, and EPICv1), adds CpG names from the arrays, then leans heavily on R packages wateRmelon and DunedinPACE to estimate age from all those clocks. So, even if your data isn't blood-based as the MADRID clock was developed for, you can still estimate ages from sequencing-based data using these clocks. Pretty nifty, huh?

## Imputation

Your dataset may very well have some missing CpGs for either the sequencing-based clock and/or the array-based clocks, which is totally okay. It's to be expected, especially if sequencing depth was on the lower end. To get around this, an imputation step has been implemented in the predictAge() function which calls the mice() function to impute missing data, using a random forest method. In my testing I found that - at least for the sequencing-based clock - mice() run with random forest then averaged across imputed datasets yielded better estimates than, say, impute.knn() did, which is the method commonly used for most array-based clocks. An advanced user can feel free to edit the code as they see fit to change the imputation method. The currently employed method just happened to be the best that I found in my testing, but may not be robust enough for all datasets. However, an astute observer will notice that the imputation step for the array-based clocks actually uses impute.knn(). I found that it worked quite well for the array clocks (who would have guessed?) and worked pretty fast.  

## Cell-type proportion estimates

I've also included a function (estimateCellCounts()) that takes your sequencing data and estimates cell-type proportions. To do this, it formats and filters data from a bsseq object then leans heavily on the R package meffil to do the calculation. There are other packages out there that do similar things (e.g., methylcc), but I found - through my own testing, using in-house datasets - that this method produced results that were more correlative with array data; correlations between
estimates from samples run on the EPIC array and the same samples that were sequenced were ~0.95.

## Stochastic epigenetic mutations (SEMs)

While the main function of this package is to estimate age from sequencing-based data, others have recently been moving to detecting stochastic epigenetic mutations (SEMs) in DNA methylation data. So, I've implemented a function to do just that here. It calculates the interquartile range (IQR) for each CpG, then goes sample by sample, CpG by CpG to identify SEMs. For more information on SEMs, please refer to Gentilini et al (2015) Aging and Markov et al (2024) GeroScience. Please note that the function is slow and can be computationally expensive. I recommend larger sample sizes (>50) and higher coverage (>10) to reduce the number of CpGs analyzed and to only keep those with substantial evidence of their estimated methylation levels.

## Word of warning

### Blood samples
It should be noted that this package - particularly the sequencing-based clock - is intended for blood biospecimens; the selected CpGs used for age prediction, along with the methylation levels used to guess sex, were tested and validated only in blood samples that I had available. A pan-tissue clock would be great to build and implement, but - as of right now - its intended use is for blood. That being said, most of the array-based clocks that are - by default - calculated in the predictAge() function are agnostic to tissue. So, at the very least, one can use this package to estimate the age of their sequencing-based samples using those pan-tissue clocks.

### Genome reference
Samples used to test and validate this clock were aligned to the human genome (hg38) using the UCSC chromsome naming scheme (e.g,. chr1, chr2, chr3, etc). If you happen to have aligned your samples to, say, hg19, those coordinates will have to be lifted over in order to properly get filtered/selected for the clock to work. I added in a function liftBuild() that does just that.  

## Questions, concerns, collaborations
Do you have your own sequencing-based methylation data that you want to build your own clock with, say in a different tissue? If so, feel free to email me at at madrid2[at]wisc.edu and we can work it out, together. I am always happy to collaborate and/or help! Also, if there's some functinoality that you think would be of great value to add, you can also let me know and I can work on implementing that, as well.

See you space cowboy...
- AM
