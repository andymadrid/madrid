# madrid: Methylation Age Determination via Read-Informed Data (MADRID) epigenetic clock

## Welcome
Hey, thanks for stopping on by! This is baby's first R package, so bear with me here.

The package is pretty straightforward; all you really need is a bsseq object to start.

At its core, all it really does is filter your bsseq object for the select few CpGs used to predict age from peripheral blood. There's a function to guess the sex of your samples, but that is very much in the beta testing stage of things. Feel free to test it out, just know that your results may vary. To increase its accuracy, more samples are needed to enhance its predictive power.

## Sequencing-based data to array-based clocks

Additional functionality has been added to the predictAge() function to now allow it to predict age from several array-based clocks, including Horvath's (original), Hannum's, Levine's (PhenoAge), Horvath's (Skin-Blood), Lin's, and DunedinPACE. Essentially, it takes a matrix of methylation values, formats them a bit, filters to only those CpGs on the arrays (27k, 450k, and EPICv1), adds CpG names from the arrays, then leans heavily on R packages wateRmelon and DunedinPACE to estimate age from all those clocks. So, even if your data isn't blood-based as the MADRID clock was developed for, you can still estiamte ages from sequencing-based data using these clocks. Pretty cool, right?

## Word of warning

### Blood samples
It should be noted that this package - particularly the sequencing-based clock - is intended for blood biospecimens; the selected CpGs used for age prediction, along with the methylation levels used to guess sex, were tested and validated only in blood samples that I had available. A pan-tissue clock would be great to build and implement, but - as of right now - its intended use is for blood. 

### Genome reference
Samples used to test and validate this clock were aligned to the human genome (hg38) using the UCSC chromsome naming scheme (e.g,. chr1, chr2, chr3, etc). If you happen to have aligned your samples to, say, hg19, those coordinates will have to be lifted over in order to properly get filtered/selected for the clock to work. Currently, there is no built-in function to do this lift over, but can be implemented if the need for it arises. 

## Questions, concerns, collaborations
Do you have your own sequencing-based methylation data that you want to build your own clock with, say in a different tissue? If so, feel free to email me at at madrid2[at]wisc.edu and we can work it out, together. I am always happy to collaborate and/or help! Also, if there's some functinoality that you think would be of great value to add, you can also let me know and I can work on implementing that, as well.

See you space cowboy...
- AM
