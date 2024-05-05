# madrid
R package/scripts for the Methylation Age Determination via Read-Informed Data (MADRID) epigenetic clock

Hey, thanks for stopping on by! This is baby's first R package, so bear with me here.

The package is pretty straightforward; all you really need is a bsseq object to start.

At its core, all it really does is filter your bsseq object for the select few CpGs used to predict age from peripheral blood. There's a function to guess the sex of your samples, but that is very much in the beta testing stage of things. It should be noted that this package is intended for blood biospecimens; the selected CpGs used for age prediction, along with the methylation levels used to guess sex, were tested and validated only in blood samples that I had available. A pan-tissue clock would be great to build and implement, but - as of right now - its intended use is for blood. 

Do you have your own sequence-based methylation data that you want to build your own clock with, say in a different tissue? If so, feel free to email me at at madrid2[at]wisc.edu and we can work it out, together. I am always happy to collaborate! Also, if there's some functinoality that you think would be of great value to add, you can also let me know and I can work on implementing that, as well.
