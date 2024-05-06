#' @import bsseq
#' @import glmnet
#' @export

predictAge <- function(bs, adult.age = 20) {

    message("[madrid]: getting methylation levels for each sample")

    meth.mat <- bsseq::getMeth(bs, type = "raw")
    rownames(meth.mat) <- paste0(seqnames(bs),":",start(bs))
    meth.mat <- t(meth.mat)

    message("[madrid]: filtering CpGs to only those predictive of age")
    data(predictiveCpGs)
    meth.mat <- meth.mat[,intersect(predictiveCpGs,colnames(meth.mat))]

    message("[madrid]: predicting age now")

    # split up into chunks if there’s too many samples to run at once
    if (nrow(meth.mat) > 20) {
        nSegments <- round(nrow(meth.mat)/20)
        nRows <- 1:nrow(meth.mat)
        x <- split(nRows, factor(sort(rank(nRows)%%nSegments)))
    }  else {
        nSegments <- 1
        nRows <- 1:nrow(meth.mat)
        x <- split(nRows, factor(sort(rank(nRows)%%nSegments)))
    }
 
    data(lambda.glmnet.Training)
    data(glmnet.Training2)
    pred_t_age <- c()
    for (ii in 1:length(x)) {
        meth.mat.sub <- meth.mat[x[[ii]],]
        set.seed(714)
        pred_t_age.sub <- predict(glmnet.Training2, meth.mat.sub, type = "response", s = lambda.glmnet.Training)		
        pred_t_age <- rbind(pred_t_age,pred_t_age.sub)
        #cat(paste0("Done with chunk ",ii,"\n"))
    }

    # take predicted transformed age and inverse to a predicted chronological age
    pred_t_age <- as.data.frame(pred_t_age)
    iAges <- c()
    for (i in 1:nrow(pred_t_age)) {
	ii <- inverse.transform( tAge = pred_t_age[i,"s1"], adult.age)
	iAges <- rbind(iAges,ii)
    }
    iAges <- as.data.frame(iAges)
    rownames(iAges) <- colnames(bs)
    colnames(iAges) <- "predictedAge"
    return(iAges)
}
