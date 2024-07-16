#' Predicts the age  of samples from a bsseq object
#' @param bs A bsseq object
#' @param adult.age Adult age for transformation. Default is 20
#' @param fastImpute Logical. Whether fast imputation should be performed. This is not advised as it leads to less accurate predictions. Default is FALSE
#' @param imp Number of imputations to be performed by mice() during imputation of missing methylation values. Default is 5
#' @param calcArrayClocks Logical. Whether array-based clocks (Horvath, Hannum, PhenoAge, SkinBlood, Lin) should be calculated. Default is TRUE.
#' @examples
#' predictAge(bs, adult.age = 20, fastImpute = FALSE, imp = 5, calcArrayClocks = TRUE)
#' @import bsseq
#' @import glmnet
#' @import mice
#' @import impute
#' @import wateRmelon
#' @import DunedinPACE
#' @export

predictAge <- function(bs, adult.age = 20, fastImpute = FALSE, imp = 5, calcArrayClocks = TRUE) {

    message("[madrid]: getting methylation levels for each sample")

    meth.mat.full <- bsseq::getMeth(bs, type = "raw")
    meth.mat <- meth.mat.full
    rownames(meth.mat) <- paste0(seqnames(bs),":",start(bs))

    message("[madrid]: filtering CpGs to only those predictive of age")
    data(predictiveCpGs)
    meth.mat <- meth.mat[match(predictiveCpGs,rownames(meth.mat)),]
    meth.mat <- t(meth.mat)

    data(crimson.standard)
    
    if (sum(is.na(meth.mat)) == 0) {
            message("[madrid]: no missing data, skipping imputation")
    } else {

        message("[madrid]: imputing missing data . . . this may take some time . . .")

    # Find columns that are missing data for all samples, impute with crimson standards
        numSamples <- nrow(meth.mat)
        for (i in 1:ncol(meth.mat)) {
            numSamplesMissing <- sum(is.na(meth.mat[,i]))
            if (numSamplesMissing == numSamples) {
                meth.mat[,i] <- crimson.standard[i]
            }
        }

    # If fast impute selected, impute remaining missing data with crimson standard
        if (fastImpute == TRUE) {
            message("[madrid]: fast imputation selected (not advised), moving along.")
            for (i in 1:ncol(meth.mat)) {
               meth.mat[is.na(meth.mat[,i]),i] <- crimson.standard[i]
            }
        } else {

    # if fast impute not selected, impute with random forest
            colnames(meth.mat) <- paste0("CpG",1:ncol(meth.mat))
            captured <- capture.output(tmpData <- mice::mice(meth.mat, m = imp, method = "rf", seed = 714))
            imputedData <- mice::complete(tmpData, 1)

        # get averages across imputed datasets
            for(m in 2:tmpData$m) {
                imputedData <- imputedData + mice::complete(tmpData, m)
            }
            imputedData <- imputedData / tmpData$m

        # bringing it all back home
            meth.mat <- imputedData
            colnames(meth.mat) <- predictiveCpGs
            for (i in 1:ncol(meth.mat)) {
               meth.mat[is.na(meth.mat[,i]),i] <- median(na.exclude(meth.mat[,i]))
            }
        }
    }

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
        pred_t_age.sub <- predict(glmnet.Training2, as.matrix(meth.mat.sub), type = "response", s = lambda.glmnet.Training)		
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
    colnames(iAges) <- "MADRID_Age"

    if ( calcArrayClocks == TRUE) {

        message("[madrid]: gathering CpGs for array-based clocks now")
        meth.mat <- bsseq::getMeth(bs, type = "raw")
        chr <- as.data.frame(seqnames(bs))
        colnames(chr) <- "chr"
        pos <- as.data.frame(start(bs))
        colnames(pos) <- "pos"
        chrPos <- cbind(chr,pos)
        beta.gr <- with(chrPos,GRanges(chr,IRanges(pos,pos+1)))
        data(arrayCpGs)
        data(arrayCpGsFull)
        arrayFull.gr <- with(arrayCpGsFull,GRanges(V1,IRanges(V2,V3)))
        overlapsFull <- suppressWarnings(as.data.frame(findOverlaps(arrayFull.gr,beta.gr)))
        arrayFull.gr.subset <- arrayCpGsFull[overlapsFull$queryHits,]
        betaFull.subset <- meth.mat[overlapsFull$subjectHits,]
        rownames(betaFull.subset) <- arrayFull.gr.subset$V4
        array.gr <- with(arrayCpGs,GRanges(V1,IRanges(V2,V3)))
        overlaps <- suppressWarnings(as.data.frame(findOverlaps(array.gr,beta.gr)))
        array.gr.subset <- arrayCpGs[overlaps$queryHits,]
        beta.subset <- meth.mat[overlaps$subjectHits,]
        rownames(beta.subset) <- array.gr.subset$V4
        message("[madrid]: calculating Dunedin methylation pace now")
        dunedin <- DunedinPACE::PACEProjector(betaFull.subset)
        message("[madrid]: imputing missing data for array-based CpGs now")
        captured <- capture.output(tmpArrayData <- impute::impute.knn(beta.subset, rng.seed = 714))
        beta.subset <- tmpArrayData$data
        message("[madrid]: calculating age from array-based clocks now")
        predAge.array <- wateRmelon::agep(beta.subset, method='all')
        predAge.array <- predAge.array[, c("horvath.horvath.age", "hannum.hannum.age", "phenoage.phenoage.age", "skinblood.skinblood.age", "lin.lin.age")]
        predAge.array$DunedinPACE <- dunedin$DunedinPACE
        colnames(predAge.array) <- c("Horvath_Age", "Hannum_Age", "PhenoAge", "SkinBlood_Age", "Lin_Age", "DunedinPACE")
        iAges <- cbind(iAges, predAge.array)
        return(iAges)

    } else {
        message("[madrid]: no ages from array-based clocks will be estimated")

        return(iAges)

    }
}
