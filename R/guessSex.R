#' Guesses the sex of samples from a bsseq object
#' @import bsseq
#' @export
#' @param bs A bsseq object, must include CpGs on sex (X and Y) chromosomes.
#' @param fastGuess Logical. Whether to guess sex using the fast method that just uses sex chromosomes, or not. The fast method is not as accurate. Default is FALSE
#' @examples
#' guessSex(bs, fastGuess = FALSE)

guessSex <- function(bs, fastGuess = FALSE) {

    if ( fastGuess == FALSE) {
        message("[madrid]: normalizing methylation values by Z score")

        bs.auto <- bsseq::chrSelectBSseq(bs,seqnames=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"))
        bs.chrX <- bsseq::chrSelectBSseq(bs,seqnames=c("chrX"))
        bs.chrY <- bsseq::chrSelectBSseq(bs,seqnames=c("chrY"))

        # make sure chrX and chrY have data 
        if (nrow(bs.chrX) == 0) {
            stop("No CpGs for chrX! Did you align to the sex chromosomes?")
        }
        if (nrow(bs.chrY) == 0) {
            stop("No CpGs for chrY! Did you align to the sex chromosomes?")
        }

        meth.auto <- bsseq::getMeth(bs.auto, type = "raw")
        meth.chrX <- bsseq::getMeth(bs.chrX, type = "raw")
        meth.chrY <- bsseq::getMeth(bs.chrY, type = "raw")
        auto.mean <- colMeans(meth.auto, na.rm=TRUE)
        auto.sd <- colSds(meth.auto, na.rm=TRUE)
        meth.chrX <- t((t(meth.chrX) - auto.mean) / auto.sd)
        meth.chrY <- t((t(meth.chrY) - auto.mean) / auto.sd)

        meth.chrX[is.na(meth.chrX)] <- 0
        meth.chrY[is.na(meth.chrY)] <- 0

        message("[madrid]: getting statistics from sex chromosomes for each sample")

        median.chrX <- suppressWarnings(colMedians(meth.chrX, na.rm = TRUE))
        median.chrY <- suppressWarnings(colMedians(meth.chrY, na.rm = TRUE))
        sexMedDiff <- median.chrY - median.chrX
        k <- kmeans(sexMedDiff, centers = c(min(sexMedDiff), max(sexMedDiff)))
        dat <- as.data.frame(cbind(median.chrX, median.chrY))
        colnames(dat) <- c("Median_Meth_chrX", "Median_Meth_chrY")

        message("[madrid]: guessing the sex for each sample")
        dat$predictedSex <- ifelse(k$cluster == which.min(k$centers), "Male", "Female")
        return(dat)
    }

    if ( fastGuess == TRUE) {

       message("[madrid]: filtering to only CpGs on sex chromosomes")

        bs.chrX <- bsseq::chrSelectBSseq(bs,seqnames=c("chrX"))
        bs.chrY <- bsseq::chrSelectBSseq(bs,seqnames=c("chrY"))

        # make sure chrX and chrY have data 
        if (nrow(bs.chrX) == 0) {
            stop("No CpGs for chrX! Did you align to the sex chromosomes?")
        }
        if (nrow(bs.chrY) == 0) {
            stop("No CpGs for chrY! Did you align to the sex chromosomes?")
        }

        meth.chrX <- bsseq::getMeth(bs.chrX, type = "raw")
        meth.chrY <- bsseq::getMeth(bs.chrY, type = "raw")

        message("[madrid]: getting medians for each sample")

        median.chrX <- suppressWarnings(colMedians(meth.chrX, na.rm = TRUE))
        median.chrY <- suppressWarnings(colMedians(meth.chrY, na.rm = TRUE))
        sexMedDiff <- median.chrY - median.chrX
        k <- kmeans(sexMedDiff, centers = c(min(sexMedDiff), max(sexMedDiff)))
        dat <- as.data.frame(cbind(median.chrX, median.chrY))
        colnames(dat) <- c("Median_Meth_chrX", "Median_Meth_chrY")

        message("[madrid]: guessing the sex for each sample")
        dat$predictedSex <- ifelse(k$cluster == which.min(k$centers), "Male", "Female")
        return(dat)
    }
}
