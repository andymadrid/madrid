#' Guesses the sex of samples from a bsseq object
#' @import bsseq
#' @export
#' @param bs A bsseq object, must include CpGs on sex (X and Y) chromosomes.
#' @param cut Cutoff for determining female/male samples. Default is 0.
#' @examples
#' guessSex(bs)

guessSex <- function(bs, cut = 0) {

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
