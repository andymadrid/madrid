#' Guesses the sex of samples from a bsseq object
#' @import bsseq
#' @export
#' @param bs A bsseq object, must include CpGs on sex (X and Y) chromosomes.
#' @examples
#' guessSex(bs)

guessSex <- function(bs) {

        message("[madrid]: filtering for only CpGs on sex chromosomes")

        bs.chrX <- bsseq::chrSelectBSseq(bs,seqnames=c("chrX"))
        bs.chrY <- bsseq::chrSelectBSseq(bs,seqnames=c("chrY"))

        # make sure chrX and chrY have data 
        if (nrow(bs.chrX) == 0) {
            stop("No CpGs for chrX! Did you align to the sex chromosomes?")
        }
        if (nrow(bs.chrY) == 0) {
            stop("No CpGs for chrY! Did you align to the sex chromosomes?")
        }

        message("[madrid]: getting statistics from sex chromosomes for each sample")

        meansX <- colMeans(bsseq::getCoverage(bs.chrX, type = "Cov"))
        meansY <- colMeans(bsseq::getCoverage(bs.chrY, type = "Cov"))
        sexRatio <- meansX / meansY
        k <- kmeans(sexRatio, centers = c(min(sexRatio), mean(sexRatio)))
        dat <- as.data.frame(cbind(meansX, meansY, sexRatio))
        colnames(dat) <- c("Mean_Cov_chrX", "Mean_Cov_chrY", "Sex_Cov_Ratio")

        message("[madrid]: guessing the sex for each sample")
        dat$predictedSex <- ifelse(k$cluster == which.min(k$centers), "Male", "Female")
        return(dat)
}
