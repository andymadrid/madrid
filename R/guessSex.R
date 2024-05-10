#' Guesses the sex of samples from a bsseq object
#' @import bsseq
#' @export
#' @param bs A bsseq object, must include CpGs on sex (X and Y) chromosomes.
#' @examples
#' guessSex(bs)

guessSex <- function(bs) {

    message("[madrid]: filtering to only CpGs on sex chromosomes")

    bs.chrX <- bsseq::chrSelectBSseq(bs,seqnames=c("chrX"))
    bs.chrY <- bsseq::chrSelectBSseq(bs,seqnames=c("chrY"))

    meth.chrX <- bsseq::getMeth(bs.chrX, type = "raw")
    meth.chrY <- bsseq::getMeth(bs.chrY, type = "raw")

    message("[madrid]: getting averages for each sample")

    means.chrX <- colMeans(na.exclude(meth.chrX))
    means.chrY <- colMeans(na.exclude(meth.chrY))

    dat <- as.data.frame(cbind(means.chrX, means.chrY))
    colnames(dat) <- c("Mean_Meth_chrX", "Mean_Meth_chrY")

    message("[madrid]: guessing the sex for each sample")
    for (i in 1:nrow(dat)) {

        if (dat[i, "Mean_Meth_chrX"] > 0.80 & dat[i, "Mean_Meth_chrX"] < 0.87) { 
            dat[i, "predictedSex"] <- "Male"
        }

        if (dat[i, "Mean_Meth_chrX"] > 0.73 & dat[i, "Mean_Meth_chrX"] < 0.80) { 
            dat[i, "predictedSex"] <- "Female"
        }

        if (dat[i, "Mean_Meth_chrX"] < 0.73 | dat[i, "Mean_Meth_chrX"] > 0.87) { 
            dat[i, "predictedSex"] <- "Unknown"

            if (dat[i, "Mean_Meth_chrX"] < 0.73) {
                dat[i, "warning"] <- "Cannot confidently guess as methylation levels were lower than those used to determine sex"
            }

            if (dat[i, "Mean_Meth_chrX"] > 0.87) {
                dat[i, "warning"] <- "Cannot confidently guess as methylation levels were higher than those used to determine sex"
            }
        }
    }
    return(dat)
}
