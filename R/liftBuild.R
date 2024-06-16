#' Convert coordinates from hg19 (or hg38) to hg38 (or hg19)
#' @param bs A bsseq object
#' @param current The current human genome build. Can be one of hg19 or hg38 
#' @param new The human genome build that you would like to convert to. Can be one of hg19 or hg38
#' @examples
#' liftBuild(bs, current = "hg38", new = "hg19")
#' @import AnnotationHub
#' @import rtracklayer
#' @import MatrixGenerics
#' @import matrixStats
#' @import GenomicRanges
#' @import dplyr
#' @export

liftBuild <- function(bs, current = c("hg19", "hg38"), new = c("hg19", "hg38")) {

    # check the current and new version of the genome you want
    if (current == "hg38" && new == "hg19") {
        direction <- "AH14108"
    }    else if (current == "hg19" && new == "hg38") {
        direction <- "AH14150"
    }  else {
        stop("Currently only can go from hg19 to hg38 (or the reverse). Consider what builds you’re trying to lift")
    }
    message(paste0("[madrid]: Will lift your bsseq object from ", current," to ", new))
    # hub <- AnnotationHub()
    # chains <- query(hub, c(current, new, "chainfile"))
    # AH14108 if hg38 to hg19
    # AH14150 if hg19 to hg38
    mcols(bs)$cpgs <- 1:length(bs)
    hgCurrent <- rowRanges(bs)
    hgCurrent$cpgs <- 1:length(hgCurrent)
    hgNew <- hgCurrent %>% rtracklayer::liftOver(AnnotationHub::AnnotationHub()[[direction]]) %>% unlist()
    bsLifted <- bs[which(hgNew$cpgs %in% hgCurrent$cpgs)]
    rowRanges(bsLifted) <- hgNew
    genome(bsLifted) <- new
    lost <- length(bs) - length(bsLifted)
    message(paste0("[madrid]: Lifting over resulted in a loss of ", lost," CpGs"))
    return(bsLifted)
}
