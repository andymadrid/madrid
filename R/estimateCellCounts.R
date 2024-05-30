#' Calculates estimated cell-type proportions from sequencing-based data, using R package meffil to do the heavy lifting
#' @param bs A bsseq object. Removing low-coverage CpGs is not required
#' @param reference Reference used for cell-type estimations. Please reference meffil.list.cell.type.references() for the available references
#' @examples
#' estimateCellCounts(bs)
#' @import bsseq
#' @import bsseq
#' @import meffil
#' @export

estimateCellCounts <- function(bs, reference = NULL) {

    # check that reference is available
    if (!is.null(reference)) {
        if (length(intersect(meffil.list.cell.type.references(),reference)) == 1) {
            message(paste0("[madrid]: will use ", reference," to estimate cell proportions"))
        } else {
            stop(paste0("Your suggested reference (", reference,") is not available.\n\tCheck available references with meffil.list.cell.type.references()"))
        }
    }
    if (is.null(reference)) {
        reference = "blood gse35069 complete"
        message("[madrid]: will use blood gse35069 complete to estimate cell proportions")
   }

    # project methylation levels onto array-based CpGs
    message("[madrid]: gathering CpGs for array-based cell-type estimation now")
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
    overlapsFull <- as.data.frame(findOverlaps(arrayFull.gr,beta.gr))
    arrayFull.gr.subset <- arrayCpGsFull[overlapsFull$queryHits,]
    betaFull.subset <- meth.mat[overlapsFull$subjectHits,]
    rownames(betaFull.subset) <- arrayFull.gr.subset$V4

    message("[madrid]: estimating cell-type proportions now")
    cellProps <- meffil.estimate.cell.counts.from.betas(betaFull.subset,  cell.type.reference=reference, verbose = FALSE)
    return(cellProps)
}
