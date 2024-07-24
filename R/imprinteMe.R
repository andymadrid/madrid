#' Overlaps DMRs with known imprinted regions
#' @param dmrs A GenomicRanges object containing the chromosome, start, and end of regions of interest (e.g., differentially methylated regions (DMRs)), with coordinates corresponding to the hg38 version of the human genome
#' @examples
#' imprintMe(dmrs)
#' @import GenomicRanges
#' @export

imprintMe <- function(dmrs) {

    # check that the object is a GenomicRanges object
    if (intersect("GRanges", class(dmrs)[1]) == "GRanges") {        
        message("[madrid]: Looking for overlaps with imprinted regions now")
    } else {
        stop("[madrid]: Your object does not seem to be a GRanges object. Please double check your input object.")
    }

    # load in known imprinted regions
    data(imprintome)
    
    # look for overlaps
    overlaps <- suppressWarnings(as.data.frame(GenomicRanges::findOverlaps(dmrs, imprintome)))

    # filter for the overlaps, then return object
    dmrs.sub <- dmrs[overlaps$queryHits, ]
    imprintome.sub <- imprintome[overlaps$subjectHits, ]    
    mcols(dmrs.sub) <- cbind(mcols(dmrs.sub), as.data.frame(imprintome.sub))
    return(dmrs.sub)
}
