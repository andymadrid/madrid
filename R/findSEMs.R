#' Calculates CpG-level statistics and identifies stochastic epigenetic mutations (SEMs) from sequencing-based data
#' @param bs A bsseq object
#' @param minCov Numeric value. The minimum number of reads to overlap a CpG sites by all samples for it to be used for downstream analysis. Note, lowering the number will likely increase the number of CpGs that are tested for SEM identification but at two big costs: 1) increase the computation time as it will be required to analyze more CpGs and 2) likely introduce false positives as less coverage generally equals less certainty of the estimation of methylation at a given CpG site. Default is 10. 
#' @param minSamples Numeric value. The minimum number of samples required to perform the SEM analysis. A larger sample size (e.g., >50) is generally required for an analysis such as this. Default is 50.
#' @param numCores Numeric value. The number of cores used to parallelize processing. Default is 1.
#' @param saveCpGStats Logical. Whether CpG-level statistics should be saved. Default is FALSE.
#' @param saveIndSEMs Logical. Whether sample-specific SEMs should be saved. Note, these files may be large in size, depending on how many SEMs are identified per sample. Default is FALSE.
#' @param saveDir Directory where outputs should be saved. Default is current working directory.
#' @param verbose Logical. Whether output of functions (e.g., number of CpGs assess per sample) should be verbose or not. Default is TRUE.
#' @examples
#' findSEMs(bs, minCov = 10, minSamples = 50, numCores = 1, saveCpGStats = FALSE, saveIndSEMs = FALSE, saveDir = "", verbose = TRUE)
#' @import bsseq
#' @import moments
#' @import parallel
#' @export
findSEMs <- function(bs, minCov = 10, minSamples = 50, numCores = 1, saveCpGStats = TRUE, saveIndSEMs = FALSE, saveDir = getwd(), verbose = TRUE) {

# register cores for parallel processing
    if (numCores > 1) {
        message("[madrid]: registering cores for parallel processing")
        availCores <- parallel::detectCores()
        if (numCores > availCores) {
            warning(paste(" Set numCores exceeds available cores. Will use ", availCores-2, "instead."))
            numCores <- availCores - 2
        }
        cl <- parallel::makeCluster(numCores)
        parallel::clusterEvalQ(cl, library(moments))
    }

    # filter CpGs based on coverage
    message("[madrid]: filtering low-coverage CpGs from further analysis")
    loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov")<minCov) == 0)
    bs <- bs[loci.idx,]

    # get matrix of methylation values
    message("[madrid]: getting methylation estimates from CpGs")
    meth.mat <- bsseq::getMeth(bs, type = "raw")
    rownames(meth.mat) <- paste0(seqnames(bs), ":", start(bs))

    # check that there are enough samples for SEM calculation
    if (ncol(meth.mat) < minSamples) {
        stop(paste0("Number of samples is less than "),minSamples,"! Cannot calculate SEMs with such few samples.")
    }

    # function to calculate first and third quantiles, interquartiles ranges, and other CpG-level statistics
    getStats <- function(cpg) {	
        c(q1 = stats::quantile(cpg, 0.25),
        q3 = stats::quantile(cpg, 0.75),
        iq = stats::IQR(cpg),
        kurtosis = moments::kurtosis(cpg),
        skewness = moments::skewness(cpg),
        mean = mean(cpg),
        min = min(cpg),
        max = max(cpg),
        sd = sd(cpg))
    }

    message("[madrid]: getting CpG-level statistics now")
    if (numCores > 1) {
        cpgStats <- parallel::parApply(cl, meth.mat, MARGIN=1, FUN=getStats)
    } else {
        cpgStats <- apply(meth.mat, MARGIN=1, FUN=getStats)
    }
    cpgStats <- t(cpgStats)
    colnames(cpgStats) <- c("q1","q3","iq","kurtosis","skewness","mean","min","max","sd")
    cpgStats <- as.data.frame(cpgStats)
    cpgStats$range <- cpgStats$max - cpgStats$min
    if (saveCpGStats == TRUE) {
        message(paste0("[madrid]: will save CpG-level statistics at ", saveDir, "cpgStats.rda"))
        save(cpgStats, file =paste0(saveDir, "/cpgStats.rda"))
    } else {
        message("[madrid]: CpG-level statistics will not be saved")
    }

    # convert meth.mat to a df for easier handling below
    meth.mat <- as.data.frame(meth.mat)

    # workhorse . . . time to find SEMs
    message("[madrid]: Finding SEMs now . . .")
    sampleSEMs <- c()

    # set up empty objects/counts for SEM data to be stored in as you move along
    vec <- 1:ncol(meth.mat)
    for (j in seq_along(vec)) {
        sample_id <- colnames(meth.mat[j])
        hyper_sems_j <- 0
        hypo_sems_j <- 0
        hypo_vector <- c()
        hyper_vector <- c()

        if (verbose ==TRUE) {
            message(paste0("\t Starting sample ", sample_id, " . . ."))
        }

        for (i in 1:nrow(meth.mat)) {

        if (verbose ==TRUE) {
            if(abs(i) %% 1000000 == 0) {
                message(paste0("\t\t Done with ", i, " CpGs . . ."))
            }
        }

        # only look at CpGs with IQRs that are not zero
        if (cpgStats[i, "iq"] > 0) {

            # look for hypoSEMs
            if ((meth.mat[i,j] < (cpgStats[i, "q1"] - 3*cpgStats[i, "iq"]))) {
                hypo_sems_j <- hypo_sems_j + 1
                if (saveIndSEMs == TRUE) {
                    hypo_vector <- rbind(hypo_vector, rownames(meth.mat[i,]))
                }
            }

            # look for hyperSEMs
             if((meth.mat[i,j] > (cpgStats[i, "q3"] + 3*cpgStats[i, "iq"]))) {
                hyper_sems_j <- hyper_sems_j + 1
                if (saveIndSEMs == TRUE) {
                    hyper_vector <- rbind(hyper_vector, rownames(meth.mat[i,]))
               }
            }
        }
    }

        # bringing it all back home (aka getting ready to deliver the goods)
        total_sems_j <- hyper_sems_j + hypo_sems_j
        total_sems_j <- as.data.frame(total_sems_j)
        rownames(total_sems_j) <- sample_id
        colnames(total_sems_j) <- "Total_SEMs"
        hyper_sems_j <- as.data.frame(hyper_sems_j)
        rownames(hyper_sems_j) <- sample_id
        colnames(hyper_sems_j) <- "Hyper_SEMs"
        hypo_sems_j <- as.data.frame(hypo_sems_j)
        rownames(hypo_sems_j) <- sample_id
        colnames(hypo_sems_j) <- "Hypo_SEMs"
        sems_j <- cbind(hyper_sems_j, hypo_sems_j, total_sems_j)
        rownames(sems_j) <- sample_id
        sampleSEMs <- rbind(sampleSEMs, sems_j)
        if (saveIndSEMs == TRUE) {
            hyper_vector <- as.data.frame(hyper_vector)
            hypo_vector <- as.data.frame(hypo_vector)
            hyper_vector$type <- "hyper"
            hypo_vector$type <- "hypo"
            colnames(hyper_vector) <- c("SEMs", "type")
            colnames(hypo_vector) <- c("SEMs", "type")
            sems <- rbind(hyper_vector, hypo_vector)
           write.table(sems, file = paste0(saveDir, "/SEMs_", sample_id), quote = F, row.names = F, sep = "\t")
        }
    }	
    sampleSEMs$log10SEMs <- log10(sampleSEMs$Total_SEMs)
    return(sampleSEMs)
}
