#' @name MBASIC.pipeline
#' @title The pipeline for fitting a MBASIC model for sequencing data.
#' @param chipfile A string vector for the ChIP files.
#' @param inputfile A string vector for the matching input files. The length must be the same as 'chipfile'.
#' @param input.suffix A string for the suffix of input files. If \code{NULL}, \code{inputfile} will be treated as the full names of the input files. Otherwise, all inputfiles with the initial \code{inputfile} and this suffix will be merged.
#' @param target A GenomicRanges object for the target intervals where the reads are mapped.
#' @param chipformat A string specifying the type of the ChIP file. Currently two file types are allowed: 'BAM' or 'BED'. Default: 'BAM'.
#' @param inputformat A string specifying the type of the input files. Currently two file types are allowed: 'BAM' or 'BED'. Default: 'BAM'.
#' @param fragLen Either a single value or a 2-column matrix of the fragment lengths for the chip and input files.  Default: 150.
#' @param pairedEnd Either a boolean value or a 2-column boolean matrix for whether each file is a paired-end data set. Currently this function only allows 'BAM' files for paired-end data. Default: FALSE.
#' @param unique A boolean value for whether only reads with distinct genomic coordinates or strands are mapped. Default: TRUE.
#' @param m.prefix A string for the prefix of the mappability files.
#' @param m.suffix A string for the suffix of the mappability files. See details for more information. Default: NULL.
#' @param gc.prefix A string for the prefix of the GC files.
#' @param gc.suffix A string for the suffix of the GC files. See details for more information. Default: NULL.
#' @param datafile The file location to save or load the data matrix. See details.
#' @param J The number of clusters to be fitted.
#' @param method The fitting algorithm, either "mbasic" (default), which runs an Expectation-Maximization algorithm, or "madbayes", which runs a K-means-like algorithm.
#' @param ... Parameters for \code{\link{MBASIC}} (if method is "mbasic") or \code{\link{MBASIC.MADBayes.full}}.
#' @details
#' This function executes three steps:\cr
#' The first step uses the \code{\link{generateReadMatrices}} function to get the ChIP and Input counts for each locus.\cr
#' The second step is to compute the covariate matrix. If any of \code{m.prefix}, \code{m.suffix}, \code{gc.prefix}, \code{gc.suffix} is NULL, then the input count matrix is directly used as the covariate matrix for MBASIC. Alternatively, it will use the \code{\link{bkng_mean}} to normalize the input count data according to the mappability and GC scores to produce the covariate matrix.\cr
#' The final step is to call the MBASIC function for model fitting.\cr
#' Because the first two steps are time consuming, we recommend in specifying a file location for '\code{datafile}. Then, when this function executes, it first checks whether \code{datafile} exists. If it exists, it will be loaded and the function will jump to the final step. If it does not exist, after the function executes the first two steps, the ChIP data matrix and the covariate matrix will be saved to this file, so that when you rerun this function you do not need to repeat the first two steps.\cr
#' @return If \code{method="mbasic"}, the return value by function \code{\link{MBASIC}} (if 'J' is scalar) or \code{\link{MBASIC.full}} (if 'J' is a vector). If \code{method="madbayes"}, the return value by function \code{\link{MBASIC.MADBayes.full}}.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' \dontrun{
#' ## This is the example in our vignette
#' target <- generateSyntheticData(dir = "syntheticData")
#' tbl <- ChIPInputMatch(dir = paste("syntheticData/", c("chip", "input"), sep = ""),suffix = ".bed", depth = 5)
#' conds <- paste(tbl$cell, tbl$factor, sep = ".")
#' MBASIC.fit <- MBASIC.pipeline(chipfile = tbl$chipfile, inputfile = tbl$inputfile, input.suffix = ".bed", target = target, format = "BED", fragLen = 150, pairedEnd = FALSE, unique = TRUE, fac = conds, struct = NULL, S = 2, J = 3, family = "lognormal", maxitr = 10, statemap = NULL)
#'}
#' @export
MBASIC.pipeline <- function(chipfile, inputfile, input.suffix, target, chipformat, inputformat, fragLen, pairedEnd, unique, m.prefix = NULL, m.suffix = NULL, gc.prefix = NULL, gc.suffix = NULL, datafile = NULL, ncores = 10, J, method = "mbasic", ...) {

  if(is.null(datafile) | (!is.null(datafile) & !file.exists(datafile))) {
    dat <- generateReadMatrices(chipfile = chipfile, inputfile = inputfile, input.suffix = input.suffix, target = target, chipformat = chipformat, inputformat = inputformat, fragLen = fragLen, pairedEnd = pairedEnd, unique = unique)
    if(!is.null(m.prefix) & !is.null(m.suffix) & !is.null(gc.prefix) & !is.null(gc.suffix)) {
      target <- averageMGC(target = target, m.prefix = m.prefix, m.suffix = m.suffix, gc.prefix = gc.prefix, gc.suffix = gc.suffix)
      Gamma <- bkng_mean(inputdat = dat$input, target = target, family = family)
    } else {
      Gamma <- t(dat$input)
    }
    if(!is.null(datafile)) {
      save(dat, Gamma, file = datafile)
    }
  } else {
    load(datafile)
  }

  if(method == "mbasic") {
    if(length(J) == 1) {
      return(MBASIC(Y = t(dat$chip), Gamma = Gamma, J = J, ...))
    } else {
      return(MBASIC.full(Y = t(dat$chip), Gamma = Gamma, J = J, ...))
    }
  } else if(method == "madbayes") {
    return(MBASIC.MADBayes.full(Y = log(t(dat$chip) + 1), Gamma = log(Gamma + 1), ...))
  }
}
