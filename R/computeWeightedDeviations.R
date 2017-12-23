#' computeWeightedDeviations
#'
#' Computes weighted deviations in chromatin accessibility
#' across sets of weights where the number of the annotation
#' meaningful.
#' 
#' @param object chromVARCounts object
#' @param weights chromVARweights object
#' @param background_peaks (optional) background peaks matrix computed using
#' \code{\link{getBackgroundPeaks}}, computed internally with default
#' paramaters if not provided
#' @param expectation (optional) expectations computed using
#' \code{\link{computeExpectations}}, computed automatically if not provided
#' @param ... additional arguments
#' @details multiprocessing using \code{\link[BiocParallel]{bplapply}}
#' @return  \code{\link{chromVARDeviations-class}}, which inherits from
#' SummarizedExperiment, and has two assays: deviations and deviation scores.
#' @import chromVAR
#' @import BiocParallel
#' @importFrom stats sd
#' @export
#' @author Caleb Lareau
#' @examples
#' # Register BiocParallel
#' BiocParallel::register(BiocParallel::SerialParam())
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' # Load mini weighted counts
#' rdsA<-paste0(system.file('rds',package='gchromVAR'),'/mini_w.rds')
#' w_se <- readRDS(rdsA)
#' 
#' # Build weights from .bed file
#' files <- list.files(system.file('extdata',package='gchromVAR'), full.names = TRUE, pattern = "*.bed$")
#' data(mini_counts, package = "chromVAR")
#' uk_se <- importBedScore(SummarizedExperiment::rowRanges(mini_counts), files)
#'
#' wdev <- computeWeightedDeviations(mini_counts, w_se)
#' ukdev <- computeWeightedDeviations(mini_counts, uk_se)
#' 
setGeneric("computeWeightedDeviations",
   function(object, weights, ...) standardGeneric("computeWeightedDeviations"))

#' @describeIn computeWeightedDeviations / weights are SummarizedExperiment
#' @export
setMethod("computeWeightedDeviations", c(object = "SummarizedExperiment",
                                         weights = "SummarizedExperiment"),
  function(object, weights, background_peaks = getBackgroundPeaks(object),
                   expectation = computeExpectations(object)) {
            
            cwdc(counts(object), getWeights(weights), background_peaks,
                  expectation, colData = colData(object), rowData = colData(weights))
          })


#' @describeIn computeWeightedDeviations object is SummarizedExperiment,
#' weights are in a Matrix
#' @export
setMethod("computeWeightedDeviations", c(object = "SummarizedExperiment", weights = "Matrix"),
  function(object, weights, background_peaks = getBackgroundPeaks(object),
                   expectation = computeExpectations(object)) {
      cwdc(counts(object), weights, background_peaks, expectation, colData = colData(object))
})

#' @describeIn computeWeightedDeviations object is SummarizedExperiment,
#' weights are in a numeric
#' @export
setMethod("computeWeightedDeviations", c(object = "SummarizedExperiment", weights = "matrix"),
  function(object, weights, background_peaks = getBackgroundPeaks(object),
                   expectation = computeExpectations(object)) {
      cwdc(counts(object), weights, background_peaks, expectation, colData = colData(object))
})

# Compute weighted deivations core
cwdc <- function(counts_mat, weights, background_peaks, expectation,
               rowData = NULL, colData = NULL) {

  if (min(getFragmentsPerPeak(counts_mat)) <= 0)
    stop("All peaks must have at least one fragment in one sample")

  stopifnot(nrow(counts_mat) == nrow(background_peaks))
  stopifnot(length(expectation) == nrow(counts_mat))

  if (is.null(colData))
    colData <- DataFrame(seq_len(ncol(counts_mat)),
                         row.names = colnames(counts_mat))[,FALSE]

  if (is.null(colnames(weights))) {
    names(weights) <- seq_along(weights)
  }

  sample_names <- colnames(counts_mat)

  results <- BiocParallel::bplapply(1:dim(weights)[2], function(i){
      gchromVAR:::cwds(as.numeric(weights[,i]), counts_mat = counts_mat,
        background_peaks = background_peaks, expectation = expectation)
  })

  z <- t(vapply(results, function(x) x[["z"]], rep(0, ncol(counts_mat))))
  dev <- t(vapply(results, function(x) x[["dev"]], rep(0, ncol(counts_mat))))
  
  if(dim(z)[1] == 1){
      z <- t(z)
      dev <- t(dev)
  }
  
  colnames(z) <- colnames(dev) <- sample_names
  
  out <- SummarizedExperiment(assays = list(deviations = dev, z = z),
                              colData = colData,
                              rowData = rowData)

  return(new("chromVARDeviations", out))
}


# Compute weighted deviations single
cwds <- function(quant, counts_mat, background_peaks, expectation = NULL, threshold = 1) {
      fragments_per_sample <- colSums(counts_mat)
  
  ### counts_mat should already be normed!
  vec <- Matrix(matrix(quant, ncol = length(quant)))
  
  observed <- as.vector(vec %*% counts_mat)
  expected <- as.vector(vec %*% expectation %*% fragments_per_sample)
  observed_deviation <- (observed - expected)/expected
  
  niterations <- ncol(background_peaks)
  sample_mat <- sparseMatrix(
        j = as.vector(background_peaks[quant != 0, seq_len(niterations)]),
        i = rep(seq_len(niterations), each = sum(quant != 0)),
        x = rep(quant[quant != 0], niterations),
        dims = c(niterations, nrow(counts_mat))
  )
  
  sampled <- as.matrix(sample_mat %*% counts_mat)
  sampled_expected <- as.matrix(sample_mat %*% expectation %*% fragments_per_sample)
  sampled_deviation <- (sampled - sampled_expected)/sampled_expected

  mean_sampled_deviation <- colMeans(sampled_deviation)
  sd_sampled_deviation <- apply(sampled_deviation, 2, sd)
  
  normdev <- (observed_deviation - mean_sampled_deviation)
  z <- normdev/sd_sampled_deviation
  return(list(z = z, dev = normdev))
}


