#' computeWeightedDeviations
#'
#' Computes weighted deviations in chromatin accessibility
#' across sets of annotations where the number of the annotation
#' meaningful.
#' 
#' @param object chromVARCounts object
#' @param annotations chromVARAnnotations object
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
#' @export
#' @author Caleb Lareau
#' @examples
#' # Register BiocParallel
#' BiocParallel::register(BiocParallel::SerialParam())
#' # Load very small example counts (already filtered)
#' data(mini_counts, package = "chromVAR")
#' # load annotation matrix; result from matchMotifs
#' data(mini_ix, package = "chromVAR")
#'
#' # computing deviations
#' dev <- computeWeightedDeviations(object = mini_counts,
#'                          annotations = mini_ix)
setGeneric("computeWeightedDeviations",
   function(object, annotations, ...) standardGeneric("computeWeightedDeviations"))

#' @describeIn computeWeightedDeviations / annotations are SummarizedExperiment
#' @export
setMethod("computeWeightedDeviations", c(object = "SummarizedExperiment",
                                         annotations = "SummarizedExperiment"),
  function(object, annotations, background_peaks = getBackgroundPeaks(object),
                   expectation = computeExpectations(object)) {
            object <- counts_check(object)
            peak_indices <- convert_to_ix_list(annotationMatches(annotations))
            compute_weighted_deviations_core(counts(object),
                                    peak_indices,
                                    background_peaks,
                                    expectation,
                                    colData = colData(object),
                                    rowData = colData(annotations))
          })


#' @describeIn computeWeightedDeviations object is SummarizedExperiment,
#' annotations are list
#' @export
setMethod("computeWeightedDeviations",
          c(object = "SummarizedExperiment", annotations = "list"),
          function(object, annotations,
                   background_peaks = getBackgroundPeaks(object),
                   expectation = computeExpectations(object)) {
            object <- counts_check(object)
            compute_weighted_deviations_core(counts(object),
                                    annotations,
                                    background_peaks,
                                    expectation,
                                    colData = colData(object))
          })

compute_weighted_deviations_core <- function(counts_mat,
                                    peak_indices,
                                    background_peaks,
                                    expectation,
                                    rowData = NULL,
                                    colData = NULL) {

  if (min(getFragmentsPerPeak(counts_mat)) <= 0)
    stop("All peaks must have at least one fragment in one sample")

  stopifnot(nrow(counts_mat) == nrow(background_peaks))
  stopifnot(length(expectation) == nrow(counts_mat))

  if (is.null(colData))
    colData <- DataFrame(seq_len(ncol(counts_mat)),
                         row.names = colnames(counts_mat))[,FALSE]

  # check that indices fall within appropriate bounds
  tmp <- unlist(peak_indices, use.names = FALSE)
  if (is.null(tmp) ||
      !(all.equal(tmp, as.integer(tmp))) ||
      max(tmp) > nrow(counts_mat) ||
      min(tmp) < 1) {
    stop("Annotations are not valid")
  }

  if (is.null(names(peak_indices))) {
    names(peak_indices) <- seq_along(peak_indices)
  }

  sample_names <- colnames(counts_mat)

  results <- bplapply(peak_indices,
                      compute_deviations_single,
                      counts_mat = counts_mat,
                      background_peaks = background_peaks,
                      expectation = expectation)

  z <- t(vapply(results, function(x) x[["z"]], rep(0, ncol(counts_mat))))
  dev <- t(vapply(results, function(x) x[["dev"]], rep(0, ncol(counts_mat))))

  colnames(z) <- colnames(dev) <- sample_names

  #rowData$fractionMatches <- vapply(results, function(x) x[["matches"]], 0)
  #rowData$fractionBackgroundOverlap <- vapply(results, function(x) x[["overlap"]], 0)
  
  out <- SummarizedExperiment(assays = list(deviations = dev, z = z),
                              colData = colData,
                              rowData = rowData)

  return(new("chromVARDeviations", out))
}

compute_weighted_deviations_single <- function(peak_set,
                                      counts_mat,
                                      background_peaks,
                                      expectation = NULL,
                                      intermediate_results = FALSE,
                                      threshold = 1) {

  if (length(peak_set) == 0) {
    return(list(z = rep(NA, ncol(counts_mat)), dev = rep(NA, ncol(counts_mat))))
  }

  fragments_per_sample <- colSums(counts_mat)

  ### counts_mat should already be normed!
  tf_count <- length(peak_set)

  if (tf_count == 1) {
    observed <- as.vector(counts_mat[peak_set, ])
    expected <- expectation[peak_set] * fragments_per_sample
    observed_deviation <- (observed - expected)/expected

    sampled <- counts_mat[background_peaks[peak_set, ], ]
    sampled_expected <- outer(expectation[background_peaks[peak_set, ]],
                              fragments_per_sample)
    sampled_deviation <- (sampled - sampled_expected)/sampled_expected
    bg_overlap <- sum(background_peaks[peak_set,] == peak_set[1]) / 
      ncol(background_peaks)
  } else {
    tf_vec <- sparseMatrix(j = peak_set,
                           i = rep(1, tf_count),
                           x = 1,
                           dims = c(1,
                                    nrow(counts_mat)))
    
    observed <- as.vector(tf_vec %*% counts_mat)
    
    expected <- as.vector(tf_vec %*% expectation %*% fragments_per_sample)
    observed_deviation <- (observed - expected)/expected
    
    niterations <- ncol(background_peaks)
    sample_mat <-
      sparseMatrix(j = 
                     as.vector(background_peaks[peak_set, 
                                                seq_len(niterations)]),
                   i = rep(seq_len(niterations), each = tf_count),
                   x = 1,
                   dims = c(niterations,
                            nrow(counts_mat)))
    
    sampled <- as.matrix(sample_mat %*% counts_mat)
    sampled_expected <- as.matrix(sample_mat %*% expectation %*% 
                                    fragments_per_sample)
    sampled_deviation <- (sampled - sampled_expected)/sampled_expected

    bg_overlap <- mean(sample_mat %*% t(tf_vec)) / tf_count
  }

  fail_filter <- which(expected < threshold)

  mean_sampled_deviation <- colMeans(sampled_deviation)
  sd_sampled_deviation <- apply(sampled_deviation, 2, sd)

  normdev <- (observed_deviation - mean_sampled_deviation)
  z <- normdev/sd_sampled_deviation

  if (length(fail_filter) > 0) {
    z[fail_filter] <- NA
    normdev[fail_filter] <- NA
  }

  if (intermediate_results) {
    out <- list(z = z, dev = normdev,
                observed = observed,
                sampled = sampled,
                expected = expected,
                sampled_expected = sampled_expected,
                observed_deviation = observed_deviation,
                sampled_deviation = sampled_deviation,
                matches = tf_count / nrow(counts_mat),
                overlap = bg_overlap)
  } else {
    out <- list(z = z, dev = normdev, matches = tf_count / nrow(counts_mat),
                overlap = bg_overlap)
  }
  return(out)
}


