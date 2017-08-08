library(chromVAR)

# Input a standard summarized experiment object
# and a numeric vector of weights

computeDeviationsQuant <- function(object = "SummarizedExperiment", quant = "numeric"){
  
  chromVAR:::getBackgroundPeaks(object)
  object <- chromVAR:::counts_check(object)
  background_peaks <- chromVAR:::getBackgroundPeaks(object)
  expectation <- chromVAR:::computeExpectations(object)
  
  compute_deviations_core_quant(counts(object), quant,  background_peaks, expectation)
}

compute_deviations_core_quant <- function(counts_mat, quant, background_peaks, expectation) {
  
  if (min(getFragmentsPerPeak(counts_mat)) <= 0)
    stop("All peaks must have at least one fragment in one sample")
  
  stopifnot(nrow(counts_mat) == nrow(background_peaks))
  stopifnot(length(expectation) == nrow(counts_mat))
  
  results <- compute_deviations_single_quant(counts_mat, quant, background_peaks, expectation)
  
}


compute_deviations_single_quant <- function(counts_mat, quant, background_peaks, expectation = NULL) {
  
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
  return(z)
}

