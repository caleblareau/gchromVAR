
context("Weighted chromVAR and import functions")

test_that("Can import bed score files", {
    
    files <- list.files(system.file('extdata',package='gchromVAR'), full.names = TRUE, pattern = "*.bed$")
    data(mini_counts, package = "chromVAR")
    w_se <- importBedScore(rowRanges(mini_counts), files)
    expect_equal(dim(w_se)[1], 1000)
    expect_equal(dim(w_se)[2], 2)
    
})

test_that("Can compute deviations and that some cells can be scored", {
    
    files <- list.files(system.file('extdata',package='gchromVAR'), full.names = TRUE, pattern = "*.bed$")
    data(mini_counts, package = "chromVAR")
    w_se <- importBedScore(rowRanges(mini_counts), files)
    ukdev <- computeWeightedDeviations(mini_counts, w_se)
    expect_gt(sum(abs(assays(ukdev)[["z"]]), na.rm = TRUE), 0)
    
})
