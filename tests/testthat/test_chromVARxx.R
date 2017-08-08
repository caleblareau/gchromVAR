
context("bagging deviations works")

test_that("Bagging deviations reduces the feature space", {
    
    # Human
    rdsA<-paste0(system.file('rds',package='chromVARxx'),'/dev_humanSimple.rds')
    object <- readRDS(rdsA)
    bagged <- bagDeviations(object, cor = 0.3, organism = "human")
    expect_equal(dim(object)[1], 1764)
    expect_equal(dim(bagged)[1], 8)
    
    # Mouse
    rdsB<-paste0(system.file('rds',package='chromVARxx'),'/dev_mouseSimple.rds')
    object <- readRDS(rdsB)
    baggedm <- bagDeviations(object, cor = 0.3, organism = "mouse")
    expect_equal(dim(object)[1], 1346)
    expect_equal(dim(baggedm)[1], 5)
})

context("Caleb's weighted chromVAR and import functions")

test_that("Can import bed score files", {
    
    files <- list.files(system.file('extdata',package='chromVARxx'), full.names = TRUE)
    data(mini_counts, package = "chromVAR")
    w_se <- importBedScore(rowRanges(mini_counts), files)
    expect_equal(dim(w_se)[1], 1000)
    expect_equal(dim(w_se)[2], 2)
    
})

test_that("Can compute deviations and that some cells can be scored", {
    
    files <- list.files(system.file('extdata',package='chromVARxx'), full.names = TRUE)
    data(mini_counts, package = "chromVAR")
    w_se <- importBedScore(rowRanges(mini_counts), files)
    ukdev <- computeWeightedDeviations(mini_counts, w_se)
    expect_gt(sum(abs(assays(ukdev)[["z"]]), na.rm = TRUE), 0)
    
})
