#' importBedScore
#'
#' Given a set of peaks (GRanges object) and a vector
#' of one or more file paths of .bed files with some
#' score, import these and summarize using a FUN
#' over them. Also can specify which column the score 
#' should be read from in the file and the default
#' value per peak if nothing is determined. 
#' 
#' @param ranges A GRanges object corresponding to the
#' peaks used to aggregrate a score over
#' @param files A character vector of files that will
#' be imported and digested for analysis
#' @param colidx The column index of the score. This
#' assumes that the first three columns specific genomic 
#' coordinates. Default is 5 for the fifth column in the .bed files
#' @param FUN Function to summarize multiple hits in the
#' .bed file over the peak. 
#' @param default.val Default value to populate the matrices
#' ahead of time. By default, 0. 
#' @import GenomicRanges 
#' @import SummarizedExperiment
#' @importFrom Matrix Matrix sparseMatrix colSums
#' @importFrom utils read.table
#' @importFrom S4Vectors queryHits subjectHits DataFrame
#' @export 
#' @author Caleb Lareau
#' @examples
#'
#' files <- list.files(system.file('extdata',package='gchromVAR'), full.names = TRUE)
#' data(mini_counts, package = "chromVAR")
#' w_se <- importBedScore(SummarizedExperiment::rowRanges(mini_counts), files)
#' 
setGeneric("importBedScore",
   function(ranges, files, colidx = 5, FUN = sum, default.val = 0) standardGeneric("importBedScore"))

#' @describeIn importBedScore 
#' @export
setMethod("importBedScore", c(ranges = "GRanges", files = "character", colidx = "ANY", FUN = "ANY", default.val = "ANY"),
          function(ranges, files, colidx = 5, FUN = sum, default.val = 0){

  la <- sapply(files, function(file){
   
    # Import Bed 
    hitdf <- read.table(file)
    stopifnot(dim(hitdf)[2] >= colidx)
    stopifnot(3 < colidx)
    gdf <- hitdf[,c(1,2,3)]
    
    #Make GRanges
    names(gdf) <- c("chr", "start", "end")
    hitG <- makeGRangesFromDataFrame(gdf)
    
    # Get overlaps / summarized statistic summarized by function
    ov <- findOverlaps(ranges, hitG)
    score <- rep(default.val, length(ranges))
    ss <- tapply(hitdf[subjectHits(ov), colidx], queryHits(ov), FUN = FUN)
    score[as.integer(names(ss))] <- unname(ss)
    score
  })
  
  # Make a new SummarizedExperiment and export
  nnames <-  gsub(".bed", "", basename(files))
  colnames(la) <- nnames
  se <- SummarizedExperiment(assays=list(weights=Matrix(la)),
                     rowRanges=ranges,
                     colData=DataFrame(file = nnames))
  return(se)
})

