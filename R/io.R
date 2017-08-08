#' importBedScore
#'
#' Given a set of peaks (GRanges object) and a vector
#' of one or more file paths of .bed files with some
#' score, import these and summarize using a FUN
#' over them 
#' 
#' @param ranges A GRanges object corresponding to the
#' peaks used to aggre 
#' @param files 
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import Matrix
#' @export 
#' @author Caleb Lareau
#' @examples
#'
#' files <- list.files(system.file('extdata',package='chromVARxx'), full.names = TRUE)
#' data(mini_counts, package = "chromVAR")
#' w_se <- importBedScore(rowRanges(mini_counts), files)
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
  
  nnames <-  gsub(".bed", "", basename(files))
  colnames(la) <- nnames
  se <- SummarizedExperiment(assays=list(weights=Matrix(la)),
                     rowRanges=ranges,
                     colData=DataFrame(file = nnames))
  return(se)
})

