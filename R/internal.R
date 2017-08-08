
getWeights <- function(object){
    if ("weights" %in% assayNames(object)){
        out <- assays(object)$weights
    } else if ("scores" %in% assayNames(object)){
        out <- assays(object)$scores
    } else {
        stop("No appropriately named assay. See Details section in man page")
    }
    
}  
