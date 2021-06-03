
#' @title Construct sets from objects
#' @description This function allows the construction of taxa sets of
#'    \code{\linkS4class{BiocSet}} type from either
#'    \code{\linkS4class{taxonomyTable}},
#'    \code{matrix}, or \code{\link[tibble]{tibble}}
#' @param obj Object to convert to standard sets to be used
#' @param ... Additional arguments to be passed
#' @export
setGeneric("const_set", function(obj, ...) standardGeneric("const_set"), signature = "obj")

#' @describeIn const_set Convert membership matrix to \code{BiocSet}
setMethod("const_set", "matrix", function(obj){
    return(0)
})

#' @describeIn const_set convert taxonomyTable to \code{BiocSet}
#' @param rank (Character). Restrict sets to certain taxonomic ranks.
setMethod("const_set", "taxonomyTable", function(obj, rank=NULL){
    if (is.null(rank)){
        message("Using all taxa ranks to construct all possible sets")
    }
    return(0)
})

#' @describeIn const_set Convert \code{data.frame} or \code{tibble} to
#'     \code{BiocSet}
setMethod("const_set", c("data.frame", "tbl_df"), function(obj){
    return(0)
})

