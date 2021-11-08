#' @title Construct sets from objects
#' @description This function allows the construction of taxa sets of
#'    \code{\linkS4class{BiocSet}} type from either
#'    \code{\linkS4class{taxonomyTable}},
#'    \code{matrix}, or \code{\link[tibble]{tibble}}
#' @param obj Object to convert to standard sets to be used
#' @param ... Additional arguments to be passed
#' @export
setGeneric("const_set", function(obj, ...) standardGeneric("const_set"), signature = "obj")

#' @describeIn const_set convert \code{phyloseq} to \code{BiocSet}
#' @param rank Character. Restrict sets to certain taxonomic ranks.
#' @importClassesFrom phyloseq phyloseq
#' @importFrom phyloseq tax_table
#' @export
setMethod("const_set", "phyloseq", function(obj, rank){
    # Avoiding no binding notes when doing R CMD CHECK for metaprogramming
    n <- set <- NULL
    table <- phyloseq::tax_table(obj)
    if (!any(colnames(table) == rank)){
        rlang::abort("Rank name not part of taxonomy table")
    }
    sets <- .const_set_taxtable(table = table, rank = rank)
    return(sets)
})

#' @describeIn const_set Convert taxonomic tables from TreeSummarizedExperiment to BiocSet
#' @param rank String. Restrict sets to certain taxonomic ranks.
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @import TreeSummarizedExperiment
#' @importFrom SummarizedExperiment rowData
#' @export
setMethod("const_set", "TreeSummarizedExperiment", function(obj, rank){
  #TODO: Replace abun with object
  table <- SummarizedExperiment::rowData(obj)
  if (!rank %in% colnames(table)){
    stop("Specified rank is not part of the taxonomic table")
  }
  sets <- .const_set_taxtable(table = table, rank = rank)
  return(sets)
})


#' @describeIn const_set Convert membership (or logical) matrix to \code{BiocSet}
setMethod("const_set", "matrix", function(obj){
  if (!all(obj %in% c(0,1))){
    warning("Coercing to data.frame type and apply that method instead")
    obj <- as.data.frame(obj)
    return(const_set(obj))
  }
  return("PLACEHOLDER")
})


#' @describeIn const_set Convert \code{data.frame} membership or \code{tibble} to
#'     \code{BiocSet}
setMethod("const_set", "data.frame", function(obj){
    return("PLACEHOLDER")
})


#' @title Internal function to construct from a tax table converted to data.frame or matrix formats
#' @param table DataFrame.The taxonomic table extracted from relevant objects
#' @param rank String. The rank to extract sets from
#' @importFrom purrr map
#' @importFrom BiocSet BiocSet es_elementset left_join_set
#' @importFrom stringr str_ends
#' @importFrom magrittr %>%
#' @importFrom dplyr pull group_by
#' @importFrom stats na.omit
#' @importFrom lobstr obj_size
#' @keywords internal
.const_set_taxtable <- function(table, rank){
  id <- which(colnames(table) == rank)
  table <- as.data.frame(as(table, "matrix")[,1:id])
  # Get all names -- have to do this because sometimes different phyla might have the same
  # genus or family.
  all_names <- apply(table, 1, function(i){
    paste(i, sep = ";_;", collapse = ";_;")
  })
  # Getting all unique names by removing NAs
  unq_names <- stats::na.omit(unique(all_names))
  unq_names <- unq_names[!stringr::str_ends(unq_names, "NA")]
  # Get all set ids
  sets <- purrr::map(unq_names, ~{
    names(all_names)[which(all_names == .x)]
  })
  names(sets) <- unq_names
  # Constructing sets
  sets <- BiocSet::BiocSet(sets)
  set_sizes <- BiocSet::es_elementset(sets) %>%
    dplyr::count(set, name = "size")
  # filter set by set sizes
  sets <- BiocSet::left_join_set(sets, set_sizes)
  if (lobstr::obj_size(sets)/1e6 > 100){
    rlang::warn("Object size is larger than 100MB")
  }
  return(sets)
}
