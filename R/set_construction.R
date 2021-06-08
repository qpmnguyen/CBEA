#' @title Construct sets from objects
#' @description This function allows the construction of taxa sets of
#'    \code{\linkS4class{BiocSet}} type from either
#'    \code{\linkS4class{taxonomyTable}},
#'    \code{matrix}, or \code{\link[tibble]{tibble}}
#' @param obj Object to convert to standard sets to be used
#' @param min_size The minimum size of a set
#' @param ... Additional arguments to be passed
#' @export
setGeneric("const_set", function(obj, min_size, ...) standardGeneric("const_set"), signature = "obj")

#' @describeIn const_set convert taxonomyTable to \code{BiocSet}
#' @param rank (Character). Restrict sets to certain taxonomic ranks.
#' @importFrom purrr map
#' @importFrom BiocSet BiocSet es_elementset mutate_set filter_set
#' @importFrom stringr str_ends
#' @importFrom magrittr %>%
#' @importFrom dplyr pull group_by
#' @export
setMethod("const_set", "taxonomyTable", function(obj, min_size, rank){
    if (!any(colnames(obj) == rank)){
        rlang::abort("Rank name not part of taxonomy table")
    }
    id <- which(colnames(obj) == rank)
    obj <- as.data.frame(as(obj, "matrix")[,1:id])
    # Get all names -- have to do this because sometimes different phyla might have the same
    # genus or family.
    all_names <- apply(obj, 1, function(i){
         paste(i, sep = ";_;", collapse = ";_;")
    })
    # Getting all unique names by removing NAs
    unq_names <- na.omit(unique(all_names))
    unq_names <- unq_names[!stringr::str_ends(unq_names, "NA")]
    # Get all set ids
    sets <- purrr::map(unq_names, ~{
        names(all_names)[which(all_names == .x)]
    })
    names(sets) <- unq_names
    # Constructing sets
    sets <- BiocSet::BiocSet(sets)
    set_sizes <- BiocSet::es_elementset(sets) %>%
        dplyr::count(set) %>% dplyr::pull(n)
    # filter set by set sizes
    sets <- BiocSet::mutate_set(sets, size = set_sizes) %>%
        BiocSet::filter_set(size >= min_size)
    if (lobstr::obj_size(sets)/1e6 > 100){
        rlang::warn("Object size is larger than 100MB")
    }
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


#' TODO: Convert set to different formats
#' This is a placeholder function
conv_set <- function(format=NULL){
    return(0)
}

#' Trim sets based on criteria
#' @param physeq (\code{phyloseqSet}). This object has taxon sets that need to be trimmed
#' @param ... Trim criteria that will be used
trim_set <- function(physeq, ...){
    return(0)
}
