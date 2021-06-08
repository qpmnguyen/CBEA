#' @title Construct sets from objects
#' @description This function allows the construction of taxa sets of
#'    \code{\linkS4class{BiocSet}} type from either
#'    \code{\linkS4class{taxonomyTable}},
#'    \code{matrix}, or \code{\link[tibble]{tibble}}
#' @param obj Object to convert to standard sets to be used
#' @param ... Additional arguments to be passed
#' @export
setGeneric("const_set", function(obj, ...) standardGeneric("const_set"), signature = "obj")

#' @describeIn const_set convert taxonomyTable to \code{BiocSet}
#' @param rank (Character). Restrict sets to certain taxonomic ranks.
#' @importFrom purrr map
#' @importFrom BiocSet BiocSet es_elementset mutate_set
#' @importFrom stringr str_ends
#' @importFrom magrittr %>%
#' @importFrom dplyr pull group_by
#' @importFrom stats na.omit
#' @export
setMethod("const_set", "taxonomyTable", function(obj, rank){
    # Avoiding no binding notes when doing R CMD CHECK for metaprogramming
    n <- set <- NULL
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
        dplyr::count(set) %>% dplyr::pull(n)
    # filter set by set sizes
    sets <- BiocSet::mutate_set(sets, size = set_sizes)
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
#' @title This is a placeholder function
#' @param format Format to convert sets to
conv_set <- function(format=NULL){
    return(0)
}

#' @title Trim sets based on user supplied criteria
#' @description A filter operation specific to phyloseqSet type objects for sets. Using the criteria that appears
#'    in the \code{set} slot of the \code{BiocSet} object, users can filter out both elements, sets,
#'    and element-sets. Usage of this function is very similar to \code{\link[dplyr]{filter}} function in the \code{dplyr}
#'    package.
#' @param physeq (\code{phyloseqSet}). This object has taxon sets that need to be trimmed
#' @param ... (Expression). Trim criteria that would be used. For example \code{size >= 5}.
#' @return A \code{phyloseqSet} object with sets and elements appropriately trimmed and filtered
#' @export
#' @importFrom dplyr pull
#' @importFrom BiocSet filter_set filter_element es_element
trim_set <- function(physeq, ...){
    # Avoid notes for no global bindings in R CMD CHECK when doing metaprogramming
    element <- . <- NULL
    set <- taxon_set(physeq)
    set <- BiocSet::filter_set(set, ...) %>%
        BiocSet::filter_element(element %in% dplyr::pull(BiocSet::es_elementset(.), element))
    taxon_set(physeq) <- set
    return(physeq)
}

