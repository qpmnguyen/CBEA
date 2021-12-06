# CONST_SETS ------------------------------------------------------------####
#' @title Construct sets from objects
#' @description This function allows the construction of taxa sets of
#'    \code{\linkS4class{BiocSet}} type from either
#'    \code{\linkS4class{taxonomyTable}},
#'    \code{matrix}, or \code{\link[tibble]{tibble}}
#' @param obj Object to convert to standard sets to be used
#' @param ... Additional arguments to be passed
#' @return A \code{\linkS4class{BiocSet}} object
#' @export
#' @import methods
#' @examples
#' data(hmp_gingival)
#' seq <- hmp_gingival$data
#' # Constructing set using taxonomic table
#' set_phylum <- const_set(obj = seq, rank = "PHYLUM")
setGeneric("const_set", function(obj, ...) standardGeneric("const_set"),
           signature = "obj")

#' @describeIn const_set Convert taxonomic tables from
#'     \code{\linkS4class{phyloseq}} to \code{\linkS4class{BiocSet}}
#'     from a certain rank.
#' @param rank Character. Restrict sets to certain taxonomic ranks.
#' @importClassesFrom phyloseq phyloseq
#' @importFrom phyloseq tax_table
#' @import methods
#' @export
setMethod("const_set", "phyloseq", function(obj, rank) {
    # Avoiding no binding notes when doing R CMD CHECK for metaprogramming
    n <- set <- NULL
    table <- phyloseq::tax_table(obj)
    if (!any(colnames(table) == rank)) {
        rlang::abort("Rank name not part of taxonomy table")
    }
    sets <- .const_set_taxtable(table = table, rank = rank)
    return(sets)
})

#' @describeIn const_set Convert taxonomic tables
#'     from \code{\linkS4class{TreeSummarizedExperiment}}
#'     to \code{\linkS4class{BiocSet}}
#'     from a certain rank.
#' @param rank String. Restrict sets to certain taxonomic ranks.
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @import TreeSummarizedExperiment
#' @importFrom SummarizedExperiment rowData
#' @import methods
#' @export
setMethod("const_set", "TreeSummarizedExperiment", function(obj, rank) {
    table <- SummarizedExperiment::rowData(obj)
    if (!rank %in% colnames(table)) {
        stop("Specified rank is not part of the taxonomic table")
    }
    sets <- .const_set_taxtable(table = table, rank = rank)
    return(sets)
})


#' @describeIn const_set Convert membership (or logical)
#'     matrix to \code{BiocSet}.
#'     This matrix is a matrix of 1 or 0s indicating membership of
#'     taxa (as row or column names)
#'     with relevant sets (as row or column names)
#' @param sets_are_rows (Logical). Indicate whether sets are
#'     rows or as columns. Affect set names.
#' @importClassesFrom BiocSet BiocSet
#' @importFrom BiocSet BiocSet
setMethod("const_set", "matrix", function(obj, sets_are_rows = FALSE) {
    if (!all(obj %in% c(0, 1))) {
        stop("Please convert element-set membership matrix to
             data.frame and specify id and member columns.
             See ?const_set for details")
    }
    if (sets_are_rows == TRUE) {
        element_names <- colnames(sets)
        set_names <- rownames(sets)
        direction <- 1
    } else {
        element_names <- rownames(sets)
        set_names <- colnames(sets)
        direction <- 2
    }
    set_list <- apply(obj, direction, function(x) {
        element_names[which(x == 1)]
    })
    if (is(set_list, "list") != TRUE) {
        stop("Sets needs to be a list")
    }
    if (length(set_list) != set_names) {
        stop("List does not have the same number of elements as possible names")
    }
    names(set_list) <- set_names
    sets <- BiocSet::BiocSet(set_list)
    return(sets)
})


#' @describeIn const_set Convert \code{data.frame} element set or
#'     \code{tibble} element set
#'     to \code{BiocSet}. We expect at least a 2-column
#'     data frame with the taxa names
#'     and the set names as the set-membership identifiers.
#' @param id String. Indicates which columns are taxa ids
#' @param member String. Indicate which columns connect taxa ids to set names
#' @importClassesFrom BiocSet BiocSet
#' @importFrom BiocSet BiocSet
#' @importFrom dplyr pull filter
#' @importFrom magrittr %>%
#' @export
setMethod("const_set", "data.frame", function(obj, id, member) {
    # if just one or zero then dispatch to the matrix format
    if (all(obj %in% c(0, 1))) {
        warning("This looks like a logical dummy variable
                matrix indicating membership. Dispatching to the matrix signature")
        obj <- as.matrix(obj)
        sets <- const_set(obj)
    }
    set_names <- obj %>%
        dplyr::pull(member) %>%
        unique()
    # convert to sym to utilize dplyr programming
    member <- rlang::sym(member)
    set_list <- vector(mode = "list", length = length(set_names))
    for (i in seq_along(set_names)) {
    set_list[[i]] <- obj %>%
        dplyr::filter({{ member }} == set_names[i]) %>%
        dplyr::pull(id)
    }
    names(set_list) <- set_names
    sets <- BiocSet::BiocSet(set_list)
    return(sets)
})


#' @title Internal function to construct from a tax table
#'     converted to data.frame or matrix formats
#' @param table DataFrame.The taxonomic table extracted from relevant objects
#' @param rank String. The rank to extract sets from
#' @importFrom purrr map
#' @importFrom BiocSet BiocSet es_elementset left_join_set
#' @importFrom stringr str_ends
#' @importFrom magrittr %>%
#' @importFrom dplyr pull group_by
#' @importFrom stats na.omit
#' @importFrom lobstr obj_size
#' @return An object of type \code{\linkS4class{BiocSet}}
#' @keywords internal
.const_set_taxtable <- function(table, rank) {
    set <- NULL # R CMD NOTE
    id <- which(colnames(table) == rank)
    table <- as.data.frame(as(table, "matrix")[, seq_len(id)])
    # Get all names -- have to do this because sometimes
    # different phyla might have the same
    # genus or family.
    all_names <- apply(table, 1, function(i) {
        paste(i, sep = ";_;", collapse = ";_;")
    })
    # Getting all unique names by removing NAs
    unq_names <- stats::na.omit(unique(all_names))
    unq_names <- unq_names[!stringr::str_ends(unq_names, "NA")]
    # Get all set ids
    sets <- purrr::map(unq_names, ~ {
        names(all_names)[which(all_names == .x)]
    })
    names(sets) <- unq_names
    # Constructing sets
    sets <- BiocSet::BiocSet(sets)
    set_sizes <- BiocSet::es_elementset(sets) %>%
        dplyr::count(set, name = "size")
    # filter set by set sizes
    sets <- BiocSet::left_join_set(sets, set_sizes)
    if (lobstr::obj_size(sets) / 1e6 > 100) {
        rlang::warn("Object size is larger than 100MB")
    }
    return(sets)
}

# UNIFY SETS --------------------------------------------------------------####
#' @title Unify sets with data containers
#' @description This function checks for consistencies in naming between
#'     the data container and the set itself.
#' @details Extracts out elements from a \code{\linkS4class{BiocSet}}
#'     and then compares it against taxa names from the data container
#'     (either \code{\linkS4class{phyloseq}} or
#'     \code{\linkS4class{TreeSummarizedExperiment}}).
#'     It then removes elements from the \code{\linkS4class{BiocSet}} that are
#'     not in the taxa_names of the data container.
#'     Throws an error if there are no elements in BiocSet that matches
#'     taxa names of the data container.
#' @export
#' @param obj Data container object containing taxonomic
#'     relative abundance data with taxa and sample names
#' @param set (BiocSet) The set being tested against
#' @param ... Additional arguments to be passed
#' @return An object of type \code{\linkS4class{BiocSet}}
#' @import methods
#' @examples
#' data(hmp_gingival)
#' set <- hmp_gingival$set
#' seq <- hmp_gingival$data
#' # adding some random elements in some set
#' set_add <- dplyr::union(set, BiocSet::BiocSet(met1 = letters[c(1:3)],
#'     met2 = letters[c(5:7)]))
#' set_new <- unify_sets(obj = seq, set = set_add)
setGeneric("unify_sets",
    function(obj, set, ...) standardGeneric("unify_sets"),
    signature = "obj"
)

#' @describeIn unify_sets Unify set and data
#'     container for \code{\linkS4class{TreeSummarizedExperiment}}
#' @importFrom BiocSet filter_element
#' @export
setMethod("unify_sets", "TreeSummarizedExperiment", function(obj, set) {
    element <- NULL # R CMD CHECK
    ref_names <- rownames(obj)
    n_set <- BiocSet::filter_elementset(set, element %in% ref_names)
    return(n_set)
})

#' @describeIn unify_sets Unify set and data
#'     container for \code{\linkS4class{phyloseq}}
#' @importFrom phyloseq taxa_names
#' @importFrom BiocSet filter_element
#' @export
setMethod("unify_sets", "phyloseq", function(obj, set) {
    element <- NULL # R CMD CHECK
    ref_names <- phyloseq::taxa_names(obj)
    n_set <- BiocSet::filter_elementset(set, element %in% ref_names)
    return(n_set)
})

#' @describeIn unify_sets Unify set and data container for \code{data.frame}
#' @param taxa_are_rows (Logical). Indicate whether taxa
#'     are stored as rows or columns of the data.frame/matrix.
#' @param tax_ids (String). Indicate column name with all the names of taxa.
#'     If NULL, rownames will be used.
#'     Only used if \code{taxa_are_rows = TRUE}.
#' @importFrom magrittr %>%
#' @importFrom dplyr pull
#' @importFrom rlang sym
#' @export
setMethod("unify_sets", "data.frame", function(obj, set,
                                               tax_ids = NULL,
                                               taxa_are_rows = TRUE) {
    element <- NULL # R CMD CHECK
    if (taxa_are_rows == TRUE) {
        if (is.null(tax_ids)) {
            ref_names <- rownames(obj)
        } else {
            tax_ids <- rlang::sym(tax_ids)
            ref_names <- obj %>% dplyr::pull({{ tax_ids }})
        }
    } else {
        ref_names <- colnames(obj)
    }
    n_set <- BiocSet::filter_elementset(set, element %in% ref_names)
    return(n_set)
})

#' @describeIn unify_sets Unify set and data container
#'     for \code{\linkS4class{matrix}}
#' @importFrom BiocSet filter_elementset
#' @export
setMethod("unify_sets", "matrix", function(obj, set, taxa_are_rows) {
    element <- NULL # R CMD CHECK
    if (taxa_are_rows == TRUE) {
        ref_names <- rownames(obj)
    } else {
        ref_names <- colnames(obj)
    }
    n_set <- BiocSet::filter_elementset(set, element %in% ref_names)
    return(n_set)
})
