# Class definition and documentation for basic functionality of the phyloseqSet class
# PhyloseqSet is an S4 object extending the existing phyloseq class, adding one more slot
# for support of taxa sets
# Part of the code adapted from the phyloseq package itself (https://github.com/joey711/phyloseq)
# Quang Nguyen

#' @keywords internal
setClassUnion("BiocSetOrNULL", c("BiocSet", "NULL"))

#' @title S4 class extending phyloseq object to also store taxon sets
#' @aliases phyloseqSet-class
#' @description This is a simple S4 object that extends \code{\link{phyloseq}} by one extra slot
#'     with a \code{taxon_set} object which is of class \code{\link{BiocSet}} class.
#' @slot otu_table A \code{\link{otu_table}} type from \code{\link{phyloseq}}
#' @slot sam_data A \code{\link{sam_data}} type from \code{\link{phyloseq}}
#' @slot tax_table A \code{\link{tax_table}} type from \code{\link{phyloseq}}
#' @slot phy_tree A \code{\link{phy_tree}} type from \code{\link{phyloseq}}
#' @slot refseq A \code{\link{refseq}} type from \code{\link{phyloseq}}
#' @slot taxon_set A \code{\link{BiocSet}} object
#' @importClassesFrom BiocSet BiocSet
#' @importClassesFrom phyloseq phyloseq
#' @name phyloseqSet-class
#' @rdname phyloseqSet-class
#' @exportClass phyloseqSet
setClass("phyloseqSet",
    contains = "phyloseq",
    slots = c(
        taxon_set = "BiocSetOrNULL"
    ),
    prototype = list(
        taxon_set = NULL
    )
)

#' @keywords internal
check_valid <- function(object){
    # Object has to be a phyloseq object
    if (is(object, "phyloseq") == FALSE){
        return("Object have to be a child of phyloseq type")
    }
    if (!is.null(object@taxon_set)){
        # Only one taxon_set object available. BiocSets supports multiple set definitions that allow for
        # overlapping or multi-sets
        if (length(object@taxon_set) > 1){
            return("There should only be one taxon set object. BiocSets supports multiple set definitions")
        }
        # Check type of taxon_sets
        if (class(object@taxon_set) != "BiocSet"){
            return("taxon_sets has to be a BiocSet")
        }
        # Check the names of taxon_sets
        names <- dplyr::pull(BiocSet::es_element(object@taxon_set))
        matches <- base::match(names, phyloseq::taxa_names(object))
        if (length(matches) != length(names)){
            return("List of element taxa does not completely match taxa_names of phyloseq object")
        }
    }
    return(TRUE)
}


setValidity("phyloseqSet", check_valid)



#' @title Construct a phyloseqSet object
#' @aliases phyloseqSet
#' @description \code{phyloseqSet} is a child class of the \code{phyloseq} object from
#'     the package of the same name. This function takes basic arguments simlar to \code{phyloseq}
#'     with a taxon_set of class \code{BiocSet} and construct an equivament \code{phyloseqSet} object.
#' @details This function can take in either the usual constructors for \code{phyloseq} object
#'     (e.g. \code{\link{otu_table}}, \code{\link{sample_data}}, etc.) or a \code{phyloseq} object itself.
#' @usage phyloseqSet(..., taxon_set)
#' @param ... This can be either the components of constructing a \code{phyloseq} class (refer to
#'     \code{\link{phyloseq}} for more information), or a \code{phyloseq} object itself. If only one element
#'     is supplied, must either be an \code{otu_table} object or a whole \code{phyloseq} object.
#' @param taxon_set \code{BiocSet} object of taxa sets
#' @rdname phyloseqSet
#' @importFrom phyloseq phyloseq
#' @importFrom methods new
#' @seealso \code{\link{phyloseq}}
#' @export
phyloseqSet <- function(..., taxon_set=NULL){
    args <- list(...)
    if (length(args) > 1){
        physeq <- do.call(phyloseq::phyloseq, args)
    } else {
        if (class(args[[1]]) != "phyloseq"){
            rlang::abort("If one argument is supplied, must be phyloseq object or an otu_table object")
        } else {
            physeq <- args[[1]]
        }
    }
    new("phyloseqSet", physeq, taxon_set=taxon_set)
}

#' @title Access taxon set object from \code{phyloseqSet} object
#' @param x \link{phyloseqSet}
setGeneric("taxon_set", function(x) standardGeneric("taxon_set"), signature = "x")
setMethod("taxon_set", "phyloseqSet", function(x) x@taxon_set)

#' @title Assign a new taxon set object to \code{x}
#' @param x \code{\link{phyloseqSet-class}}
#' @param value \code{\link{BiocSet}} class of taxa sets
setGeneric("taxon_set<-", function(x, value) standardGeneric("taxon_set<-"), signature = "x")
setMethod("taxon_set<-", "phyloseqSet", function(x, value){
    x@taxon_set <- value
    validObject(x)
    x
})

# TODO: Please remove placeholder objects here please Quang you're my only hope
setMethod("show", "phyloseqSet", function(object){
    # callNextMethod calls the same inherited method
    callNextMethod()
    if (!is.null(taxon_set(object))){
        cat(paste("taxon_set()   Taxon Set Tables:  [ ",
                  nrow(BiocSet::es_element(taxon_set(object))),
                  " taxa by ",
                  nrow(BiocSet::es_set(taxon_set(object))),
                  " taxa sets ]", sep = ""), fill = TRUE)
    }
})


