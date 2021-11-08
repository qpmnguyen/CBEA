#' @title Enrichment analysis using competitive compositional balances (CBEA)
#' @param obj The element of class \code{TreeSummarizedExperiment}, \code{phyloseq}, \code{data.frame},
#'     or \code{matrix}
#' @param set \code{BiocSet}. Sets to be tested for enrichment in the \code{BiocSet} format.
#'     Taxa names must be in the same format as elements in the set.
#' @param output String. The form of the output of the model. Has to be either \code{zscore},
#'     \code{cdf}, \code{zscore}, \code{pval}, \code{sig}
#' @param distr String. The choice of distribution for the null.
#' @param adj Logical. Whether correlation adjustment procedure is utilized.
#' @param thresh Numeric. Threshold for significant p-values if \code{sig} is the output.
#' @param init Named List. Initialization parameters for estimating the null distribution. Default is NULL.
#' @param raw Logical. Whether scores are returned as raw (no parameter estimation step). Default is FALSE.
#' @param ... Named List. Additional arguments to be passed to \code{fitdistr} and \code{normmixEM}
#'
#' @return \code{R}    An \code{n} by \code{m} matrix of enrichment scores at the sample level
#' @name cbea
NULL


#' @rdname cbea
setGeneric("cbea", function(obj, set,
                            output = c("cdf", "zscore", "pval", "sig"),
                            distr = c("mnorm", "norm"),
                            adj = FALSE, thresh = 0.05,
                            init = NULL, raw = FALSE, ...) standardGeneric("cbea"), signature = "obj")

#' @rdname cbea
#' @importClassesFrom phyloseq phyloseq
#' @importFrom phyloseq otu_table taxa_are_rows
setMethod("cbea", "phyloseq", function(obj, set,
                                       output = c("cdf", "zscore", "pval", "sig"),
                                       distr = c("mnorm", "norm"),
                                       adj = FALSE, thresh = 0.05,
                                       init = NULL, raw = FALSE, ...){
    tab <- phyloseq::otu_table(obj)
    tab <- as(tab, "matrix")

    if (phyloseq::taxa_are_rows(obj) == TRUE){
        tab <- t(tab)
    }

    if (length(which(tab == 0)) > 0){
        warning("Taxonomic count table contains zeros, which would invalidate the log-ratio transform. Adding a pseudocount of 1...")
        tab <- tab + 1
    }
    set_list <- as(set, "list")
    model <- .cbea(ab_tab = tab, set_list = set_list, output = output,
                   distr = distr, adj = adj, thresh = thresh, init = init,
                   raw = raw, ...)
    return(model)
})

#' @rdname cbea
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @import TreeSummarizedExperiment
setMethod("cbea", "TreeSummarizedExperiment", function(obj, set,
                                                       output = c("cdf", "zscore", "pval", "sig"),
                                                       distr = c("mnorm", "norm"),
                                                       adj = FALSE, thresh = 0.05,
                                                       init = NULL, raw = FALSE, ...){
    tab <- SummarizedExperiment::assays(obj)[[1]]
    # TreeSummarizedExperiment data sets are always transposed
    tab <- as(tab, "matrix")
    tab <- t(tab)
    if (length(which(tab == 0)) > 0){
        warning("Taxonomic count table contains zeros, which would invalidate the log-ratio transform. Adding a pseudocount of 1...")
        tab <- tab + 1
    }
    set_list <- as(set, "list")
    model <- .cbea(ab_tab = tab, set_list = set_list, output = output,
                   distr = distr, adj = adj, thresh = thresh, init = init,
                   raw = raw, ...)
    return(model)
})








