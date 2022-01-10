#' @title Enrichment analysis using competitive compositional balances (CBEA)
#' @description \code{cbea} is used compute enrichment scores
#'     per sample for pre-defined sets using the CBEA
#'     (Competitive Balances for Enrichment Analysis).
#' @details This function support different formats of the OTU table, however
#'     for best results please use either \code{\linkS4class{TreeSummarizedExperiment}} or
#'     \code{\linkS4class{phyloseq}}. If use \code{data.frame} or \code{matrix}, additional
#'     arguments specifying if \code{taxa_are_rows} should be true or not. \cr
#'     The \code{output} specifies what type of values will be returned in the final matrix.
#'     The options \code{pval} or \code{sig} returns either unadjusted p-values or dummy variables
#'     indicating whether a set is significantly enriched in that sample (based on unadjusted
#'     p-values thresholded at \code{thresh}).
#'     The option \code{raw} returns raw scores computed for each set without any distribution fitting
#'     or inference procedure. Users can use this option to examine the distribution of CBEA scores
#'     under the null.
#' @param obj The element of class \code{TreeSummarizedExperiment},
#'     \code{phyloseq}, \code{data.frame}, or \code{matrix}
#' @param set \code{BiocSet}. Sets to be tested for
#'     enrichment in the \code{BiocSet}
#'     format. Taxa names must be in the same format
#'     as elements in the set.
#' @param output (String). The form of the output of the model.
#'     Has to be either \code{zscore},
#'     \code{cdf}, \code{raw}, \code{pval}, or \code{sig}
#' @param distr (String). The choice of distribution for the null. Can be either \code{mnorm}
#'     (2 component mixture normal) or \code{norm} (Normal distribution).
#' @param adj (Logical). Whether correlation adjustment procedure is utilized. Defaults to FALSE.
#' @param n_boot (Numeric). Add bootstrap resamples to both the permuted and unpermuted data set.
#'     This might help with stabilizing the distribution fitting procedure, especially if the sample
#'     size is low. CBEA will throw a warning if the object size in memory of the resampled data set
#'     is larger than 250 MiB. Defaults to 1.
#' @param thresh (Numeric). Threshold for significant p-values if \code{sig}
#'     is the output. Defaults to 0.05
#' @param init (Named List). Initialization parameters
#'     for estimating the null distribution.
#'     Default is NULL.
#' @param control (Named List). Additional arguments to be passed to
#'     \code{fitdistr} and \code{normmixEM}. Defaults to NULL.
#' @param ... Additional arguments not used at the moment.
#' @return \code{R}    An \code{n} by \code{m} matrix of enrichment scores at
#'     the sample level
#' @name cbea
#' @examples
#' data(hmp_gingival)
#' seq <- hmp_gingival$data
#' set <- hmp_gingival$set
#' mod <- cbea(obj = seq, set = set, output = "zscore", distr = "norm",
#'     adj = TRUE, thresh = 0.05)
NULL


#' @rdname cbea
#' @export
#' @import methods
setGeneric("cbea", function(obj, set,
                            output,
                            distr,
                            adj = FALSE,
                            n_boot = 1,
                            thresh = 0.05,
                            init = NULL,
                            control = NULL, ...) {
  standardGeneric("cbea")
}, signature = "obj")

#' @rdname cbea
#' @importClassesFrom phyloseq phyloseq
#' @importFrom phyloseq otu_table taxa_are_rows
#' @import methods
#' @export
setMethod("cbea", "phyloseq", function(obj, set,
                                       output,
                                       distr,
                                       adj = FALSE,
                                       n_boot = 1,
                                       thresh = 0.05,
                                       init = NULL,
                                       control = NULL, ...) {
    output <- match.arg(output, choices = c("cdf", "zscore", "pval", "sig", "raw"))
    distr <- match.arg(distr, choices = c("mnorm", "norm"))
    tab <- phyloseq::otu_table(obj)
    tab <- as(tab, "matrix")

    if (phyloseq::taxa_are_rows(obj) == TRUE) {
        tab <- t(tab)
    }

    if (length(which(tab == 0)) > 0) {
        warning("Taxonomic count table contains zeros,
            which would invalidate the log-ratio transform.
            Adding a pseudocount of 1...")
        tab <- tab + 1
    }
    set_list <- as(set, "list")
    model <- .cbea(
        ab_tab = tab, set_list = set_list, output = output,
        distr = distr, adj = adj, n_boot = n_boot,
        thresh = thresh, init = init,
        control = control, ...
    )
    return(model)
})

#' @rdname cbea
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @import TreeSummarizedExperiment
#' @import methods
#' @export
setMethod("cbea", "TreeSummarizedExperiment", function(obj, set,
                                                       output,
                                                       distr,
                                                       adj = FALSE,
                                                       n_boot = 1,
                                                       thresh = 0.05,
                                                       init = NULL,
                                                       control = NULL, ...) {
    output <- match.arg(output, choices = c("cdf", "zscore", "pval", "sig", "raw"))
    distr <- match.arg(distr, choices = c("mnorm", "norm"))
    tab <- SummarizedExperiment::assays(obj)[[1]]
    # TreeSummarizedExperiment data sets are always transposed
    tab <- as(tab, "matrix")
    tab <- t(tab)
    if (length(which(tab == 0)) > 0) {
        warning("Taxonomic count table contains zeros,
            which would invalidate the log-ratio transform.
            Adding a pseudocount of 1...")
        tab <- tab + 1
    }
    set_list <- as(set, "list")
    model <- .cbea(
        ab_tab = tab, set_list = set_list, output = output,
        distr = distr, adj = adj, n_boot = n_boot,
        thresh = thresh, init = init,
        raw = raw, control = control, ...
    )
    return(model)
})

#' @rdname cbea
#' @import methods
#' @export
#' @param taxa_are_rows (Logical). Indicate whether the data frame
#'     or matrix has taxa as rows
setMethod("cbea", "data.frame", function(obj, set,
                                         output,
                                         distr,
                                         adj = FALSE,
                                         thresh = 0.05,
                                         init = NULL,
                                         control = NULL,
                                         taxa_are_rows = FALSE, ...) {
    output <- match.arg(output, choices = c("cdf", "zscore", "pval", "sig", "raw"))
    distr <- match.arg(distr, choices = c("mnorm", "norm"))
    if (taxa_are_rows == TRUE) {
        check_numeric <- vapply(tab, class)
        idx <- which(check_numeric != "numeric")
        if (length(idx) >= 2) {
            stop("There are more than two non-numeric columns,
                please remove columns un-associated with taxa identifiers")
        } else if (length(idx) == 1) {
            rownames(tab) <- tab[, idx]
            tab <- tab[, -idx]
        }
        tab <- as(obj, "matrix")
        tab <- t(tab)
    }

    if (length(which(tab == 0)) > 0) {
        warning("Taxonomic count table contains zeros,
            which would invalidate the log-ratio transform.
            Adding a pseudocount of 1...")
        tab <- tab + 1
    }

    set_list <- as(set, "list")
    model <- .cbea(
        ab_tab = tab, set_list = set_list,
        output = output, distr = distr, adj = adj,
        thresh = thresh, init = init, control = control, ...
    )
    return(model)
})

#' @rdname cbea
#' @import methods
#' @export
setMethod("cbea", "matrix", function(obj, set,
                                     output,
                                     distr,
                                     adj = FALSE,
                                     thresh = 0.05,
                                     init = NULL,
                                     control = NULL,
                                     taxa_are_rows = FALSE, ...) {
    output <- match.arg(output, choices = c("cdf", "zscore", "pval", "sig", "raw"))
    distr <- match.arg(distr, choices = c("mnorm", "norm"))
    if (taxa_are_rows == TRUE) {
        tab <- t(tab)
    }
    if (length(which(tab == 0)) > 0) {
        warning("Taxonomic count table contains zeros,
            which would invalidate the log-ratio transform.
            Adding a pseudocount of 1...")
        tab <- tab + 1
    }

    set_list <- as(set, "list")
    model <- .cbea(
        ab_tab = tab, set_list = set_list,
        output = output, distr = distr, adj = adj,
        thresh = thresh, init = init, control = control, ...
    )
    return(model)
})
