#' @title Enrichment analysis using competitive compositional balances (CBEA)
#' @description \code{cbea} is used compute enrichment scores
#'     per sample for pre-defined sets using the CBEA
#'     (Competitive Balances for Enrichment Analysis).
#' @details This function support different formats of the OTU table, however
#'     for best results please use \code{\linkS4class{TreeSummarizedExperiment}}.
#'     \code{phyloseq} is supported, however \code{CBEA} will not explicitly import
#'     \code{phyloseq} package and will require users to install them separately.
#'     If use \code{data.frame} or \code{matrix}, users should specify whether taxa are rows using the
#'     \code{taxa_are_rows} option. Additionally, for \code{data.frame}, users can specify
#'     metadata columns to be kept via the \code{id_col} argument. \cr
#'     The \code{output} argument specifies what type of values will be returned in the final matrix.
#'     The options \code{pval} or \code{sig} returns either unadjusted p-values or dummy variables
#'     indicating whether a set is significantly enriched in that sample (based on unadjusted
#'     p-values thresholded at \code{thresh}).
#'     The option \code{raw} returns raw scores computed for each set without any distribution fitting
#'     or inference procedure. Users can use this option to examine the distribution of CBEA scores
#'     under the null.
#' @param obj The element of class \code{TreeSummarizedExperiment},
#'     \code{data.frame}, or \code{matrix}. \code{phyloseq} is not supported
#'     due to conflicting dependencies and \code{TreeSummarizedExperiment} is much
#'     more compact.
#' @param set \code{BiocSet}. Sets to be tested for
#'     enrichment in the \code{BiocSet}
#'     format. Taxa names must be in the same format
#'     as elements in the set.
#' @param output (String). The form of the output of the model.
#'     Has to be either \code{zscore},
#'     \code{cdf}, \code{raw}, \code{pval}, or \code{sig}
#' @param distr (String). The choice of distribution for the null. Can be either \code{mnorm}
#'     (2 component mixture normal), \code{norm} (Normal distribution), or NULL if \code{parametric}
#'     is \code{TRUE}.
#' @param adj (Logical). Whether correlation adjustment procedure is utilized. Defaults to FALSE.
#' @param n_perm (Numeric). Add bootstrap resamples to both the permuted and unpermuted data set.
#'     This might help with stabilizing the distribution fitting procedure, especially if the sample
#'     size is low. Defaults to 1.
#' @param parametric (Logical). Indicate whether a parametric distribution will be fitted to estimate
#'     z-scores, CDF values, and p-values. Defaults to \code{TRUE}
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
#' # n_perm = 10 to reduce runtime
#' mod <- cbea(obj = seq, set = set, output = "zscore",
#'     abund_values = "16SrRNA",
#'     distr = "norm", parametric = TRUE,
#'     adj = TRUE, thresh = 0.05, n_perm = 10)
NULL


#' @rdname cbea
#' @importClassesFrom BiocSet BiocSet
#' @export
setGeneric("cbea", function(obj, set,
                            output,
                            distr = NULL,
                            adj = FALSE,
                            n_perm = 100,
                            parametric = TRUE,
                            thresh = 0.05,
                            init = NULL,
                            control = NULL, ...) {
  standardGeneric("cbea")
}, signature = "obj")



#' @rdname cbea
#' @param abund_values (Character). Character value for selecting the \code{assay} to be
#'     the input to \code{cbea}
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom SummarizedExperiment assay assays
#' @importFrom methods as
#' @export
setMethod("cbea", "TreeSummarizedExperiment", function(obj, set,
                                                       output,
                                                       distr = NULL,
                                                       abund_values,
                                                       adj = FALSE,
                                                       n_perm = 100,
                                                       parametric = TRUE,
                                                       thresh = 0.05,
                                                       init = NULL,
                                                       control = NULL, ...) {
    # Validate inputs ####
    check_args()
    call <- match.call()
    # generate table
    if (!abund_values %in% names(assays(obj))){
        stop("abund_values must be part of the assays of the object.
             Check names(assays(obj)) for details")
    }
    tab <- assay(obj, abund_values)
    # TreeSummarizedExperiment data sets are always transposed
    tab <- as(tab, "matrix")
    tab <- t(tab)
    if (length(which(tab == 0)) > 0) {
        warning("Taxonomic count table contains zeros,
            which would invalidate the log-ratio transform.
            Adding a pseudocount of 1e-5...")
        tab <- tab + 1e-5
    }
    set_list <- as(set, "list")
    model <- .cbea(
        ab_tab = tab, set_list = set_list, output = output,
        distr = distr, adj = adj, n_perm = n_perm, parametric = parametric,
        thresh = thresh, init = init,
        raw = raw, control = control, ...
    )
    
    out <- new_CBEAout(model, sample_ids = rownames(tab), 
                       output = output, 
                       parametric, 
                       distr) 
    return(out)
})

#' @rdname cbea
#' @param taxa_are_rows (Logical). Indicate whether the data frame
#'     or matrix has taxa as rows
#' @param id_col (Character Vector). Vector of character to indicate metadata
#'     columns to keep (for example, \code{sample_id})
#' @importFrom dplyr select
#' @importFrom methods as is
#' @export
setMethod("cbea", "data.frame", function(obj, set,
                                         taxa_are_rows = FALSE,
                                         id_col = NULL,
                                         output,
                                         distr = NULL,
                                         adj = FALSE,
                                         n_perm = 100,
                                         parametric = TRUE,
                                         thresh = 0.05,
                                         init = NULL,
                                         control = NULL,
                                         ...) {
    # Validate inputs ####
    check_args()

    if (!is.null(id_col)){
        tab <- dplyr::select(obj, -c({{id_col}}))
    } else {
        tab <- obj
    }

    check_numeric <- vapply(tab, FUN = function(x) is(x, "numeric"), FUN.VALUE = TRUE)
    idx <- which(check_numeric == FALSE)
    if (length(idx) >= 1){
        stop("There are more than two non-numeric columns after removing
             row identifiers. Please explicitly remove these additional non-numeric columns")
    }
    tab <- as(tab, "matrix")
    # transformations and issues with data.frames
    if (taxa_are_rows == TRUE){
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
        output = output, distr = distr, adj = adj, n_perm = n_perm,
        parametric = parametric,
        thresh = thresh, init = init, control = control, ...
    )
    return(model)
})

#' @rdname cbea
#' @importFrom methods as
#' @export
setMethod("cbea", "matrix", function(obj, set,
                                     taxa_are_rows = FALSE,
                                     output,
                                     distr = NULL,
                                     adj = FALSE,
                                     n_perm = 100,
                                     parametric = TRUE,
                                     thresh = 0.05,
                                     init = NULL,
                                     control = NULL,
                                     ...) {
    # Validate inputs ####
    check_args()

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
        n_perm = n_perm, parametric = parametric,
        thresh = thresh, init = init, control = control, ...
    )
    return(model)
})
