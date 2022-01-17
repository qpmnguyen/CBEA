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
#' mod <- cbea(obj = seq, set = set, output = "zscore", distr = "norm", parametric = TRUE,
#'     adj = TRUE, thresh = 0.05)
NULL


#' @rdname cbea
#' @export
#' @import methods
setGeneric("cbea", function(obj, set,
                            output,
                            distr,
                            adj = FALSE,
                            n_perm = 1,
                            parametric = TRUE,
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
                                       distr = NULL,
                                       adj = NULL,
                                       n_perm = 1,
                                       parametric = TRUE,
                                       thresh = 0.05,
                                       init = NULL,
                                       control = NULL, ...) {
    # Validate inputs ####
    output <- match.arg(output, choices = c("cdf", "zscore", "pval", "sig", "raw"))
    distr <- match.arg(distr, choices = c("mnorm", "norm", NULL))
    
    args <- match.call()
    print(args)

    # wrangle data into the correct format
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
        distr = distr, adj = adj, n_perm = n_perm, parametric = parametric,
        thresh = thresh, init = init,
        control = control, ...
    )
    return(model)
})

#' @rdname cbea
#' @importClassesFrom TreeSummarizedExperiment TreeSummarizedExperiment
#' @importFrom SummarizedExperiment assays
#' @import methods
#' @export
setMethod("cbea", "TreeSummarizedExperiment", function(obj, set,
                                                       output,
                                                       distr,
                                                       adj = FALSE,
                                                       n_perm = 1,
                                                       parametric = TRUE,
                                                       thresh = 0.05,
                                                       init = NULL,
                                                       control = NULL, ...) {
    # Validate inputs ####
    output <- match.arg(output, choices = c("cdf", "zscore", "pval", "sig", "raw"))
    distr <- match.arg(distr, choices = c("mnorm", "norm", NULL))
    # handling instances if distr was missing
    if (missing(distr) | is.null(distr)){
        if (parametric == TRUE){
            message("Distribution was not specified, returning raw scores")
            output <- "raw"
        } else {
            distr <- NULL
        }
    }
    # handling instances if adj was missing
    if (missing(adj) | is.null(adj)){
        if (parametric == TRUE){
            message("Correlation adjustment was not specified, defaulting to FALSE")
            adj <- FALSE
        } else {
            adj <- NULL
        }
    }

    # if parametric is false cannot get either cdf values or z-scores
    if (parametric == FALSE){
        if (output %in% c("zscore", "cdf")){
            stop("Output cannot be either z-scores or CDF values if no parametric fit was performed")
        }
    }

    # generate table
    tab <- assays(obj)[[1]]
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
        distr = distr, adj = adj, n_perm = n_perm, parametric = parametric,
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
                                         n_perm = 1,
                                         parametric = TRUE,
                                         thresh = 0.05,
                                         init = NULL,
                                         control = NULL,
                                         taxa_are_rows = FALSE, ...) {
    # Validate inputs ####
    output <- match.arg(output, choices = c("cdf", "zscore", "pval", "sig", "raw"))
    distr <- match.arg(distr, choices = c("mnorm", "norm", NULL))
    # handling instances if distr was missing
    if (missing(distr) | is.null(distr)){
        if (parametric == TRUE){
            message("Distribution was not specified, returning raw scores")
            output <- "raw"
        } else {
            distr <- NULL
        }
    }
    # handling instances if adj was missing
    if (missing(adj) | is.null(adj)){
        if (parametric == TRUE){
            message("Correlation adjustment was not specified, defaulting to FALSE")
            adj <- FALSE
        } else {
            adj <- NULL
        }
    }

    # if parametric is false cannot get either cdf values or z-scores
    if (parametric == FALSE){
        if (output %in% c("zscore", "cdf")){
            stop("Output cannot be either z-scores or CDF values if no parametric fit was performed")
        }
    }
    if (taxa_are_rows == TRUE) {
        check_numeric <- vapply(tab, FUN = function(x) is(x, "numeric"), FUN.VALUE = TRUE)
        idx <- which(check_numeric == FALSE)
        if (length(numeric) >= 2) {
            stop("There are more than two non-numeric columns,
                please remove columns un-associated with taxa identifiers")
        } else if (length(numeric) == 1) {
            warning("Removing the one non-numeric column, assuming that it's the sample/taxa identifiers")
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
        output = output, distr = distr, adj = adj, n_perm = n_perm,
        parametric = parametric,
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
                                     n_perm = 1,
                                     parametric = TRUE,
                                     thresh = 0.05,
                                     init = NULL,
                                     control = NULL,
                                     taxa_are_rows = FALSE, ...) {
    # Validate inputs ####
    output <- match.arg(output, choices = c("cdf", "zscore", "pval", "sig", "raw"))
    distr <- match.arg(distr, choices = c("mnorm", "norm", NULL))
    # handling instances if distr was missing
    if (missing(distr) | is.null(distr)){
        if (parametric == TRUE){
            message("Distribution was not specified, returning raw scores")
            output <- "raw"
        } else {
            distr <- NULL
        }
    }
    # handling instances if adj was missing
    if (missing(adj) | is.null(adj)){
        if (parametric == TRUE){
            message("Correlation adjustment was not specified, defaulting to FALSE")
            adj <- FALSE
        } else {
            adj <- NULL
        }
    }

    # if parametric is false cannot get either cdf values or z-scores
    if (parametric == FALSE){
        if (output %in% c("zscore", "cdf")){
            stop("Output cannot be either z-scores or CDF values if no parametric fit was performed")
        }
    }
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