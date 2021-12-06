#' @title Enrichment analysis using competitive compositional balances (CBEA)
#' @param obj The element of class \code{TreeSummarizedExperiment}, \code{phyloseq}, \code{data.frame},
#'     or \code{matrix}
#' @param set \code{BiocSet}. Sets to be tested for enrichment in the \code{BiocSet} format.
#'     Taxa names must be in the same format as elements in the set.
#' @param output (String). The form of the output of the model. Has to be either \code{zscore},
#'     \code{cdf}, \code{zscore}, \code{pval}, \code{sig}
#' @param distr (String). The choice of distribution for the null.
#' @param adj (Logical). Whether correlation adjustment procedure is utilized.
#' @param thresh (Numeric). Threshold for significant p-values if \code{sig} is the output.
#' @param init (Named List). Initialization parameters for estimating the null distribution. Default is NULL.
#' @param raw (Logical). Whether scores are returned as raw (no parameter estimation step). Default is FALSE.
#' @param control (Named List). Additional arguments to be passed to \code{fitdistr} and \code{normmixEM}
#' @param ... Additional arguments not used at the moment.
#' @return \code{R}    An \code{n} by \code{m} matrix of enrichment scores at the sample level
#' @name cbea
NULL


#' @rdname cbea
#' @export
#' @import methods
setGeneric("cbea", function(obj, set,
                            output,
                            distr,
                            adj = FALSE, thresh = 0.05,
                            init = NULL, raw = FALSE, control = NULL, ...) {
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
                                       adj = FALSE, thresh = 0.05,
                                       init = NULL, raw = FALSE, control = NULL, ...) {
  output <- match.arg(output, choices = c("cdf", "zscore", "pval", "sig"))
  distr <- match.arg(distr, choices = c("mnorm", "norm"))
  tab <- phyloseq::otu_table(obj)
  tab <- as(tab, "matrix")

  if (phyloseq::taxa_are_rows(obj) == TRUE) {
    tab <- t(tab)
  }

  if (length(which(tab == 0)) > 0) {
    warning("Taxonomic count table contains zeros, which would invalidate the log-ratio transform. Adding a pseudocount of 1...")
    tab <- tab + 1
  }
  set_list <- as(set, "list")
  model <- .cbea(
    ab_tab = tab, set_list = set_list, output = output,
    distr = distr, adj = adj, thresh = thresh, init = init,
    raw = raw, control = control, ...
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
                                                       adj = FALSE, thresh = 0.05,
                                                       init = NULL, raw = FALSE, control = NULL, ...) {
  output <- match.arg(output, choices = c("cdf", "zscore", "pval", "sig"))
  distr <- match.arg(distr, choices = c("mnorm", "norm"))
  tab <- SummarizedExperiment::assays(obj)[[1]]
  # TreeSummarizedExperiment data sets are always transposed
  tab <- as(tab, "matrix")
  tab <- t(tab)
  if (length(which(tab == 0)) > 0) {
    warning("Taxonomic count table contains zeros,
                which would invalidate the log-ratio transform. Adding a pseudocount of 1...")
    tab <- tab + 1
  }
  set_list <- as(set, "list")
  model <- .cbea(
    ab_tab = tab, set_list = set_list, output = output,
    distr = distr, adj = adj, thresh = thresh, init = init,
    raw = raw, control = control, ...
  )
  return(model)
})

#' @rdname cbea
#' @import methods
#' @export
#' @param taxa_are_rows (Logical). Indicate whether the data frame or matrix has taxa as rows
setMethod("cbea", "data.frame", function(obj, set,
                                         output,
                                         distr,
                                         adj = FALSE, thresh = 0.05,
                                         init = NULL, raw = FALSE, control = NULL, taxa_are_rows = FALSE, ...) {
  output <- match.arg(output, choices = c("cdf", "zscore", "pval", "sig"))
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
    ab_tab = tab, set_list = set_list, output = output, distr = distr, adj = adj,
    thresh = thresh, init = init, raw = raw, control = control, ...
  )
  return(model)
})

#' @rdname cbea
#' @import methods
#' @export
setMethod("cbea", "matrix", function(obj, set,
                                     output,
                                     distr,
                                     adj = FALSE, thresh = 0.05,
                                     init = NULL, raw = FALSE, control = NULL, taxa_are_rows = FALSE, ...) {
  output <- match.arg(output, choices = c("cdf", "zscore", "pval", "sig"))
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
    ab_tab = tab, set_list = set_list, output = output, distr = distr, adj = adj,
    thresh = thresh, init = init, raw = raw, control = control, ...
  )
  return(model)
})
