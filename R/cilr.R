
#' Main cILR function to perform enrichment analysis using cILR method
#' @param otu_table Named \code{n} by \code{p} matrix. This is the OTU/ASV/Strain table where taxa are columns.
#' @param tax_sets List of length \code{m}. This is a list of set membership by column names.
#' @param output String. The form of the output of the model.
#' @param distr String. The choice of distribution for the null.
#' @param adj Logical. Whether correlation adjustment procedure is utilized.
#' @param nperm Integer. The number of permutations/bootstrap re-sampling. Recommended not to exceed 5
#'     due to memory issues.
#' @param thresh Numeric. Threshold for significant p-values if \code{sig} is the output.
#' @param init Named List. Initialization parameters for estimating the null distribution. Default is NULL.
#' @param raw Logical. Whether scores are returned as raw (no parameter estimation step). Default is FALSE.
#' @param min_size Integer. The lower bound of set size. Sets of size lower than min_size will be filtered out.
#'     Default is 2 (no singleton sets).
#' @param ... Named List. Additional arguments to be passed to \code{fitdistr} and \code{normmixEM}
#'
#' @return \code{R}    An \code{n} by \code{m} matrix of enrichment scores at the sample level
#'
#' @import rlang GSEABase
#' @export
cilr <- function(otu_table, tax_sets,
                 output = c("cdf", "zscore", "pval", "sig"),
                 distr = c("mnorm", "norm"),
                 adj = TRUE,
                 nperm = 5,
                 thresh = 0.05,
                 init = NULL,
                 min_size = 2,
                 raw = FALSE, ...){
    # check arguments
    output <- match.arg(output)
    distr <- match.arg(distr)
    # check input
    if (!is.matrix(otu_table)){
        rlang::warn("Coercing OTU/ASV table into matrix format")
        X <- as.matrix(tab)
    }

    if (class(tax_sets) != "GeneSetCollection"){
        rlang::abort("A has to be of the GeneSet type")
    }

    p <- ncol(X) # number of features
    n <- nrow(X) # number of samples

    return(0)
}






