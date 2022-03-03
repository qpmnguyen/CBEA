# CREATING CBEAout ####

#' @title Creating an output object of type CBEAout
#' @description This function takes a list of lists from each object and
#'  turns it into a CBEAout type object
#' @param out A list containing scores for each set
#' @param call A list containing all important arguments for printing
#' @return A new CBEAout object (which is a cleaner list of lists)
new_CBEAout <- function(out, call){
    req_values <- c("parametric", "distr", "output", "n_perm", "adj", "sample_ids")
    if (length(intersect(names(call), req_values)) != length(req_values)){
        stop("Require the 'call' argument to have all important elements of a function call")
    }

    # first process the score_matrix
    score_matrix <- lapply(out, function(x) x$score)
    if (call$output == "raw" | call$parametric == FALSE){
        parameters <- NULL
        diagnostic <- NULL
        fit_comparison <- NULL
    }  else {
        # second, get parameters
        parameters <- lapply(out, function(x) unlist(x$diagnostic$parameters, use.names = TRUE))

        # third, get diagnostic values
        diagnostic <- lapply(seq_along(out), function(x) {
            out[[x]]$diagnostic[-c(which(names(out[[x]]$diagnostic) %in% c("distr", "parameters")))]
        })
        names(diagnostic) <- names(out)
        # fourth, get fit comparisons
        fit_comparison <- lapply(out, function(x){
            x$lmoment
        })
        names(fit_comparison)
    }
    object <- list(
        R = score_matrix,
        parameters = parameters,
        diagnostic = diagnostic,
        fit_comparison = fit_comparison,
        call_param = call
    )
    class(object) <- "CBEAout"
    return(object)
}


# METHODS FOR CBEAout ####

#' @title Print dispatch for CBEAout objects
#' @param x The \code{CBEAout} object
#' @param ... Undefined arguments, keeping consistency for generics
#' @importFrom glue glue
#' @return Text for printing
#' @export
print.CBEAout <- function(x, ...){
    if (!is.null(x$call_param$distr)){
        distr <- switch(x$call_param$distr,
                        norm = "Gaussian",
                        mnorm = "2-component Gaussian Mixture",
                        lst = "Location-scale Student's t")
    } else {
        distr <- "No"
    }
    output <- switch(x$call_param$output,
                     raw = "Raw scores (raw)",
                     cdf = "CDF values (cdf)",
                     zscore = "Z-scores (zscore)",
                     sig = "Dummy variable for significance at threshold (sig)",
                     pval = "Unadjusted p-values (pval)")
    par <- ifelse(x$call_param$parametric, "Parametric", "Non-parametric")
    cat(
        glue("CBEA output of class 'CBEAout' with {nsamp} samples and {nset} sets",
             nsamp = length(x$R[[1]]), nset = length(x$R)), "\n",
        glue("Fit type: {par} with {distr} Distribution"), "\n",
        glue("Number of permutations: {n_perm}", n_perm = x$call_param$n_perm), "\n",
        glue("Output type: {output}"), "\n"
    )
}
