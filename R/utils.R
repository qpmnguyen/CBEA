#' @title This function handles the ability to merge supplied and defaults
#' @param defaults (List). Default options
#' @param supplied (List). Supplied options
#' @return A merged list
#' @keywords internal
merge_lists <- function(defaults, supplied) {
    similar_idx <- which(names(defaults) %in% names(supplied))
    if (length(similar_idx) == 0) {
        merged <- c(defaults, supplied)
    } else {
        merged <- c(defaults[-similar_idx], supplied)
    }
    return(merged)
}

#' @title The Two Component Mixture Normal Distribution
#' @param q,x (Vector). Values to calculate distributional values of.
#' @param mu (Vector). A two value vector of mean values.
#' @param sigma (Vector). A two value vector of component-wise variances
#' @param lambda (Vector). A two value vector of component mixing coefficients
#' @param log (Boolean). Whether returning probabilities are in log format
#' @param verbose (Boolean). Whether to return component values.
#' @return A numeric value representing the probability
#'     density value of a two-component mixture distribution
#' @describeIn pmnorm Cumulative Distribution Function
#' @export
#' @importFrom stats pnorm
pmnorm <- function(q, mu, sigma, lambda, log = FALSE, verbose = FALSE) {
    q <- as.vector(q)
    n_components <- length(sigma)
    if (verbose == TRUE) {
        message(paste(n_components, "components!", "\n"))
    }
    comp <- vector(mode = "list", length = n_components)
    for (i in seq_len(n_components)) {
        comp[[i]] <- lambda[i] * pnorm(q, mu[i], sigma[i], log.p = log)
    }
    return(Reduce("+", comp))
}

#' @describeIn pmnorm Probability Density Function
#' @export
#' @importFrom stats dnorm
dmnorm <- function(x, mu, sigma, lambda, log = FALSE, verbose = FALSE) {
    x <- as.vector(x)
    n_components <- length(sigma)
    if (verbose == TRUE) {
        message(paste(n_components, "components!", "\n"))
    }
    comp <- vector(mode = "list", length = n_components)
    for (i in seq_len(n_components)) {
        comp[[i]] <- lambda[i] * dnorm(x, mu[i], sigma[i], log = log)
    }
    return(Reduce("+", comp))
}

#' @title Defintions for location-scale t distribution
#' @description Internal functions for defining the t-distribution in terms
#'     of location-scale.
#' @param x,q The data vector
#' @param df Degrees of freedom
#' @param mu The location parameter
#' @param sigma The scale parameter
#' @param log Indicate whether probabilities are return as log
#' @importFrom stats dt
#' @describeIn dlst Probability Density Function
#' @keywords internal
#' @export
#' @examples
#'     val <- rnorm(10)
#'     dlst(val, df = 1, mu = 0, sigma = 1)
dlst <- function(x, df=1, mu=0, sigma=1, log = FALSE){
    prob <- (log(1) - log(sigma)) + dt((x - mu)/sigma, df, log = TRUE)
    if (log == TRUE){
        return(prob)
    } else {
        return(exp(prob))
    }
}

#' @describeIn dlst Cumulative distribution function
#' @importFrom stats pt
#' @export
#' @examples
#'     val <- rnorm(10)
#'     plst(q = val, df = 1, mu = 0, sigma = 1)
plst <- function(q, df=1, mu=0, sigma=1, log = FALSE){
    prob <- pt((q - mu)/sigma, df, log.p = log)
    return(prob)
}



#' @title Get the overall standard deviation of a
#'     two component mixture distribution
#' @param mu (Vector). A two value vector of mean values.
#' @param sigma (Vector). A two value vector of component-wise variances
#' @param lambda (Vector). A two value vector of component mixing coefficients
#' @param mean (Numeric Value). The overall mean.
#' @return A numeric value representing the overall standard deviation
#' @keywords internal
get_sd <- function(sigma, mu, mean, lambda) {
    sqrt(sum((sigma + mu - mean) * lambda))
}

#' @title Get the overall mean of a two component mixture distribution
#' @param mu (Vector). A two value vector of mean values.
#' @param lambda (Vector). A two value vector of component mixing coefficients
#' @return A numeric value representing the overall mean
#' @keywords internal
get_mean <- function(mu, lambda) {
    sum(lambda * mu)
}

#' @title Checking arguments of the function
#' @description This function extracts the parent environment (when called under the cbea function)
#'     and then check all the arguments.
#' @return None
#' @keywords internal
check_args <- function(){
    env <- parent.frame()

    # check arguments for control
    if (!is.null(env$control)){
        if ("fix_comp" %in% names(env$control) & (env$distr != "mnorm"| env$parametric == FALSE)) {
            stop("If fix_comp is part of control then distribution
                 has to be 'mnorm' and the fit is parametric")
        }
    }

    # output
    if (!env$output %in% c("cdf", "zscore", "pval", "sig", "raw")){
        stop("Output has to be of options 'cdf', 'zscore', 'pval', 'sig', 'raw'")
    }




    # first, check if distr is null
    if (is.null(env$distr)){
        if (env$parametric == TRUE){
            stop("Distribution needs to be specified if parametric fit is desired")
        }
    } else {
        if (!env$distr %in% c("norm", "mnorm", "lst")){
            stop("Distribution choices has to either be 'norm', 'mnorm', 'lst',
             or 'NULL' (equivalent to 'parametric' = FALSE)")
        }
    }

    # handling if adj argument is null
    if (is.null(env$adj)){
        if (env$parametric == TRUE){
            stop("Correlation adjustment option needs to be specified if parametric fit is desired")
        }
    }
    # if parametric is false cannot get either cdf values or z-scores
    if (env$parametric == FALSE){
        if (env$output %in% c("zscore", "cdf")){
            stop("Output cannot be either z-scores or CDF values if no parametric fit was performed")
        }
        # if parametric fit is false then needs to perform more permutations
        if (env$n_perm < 100){
            message("For non-parametric fits, the number of permutations should be higher (Rec: 100)")
        }
    }

}


