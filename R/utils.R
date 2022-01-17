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
#' @importFrom stats pnorm
pmnorm <- function(q, mu, sigma, lambda, log = FALSE, verbose = FALSE) {
    q <- as.vector(q)
    n_components <- length(sigma)
    if (verbose == TRUE) {
        message(paste(n_components, "components!", "\n"))
    }
    comp <- vector(mode = "list", length = n_components)
    for (i in seq_len(n_components)) {
        comp[[i]] <- lambda[i] * stats::pnorm(q, mu[i], sigma[i], log.p = log)
    }
    return(Reduce("+", comp))
}

#' @describeIn pmnorm Probability Density Function
#' @importFrom stats dnorm
dmnorm <- function(x, mu, sigma, lambda, log = FALSE, verbose = FALSE) {
    x <- as.vector(x)
    n_components <- length(sigma)
    if (verbose == TRUE) {
        message(paste(n_components, "components!", "\n"))
    }
    comp <- vector(mode = "list", length = n_components)
    for (i in seq_len(n_components)) {
        comp[[i]] <- lambda[i] * stats::dnorm(x, mu[i], sigma[i], log = log)
    }
    return(Reduce("+", comp))
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
#' @param args (List) A list of parameters
check_args <- function(args){
    arg_names <- c("obj", "set", "output", 
                   "distr", "adj", "n_perm", "parametric", "init", "control")
    if (intersect(names(args), arg_names) != length(arg_names)){
        stop("Please specify all required inputs")
    } 
    
    # first, check if distr is null 
    if (is.null(args$distr)){
        if (args$parametric == TRUE){
            stop("Distribution needs to be specified if parametric fit is desired")
        } 
    }
    # handling if adj argument is null 
    if (is.null(args$adj)){
        if (args$parametric == TRUE){
            stop("Correlation adjustment option needs to be specified if parametric fit is desired")
        }
    }
    # if parametric is false cannot get either cdf values or z-scores
    if (args$parametric == FALSE){
        if (args$output %in% c("zscore", "cdf")){
            stop("Output cannot be either z-scores or CDF values if no parametric fit was performed")
        }
        # if parametric fit is false then needs to perform more permutations 
        if (args$n_perm <= 200){
            message("For non-parametric fits, the number of permutations should be higher (Rec: 200)")
        }
    }
    
}

