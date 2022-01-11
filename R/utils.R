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

