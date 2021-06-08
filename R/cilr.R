
#' Main cILR function to perform enrichment analysis using cILR method
#' @param ab_tab Named \code{n} by \code{p} matrix. This is the OTU/ASV/Strain table where taxa are columns.
#' @param set_list List of length \code{m}. This is a list of set membership by column names.
#' @param output String. The form of the output of the model.
#' @param distr String. The choice of distribution for the null.
#' @param adj Logical. Whether correlation adjustment procedure is utilized.
#' @param thresh Numeric. Threshold for significant p-values if \code{sig} is the output.
#' @param init Named List. Initialization parameters for estimating the null distribution. Default is NULL.
#' @param raw Logical. Whether scores are returned as raw (no parameter estimation step). Default is FALSE.
#' @param ... Named List. Additional arguments to be passed to \code{fitdistr} and \code{normmixEM}
#'
#' @return \code{R}    An \code{n} by \code{m} matrix of enrichment scores at the sample level
#'
#' @export
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfc
#' @importFrom rlang warn abort
cilr <- function(ab_tab, set_list,
                 output = c("cdf", "zscore", "pval", "sig"),
                 distr = c("mnorm", "norm"),
                 adj = TRUE,
                 thresh = 0.05,
                 init = NULL,
                 raw = FALSE, ...){
    # check arguments
    output <- match.arg(output)
    distr <- match.arg(distr)
    # check input
    if (!is.matrix(ab_tab)){
        rlang::warn("Coercing OTU/ASV table into matrix format")
        ab_tab <- as.matrix(ab_tab)
    }

    if (class(set_list) != "list"){
        rlang::warn("Set has to be converted to list format first")
        set_list <- as(set_list, "list")
    }

    p <- ncol(ab_tab) # number of features
    n <- nrow(ab_tab) # number of samples
    ab_perm <- ab_tab[,sample(1:p, size = p, replace = FALSE)]
    # 1. generate raw scores per set
    R <- purrr::map_dfc(set_list, ~{
        index <- which(colnames(ab_tab) %in% .x)
        if (raw == FALSE){
            raw_scores <- getScore(ab_tab, index)
            perm_scores <- getScore(ab_perm, index)
            if (adj == TRUE){
                # estimate perm and unperm distributions
                perm_distr <- estimate_distr(perm_scores, distr = distr, init = init, ...)
                unperm_distr <- estimate_distr(raw_scores, distr = distr, init = init, ...)
                # combine them into one final distribution
                final_distr <- combine_distr(perm_distr, unperm_distr, distr = distr)
            } else {
                # only estimate the distribution from the permuted data
                final_distr <- estimate_distr(perm_scores, distr = distr, init = init, ...)
            }
            # have to scale scores if raw is false (why else estimate the distribution)
            scores <- scale_score(raw_scores, method = output, param = final_distr, thresh = thresh)
        } else {
            # if raw is true then just get the score
            scores <- getScore(ab_tab, index)
        }
        return(scores)
    })
    colnames(R) <- names(set_list)
    return(R)
}
#' @title Get cILR scores for a given matrix and a vector of column indices
#' @param X (Matrix). OTU table of matrix format where taxa are columns and samples are rows
#' @param idx (Integer vector). Vector of integers indicating the column ids of taxa in a set
#' @export
get_score <- function(X, idx){
    # check type
    if (is.matrix(X) == F){
        message("Coercing X to matrix")
        X <- as.matrix(X)
    }
    if (is.vector(idx) == F){
        message("Coercing idx to a numeric vector")
        idx <- as.vector(idx)
    }
    # get total columns and define size ans scale funciton
    p <- ncol(X)
    size <- sum(idx == 1)
    scale <- sqrt(size * (p - size)/p)
    # calculate geometric mean
    num <- gmeanRow(X[,idx == 1])
    denom <- gmeanRow(X[,idx == 0])
    # return the ilr like statistic
    return(scale*(log(num/denom)))
}

#' @title Estimate distribution parameters from data
#' @description This function takes a numeric vector input and attempts to find the most optimal solution
#'     for the parameters of the distribution of choice. Right now only \code{norm} and \code{mnorm} distributions
#'     are supported.
#' @details The package \code{\link[fitdistrplus]{fitdistrplus}} is used to estimate parameters of the normal distribution
#'     while the package \code{\link[mixtools]{normalmixEM}} is used to estimate parameters of the mixture normal
#'     distribution. So far we suggest only estimating two components for the mixture normal distribution.
#'     For default options, we use mostly defaults from the packages themselves. The only difference was the mixture normal distribution
#'     where the convergence parameters were loosened and requiring more iterations to converge.
#' @param data (Numeric Vector). A vector of numbers that can be inputted to estimate the
#'     parameters of the distributional forms.
#' @param distr (String). The distribution to be fitted. Right now only \code{norm} or \code{mnorm} is supported
#' @param init (List). Initialization parameters for each distribution. For mixtures, each named
#'     element in the list should be a vector with length equal to the number of components
#' @param ... Other paremteres passed to fitdistrplus or normalmixEM
#' @keywords internal
#' @importFrom fitdistrplus fitdist
#' @importFrom mixtools normalmixEM
#' @importFrom rlist list.append
estimate_distr <- function(data, distr = c("mnorm", "norm"), init, ...){
    # matching argument
    distr <- match.arg(distr)
    # Default initialization is all zeroes
    if (missing(init) | is.null(init)){
        if (distr == "norm"){
            init <- list(mu = NULL, sd = NULL)
        } else if (distr == "mnorm"){
            init <- list(lambda = NULL, mu = NULL, sigma = NULL)
        }
    }

    # supplied and defaults for additional parameters
    supplied <- list(...)
    if (distr == "mnorm"){
        defaults <- list(maxrestarts=1000, epsilon = 1e-06, maxit= 1e5, arbmean = TRUE, k = 2)
    } else if (distr == "norm") {
        # so far there have no opinions
        defaults <- list(method = "mle")
    }
    params <- merge_lists(defaults = defaults, supplied = supplied)

    # Try catch would return NULL if the procedure errored out
    dist <- tryCatch(
        error = function(cnd) {
            message("There were errors with the fitting procedure, adjusting the fitting parameters might help")
            message(cnd)
            return(NULL)
        },
        {
            if (distr %in% c("norm")){
                params <- rlist::list.append(start = init, distr = distr, data = data)
                fit <- do.call(fitdistrplus::fitdist, params)
                list(mean = fit$estimate[['mean']], sd = fit$estimate[['sd']])
            } else if (distr == "mnorm"){
                params <- rlist::list.append(params, x = data)
                params <- c(init, params)
                fit <- do.call(normalmixEM, params)
                list(mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda)
            }
        }
    )
    return(dist)
}

#' @title Combining two distributions
#' @description Pass along handling of combining distributions to avoid clogging up the main function
#' @param perm (List). A list of parameters for permuted distribution
#' @param unperm (List). A list of parameters for the unpermuted distribtion
#' @param distr (String). Distribution of choice
#' @return A list of the combined distribution form based on the initial distribution of choice
#' @keywords internal
#' @importFrom rlang abort
combine_distr <- function(perm, unperm, distr){
    if (length(intersect(names(perm), names(unperm))) != length(perm)){
        rlang::abort("The two distributions to combine have to have the same parameters")
    }
    if (distr == "norm"){
        if (length(intersect(names(perm), c("mean", "sd"))) != 2){
            rlang::abort("Normal requires both mean and standard deviation")
        }
        final <- list(mean = perm$mean, sd = unperm$sd)
    } else if (distr == "mnorm"){
        if (length(intersect(names(perm), c("mu", "lambda", "sigma"))) != 3){
            rlang::abort("Mixture normal requires mu, lambda and sigma arguments")
        }
        check_perm <- all(vapply(perm, FUN = length, FUN.VALUE = 1) == 2)
        check_unperm <- all(vapply(unperm, FUN = length, FUN.VALUE = 1) == 2)
        if (check_perm == FALSE | check_unperm == FALSE){
            rlang::abort("Each named parameter much have values for each of the two components. Number of components restricted to 2")
        }
        final <- get_adj_mnorm(perm = perm, unperm = unperm)
    }
    return(combine_distr)
}

#' @title Function to perform the adjustment for the mixture normal distribution
#' @param perm (List). Parameter values of the distribution of scores computed on permuted data
#' @param unperm (List). Parameter values of the distribution of scores computed on unpermuted data
#' @return A List of parameters for the adjusted mixture normal.
#' @keywords internal
#' @importFrom stats optim
get_adj_mnorm <- function(perm, unperm, verbose = FALSE){
    # get the overall mean first
    perm_mean <- get_mean(mu = perm$mu, lambda = perm$lambda)
    unperm_mean <- get_mean(mu = unperm$mu, lambda = unperm$lambda)

    # get the overall standard deviation
    unperm_sd <- get_sd(sigma = unperm$sigma, mu = unperm$mu,
                        mean = unperm_mean, lambda = unperm$lambda)
    perm_sd <- get_sd(sigma = perm$sigma, mu = perm$mu, mean = perm_mean,
                      lambda = perm$lambda)
    # define objective function
    obj_function <- function(sigma, mu, lambda, mean, sd){
        s1 <- sigma[1]
        s2 <- sigma[2]
        m1 <- mu[1]
        m2 <- mu[2]
        l1 <- lambda[1]
        l2 <- lambda[2]
        estimate_sd <- sqrt((s1 + m1 - mean)*l1 + (s2 + m2 - mean)*l2)

        # estimate the objective funciton which is an l2 norm
        obj <- sqrt(sum((estimate_sd - sd)^2))
        return(obj)
    }
    opt <- stats::optim(par = c(0.01,0.01),
                 obj_function,
                 method = "L-BFGS-B",
                 lower = c(1e-5,1e-5),
                 upper = c(Inf, Inf),
                 mu = perm$mu, lambda = perm$lambda, mean = perm_mean, sd = unperm_sd)
    estim_sd <- get_sd(sigma = opt$par, lambda = perm$lambda, mu = perm$mu, mean = perm_mean)
    if (verbose == TRUE){
        print(paste("Total sd is", unperm_sd,"and estimated sd is", estim_sd))
    }
    param <- list(mu = perm$mu, sigma = opt$par, lambda = perm$lambda)
    return(param)
}


#' @title This function handles the ability to merge supplied and defaults
#' @param defaults (List). Default options
#' @param supplied (List). Supplied options
#' @return A merged list
#' @keywords internal
merge_lists <- function(defaults, supplied){
    similar_idx <- which(names(defaults) %in% names(supplied))
    if (length(similar_idx) == 0){
        merged <- c(defaults, supplied)
    }
    else {
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
#' @describeIn pmnorm Cumulative Distribution Function
#' @importFrom stats pnorm
#' @export
pmnorm <- function(q, mu, sigma, lambda, log = FALSE, verbose = FALSE){
    q <- as.vector(q)
    n_components <- length(sigma)
    if (verbose == TRUE){
        cat(paste(n_components, "components!", "\n"))
    }

    comp <- vector(mode = "list", length = n_components)
    for (i in 1:n_components){
        comp[[i]] <- lambda[i] * stats::pnorm(q, mu[i], sigma[i], log.p = log)
    }
    return(Reduce("+", comp))
}

#' @describeIn pmnorm Probability Density Function
#' @importFrom stats dnorm
#' @export
dmnorm <- function(x, mu, sigma, lambda, log=FALSE, verbose = FALSE){
    x <- as.vector(x)
    n_components <- length(sigma)
    if (verbose == TRUE){
        cat(paste(n_components, "components!", "\n"))
    }
    comp <- vector(mode = "list", length = n_components)
    for (i in 1:n_components){
        comp[[i]] <- lambda[i] * stats::dnorm(x, mu[i], sigma[i], log = log)
    }
    return(Reduce("+", comp))
}

#' @title Get the overall standard deviation of a two component mixture distribution
#' @param mu (Vector). A two value vector of mean values.
#' @param sigma (Vector). A two value vector of component-wise variances
#' @param lambda (Vector). A two value vector of component mixing coefficients
#' @param mean (Numeric Value). The overall mean.
#' @return A numeric value representing the overall standard deviation
#' @keywords internal
get_sd <- function(sigma, mu, mean, lambda){
    sqrt(sum((sigma + mu - mean)*lambda))
}

#' @title Get the overall mean of a two component mixture distribution
#' @param mu (Vector). A two value vector of mean values.
#' @param lambda (Vector). A two value vector of component mixing coefficients
#' @return A numeric value representing the overall mean
#' @keywords internal
get_mean <- function(mu, lambda){
    sum(lambda * mu)
}

