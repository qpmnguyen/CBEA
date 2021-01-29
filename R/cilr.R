# This script contains the main cilr functions
# Quang Nguyen
# Jan 27th 2021



#' @title Applying the cILR method to a data.frame object and pre-defined sets
#' @description This function performs the cILR aggregation algorithm as defined in the Nguyen et al. paper. This
#'    function takes in two inputs: an otu table as a \code{matrix} object where rows are samples and columns are
#'    taxa; a named \code{list} where each list element is a numeric vector of indices for the columns of interest
#' @param data \code{matrix}. The otu or asv table in data.frame format
#' @param sets \code{list}. The sets exists as named list with set names as strings and column index as integer
#' @param distr  \code{string}. What is the distribution, required only when \code{estim} is TRUE,
#'     can only be \code{norm} or \code{mnorm}
#' @param output \code{string}. What is the preferred output, can be either
#'     \code{zscores}, \code{cdf}, \code{pvalues}, or \code{sig}
#' @param min_size \code{numeric}. What is the minimum size of a set
#' @param n_perm \code{numeric}. Number of resamples for the data set for good estimation, defaults to 5.
#' @param adj \code{logical}. Whether correlation adjustment should be done.
#' @param thresh \code{numeric}. Threshold to consider rejecting the null hypothesis based on p-values.
#'     Only relevant when output is \code{sig}. Defaults to 0.05
#' @param p_adj \code{logical}. Whether to adjust p-values using the BH procedure. Defaults to FALSE
#' @param init \code{list}. List of initialization parameters for the fitting portion, if not specified, defaults to
#'     default options in \code{mixtools} and \code{fitdistrplus} packages
#' @param ... Additional parameters passed to fitdistrplus or mixtools functions
#' @return A \code{data.frame} object with rows equals the number of samples and
#'     columns specifies the number of samples defined in the set
#' @export
cilr <- function(data, sets, output = NULL, distr=c("mnorm", "norm"), min_size = 2, nperm = 5,
                 adj = TRUE, thresh=NULL, init = NULL, ...){
    # first, let's check for types
    if (!"matrix" %in% class(data)) rlang::warn(message = "Coercing data to matrix format")
    if (!"list" %in% class(sets)) rlang::stop(message = "Set definition must be a named list of numeric indices")
    if (is.null(names(test))) rlang::stop(message = "Set definition must be a named list of numeric indices")
    typecheck <- sapply(sets, typeof) == "double"
    if (all(typecheck, TRUE)) rlang::stop(message = "Set definition must be a named list of integer indices")

    # Drop all sets with min_size
    sets <- sets[sapply(sets, length) >= min_size]

    # Save ncol and nrows
    p <- ncol(data)
    n <- nrow(data)

    # Initialize the final matrix
    R <- matrix(nrow = nrow(X), ncol = length(sets))
    # Generate permuted data
    perm <- map_df(seq_len(n_perm), ~{
        data[,sample(seq_len(p), replace = FALSE)]
    })
    if (adj == TRUE){
        unperm <- map_df(seq_len(n_perm), ~{
            data[sample(seq_len(n), replace = T),]
        })
    }

    # Loop
    for (i in seq_along(sets)){
        score <- get_score(data, sets[[i]])
        perm_score <- get_score(perm, sets[[i]])
        perm_dist <- estimate_distr(perm_score, distr = distr, init = init, ...)
        if (adj == TRUE){
            unperm_score <- get_score(unperm, sets[[i]])
            unperm_dist <- estimate_distr(unperm_score, distr = distr, init = init, ...)
            if (distr == "norm"){
                final_distr <- list(mean = perm_dist$mean, sd = unperm_dist$sd)
            } else if (distr == "mnorm"){
                final_distr <- get_adj_mnorm(perm = perm_dist, unperm = unperm_dist)
            }
        } else {
            final_distr <- perm_dist
        }
        final_score <- transform_scores(score, method = output, param = final_distr,
                                  thresh = thresh)
        R[,i] <- score
    }

    rownames(R) <- rownames(X)
    colnames(R) <- names(sets)

    R <- as.data.frame(R)
    return(R)
}

#' @title Internal function calculate scores based on a certain index
#' @param data \code{matrix}. The otu table where taxa are columns and samples are rows
#' @param set_idx \code{vector}. Numeric vector of index of columns belonging to a certain set
get_score <- function(data, set_idx){
    # get total columns and define size ans scale funciton
    p <- get("p", parent.frame()) # so we don't have to compute p twice
    size <- length(set_idx)
    scale <- sqrt(size * (p - size)/p)
    # calculate geometric mean
    num <- compositions::geometricmeanRow(data[,set_idx])
    denom <- compositions::geometricmeanRow(data[,-set_idx])
    score <- scale * log(num/denom)
    return(score)
}

#' @title Function to throw error in the distribution fitting procedure
#' @param cond Condition of the function
#' @return Throw warning and return null
throw_error <- function(cond){
    rlang::warn("There were errors with the fitting procedure, adjusting the fitting parameters might help")
    rlang::inform(cond)
    return(NULL)
}


#' @title Estimating distribution given a vector of scores
#' @param data \code{Matrix} otu table in matrix format where taxa are columns
#' @param distr \code{string} the distribution of the function
#' @param init \code{list} List of starting parameters
#' @param ... Additional arguments to normix and fitdistrplus
estimate_distr <- function(data, distr = c("mnorm", "norm"), init, ...){
    # mnorm have suggested defaults
    def_mnorm <- list(
        maxrestarts=1000,
        epsilon = 1e-06,
        maxit= 1e5
    )
    # matching argument
    distr <- match.arg(distr)

    if (missing(init) | is.null(init)){
        if (distr == "norm"){
            init <- NULL
        } else if (distr == "mnorm"){
            init <- list(lambda = NULL, mu = NULL, sigma = NULL)
        }
    }
    dist <- tryCatch({
        if (distr == "norm"){
            fit <- fitdistrplus::fitdist(data, distr = distr, method = "mle", start = init)
            list(mean = fit$estimate[['mean']], sd = fit$estimate[['sd']])
        } else if (distr == "mnorm"){
            params <- rlist::list.append(init, x = data)
            params <- c(params, list(...)) # grabbing arguments to put in normalmixEM
            fit <- do.call(normalmixEM, params)
            list(mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda)
        }
    }, error = throw_error)
    # if you cannot estimate the distribution of dist
    if (is.null(dist)){
        return(NULL)
    } else {
        return(dist)
    }
}





