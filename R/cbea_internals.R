#' @title Internal cbea function
#' @description See main function \code{cbea} documentation for more details.
#' @param ab_tab Named \code{n} by \code{p} matrix.
#'     This is the OTU/ASV/Strain table where taxa are columns.
#' @param set_list List of length \code{m}.
#'     This is a list of set membership by column names.
#' @param output See documentation \code{\link{cbea}}
#' @param distr See documentation \code{\link{cbea}}
#' @param adj See documentation \code{\link{cbea}}
#' @param thresh See documentation \code{\link{cbea}}
#' @param init See documentation \code{\link{cbea}}
#' @param control See documentation \code{\link{cbea}}
#' @param ... See documentation \code{\link{cbea}}
#' @importFrom magrittr %>%
#' @importFrom purrr map_dfc
#' @importFrom rlang warn abort
#' @return A \code{data.frame} of size \code{n} by \code{m}.
#'     \code{n} is the total number of samples and \code{m}
#'     is the total number of sets with elements represented in the data.
#' @keywords internal
.cbea <- function(ab_tab,
                  set_list,
                  output, distr,
                  adj = TRUE,
                  thresh = 0.05,
                  init = NULL,
                  control = NULL, ...) {
    p <- ncol(ab_tab) # number of features
    n <- nrow(ab_tab) # number of samples
    ab_perm <- ab_tab[, sample(seq_len(p), size = p, replace = FALSE)]

    # Loop through each set
    R <- purrr::map_dfc(set_list, ~ {
        # first, retrieve indices for the set of interest
        index <- which(colnames(ab_tab) %in% .x)
        if (output != "raw") {
            raw_scores <- get_score(ab_tab, index)
            perm_scores <- get_score(ab_perm, index)
            if (adj == TRUE) {
            # estimate perm and unperm distributions
                perm_distr <- estimate_distr(perm_scores,
                                distr = distr, init = init,
                                args_list = control)
                unperm_distr <- estimate_distr(raw_scores,
                                distr = distr, init = init,
                                args_list = control)
            # combine them into one final distribution
                final_distr <- combine_distr(perm_distr,
                                             unperm_distr, distr = distr)
            } else {
                # only estimate the distribution from the permuted data
                final_distr <- estimate_distr(perm_scores, distr = distr,
                                                  init = init,
                                                  args_list = control)
            }
            # have to scale scores if raw is false
            scores <- scale_scores(raw_scores, method = output,
                                     param = final_distr,
                                     thresh = thresh, distr = distr)
        } else {
        # if raw is true then just get the score
            scores <- get_score(ab_tab, index)
        }
        return(scores)
    })
    colnames(R) <- names(set_list)
    R <- tibble::add_column(R, sample_id = rownames(ab_tab), .before = 1)
    return(R)
}


#' @title Get CBEA scores for a given matrix and a vector of column indices
#' @param X (Matrix). OTU table of matrix format where taxa are
#'    columns and samples are rows
#' @param idx (Integer vector). Vector of integers indicating
#'    the column ids of taxa in a set
#' @return A matrix of size \code{n} by \code{1} where \code{n}
#'    is the total number of samples
#' @examples
#' data(hmp_gingival)
#' seq <- hmp_gingival$data
#' seq_matrix <- as(phyloseq::otu_table(seq), "matrix")
#' seq_matrix <- t(seq_matrix) + 1
#' rand_set <- sample(seq_len(ncol(seq_matrix)), size = 10)
#' scores <- get_score(X = seq_matrix, idx = rand_set)
#' @export
get_score <- function(X, idx) {
  # check type
    if (is.matrix(X) == FALSE) {
        message("Coercing X to matrix")
        X <- as.matrix(X)
    }
    if (is.vector(idx) == FALSE) {
        message("Coercing idx to a numeric vector")
        idx <- as.vector(idx)
    }
    # get total columns and define size ans scale funciton
    p <- ncol(X)
    size <- length(idx)
    scale <- sqrt(size * (p - size) / p)
    # calculate geometric mean
    num <- gmeanRow(as.matrix(X[, idx]))
    denom <- gmeanRow(as.matrix(X[, -idx]))
    # return the ilr like statistic
    return(scale * (log(num / denom)))
}

#' @title Estimate distribution parameters from data
#' @description This function takes a numeric vector input and
#'     attempts to find the most optimal solution
#'     for the parameters of the distribution of choice.
#'     Right now only \code{norm} and \code{mnorm} distributions
#'     are supported.
#' @details The package \code{\link[fitdistrplus]{fitdistrplus}} is used to
#'     estimate parameters of the normal distribution
#'     while the package \code{\link[mixtools]{normalmixEM}} is used to
#'     estimate parameters of the mixture normal distribution.
#'     So far we suggest only estimating two components for the mixture
#'     normal distribution.
#'     For default options, we use mostly defaults from the packages themselves.
#'     The only difference was the mixture normal distribution
#'     where the convergence parameters were loosened and requiring
#'     more iterations to converge.
#' @param data (Numeric Vector). A vector of numbers that can be inputted
#'     to estimate the parameters of the distributional forms.
#' @param distr (String). The distribution to be fitted.
#'     Right now only \code{norm} or \code{mnorm} is supported
#' @param init (List). Initialization parameters for each distribution.
#'     For mixtures, each named element in the list should be a vector
#'     with length equal to the number of components
#' @param args_list (List). Named list of additional arguments
#'     passed onto fitdist and normalmixEM
#' @param ... Other paremteres passed to fitdistrplus or normalmixEM
#' @return A named list with all the parameter names and values
#' @keywords internal
#' @importFrom fitdistrplus fitdist
#' @importFrom mixtools normalmixEM
#' @importFrom rlist list.append
#' @importFrom stats na.omit
#' @importFrom rlang abort
#' @importFrom utils capture.output
estimate_distr <- function(data, distr = c("mnorm", "norm"),
                           init = NULL, args_list = NULL) {
    distr <- match.arg(distr)
    # check if data has NAs
    if (any(is.na(data))) {
        message("There are NAs in the data vector, omitting NA values")
        org_length <- length(data)
        data <- stats::na.omit(data)
        if (length(data) < 0.5 * org_length) {
            rlang::abort("More than 50% of the data is NA,
                       aborting distribution fitting")
        }
    }
    # Default initialization is all zeroes
    if (missing(init) | is.null(init)) {
        if (distr == "norm") {
            init <- list(mean = 0, sd = 1)
        } else if (distr == "mnorm") {
            init <- list(lambda = NULL, mu = NULL, sigma = NULL)
        }
    }
    # supplied and defaults for additional parameters
    if (distr == "mnorm") {
        defaults <- list(maxrestarts = 1000, epsilon = 1e-06,
                     maxit = 1e5, arbmean = TRUE, k = 2)
    } else if (distr == "norm") {
    # so far there have no opinions
        defaults <- list(method = "mle", fix.arg = NULL,
                     discrete = FALSE, keepdata = FALSE, keepdata.nb = 100)
    }
    params <- merge_lists(defaults = defaults, supplied = args_list)
    # Try catch would return NULL if the procedure errored out
    dist <- tryCatch(
        error = function(cnd) {
              message("Adjusting the fitting parameters might
                      help with the fitting")
              message(cnd)
              return(NULL)
        }, {
            if (distr %in% c("norm")) {
                params <- rlist::list.append(params, start = init,
                                         distr = distr, data = data)
                fit <- do.call(fitdistrplus::fitdist, params)
                list(mean = fit$estimate[["mean"]], sd = fit$estimate[["sd"]])
            } else if (distr == "mnorm") {
                params <- rlist::list.append(params, x = data)
                params <- c(init, params)
                fit <- do.call(normalmixEM, params)
                list(mu = fit$mu, sigma = fit$sigma, lambda = fit$lambda)
            }
        }
    )
    return(dist)
}

#' @title Scaling scores based on estimated null distribution
#' @param scores (Numeric Vector). Raw CBEA scores
#'    generated without permutations
#' @param method (String). The final form that the user want to return.
#'    Options include \code{cdf}, \code{zscore}, \code{pval} and \code{sig}.
#' @param param (List). The parameters of the estimated null distribution.
#'    Names must match distribution
#' @param thresh (Numeric). The threshold to decide whether a set
#'    is significantly enriched. Only available if \code{method} is \code{sig}
#' @return A vector of size \code{n} where \code{n} is the sample size
#' @keywords internal
#' @importFrom rlist list.append
#' @importFrom rlang abort
scale_scores <- function(scores, method = c("cdf", "zscore", "pval", "sig"),
                         param, distr, thresh = 0.05) {
    # detect if parameter length is concordant with distribution type
    if (distr == "norm") {
        f <- "pnorm"
        if (length(intersect(names(param), c("mean", "sd"))) != 2) {
            rlang::abort("Normal requires both mean and standard deviation")
        }
    } else {
        f <- "pmnorm"
        if (length(intersect(names(param), c("mu", "lambda", "sigma"))) != 3) {
            rlang::abort("Mixture normal requires mu,
                         lambda and sigma arguments")
        }
        check_param <- all(vapply(param, FUN = length, FUN.VALUE = 1) == 2)
        if (check_param == FALSE) {
            rlang::abort("Each named parameter much have values for
                       each of the two components.
                       Number of components restricted to 2")
        }
    }
    if (method %in% c("cdf", "sig", "pval")) {
        param <- rlist::list.append(q = as.vector(scores), param)
        scale <- do.call(f, param)
        if (sum(is.na(scale)) >= 1) {
            cat("There are NA values here")
            if (file.exists("output") == FALSE) {
                dir.create("output")
            }
            saveRDS(param, file = "output/null_values.rds")
        }
        if (method == "pval") {
            scale <- 1 - scale
        } else if (method == "sig") {
            scale <- 1 - scale
            scale <- ifelse(scale <= thresh, 1, 0)
        }
    } else if (method == "zscore") {
        if (f == "pmnorm") {
            mean <- get_mean(mu = param$mu, lambda = param$lambda)
            sd <- get_sd(sigma = param$sigma, mu = param$mu,
                       mean = mean, lambda = param$lambda)
        } else if (f == "pnorm") {
            mean <- param$mean
            sd <- param$sd
        }
        scale <- (scores - mean) * 1 / sd
    }
    return(scale)
}


#' @title Combining two distributions
#' @description Pass along handling of combining distributions to avoid
#'    clogging up the main function
#' @param perm (List). A list of parameters for permuted distribution
#' @param unperm (List). A list of parameters for the unpermuted distribtion
#' @param distr (String). Distribution of choice
#' @return A list of the combined distribution form based on
#'    the initial distribution of choice
#' @keywords internal
#' @importFrom rlang abort
combine_distr <- function(perm, unperm, distr) {
  # If unable to estimate anything, then return NULL final distribution
    if (is.null(perm) | is.null(unperm)) {
        return(NULL)
    }
    if (length(intersect(names(perm), names(unperm))) != length(perm)) {
        rlang::abort("The two distributions to combine
                     have to have the same parameters")
    }
    if (distr == "norm") {
        if (length(intersect(names(perm), c("mean", "sd"))) != 2) {
            rlang::abort("Normal requires both mean and standard deviation")
        }
        final <- list(mean = perm$mean, sd = unperm$sd)
    } else if (distr == "mnorm") {
        if (length(intersect(names(perm), c("mu", "lambda", "sigma"))) != 3) {
            rlang::abort("Mixture normal requires mu,
                         lambda and sigma arguments")
        }
        check_perm <- all(vapply(perm, FUN = length, FUN.VALUE = 1) == 2)
        check_unperm <- all(vapply(unperm, FUN = length, FUN.VALUE = 1) == 2)
        if (check_perm == FALSE | check_unperm == FALSE) {
              rlang::abort("Each named parameter much have values for
                           each of the two components.
                           Number of components restricted to 2")
        }
        final <- get_adj_mnorm(perm = perm, unperm = unperm)
    }
    return(final)
}

#' @title Function to perform the adjustment for the mixture normal distribution
#' @param perm (List). Parameter values of the distribution of scores
#;    computed on permuted data
#' @param unperm (List). Parameter values of the distribution of scores
#'    computed on unpermuted data
#' @return A List of parameters for the adjusted mixture normal.
#' @keywords internal
#' @importFrom stats optim
get_adj_mnorm <- function(perm, unperm, verbose = FALSE) {
    # get the overall mean first
    perm_mean <- get_mean(mu = perm$mu, lambda = perm$lambda)
    unperm_mean <- get_mean(mu = unperm$mu, lambda = unperm$lambda)

    # get the overall standard deviation
    unperm_sd <- get_sd(
        sigma = unperm$sigma, mu = unperm$mu,
        mean = unperm_mean, lambda = unperm$lambda
    )
    perm_sd <- get_sd(
        sigma = perm$sigma, mu = perm$mu, mean = perm_mean,
        lambda = perm$lambda
    )
    # define objective function
    obj_function <- function(sigma, mu, lambda, mean, sd) {
        s1 <- sigma[1]
        s2 <- sigma[2]
        m1 <- mu[1]
        m2 <- mu[2]
        l1 <- lambda[1]
        l2 <- lambda[2]
        estimate_sd <- sqrt((s1 + m1 - mean) * l1 + (s2 + m2 - mean) * l2)

        # estimate the objective funciton which is an l2 norm
        obj <- sqrt(sum((estimate_sd - sd)^2))
        return(obj)
    }
    opt <- stats::optim(
        par = c(0.01, 0.01),
        obj_function,
        method = "L-BFGS-B",
        lower = c(1e-5, 1e-5),
        upper = c(Inf, Inf),
        mu = perm$mu, lambda = perm$lambda, mean = perm_mean, sd = unperm_sd
    )
    estim_sd <- get_sd(sigma = opt$par, lambda = perm$lambda,
                        mu = perm$mu, mean = perm_mean)
    if (verbose == TRUE) {
        print(paste("Total sd is", unperm_sd, "and estimated sd is", estim_sd))
    }
    param <- list(mu = perm$mu, sigma = opt$par, lambda = perm$lambda)
    return(param)
}
