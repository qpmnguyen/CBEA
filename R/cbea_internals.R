#' @title Internal cbea function
#' @description See main function \code{cbea} documentation for more details.
#' @param ab_tab (Matrix). Named \code{n} by \code{p} matrix.
#'     This is the OTU/ASV/Strain table where taxa are columns.
#' @param set_list (List). List of length \code{m}.
#'     This is a list of set membership by column names.
#' @param output See documentation \code{\link{cbea}}
#' @param distr See documentation \code{\link{cbea}}
#' @param adj See documentation \code{\link{cbea}}
#' @param n_perm See documentation \code{\link{cbea}}
#' @param parametric See documentation \code{\link{cbea}}
#' @param thresh See documentation \code{\link{cbea}}
#' @param init See documentation \code{\link{cbea}}
#' @param control See documentation \code{\link{cbea}}
#' @param parallel_backend See documentation \code{\link{cbea}}
#' @param ... See documentation \code{\link{cbea}}
#' @importFrom magrittr %>%
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom tibble as_tibble add_column
#' @importFrom stats rpois
#' @return A \code{data.frame} of size \code{n} by \code{m}.
#'     \code{n} is the total number of samples and \code{m}
#'     is the total number of sets with elements represented in the data.
#' @keywords internal
.cbea <- function(ab_tab,
                  set_list,
                  output, distr,
                  adj = TRUE,
                  n_perm = 100,
                  parametric = TRUE,
                  thresh = 0.05,
                  init = NULL,
                  control = NULL,
                  parallel_backend = NULL,
                  ...) {
    if (is.null(parallel_backend)){
        seed <- rpois(1, lambda = 1e4)
        parallel_backend <- SerialParam(RNGseed = seed)
    }
    R <- bplapply(set_list, fit_scores, ab_tab = ab_tab, adj = adj, distr = distr,
                  output = output, n_perm = n_perm,
                  parametric = parametric, init = init, control = control,
                  thresh = thresh, BPPARAM = parallel_backend)
    R <- as_tibble(R)
    colnames(R) <- names(set_list)
    R <- add_column(R, sample_id = rownames(ab_tab), .before = 1)
    return(R)
}

#' @title Function to compute CBEA scores for each set
#' @param ab_tab (Matrix). Named \code{n} by \code{p} matrix.
#'     This is the OTU/ASV/Strain table where taxa are columns.
#' @param index_vec (Character Vector). A character vector indicating the
#'     elements of the set of interest
#' @param adj (Logical). See documentation \code{\link{cbea}}
#' @param distr (Character). See documentation \code{\link{cbea}}
#' @param output (Character). See documentation \code{\link{cbea}}
#' @param init (List). See documentation \code{\link{cbea}}
#' @param control (List). See documentation \code{\link{cbea}}
#' @param n_perm (Numeric). The total number of permutations.
#' @param parametric (Logical). See documentation \code{\link{cbea}}
#' @param thresh (Numeric). See documentation \code{\link{cbea}}
#' @keywords internal
#' @importFrom stats quantile
fit_scores <- function(index_vec, ab_tab, adj, distr, output, n_perm, parametric, thresh,
                       init, control){
    if (!is.null(control) & ("fix_comp" %in% names(control))){
        fix_comp <- control$fix_comp
        control <- control[-which(names(control) == "fix_comp")]
        if (length(control) == 0){
            control <- NULL
        }
    } else {
        fix_comp <- "none"
    }
    # get the true index
    true_index <- which(colnames(ab_tab) %in% index_vec)
    # true_index is different than index_vec bc index_vec can have more elements not present
    # in the data set (taxa not present here but still part of the category)
    raw_scores <- get_raw_score(ab_tab, true_index)
    if (output == "raw"){
        scores <- raw_scores
    } else {
        p <- ncol(ab_tab) # total number of taxa
        # first, get a list of indices to shuffle
        # get random indices the same size as the true index
        perm_index_list <- replicate(n_perm,
                                     sample(seq_len(p), size = length(true_index),
                                            replace = FALSE),
                                     simplify = FALSE)
        # get perm scores fromt he perm_index_list
        perm_scores <- lapply(perm_index_list, function(x) get_raw_score(ab_tab, x))
        perm_scores <- unname(do.call(c, perm_scores))

        # if parametric is true estimate distribution
        if (parametric == TRUE) {
            perm_distr <- estimate_distr(perm_scores,
                                         distr = distr,
                                         init = init,
                                         args_list = control)
            if (adj == TRUE) {
                unperm_distr <- estimate_distr(raw_scores, distr = distr,
                                               init = init, args_list = control)
                final_distr <- combine_distr(perm_distr, unperm_distr, distr = distr,
                                             fix_comp = fix_comp)
            } else {
                final_distr <- perm_distr
            }
            scores <- scale_scores(raw_scores,
                                   method = output,
                                   param = final_distr,
                                   thresh = thresh, distr = distr)
        } else {
            if (output == "sig"){
                scores <- ifelse(raw_scores >= quantile(perm_scores, 1 - thresh),1,0)
            } else if (output == "pval") {
                scores <- vapply(raw_scores, function(x) sum(perm_scores >= x)/length(perm_scores), FUN.VALUE = 0.1)
            }
        }
    }
    return(scores)
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
#' seq_matrix <- SummarizedExperiment::assays(seq)[[1]]
#' seq_matrix <- t(seq_matrix) + 1
#' rand_set <- sample(seq_len(ncol(seq_matrix)), size = 10)
#' scores <- get_raw_score(X = seq_matrix, idx = rand_set)
#' @export
get_raw_score <- function(X, idx) {
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
#' @importFrom stats na.omit
#' @importFrom extraDistr dlst plst qlst rlst
estimate_distr <- function(data, distr,
                           init = NULL, args_list = NULL) {
    distr <- match.arg(distr, c("norm", "mnorm", "lst"))
    # check if data has NAs
    if (any(is.na(data))) {
        message("There are NAs in the data vector, omitting NA values")
        org_length <- length(data)
        data <- na.omit(data)
        if (length(data) < 0.5 * org_length) {
            stop("More than 50% of the data is NA,
                       aborting distribution fitting")
        }
    }
    # Default initialization is all zeroes
    if (missing(init) | is.null(init)) {
        if (distr == "norm") {
            init <- list(mean = 0, sd = 1)
        } else if (distr == "mnorm") {
            init <- list(lambda = NULL, mu = NULL, sigma = NULL)
        } else if (distr == "lst") {
            init <- list(df = 1, mu = 0, sigma = 1)
        }
    }
    # supplied and defaults for additional parameters
    if (distr == "mnorm") {
        defaults <- list(
            maxrestarts = 1000, epsilon = 1e-06,
            maxit = 1e5, arbmean = TRUE, k = 2
        )
    } else if (distr %in% c("norm", "lst")) {
        # so far there have no opinions
        defaults <- list(
            method = "mle", fix.arg = NULL,
            discrete = FALSE, keepdata = FALSE, keepdata.nb = 100
        )
    }
    params <- merge_lists(defaults = defaults, supplied = args_list)
    # Try catch would return NULL if the procedure errored out
    dist <- tryCatch(
        error = function(cnd) {
            message("Adjusting the fitting parameters might
                      help with the fitting")
            message(cnd)
            return(NULL)
        },
        {
            if (distr %in% c("norm", "lst")) {
                print(exists("dlst", mode = "function"))
                params <- c(params, list(start = init, distr = distr, data = data))
                fit <- do.call(fitdist, params)
                if (distr == "norm"){
                    list(mean = fit$estimate[["mean"]], sd = fit$estimate[["sd"]])
                } else if (distr == "lst"){
                    list(df = fit$estimate[["df"]], mu = fit$estimate[["mu"]],
                         sigma = fit$estimate[["sigma"]])
                }
            } else if (distr == "mnorm") {
                params <- c(params, list(x = data))
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
#'    Names must match distribution.
#' @param thresh (Numeric). The threshold to decide whether a set
#'    is significantly enriched. Only available if \code{method} is \code{sig}
#' @return A vector of size \code{n} where \code{n} is the sample size
#' @keywords internal
scale_scores <- function(scores, method,
                         param, distr, thresh = 0.05) {
    distr <- match.arg(distr, c("norm", "mnorm", "lst"))
    method <- match.arg(method, c("cdf", "zscore", "pval", "sig"))
    # detect if parameter length is concordant with distribution type
    if (distr == "norm") {
        f <- "pnorm"
        if (length(intersect(names(param), c("mean", "sd"))) != 2) {
            stop("Normal requires both mean and standard deviation")
        }
    } else if (distr == "mnorm"){
        f <- "pmnorm"
        if (length(intersect(names(param), c("mu", "lambda", "sigma"))) != 3) {
            stop("Mixture normal requires mu,
                         lambda and sigma arguments")
        }
        check_param <- all(vapply(param, FUN = length, FUN.VALUE = 1) == 2)
        if (check_param == FALSE) {
            stop("Each named parameter much have values for
                       each of the two components.
                       Number of components restricted to 2")
        }
    }  else if (distr == "lst"){
        f <- "plst"
        if (length(intersect(names(param), c("mu", "df", "sigma"))) != 3) {
            stop("Location-scale t-distribution requires mu, sigma, and df arguments")
        }
    }
    # compute values
    if (method %in% c("cdf", "sig", "pval")) {
        param <- c(param, list(q = as.vector(scores)))
        scale <- do.call(f, param)
        if (sum(is.na(scale)) >= 1) {
            message("There are NA values here")
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
            sd <- get_sd(
                sigma = param$sigma, mu = param$mu,
                mean = mean, lambda = param$lambda
            )
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
combine_distr <- function(perm, unperm, distr, ...) {
    distr <- match.arg(distr, c("norm", "mnorm", "lst"))
    # If unable to estimate anything, then return NULL final distribution
    if (is.null(perm) | is.null(unperm)) {
        return(NULL)
    }
    if (length(intersect(names(perm), names(unperm))) != length(perm)) {
        stop("The two distributions to combine
                     have to have the same parameters")
    }
    if (distr == "norm") {
        if (length(intersect(names(perm), c("mean", "sd"))) != 2) {
            stop("Normal requires both mean and standard deviation")
        }
        final <- list(mean = perm$mean, sd = unperm$sd)
    } else if (distr == "mnorm") {
        if (length(intersect(names(perm), c("mu", "lambda", "sigma"))) != 3) {
            stop("Mixture normal requires mu,
                         lambda and sigma arguments")
        }
        check_perm <- all(vapply(perm, FUN = length, FUN.VALUE = 1) == 2)
        check_unperm <- all(vapply(unperm, FUN = length, FUN.VALUE = 1) == 2)
        if (check_perm == FALSE | check_unperm == FALSE) {
            stop("Each named parameter much have values for
                           each of the two components.
                           Number of components restricted to 2")
        }
        final <- get_adj_mnorm(perm = perm, unperm = unperm, ...)
    } else if (distr == "lst"){
        if (length(intersect(names(perm), c("df", "mu", "sigma"))) != 3){
            stop("Location-scale t-distribution requires mu, df and sigma arguments")
        }
        final <- list(mu = perm$mu, sigma = unperm$sigma, df = unperm$df)
    }
    return(final)
}

#' @title Function to perform the adjustment for the mixture normal distribution
#' @param perm (List). Parameter values of the distribution of scores
# ;    computed on permuted data
#' @param unperm (List). Parameter values of the distribution of scores
#'    computed on unpermuted data
#' @param fix_comp (Character). Which component to keep
#' @return A List of parameters for the adjusted mixture normal.
#' @keywords internal
#' @importFrom stats optim
get_adj_mnorm <- function(perm, unperm, verbose = FALSE, fix_comp = "none") {
    fix_comp <- match.arg(fix_comp, c("large", "small", "none"))
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
    opt_params <- list(method = "L-BFGS-B", mean = perm_mean,
                       sd = unperm_sd)
    if (fix_comp == "large"){
        c_idx <- which(perm$lambda == max(perm$lambda))
    } else if (fix_comp == "small"){
        c_idx <- which(perm$lambda == min(perm$lambda))
    }
    if (fix_comp %in% c("large", "small")){
        opt_params <- c(opt_params, list(par = 0.01,
                                         fn = obj_function_single,
                                         lower = .Machine$double.eps,
                                         upper = Inf,
                                         m1 = perm$mu[c_idx],
                                         m2 = perm$mu[-c_idx],
                                         l1 = perm$lambda[c_idx],
                                         l2 = perm$lambda[-c_idx],
                                         s1 = unperm$sigma[c_idx]))
    } else {
        opt_params <- c(opt_params, list(
            par = c(0.01, 0.01),
            fn = obj_function_double,
            lower = c(.Machine$double.eps, .Machine$double.eps),
            upper = c(Inf, Inf),
            mu = perm$mu, lambda = perm$lambda
        ))
    }
    opt <- do.call(optim, opt_params)
    if (sum(opt$par == 0.01) >= 1 | sum(opt$par == .Machine$double.eps) >= 1){
        warning("Could not find an optimized solution for the adjustment. At least one of the estimated values
                are equal to initial values or the machine limit")
    }
    # calculate the estimated sd
    if (fix_comp %in% c("large", "small")){
        est_sig <- rep(0, length(perm$lambda))
        est_sig[c_idx] <- unperm$sigma[c_idx]
        est_sig[-c_idx] <- opt$par
    } else {
        est_sig <- opt$par
    }
    estim_sd <- get_sd(
        sigma = est_sig, lambda = perm$lambda,
        mu = perm$mu, mean = perm_mean
    )
    if (verbose == TRUE) {
        print(paste("Total sd is", unperm_sd, "and estimated sd is", estim_sd))
    }
    param <- list(mu = perm$mu, sigma = est_sig, lambda = perm$lambda)
    return(param)
}

#' @keywords internal
obj_function_single <- function(s1, m1, l1, s2, l2, m2, mean, sd) {
    estimate_sd <- sqrt((s1 + m1 - mean) * l1 + (s2 + m2 - mean) * l2)
    # estimate the objective funciton which is an l2 norm
    obj <- sqrt(sum((estimate_sd - sd)^2))
    return(obj)
}

#' @keywords internal
obj_function_double <- function(sigma, mu, lambda, mean, sd) {
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
