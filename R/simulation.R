# Simulations to try to fit cILR functions
# Will be superceded by microbesim package once it is finished.

#' @title Shorthand to create parameters
#' @param params List of parameters
#' @keywords internal
#' @importFrom dplyr mutate group_by
#' @importFrom tidyr nest
#' @importFrom purrr cross_df
#' @importFrom magrittr %>%
create_parameters <- function(params){
    par <- purrr::cross_df(params)
    par <- par %>% dplyr::mutate(id = seq(nrow(par))) %>% dplyr::group_by(id) %>% tidyr::nest()
    par <- par %>% dplyr::transmute(param = data)
    return(par)
}


#' @title Simulating microbiome relative abundance data
#'     according to zero inflated negative binomial distribution
#' @description Function to simulate microbiome relative abundance data. Each taxa is simulated
#'     using a gaussian copula with a zero-inflated negative binomial marginals. Default parameters
#'     taken from fitting zero-inflated data to stool samples from the Human Microbiome Project.
#'     Note: n_inflate * n_sets must be equal to n_taxa.
#' @param n_samp Number of samples
#' @param spar Additive sparsity
#' @param b_rho Baseline inter-taxa correlation
#' @param s_rho Correlation within the set
#' @param n_tax Number of taxa
#' @param n_inflate Number of taxa with inflated counts
#' @param n_sets Number of sets to simulate
#' @param prop_set_inflate Number of proportions of sets that are inflated
#' @param prop_inflate The number of differentially abundant taxa within a set that is inflated
#' @param samp_prop The proportion of samples have an inflated taxa
#' @param eff_size Effect size that is a multiplier to the base mu of a sample
#' @param method Method of simulation. Can be "compensation", "non-compensation" or "normal"
#' @param parameters Path to rds file containing estimated values
#' @param vary_params Whether we stochastically draw parameter values for each feature
#' @importFrom MASS mvrnorm
#' @importFrom purrr map_dfc
#' @importFrom stats qnbinom rbinom
#' @export
simulate_data <- function(n_samp=1000, n_tax = 300, spar = 0.1, s_rho = 0.3, eff_size = 3,
                          b_rho = 0, n_inflate = 50, n_sets = 1,
                          prop_set_inflate = 1, prop_inflate = 1, samp_prop = 0.5,
                          method = "compensation", vary_params = TRUE, parameters = NULL){

    # generate the the diagnonal matrix
    sigma <- diag(n_tax)
    sigma[sigma == 0] <- b_rho
    print(dim(sigma))
    set_sigma <- sigma[1:n_inflate, 1:n_inflate]
    set_sigma[set_sigma != 1] <- s_rho
    sigma[1:n_inflate,1:n_inflate] <- set_sigma

    # First create mvrnorm variables with correlation set by sigma
    margins <- pnorm(MASS::mvrnorm(n = n_samp, mu = rep(0, n_tax), Sigma = sigma))
    # Second, set marginals
    # default marginals for negative binomial is from size of 0.595 and mu of 6.646 from HMP data
    true_size <- round(n_inflate * prop_inflate, 0)

    # Create parameters based on situation
    if (vary_params == T){
        if (is.null(parameters)){
            message("Randomly sample means from 1 to 10 and sizes from 1 to 5")
            means <- runif(n_tax, 3,5)
            sizes <- runif(n_tax, 1,3)
        } else {
            estimated <- readRDS(file = parameters)
            means <- sample(estimated$mean, size = n_tax, replace = T)
            sizes <- sample(estimated$size, size = n_tax, replace = T)
        }
    } else {
        message("Setting mean to be constant at 3.06 and size at 1.67 estimated from HMP data...")
        means <- rep(3.05,n_tax)
        sizes <- rep(1.67, n_tax)
    }


    # first n_samp * samp_prop samples will always be inflated
    inf_size <- round(n_samp * samp_prop,0)

    # Number of inf_taxa
    inf_tax <- round(n_inflate * n_sets * prop_set_inflate, 0)

    # inflating samples
    suppressMessages(
        inf_samples <- purrr::map_dfc(seq(n_tax), .f = function(.x){
            if(method == "compensation"){
                prop <- 1/(eff_size + 1)
                # randomly select first prop of those selected to be upregulated
                idx <- sample(seq(inf_tax), size = round(prop * inf_tax, 0), replace = F)
                remainder <- seq(inf_tax)[-idx]
                a <- sum(means[idx])
                b <- sum(means[remainder])
                # randomly select second prop of those selected to be downregulated
                if (.x %in% idx){
                    result <- stats::qnbinom(p = margins[seq(inf_size),.x], size = sizes[.x],
                                      mu = means[.x]*eff_size)
                } else if (.x %in% remainder){
                    result <- stats::qnbinom(p = margins[seq(inf_size),.x], size = sizes[.x],
                                      mu = means[.x]*((a/b)*(1-eff_size) + 1))
                } else {
                    result <- stats::qnbinom(p = margins[seq(inf_size), .x], size = sizes[.x],
                                      mu = means[.x])
                }
            } else if (method == "no_compensation"){
                # randomly select first prop of those selected to be upregulated
                idx <- sample(seq(inf_tax), size = inf_tax/2, replace = F)
                remainder <- seq(inf_tax)[-idx]
                # randomly select second prop of those selected to be downregulated
                if (.x %in% idx){
                    result <- stats::qnbinom(p = margins[seq(inf_size),.x], size = sizes[.x],
                                      mu = means[.x]*eff_size)
                } else if (.x %in% remainder){
                    result <- stats::qnbinom(p = margins[seq(inf_size),.x], size = sizes[.x],
                                      mu = means[.x]/eff_size)
                } else {
                    result <- stats::qnbinom(p = margins[seq(inf_size), .x], size = sizes[.x],
                                      mu = means[.x])
                }

            } else if (method == "normal"){
                if (.x %in% seq(inf_tax)){
                    result <- stats::qnbinom(p = margins[seq(inf_size),.x], size = sizes[.x],
                                      mu = means[.x]*eff_size)
                } else {
                    result <- stats::qnbinom(p = margins[seq(inf_size), .x], size = sizes[.x],
                                      mu = means[.x])
                }
            }
            return(result)
        })
    )
    # not inflated samples
    suppressMessages(
        notinf_samples <- purrr::map_dfc(seq(n_tax), .f = function(.x){
            result <- qnbinom(p = margins[-seq(inf_size),.x], size = sizes[.x],
                              mu = means[.x])
            return(result)
        })
    )
    message("Completed loop!")
    if (inf_size == n_samp){
        abundance <- inf_samples
    } else {
        abundance <- rbind(inf_samples, notinf_samples)
    }
    abundance <- as.matrix(abundance)
    if (!is.null(spar)){
        zeroes <- stats::rbinom(length(abundance), size = 1, prob = 1 - spar)
        abundance <- abundance * zeroes
    }
    label <- c(rep(1, inf_size), rep(0, n_samp - inf_size))

    colnames(abundance) <- paste0("Tax", seq(n_tax))
    rownames(abundance) <- paste0("Samp", seq(n_samp))

    if (n_sets > 1){
        A <- diag(n_sets)
        vec <- as.matrix(rep(1, n_inflate))
        A <- kronecker(A,vec)

    } else {
        #TODO if n_sets * n_inflate != n_tax, then there are issues.
        message("Only one set!")
        A <- rep(0,n_tax)
        A[1:n_inflate] <- 1
        A <- as.matrix(A)
    }
    colnames(A) <- paste0("Set",1:n_sets,"ss")
    rownames(A) <- colnames(abundance)
    sets <- apply(A, 2, function(x) {
        rownames(A)[which(x == 1)]
    }, simplify = FALSE)
    names(sets) <- colnames(A)
    sets_inf <- rep(0, n_sets)
    sets_inf[seq(round(n_sets * prop_set_inflate,0))] <- 1
    output <- list(X = as.data.frame(abundance), sets = sets, label = label, sets_inf = sets_inf)
    return(output)
}
