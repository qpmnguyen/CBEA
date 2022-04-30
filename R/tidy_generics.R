#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance


#' @title Tidy a CBEAout object
#' @description This function takes in a CBEA type object and collects all values across all
#'     sets and samples that were evaluated.
#' @param x A `CBEAout` object.
#' @param ... Unused, included for generic consistency only.
#' @return A tidy [tibble::tibble()] summarizing scores per sample per set.
#'
#' @examples
#' # load the data
#' data(hmp_gingival)
#' mod <- cbea(hmp_gingival$data, hmp_gingival$set, abund_values = "16SrRNA",
#'     output = "sig", distr = "norm", adj = FALSE, n_perm = 5, parametric = TRUE)
#' tidy(mod)
#' @importFrom tibble as_tibble add_column
#' @export
tidy.CBEAout <- function(x, ...){
    output <- as_tibble(x$R)
    colnames(output) <- gsub(" ", "_", colnames(output))
    output <- add_column(output, sample_ids = x$call_param$sample_ids, .before = 1)
    return(output)
}

#' @title Glance at \code{CBEAout} object
#' @description This function cleans up all diagnostics of the \code{cbea} method
#'     (from the \code{CBEAout} object) into a nice [tibble::tibble()]
#' @param x An object of type \code{CBEAout}
#' @param statistic What type of diagnostic to return. Users can choose to return
#'     \code{fit_diagnostic} which returns goodness of fit statistics for the
#'     different fitted distributions (e.g. log likelihoods) while \code{fit_comparison}
#'     returns comparisons across different distributions and raw values (and data) across the
#'     4 l-moments.
#' @param ... Unused, kept for consistency with generics
#' @return A [tibble::tibble()] summarizing diagnostic fits per set (as row)
#' @importFrom tibble tibble as_tibble rownames_to_column
#' @importFrom tidyr pivot_longer unnest nest
#' @importFrom dplyr mutate select rename left_join everything
#' @importFrom magrittr %>%
#' @examples
#' # load the data
#' data(hmp_gingival)
#' mod <- cbea(hmp_gingival$data, hmp_gingival$set, abund_values = "16SrRNA",
#'     output = "sig", distr = "norm", adj = FALSE, n_perm = 5, parametric = TRUE)
#' glance(mod, "fit_diagnostic")
#' @export
glance.CBEAout <- function(x, statistic, ...){
    # Assign in R CMD CHECK
    fit_comparison <- NULL
    value <- NULL

    statistic <- match.arg(statistic, c("fit_diagnostic", "fit_comparison"))
    if (x$call_param$output == "raw"){
        message("There are no diagnostics if raw scores are returned")
        output <- NA
    } else {
        names(x$diagnostic) <- gsub(" ", "_", names(x$diagnostic))
        names(x$parameters) <- gsub(" ", "_", names(x$parameters))
        names(x$fit_comparison) <- gsub(" ", "_", names(x$fit_comparison))
        
        parameters <- tibble(set_ids = names(x$parameters),
                             final_param = x$parameters)
        if (statistic == "fit_diagnostic"){
            diagnostic <- as_tibble(x$diagnostic) %>%
                pivot_longer(everything()) %>%
                mutate(type = names(value)) %>%
                mutate(value = lapply(value, as_tibble)) %>%
                unnest(value) %>%
                rename("set_ids" = "name")
            output <- left_join(parameters, diagnostic, by = "set_ids")
        } else if (statistic == "fit_comparison"){
            fit_parameters <- tibble(set_ids = names(x$fit_comparison),
                                     fit_comparison = lapply(x$fit_comparison, function(x) {
                                         as.data.frame(x) %>%
                                             rownames_to_column(var = "distr") %>%
                                             as_tibble()
                                     })) %>% unnest(fit_comparison)
            output <- left_join(parameters, fit_parameters, by = "set_ids")
        }
    }
    return(output)
}






