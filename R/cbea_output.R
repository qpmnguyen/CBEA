# CBEAout

#' @title Creating an output object of type CBEAout  
#' @description This function takes a list of lists from each object and 
#'  turns it into a CBEAout type object
#' @param out A list containing scores for each set
#' @param sample_ids A vector of sample ids from the original data frame 
#' @param call A list containing all important arguments for printing 
#' @importFrom tibble add_column as_tibble rownames_to_column
#' @importFrom tidyr unnest
#' @importFrom dplyr bind_rows arrange
#' @importFrom magrittr %>%
new_CBEAout <- function(out, sample_ids, call){
    req_values <- c("parametric", "distr", "output", "n_perm", "adj")
    if (length(intersect(names(call), req_values)) != length(req_values)){
        stop("Require the 'call' argument to have all important elements of a function call")
    }
    # first process the score_matrix
    score_matrix <- data.frame(sapply(out, function(x) x$score))
    score_matrix <- as_tibble(add_column(score_matrix, sample_id = sample_ids, .before = 1))
    
    if (output == "raw"){
        parameters <- NULL
        diagnostic <- NULL
        fit_comparison <- NULL
    }  else {
        # second, get parameters
        parameters <- lapply(out, function(x) unlist(x$diagnostic$parameters, use.names = TRUE))
        
        # third, get diagnostic values
        diagnostic <- lapply(seq_along(out), function(x) {
            as_tibble(out[[x]]$diagnostic[-c(which(names(out[[x]]$diagnostic) %in% c("distr", "parameters")))]) %>% 
                unnest(everything()) %>% add_column(set = names(out)[[x]], .before = 1) %>% 
                add_column(eval = c("loglik", "D"), .before = 1)
            
        })
        diagnostic <- do.call(bind_rows, diagnostic) %>% arrange(eval)
        
        # fourth, get fit comparisons 
        fit_comparison <- lapply(seq_along(out), function(x){
            as.data.frame(out[[x]]$lmoment) %>% rownames_to_column(var = "type") %>% 
                add_column(set = names(out)[[x]], .before = 1) %>% as_tibble()
        })
        fit_comparison <- do.call(bind_rows, fit_comparison)
        
    }
    object <- list(
        R = score_matrix, 
        parameters = parameters, 
        diagnostic = diagnostic, 
        fit_comparison = fit_comparison
    )
    class(object) <- "CBEAout"
    return(object)
}


validate_CBEAout <- function(obj){
    
}