# Methods for cbea_output

#' @title Print dispatch for CBEAout objects 
#' @importFrom glue glue 
#' @export
print.CBEAout <- function(x, ...){
    cat(
        glue("CBEA output of class 'CBEAout' with {nsamp} samples and {nset} sets", 
             nsamp = nrow(x$R), nset = ncol(x$R)), "\n", 
        glue("Fit type: {par} with {distr}", par = ifelse(paramateric == TRUE,)), "\n", 
        glue("Number of permutations: 200"), "\n",
        glue("Output type: cdf"), "\n",
        glue("Parameters:")
    )
}