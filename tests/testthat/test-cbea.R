requireNamespace("tidyverse", quietly = TRUE)

# reference gmean and gmean rows
ref_gmean <- function(vec){
    return(exp(mean(log(unclass(vec)))))
}

ref_gmeanrow <- function(matrix){
    return(apply(matrix, 1, ref_gmean))
}

ref_ilr <- function(matrix, index){
    p <- ncol(matrix)
    size <- length(index)
    scale <- sqrt(size * (p - size)/p)
    num <- ref_gmeanrow(matrix[,index])
    denom <- ref_gmeanrow(matrix[,-index])
    return(scale * log(num/denom))
}

# generate a fake data set
df_pos <- matrix(rpois(1000,3), 100, 10) %>% replace(which(. == 0), 1)
df_neg <- replace(df_pos, sample(seq_len(length(df_pos)), size = 3), -2)
df_zero <- replace(df_pos, sample(seq_len(length(df_pos)), size = 3), 0)
index <- c(1:3)

###### Testing for get_score function ###################################
# Testing correctness
test_that("Testing correctness for get_score",{
    expect_equal(get_score(df_pos, index), ref_ilr(df_pos, index))
})

test_that("get_score should return error if there are zeroes or negatives", {
    expect_error(get_score(df_zero, index))
    expect_error(get_score(df_neg, index))
})

# Testing input output
test_that("Testing get_score returns a vector",{
    expect_vector(get_score(df_pos, index))
})

test_that("Testing if warning is displayed if data.frame is included", {
    expect_message(get_score(as.data.frame(df_pos), index))
})

##### Testing for the estimate_distr function ####
# Testing for input output first
scores <- get_score(df_pos, index)
test_that("Making sure estimate_distr returns a list",{
    expect_type(estimate_distr(scores, distr = "norm"), "list")
    expect_type(estimate_distr(scores, distr = "mnorm"), "list")
})

test_that("Making sure the names test_that are using is correct", {
    expect_named(estimate_distr(scores, distr = "norm"), c("mean", "sd"))
    expect_named(estimate_distr(scores, distr = "mnorm"), c("mu", "sigma", "lambda"))
})

test_that("Return error if have too much NAs", {
    expect_error(estimate_distr(rep(NA,100)))
    expect_error(estimate_distr(c(rnorm(10), rep(NA, 90))))
})

##### Testing for the main CBEA function #####




