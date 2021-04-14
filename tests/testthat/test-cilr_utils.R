library(teaR)
library(testthat)


sim <- sapply(1:100, function(x) rpois(1000, lambda = sample(5:10, size = 1,replace = T)))
sim_pcount <- sim[sim == 0] <- 1
index <- 1:10
base_gmean <- function(sim, index){
    p <- dim(sim)[2]
    size <- length(index)
    up <- apply(sim[,index], 1, function(x) exp(mean(log(x))))
    down <- apply(sim[,-index],1, function(x) exp(mean(log(x))))
    scale <- sqrt(size * (p - size)/p)
    score <- scale * log(up/down)
    return(score)
}



get_score_wrapper <- function(sim, index){
    p <- ncol(sim)
    get_score(sim,index)
}

test_that("get_score returns the appropriate vector",{
    # first let's test for specific
    expect_equivalent(get_score_wrapper(sim_pcount, index), base_gmean(sim_pcount, index))

})


