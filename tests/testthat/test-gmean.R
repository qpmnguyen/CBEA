library(magrittr)

df_pos <- matrix(rpois(100,3), 10, 10) %>% replace(which(. == 0), 1)
df_neg <- replace(df_pos, sample(seq_len(length(df_pos)), size = 3), -2)
df_zero <- replace(df_pos, sample(seq_len(length(df_pos)), size = 3), 0)

ref_gmean <- function(vec){
    return(exp(mean(log(unclass(vec)))))
}

ref_gmeanrow <- function(matrix){
    return(apply(matrix, 1, ref_gmean))
}

test_that("gmeanRow handles negatives", {
    expect_error(gmeanRow(df_neg), "X has 0 or negative values")
    expect_error(gmeanRow(df_zero), "X has 0 or negative values")
})

test_that("gmeanRow accepts and returns correct output formats", {
    expect_vector(gmeanRow(df_pos), ptype = numeric(), size = 10)
    expect_error(gmeanRow(c(1,2,3)))
    expect_error(gmeanRow(as.data.frame(df_pos)))
})

test_that("gmeanRow and gmean returns the correct values", {
    expect_equal(gmean(df_pos[1,]), ref_gmean(df_pos[1,]))
    expect_equal(gmeanRow(df_pos), ref_gmeanrow(df_pos))
})
