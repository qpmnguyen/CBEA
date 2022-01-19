data("hmp_gingival")
library(mia)
library(SummarizedExperiment)
library(glue)
library(purrr)

physeq <- hmp_gingival$data
sets <- hmp_gingival$set

SummarizedExperiment::assay(physeq) <- SummarizedExperiment::assay(physeq, "16SrRNA") + 1

test_that("Basic outputs", {
    settings <- purrr::cross_df(list(
        distr = c("mnorm", "norm"),
        adj = c(TRUE, FALSE),
        output = c("sig", "pval", "cdf", "zscore"),
        parametric = c(TRUE, FALSE)
    ))
    settings <- settings %>% dplyr::anti_join(
        purrr::cross_df(list(
            output = c("zscore", "cdf"),
            parametric = FALSE
        ))
    )
    purrr::pmap(settings, function(distr, adj, output, parametric){
        print(glue::glue("Currently at {distr}, {adj}, {output}, {parametric}"))
        expect_s3_class(cbea(obj = physeq, abund_values = "16SrRNA",
                             set = sets, output = output,
                           distr = distr, adj = adj, n_perm = 10,
                           parametric = parametric), class = "tbl_df")
        return(NULL)
    })
})

test_that("Expected errors and warnings combining different options", {
    # parametric is false but output is either cdf or zscore
    expect_error(cbea(obj = physeq, set = sets, output = "cdf", distr = "norm",
                      adj = TRUE, n_perm = 1,parametric = FALSE))
    expect_error(cbea(obj = physeq, set = sets, output = "zscore", distr = "norm",
                      adj = TRUE, n_perm = 1, parametric = FALSE))
    # parametric is true but no adj argument was performed
    expect_error(cbea(obj = physeq, set = sets, output = "zscore", distr = "norm",
                        n_perm = 1, parametric = TRUE))
    # parametric is true but no distr argument was performed
})
