data("hmp_gingival")
library(mia)
library(SummarizedExperiment)
library(glue)
library(purrr)

physeq <- hmp_gingival$data
sets <- hmp_gingival$set

SummarizedExperiment::assay(physeq) <- SummarizedExperiment::assay(physeq, "16SrRNA") + 1

test_that("Testing parametric fits w/o adjustment", {
    expect_s3_class(cbea(obj = physeq, abund_values = "16SrRNA",
                         set = sets, output = "sig",
                         distr = "norm", adj = FALSE, n_perm = 1,
                         parametric = TRUE), class = "CBEAout")
    expect_s3_class(cbea(obj = physeq, abund_values = "16SrRNA",
                         set = sets, output = "cdf",
                         distr = "norm", adj = FALSE, n_perm = 1,
                         parametric = TRUE), class = "CBEAout")
    expect_s3_class(cbea(obj = physeq, abund_values = "16SrRNA",
                         set = sets, output = "sig",
                         distr = "mnorm", adj = FALSE, n_perm = 1,
                         parametric = TRUE), class = "CBEAout")
    expect_s3_class(cbea(obj = physeq, abund_values = "16SrRNA",
                         set = sets, output = "cdf",
                         distr = "mnorm", adj = FALSE, n_perm = 1,
                         parametric = TRUE), class = "CBEAout")
})

test_that("Testing non-parametric fits", {
    expect_s3_class(cbea(obj = physeq, abund_values = "16SrRNA",
                         set = sets, output = "sig", n_perm = 1,
                         parametric = FALSE), class = "CBEAout")

    expect_s3_class(cbea(obj = physeq, abund_values = "16SrRNA",
                         set = sets, output = "pval", n_perm = 1,
                         parametric = FALSE), class = "CBEAout")
})

test_that("Testing parametric fits w/ adjustment", {
    expect_s3_class(cbea(obj = physeq, abund_values = "16SrRNA",
                         set = sets, output = "sig",
                         distr = "norm", adj = TRUE, n_perm = 1,
                         parametric = TRUE), class = "CBEAout")
    expect_s3_class(cbea(obj = physeq, abund_values = "16SrRNA",
                         set = sets, output = "cdf",
                         distr = "norm", adj = TRUE, n_perm = 1,
                         parametric = TRUE), class = "CBEAout")
    expect_s3_class(cbea(obj = physeq, abund_values = "16SrRNA",
                         set = sets, output = "sig",
                         distr = "mnorm", adj = TRUE, n_perm = 1,
                         parametric = TRUE), class = "CBEAout")
    expect_s3_class(cbea(obj = physeq, abund_values = "16SrRNA",
                         set = sets, output = "cdf",
                         distr = "mnorm", adj = TRUE, n_perm = 1,
                         parametric = TRUE), class = "CBEAout")

})

test_that("Testing producing only raw outputs", {
    expect_s3_class(cbea(obj = physeq, abund_values = "16SrRNA", set = sets, 
         output = "raw", parametric = FALSE, adj = FALSE, n_perm = 1))
})


test_that("Expected errors and warnings combining different options", {
    # parametric is false but output is either cdf or zscore
    expect_error(cbea(obj = physeq, abund_values = "16SrRNA",
                      set = sets, output = "cdf", distr = "norm",
                      adj = TRUE, n_perm = 1,parametric = FALSE))
    expect_error(cbea(obj = physeq, abund_values = "16SrRNA",
                      set = sets, output = "zscore", distr = "norm",
                      adj = TRUE, n_perm = 1, parametric = FALSE))
    # parametric is true but no distr argument was performed
    expect_error(cbea(obj = physeq, abund_values = "16SrRNA",
                      set = sets, output = "zscore", adj = TRUE,
                      n_perm = 1, parametric = TRUE))
    # n_perm < 100 if parametric is FALSE
    expect_message(cbea(obj = physeq, abund_values = "16SrRNA",
                        set = sets, output = "sig", adj = TRUE,
                      n_perm = 1, parametric = FALSE))
})
