# TODO: Test input output type
# TODO: Test standard errors
# TODO: Come up with simple data sets and test for correctness
data(hmp_gingival)
requireNamespace("purrr", quietly = TRUE)
library(magrittr)
library(phyloseq)
library(BiocSet)
library(TreeSummarizedExperiment)

seq <- hmp_gingival$data
set <- hmp_gingival$set
# for some reason makeTreeSummarizedExperimentFromPhyloseq(seq) doesn't work
seq_ts <- TreeSummarizedExperiment::TreeSummarizedExperiment(
  assays = otu_table(seq),
  rowData = tax_table(seq),
  colData = sample_data(seq)
)

seq_matrix <- as(otu_table(seq), "matrix")
seq_df <- as.data.frame(seq_matrix)
objects <- list(seq, seq_ts, seq_df, seq_matrix)
names(objects) <- c("phyloseq", "ts", "df", "matrix")


### test unify sets ------------------------------- ####
test_that("unify_sets remove redundant elements from a set that is not part of the data", {
  set_add <- dplyr::union(set, BiocSet::BiocSet(met1 = letters[c(1:3)], met2 = letters[c(5:7)]))
  purrr::imap(objects, ~ {
    if (.y %in% c("phyloseq", "ts")) {
      set_new <- unify_sets(obj = .x, set = set_add)
    } else {
      set_new <- unify_sets(obj = .x, set = set_add, taxa_are_rows = TRUE)
    }
    print(paste("Currently at", .y))
    expect_s4_class(set_new, "BiocSet")
    expect_identical(object = set_new, expected = set)
  })
})
