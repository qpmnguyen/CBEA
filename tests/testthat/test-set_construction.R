#TODO: Test input output type
#TODO: Test standard errors
#TODO: Come up with simple data sets and test for correctness
library(phyloseq)
data(GlobalPatterns)
requireNamespace("tidyverse", quietly = TRUE)
# preprocess data
GlobalPatterns <- GlobalPatterns %>%
    filter_taxa(function(x) (sum(x > 0)/length(x)) > 0.1, TRUE ) %>%
    filter_taxa(function(x) mean(x) >= 1e-5, TRUE) %>%
    transform_sample_counts(function(OTU) OTU + 1)
GlobalPatterns %>%
    otu_table() %>% as(.,"matrix") %>% t() %>% .[,1:30]
# generate set
# phylum_set <- const_set(tax_table(GlobalPatterns), rank = "Phylum")
# physeq <- phyloseqSet(GlobalPatterns, taxon_set = phylum_set)
# physeq <- trim_set(physeq, size >= 5)
