#TODO: Test input output type
#TODO: Test standard errors
#TODO: Come up with simple data sets and test for correctness
library(tidyverse)
library(phyloseq)
data(GlobalPatterns)

# preprocess data
GlobalPatterns <- GlobalPatterns %>%
    filter_taxa(function(x) (sum(x > 0)/length(x)) > 0.1, TRUE ) %>%
    filter_taxa(function(x) mean(x) >= 1e-5, TRUE) %>%
    transform_sample_counts(function(OTU) OTU + 1) %>%
    transform_sample_counts(function(OTU) OTU/sum(OTU))

# generate set
phylum_set <- const_set(tax_table(GlobalPatterns), rank = "Phylum")
physeq <- phyloseqSet(GlobalPatterns, taxon_set = phylum_set)
physeq <- trim_set(physeq, size >= 5)

table <- otu_table(physeq)
table <- t(as(table, "matrix"))
sets <- as(taxon_set(physeq), "list")

devtools::load_all()
cilr(ab_tab = table, set_list = sets, output = "cdf", distr = "norm", adj = TRUE, thresh = 0.05, raw = FALSE)

index <- which(colnames(table) %in% sets[[1]])

distr <- estimate_distr(raw_scores, distr = "norm")
ab_perm <- table[,sample(1:p, size = p, replace = FALSE)]

raw_scores <- get_score(table, index)
perm_scores <- get_score(ab_perm, index)

perm_distr <- estimate_distr(perm_scores, distr = "norm", init = NULL)
unperm_distr <- estimate_distr(perm_scores, distr = "norm", init = NULL)

combine_distr(perm_distr, unperm_distr, distr = "norm")
