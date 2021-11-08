## code to prepare `DATASET` dataset goes here
library(HMP16SData)
library(tidyverse)
library(phyloseq)
library(TreeSummarizedExperiment)
library(BiocSet)
library(mia)

# basic filtering
process_16S <- function(physeq){
    physeq <- physeq %>% filter_taxa(function(x) (sum(x == 0)/length(x)) < 0.9, TRUE)
    physeq <- prune_samples(sample_sums(physeq) >= 1000, physeq)
    return(physeq)
}

# create sets
generate_sets <- function(physeq, metadata){
    names <- taxa_names(physeq)
    sets <- unique(metadata$Meth)
    set_list <- vector(length = length(sets), mode = "list")
    for (i in seq_along(sets)){
        g_names <- metadata %>% filter(Meth == sets[i]) %>% pull(Genera)
        idx <- which(as.vector(tax_table(physeq)[,"GENUS"]) %in% g_names)
        set_list[[i]] <- taxa_names(physeq)[idx]
    }
    names(set_list) <- sets
    set <- BiocSet::BiocSet(set_list)
    return(set)
}

supra <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Supragingival Plaque" & VISITNO == 1) %>%
    as_phyloseq()

sub <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Subgingival Plaque" & VISITNO == 1) %>%
    as_phyloseq()

merged <- merge_phyloseq(supra, sub)
merged <- merged %>% subset_samples(!duplicated(RSID)) %>% process_16S()
merged_mia <- mia::makeTreeSummarizedExperimentFromPhyloseq(merged)
metadata <- read_tsv("https://raw.githubusercontent.com/mcalgaro93/sc2meta/master/data/genera_methabolism.tsv")

set <- generate_sets(merged, metadata)

hmp_gingival_physeq <- list(data = merged, set=set)
hmp_gingival_treesum <- list(data = merged_mia, set=set)

usethis::use_data(hmp_gingival_physeq, overwrite = TRUE)
usethis::use_data(hmp_gingival_treesum, overwrite = TRUE)
