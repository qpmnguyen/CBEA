## code to prepare `DATASET` dataset goes here
library(HMP16SData)
library(tidyverse)
library(TreeSummarizedExperiment)
library(BiocSet)

# create sets
generate_sets <- function(physeq, metadata){
    names <- rownames(physeq)
    sets <- unique(metadata$Meth)
    set_list <- vector(length = length(sets), mode = "list")
    for (i in seq_along(sets)){
        g_names <- metadata %>% filter(Meth == sets[i]) %>% pull(Genera)
        idx <- which(as.vector(rowData(physeq)$GENUS) %in% g_names)
        set_list[[i]] <- names[idx]
    }
    names(set_list) <- sets
    set <- BiocSet::BiocSet(set_list)
    return(set)
}

supra <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Supragingival Plaque" & VISITNO == 1)
supra <- as(supra, "TreeSummarizedExperiment")

sub <- V35() %>% subset(select = HMP_BODY_SUBSITE == "Subgingival Plaque" & VISITNO == 1) %>%
    as(., "TreeSummarizedExperiment")

merged <- cbind(supra, sub)
filt_tax <- apply(assays(merged)$`16SrRNA`, 1, function(x) sum(x == 0)/length(x) <= 0.95)
filt_samples <- apply(assays(merged)$`16SrRNA`, 2, function(x) sum(x) >= 1000)
merged <- merged[filt_tax,filt_samples]

metadata <- read_tsv("https://raw.githubusercontent.com/mcalgaro93/sc2meta/master/data/genera_methabolism.tsv")

set <- generate_sets(merged, metadata)

hmp_gingival <- list(data = merged, set=set)
usethis::use_data(hmp_gingival, overwrite = TRUE)
