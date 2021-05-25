library(phyloseq)

data(enterotype)

# Assign random GENUS rank to enterotype
tax_table(enterotype) |> class()
