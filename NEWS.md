# CBEA 0.99.2  
NEW FEATURES  

- Added an option (`parametric`) to specify whether the null is estimated parametrically or
via pure permutation. To support this, an option (`nboot`) was also added to specify the number of permutations. A warning will be added if `parametric` is `FALSE` but `n_boot` is small (< 200)    
- Added option (`parallel_backend`) to specify the parallel backend of the loop using `BiocParallel`   
- Argument `control` now allow for a special slot titled `fix_comp` that can specify which component of the two-component mixture distribution to fix during the adjustment process.    

SIGNIFICANT USER-VISIBLE CHANGES  

- Combined `raw` argument with `output` argument for the CBEA function to specify
returning raw CBEA scores (without any distribution fitting and transformation).  
- New and revamped vignettes  
- Significant reduction in dependencies, including removing native support for `phyloseq`  

BUG FIXES  



# CBEA 0.99.1  
NEW FEATURES     

None  

SIGNIFICANT USER-VISIBLE CHANGES    

None  

BUG FIXES    

* Removed .Rproj file to conform with Bioconductor error    

# CBEA 0.99.0

NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.  
* Added a complete functionality to perform CBEA from scratch with bundled data set  

SIGNIFICANT USER-VISIBLE CHANGES

None

BUG FIXES

None 
