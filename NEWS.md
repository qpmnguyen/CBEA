# CBEA 1.0.1   
NEW FEATURES  
SIGNIFICANT USER-VISIBLE CHANGES  
BUG FIXES  

- Fixed an issue where `perm_scores` was not found when `output` are raw scores.  
- Fixed an issue where a permutation warning was returned if `parametric` is FALSE but `output` are raw scores (hence no need to specify the number of permutations)  



# CBEA 1.0.0  
NEW FEATURES  

- Bioconductor release version 3.15  

SIGNIFICANT USER-VISIBLE CHANGES  

None  

BUG FIXES  

None

# CBEA 0.99.3  
NEW FEATURES 

- Created an output type object (`CBEAout`). This is an S3 type object that is essentially a list that incorporates the final score matrix as well as other diagnostic details.  

SIGNIFICANT USER-VISIBLE CHANGES

- Due to the new feature above, now instead of getting a tibble, users would have to extract the scores out either using the provided function or use a custom approach based on the list-of-list format of the `CBEAout` objects.    
- Implemented `tidy` and `glance` methods to deal with `CBEAout` objects  

BUG FIXES

None

# CBEA 0.99.2  
NEW FEATURES  

- Added an option (`parametric`) to specify whether the null is estimated parametrically or
via pure permutation. To support this, an option (`n_perm`) was also added to specify the number of permutations. A warning will be added if `parametric` is `FALSE` but `n_boot` is small (< 100)    
- Added option (`parallel_backend`) to specify the parallel backend of the loop using `BiocParallel`  
- Argument `control` now allow for a special slot titled `fix_comp` that can specify which component of the two-component mixture distribution to fix during the adjustment process.    

SIGNIFICANT USER-VISIBLE CHANGES  

- Combined `raw` argument with `output` argument for the CBEA function to specify
returning raw CBEA scores (without any distribution fitting and transformation).  
- New and revamped vignettes  
- Significant reduction in dependencies, including removing native support for `phyloseq`  

BUG FIXES  

None

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
