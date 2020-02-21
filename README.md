# CuPhyR

This repository contains the R package CuPhyR (Custom Phyloseq R package). 
The package gives customized shortcuts for tidyverse and phyloseq for interacting with 
phyloseq class objects generated with the [Phyloseq package by Paul J. McMurdie](https://joey711.github.io/phyloseq/).
Since it is tested and developed with ASV data from the [dada2 package by Bejamin Callahan](https://benjjneb.github.io/dada2/index.html), 
the function names are referring to ASVs instead of 'taxa' (phyloseq terminology). 

The package is in early development. Use at your own risk, feedback welcome!

# Installation



### Packages required

* [phyloseq](https://joey711.github.io/phyloseq/)
* [dplyr](https://dplyr.tidyverse.org)
* [magrittr](https://magrittr.tidyverse.org/)

### Packages recommended
* [ape](https://cran.r-project.org/web/packages/ape/index.html)
* [data.table](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html)
