# CuPhyR

This repository contains the R package CuPhyR (Custom Phyloseq R package). 
The package gives customized shortcuts for tidyverse and phyloseq for interacting with 
phyloseq class objects generated with the [Phyloseq package by Paul J. McMurdie](https://joey711.github.io/phyloseq/).
Since it is tested and developed with ASV data from the [dada2 package by Bejamin Callahan](https://benjjneb.github.io/dada2/index.html), 
the function names are referring to ASVs instead of 'taxa' (phyloseq terminology). 

The package is in early development. Use at your own risk, feedback welcome!

# Installation
Install the package from the GitHub repository using devtools:  
```
library(devtools)
devtools::install_github("simeross/cuphyr")
```


# Info
## Functions contained in this package
#### *abundant_tax_physeq* - Return top n most abundand taxonomic terms at a given level
This function extracts the top n most abundant taxonomic terms at a given level.
It uses dplyr to calculate the most abundant elements.
Output_format will take "ps" or "tops" and either return the phyloseq object or a character vector of the n most common terms.
Ignores NA's by default but can show them (see parameters).

#### *list_subset_ASVs* - Returns a vector of ASVs belonging to groups at a taxonomic level
A function to list ASVs belonging to a set of taxonomic terms

#### *make_ranked_sums* - Rank ASV abundance sums for a given phyloseq object
A function that ranks all samples in a phyloseq object by the sum of reads assigned to all ASVs 
in that sample (Abundance sum/Total ASV counts)

#### *root_tree_in_outgroup* - Root phylogenetic tree of a phyloseq object
This funtion defines the leaf with the longest path as the root of the phylogenetic tree.
This makes results reproducible by avoiding the behaviour of some functions that would otherwise pick a
random leaf as the root of an unrooted phylogenetic tree.
Based on answers on [joey711/phyloseq/issues/597](https://github.com/joey711/phyloseq/issues/597). 
The function requires the packages 'ape' and 'data.table' to be installed.

#### *summarise_physeq* - Summarize a phyloseq object
Summarize the total number of ASVs, number of ASVs in a subset and percentage
of reads this ASV subset contains (of total) of a given phyloseq object.

# Dependencies
These packages have to be installed for CuPhyR to work properly. Except for phyloseq, 
all packages are distributed through CRAN and should be installed automatically.
### Packages required
* [phyloseq](https://joey711.github.io/phyloseq/) - Check the page for installation instructions
* [dplyr](https://dplyr.tidyverse.org)
* [magrittr](https://magrittr.tidyverse.org/)

### Packages recommended
* [ape](https://cran.r-project.org/web/packages/ape/index.html)
* [data.table](https://cran.r-project.org/web/packages/data.table/vignettes/datatable-intro.html)
