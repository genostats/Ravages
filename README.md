# Ravages
Ravages is an R package to perform rare variant association tests (RVAT) and genetic simulations in the whole genome.

## Description
Ravages enables to perform rare variant association tests (burden and SKAT) with a multi-category, binary, and continuous phenotypes.  
These association tests can be applied genome-wide using the RAVA-FIRST approach based on the CADD regions. 
Ravages also enables to perform genetic simulations based on real data to mimic allele frequency spectrum and linkage disequilibrium pattern observed on these data.

## Installation
```
library(devtools)
devtools::install_github("genostats/Ravages")
```

## Usage
Please look at the <a href="doc/Ravages_vignette.pdf">vignette on RVAT</a> for an example of how to perform both burden and SKAT rare variant association tests.
Please look at the <a href="doc/Ravages_Simulations_vignette.pdf">vignette on simulations</a> for an example of how to perform genetic simulations with rare variants and how to compute power of RVAT using these simulations.
Please look at the <a href="doc/Ravages.pdf">Ravages user manual</a> for details about individual functions from Ravages. 
