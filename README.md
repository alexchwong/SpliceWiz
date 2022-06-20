# SpliceWiz
SpliceWiz is an R package for exploring differential alternative splicing events in splice-aware alignment BAM files.

## Installation 

### Enabling OpenMP multi-threading (for MacOS users)

OpenMP is installed by default on Windows and Linux systems. For MacOS, OpenMP
is not officially supported. To install SpliceWiz with OpenMP support, first
install the `libomp` libraries via brew:

```
brew install libomp
```

### On R (version >= 4.2.0) using Bioconductor version 3.1.6 (current devel)

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SpliceWiz")
```

### On R (version >= 4.2.0) using Bioconductor version 3.1.5

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")
BiocManager::valid()              # checks for out of date packages

library("devtools")
install_github("alexchwong/SpliceWiz", dependencies=TRUE, build_vignettes=TRUE)
```

### On R (version < 4.2.0, >= 4.1.0)

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::valid()              # checks for out of date packages

library("devtools")
install_github("alexchwong/ompBAM")
install_github("alexchwong/SpliceWiz", dependencies=TRUE, build_vignettes=TRUE)
```

### On R (version < 4.1.0)

Please note that vignettes do not work with Bioconductor 3.13 and below (as the BAM files are not accessible from old versions of ExperimentHub

```
library("devtools")
install_github("alexchwong/ompBAM")
install_github("alexchwong/NxtIRFdata")
install_github("alexchwong/SpliceWiz", dependencies=TRUE, build_vignettes=FALSE)
```

## Viewing the vignette

After installing SpliceWiz as detailed above, run the following to view the
SpliceWiz quick-start vignette

```
vignette("SW_QuickStart", package = "SpliceWiz")
```

