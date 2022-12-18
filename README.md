# SpliceWiz
SpliceWiz is an R package for exploring differential alternative splicing events in splice-aware alignment BAM files.

## Publication

Check out our (latest) pre-print publication for SpliceWiz:

[SpliceWiz: easy, optimized, and accurate alternative splicing analysis in R](https://www.biorxiv.org/content/10.1101/2022.07.05.498887v1)

## Documentation

[Bioconductor Landing Page](https://bioconductor.org/packages/devel/bioc/html/SpliceWiz.html)

[QuickStart Vignette - html](https://bioconductor.org/packages/devel/bioc/vignettes/SpliceWiz/inst/doc/SW_QuickStart.html)

[Reference Manual - PDF](https://bioconductor.org/packages/devel/bioc/manuals/SpliceWiz/man/SpliceWiz.pdf) 

## Installation 

### Enabling OpenMP multi-threading (for MacOS users)

OpenMP is installed by default on Windows and Linux systems. For MacOS, OpenMP
is not officially supported. To install SpliceWiz with OpenMP support, first
install the `libomp` libraries via brew:

```
brew install libomp
```

### On R (version >= 4.2.0) using Bioconductor version 3.16 (current release)

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::valid()              # checks for out of date packages

BiocManager::install("SpliceWiz")
```

### On older versions of Bioconductor (3.15 or earlier)

```
library("devtools")
install_github("alexchwong/ompBAM")
install_github("alexchwong/SpliceWiz", dependencies=TRUE)
```

Note that prior to version 3.14 or earlier, you may need to retrieve
Mappability files from https://github.com/alexchwong/SpliceWizResources

### Enabling htslib mode

We recently implemented an htslib version of SpliceWiz's `processBAM` and
`BAM2COV` functions, (namely `processBAM_hts()` and `BAM2COV_hts()`). To
install an htslib-dependent version of SpliceWiz, you will need to clone
the SpliceWiz repo and overwrite some files with htslib source files. 

First, we recommend you install SpliceWiz as normal from the above instructions,
if only to ensure you have the required R package dependencies. Then, from the
command line:

```
# Navigate to a directory of your choice, e.g.:
cd /path/to

git clone https://github.com/alexchwong/SpliceWiz
cd SpliceWiz
cp -r inst/htslib_version/* .
```

Then, from R:

```
# Roxygenize for documentation on processBAM_hts() and BAM2COV_hts()
devtools::document("/path/to/SpliceWiz")

# Install the package
install.packages(
    "/path/to/SpliceWiz", 
    repos = NULL, 
    type = "source"
)
```

