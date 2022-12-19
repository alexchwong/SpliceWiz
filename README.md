# SpliceWiz

SpliceWiz is an R package for exploring differential alternative splicing events in splice-aware alignment BAM files.

# Table of Contents

1. [Documentation](#doco)
    1. [Bioconductor Release](#docorelease)
    2. [Bioconductor Devel](#docodevel)
2. [Installation](#inst)
    1. [Enabling OpenMP (MacOS)](#ompmac)
    2. [Installation for Bioconductor (current release)](#instrelease)
    3. [Installation for Bioconductor (devel)](#instdevel)
    4. [Installation for older Bioconductor versions](#instlegacy)
3. [Publication](#pub)

# Documentation <a name="doco"></a>

## Bioconductor Release (Bioc 3.16 / R 4.2) <a name="docorelease"></a>

[Bioconductor Landing Page (Release - 3.16)](https://bioconductor.org/packages/release/bioc/html/SpliceWiz.html)

[QuickStart Vignette (Release - 3.16)](https://bioconductor.org/packages/release/bioc/vignettes/SpliceWiz/inst/doc/SW_QuickStart.html)

[Reference Manual (Release - 3.16)](https://bioconductor.org/packages/release/bioc/manuals/SpliceWiz/man/SpliceWiz.pdf) 

## Bioconductor Devel (Bioc 3.17 / R 4.3) <a name="docodevel"></a>

[Bioconductor Landing Page (Devel - 3.17)](https://bioconductor.org/packages/devel/bioc/html/SpliceWiz.html)

[QuickStart Vignette (Devel - 3.17)](https://bioconductor.org/packages/devel/bioc/vignettes/SpliceWiz/inst/doc/SW_QuickStart.html)

[Reference Manual (Devel - 3.17)](https://bioconductor.org/packages/devel/bioc/manuals/SpliceWiz/man/SpliceWiz.pdf) 

# Installation (Release - Bioc 3.16 / R 4.2) <a name="inst"></a>

## Enabling OpenMP multi-threading (for MacOS users) <a name="ompmac"></a>

OpenMP is installed by default on Windows and Linux systems. For MacOS, OpenMP
is not officially supported. To install SpliceWiz with OpenMP support, first
install the `libomp` libraries via brew:

```
brew install libomp
```

## On R (version >= 4.2.0) using Bioconductor version 3.16 (current release) <a name="instrelease"></a>

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::valid()              # checks for out of date packages

BiocManager::install("SpliceWiz")
```

## On R (devel version >= 4.3) using Bioconductor version 3.17 (devel) <a name="instdevel"></a>

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::valid()              # checks for out of date packages

BiocManager::install("SpliceWiz")
```

## On older versions of Bioconductor (3.15 or earlier) <a name="instlegacy"></a>

```
library("devtools")
install_github("alexchwong/ompBAM")
install_github("alexchwong/SpliceWiz", dependencies=TRUE)
```

Note that prior to version 3.14 or earlier, you may need to retrieve
Mappability files from https://github.com/alexchwong/SpliceWizResources

# Publication <a name="pub"></a>

Check out our (latest) pre-print publication for SpliceWiz:

[SpliceWiz: easy, optimized, and accurate alternative splicing analysis in R](https://www.biorxiv.org/content/10.1101/2022.07.05.498887v1)