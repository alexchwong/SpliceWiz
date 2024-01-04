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
    4. [Installing via GitHub](#instlegacy)
3. [Publication](#pub)

# Documentation <a name="doco"></a>

## Bioconductor Release (Bioc 3.18 / R 4.3) <a name="docorelease"></a>

[Bioconductor Landing Page (Release)](https://bioconductor.org/packages/release/bioc/html/SpliceWiz.html)

[QuickStart Vignette (Release)](https://bioconductor.org/packages/release/bioc/vignettes/SpliceWiz/inst/doc/SW_QuickStart.html)

[Reference Manual (Release)](https://bioconductor.org/packages/release/bioc/manuals/SpliceWiz/man/SpliceWiz.pdf) 

## Bioconductor Devel (future Bioc 3.19 / R 4.3) <a name="docodevel"></a>

[Bioconductor Landing Page (Devel)](https://bioconductor.org/packages/devel/bioc/html/SpliceWiz.html)

[QuickStart Vignette (Devel)](https://bioconductor.org/packages/devel/bioc/vignettes/SpliceWiz/inst/doc/SW_QuickStart.html)

[Reference Manual (Devel)](https://bioconductor.org/packages/devel/bioc/manuals/SpliceWiz/man/SpliceWiz.pdf) 

# Installation (Release - Bioc 3.18 / R 4.3) <a name="inst"></a>

## Enabling OpenMP multi-threading (for MacOS users) <a name="ompmac"></a>

OpenMP is installed by default on Windows and Linux systems. For MacOS, OpenMP
is not officially supported. To install SpliceWiz with OpenMP support, first
install the `libomp` libraries via brew:

```
brew install libomp
```

## On R (version >= 4.3) using Bioconductor version 3.18 (current release) <a name="instrelease"></a>

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::valid()              # checks for out of date packages

BiocManager::install("SpliceWiz")
```

## On R-devel using Bioconductor devel (future 3.19) <a name="instdevel"></a>

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::valid()              # checks for out of date packages

BiocManager::install("SpliceWiz")
```

## Installing via GitHub <a name="instlegacy"></a>

Reasons for installing via GitHub:
* Using the latest development version on current release of Bioconductor -
means you don't have to install Bioconductor devel
* You are using Bioconductor version 3.16 or earlier

```
library("devtools")
install_github("alexchwong/ompBAM")

# To install the latest devel version, install from the "main" branch
install_github("alexchwong/SpliceWiz", "main", dependencies=TRUE)
```

Note that prior to Bioconductor versions 3.14 or earlier, you may need to retrieve
Mappability files from https://github.com/alexchwong/SpliceWizResources



# Publication <a name="pub"></a>

SpliceWiz is now published!

[SpliceWiz: interactive analysis and visualization of alternative splicing in R](https://academic.oup.com/bib/article/25/1/bbad468/7502685)
