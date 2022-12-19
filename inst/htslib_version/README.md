# Enabling htslib mode

We recently implemented an htslib version of SpliceWiz's `processBAM` and
`BAM2COV` functions, (namely `processBAM_hts()` and `BAM2COV_hts()`). To
install an htslib-dependent version of SpliceWiz, you will need to clone
the SpliceWiz repo and overwrite some files with htslib source files. 

First, we recommend you install SpliceWiz as normal,
if only to ensure you have the required R package dependencies.

From R:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::valid()              # checks for out of date packages

BiocManager::install("SpliceWiz")
```

Then, from the command line, clone the repository and copy the contents of this
directory to overwrite corresponding files from the root directory:

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