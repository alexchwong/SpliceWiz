#' SpliceWiz: a command line interface for NxtIRF - IRFinder-based
#' differential Alternative Splicing and Intron Retention analysis
#'
#' NxtIRF is a computationally efficient and user friendly workflow that
#' analyses aligned short-read RNA sequencing for differential
#' intron retention and alternative splicing. It utilises an improved
#' IRFinder-based OpenMP/C++ algorithm. A streamlined downstream analysis
#' pipeline allows for GLM-based differential IR and splicing analysis, suited
#' for large datasets of up to hundreds of samples. Additionally NxtIRF provides
#' a novel visualisation of per-nucleotide mean and variations of
#' alignment coverage across splice and IR events, grouped by user-defined
#' experimental conditions.
#'
#' @details
#' [IRFinder](https://doi.org/10.1186/s13059-017-1184-4) is a well-established
#' bioinformatic tool that measures intron
#' retention (IR) in annotated and novel retained introns in short-read RNA
#' sequencing samples. It is a computationally-efficient algorithm that measures
#' alignment coverage across introns, accounting for regions of low-mappable
#' intronic regions. Unlike other algorithms that measure exon-intron spanning
#' reads, IRFinder considers the alignment coverage across the whole intron,
#' allowing it to distinguish between full-length and partial IR. This
#' distinction is important as partial IR is often confounded with
#' novel alternate splice site usage, alternate transcription start site and
#' intronic polyadenylation events.
#'
#' NxtIRF is a R/Bioconductor package that provides a user-friendly workflow
#' using the IRFinder algorithm to perform both IR and
#' alternative splicing analysis in large datasets. By incorporating the core
#' C++ based IRFinder algorithm using Rcpp, NxtIRF is multi-platform and further
#' improves computational efficiency using OpenMP-based multi-threading.
#' Besides analysing IR, NxtIRF analyses other forms of
#' alternative splicing events that depend on alternate splice site selection,
#' including skipped exons, mutually exclusive exons, alternate 5'- and 3'-
#' splice sites, alternate first exons and alternate last exons.
#'
#' Downstream, NxtIRF provides functions to collate individual NxtIRF/IRFinder
#' outputs of multiple samples in an experiment / dataset, and assembles
#' these into a specialised \linkS4class{NxtSE} object that
#' inherits the SummarizedExperiment class. Users can easily define experimental
#' conditions, perform differential analysis and filter out lowly-expressed
#' splice events.
#'
#' Finally, NxtIRF provides visualisation tools to illustrate alternative
#' splicing using coverage plots, including a novel method to normalise RNA-seq
#' coverage grouped by experimental condition. This approach accounts for
#' variations introduced by sequenced library size and gene expression. NxtIRF
#' efficiently computes and visualises means and variations in per-nucleotide
#' coverage depth across alternate exons in genomic loci.
#'
#' SpliceWiz is the command line interface for R/Bioconductor. NxtIRF (coming
#' soon) will feature an interactive graphical user interface with additional
#' functions.
#'
#' **Features include:**
#' * Reference generation from user-supplied local and web resources, as well as
#'   connectivity to the AnnotationHub repository for Ensembl-based genomes
#'   and gene annotations;
#' * OpenMP and BiocParallel-based multi-threaded support to process short-read
#'   BAM files using the IRFinder algorithm written in native C++;
#' * Stores alignment coverage using the *COV* format, which is a binary
#'   compressed and indexed format for rapid recall of RNA-seq coverage. In
#'   contrast to the *BigWig* format, *COV* files store coverage of unstranded
#'   as well as stranded alignment coverage, and is much more space-efficient,
#'   allowing for better portability;
#' * Memory-efficient collation of hundreds of samples using on-disk memory
#'   approaches and H5-based assay storage;
#' * Streamlined user-friendly functions to construct multi-factor complex
#'   experimental designs, and perform differential IR and alternative splicing
#'   analysis using well-established statistical methods including limma and
#'   DESeq2;
#' * Advanced RNA-seq coverage visualisation, including the ability to
#'   combine RNA-seq coverage of multiple samples using advanced library
#'   normalisation methods across samples grouped by conditions;
#'
#' The main functions are:
#'
#' * [BuildReference] - Prepares genome and gene annotation
#'   references from FASTA and GTF files, and synthesises the NxtIRF reference
#'   for the IRFinder engine and NxtIRF-based downstream analysis.
#' * [STAR-methods] - (Optional) Provides wrapper functions to build the STAR
#'   genome reference and alignment of short-read FASTQ raw sequencing files.
#'   This functionality is only available on systems with STAR installed.
#' * [IRFinder] - OpenMP/C++ based IRFinder algorithm to analyse
#'   single or multiple BAM files using the NxtIRF/IRFinder reference.
#' * [CollateData] - Collates an experiment based on multiple IRFinder outputs
#'   for individual samples, into one unified H5-based data structure.
#' * [MakeSE] - Constructs a \linkS4class{NxtSE} (H5-based
#'   SummarizedExperiment) object, specialised to house measurements of retained
#'   introns and junction counts of alternative splice events.
#' * [apply_filters] - Use default or custom filters to remove alternative
#'   splicing or IR events pertaining to low-abundance genes and transcripts.
#' * [ASE-methods] - one-step method to perform differential alternate splice
#'   event (ASE) analysis on a NxtSE object using limma or DESeq2.
#' * [make_plot_data]: Functions that compile individual and group-mean percent
#'   spliced in (PSI) values of IR and alternative splice events; useful to
#'   produce scatter plots or heatmaps.
#' * [Plot_Coverage]: Generate RNA-seq coverage plots of
#'   individual samples or across samples grouped by user-specified conditions
#'
#' See the
#' [NxtIRF vignette](../doc/NxtIRF.html)
#' for worked examples on how to use NxtIRF
#'
#' @author Alex Wong
#'
#' @docType package
#' @name SpliceWiz-package
#' @aliases SpliceWiz-package
#' @keywords package
#' @references
#' Middleton R, Gao D, Thomas A, Singh B, Au A, Wong JJ, Bomane A, Cosson B,
#' Eyras E, Rasko JE, Ritchie W.
#' IRFinder: assessing the impact of intron retention on mammalian gene
#' expression. Genome Biol. 2017 Mar 15;18(1):51.
#' \url{https://doi.org/10.1186/s13059-017-1184-4}
#' @md
NULL
