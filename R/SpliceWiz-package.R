#' SpliceWiz: efficient and precise alternative splicing analysis in R
#'
#' SpliceWiz is a computationally efficient and user friendly workflow that
#' analyses aligned short-read RNA sequencing for differential
#' intron retention and alternative splicing. 
#'
#' @details
#' SpliceWiz uses isoform-specific alignments to quantify percent-spliced-in 
#' ratios (i.e. ratio of the "included" isoform, as a proportion of "included" 
#' and "excluded" isoforms). For intron retention (IR), the abundance of the 
#' intron-retaining transcript (included isoform) is quantified using the 
#' trimmed-mean depth of intron coverage with reads, whereas the spliced
#' transcript (excluded isoform) is measured as the splicing of the intron as 
#' well as that of overlapping introns (since splicing of any overlapping intron
#' implies the intron of interest is not retained). For other forms of 
#' alternative splicing, junction reads (reads aligned across splice junctions) 
#' are used to quantify included and excluded isoforms.
#' 
#' SpliceWiz processes BAM files (aligned RNA sequencing) using 
#' [ompBAM::ompBAM-package]. ompBAM
#' is a C++ library that allows R packages (via the Rcpp framework) to 
#' efficiently read BAM files using OpenMP-based multi-threading. SpliceWiz
#' processes BAM files via the [processBAM] function, using a splicing and
#' intron reference built from any given genome / gene annotation resource
#' using the [buildRef] function. [processBAM] generates two outputs per
#' BAM file: a `txt.gz` file which is a gzip-compressed text file with multiple
#' tables, containing information including junction read counts and intron
#' retention metrics. This output is very similar to that of 
#' [IRFinder](https://github.com/williamritchie/IRFinder), as the analysis
#' steps of SpliceWiz's BAM processing was built on an improved version of
#' IRFinder's source code (version 1.3.1). Additionally, [processBAM] outputs
#' a COV file, which is a binary bgzf-compressed file that contains
#' strand-specific coverage data.
#'
#' Once individual files have been analysed, SpliceWiz compiles a dataset using
#' these individual outputs, using [collateData]. This function unifies 
#' junctions detected across the dataset, and generates included / excluded
#' counts of all putative IR events and annotated alternative splicing events
#' (ASEs). This dataset is exported as a collection of files including an
#' H5 database. The data is later imported into the R session using the
#' [makeSE] function, as a \linkS4class{NxtSE} object.
#' 
#' The \linkS4class{NxtSE} object is a specialized 
#' \linkS4class{SummarizedExperiment} object tailored for use in SpliceWiz.
#' Annotation of rows provide information about ASEs via [rowData], while
#' columns allows users to provide annotations via [colData].
#'
#' SpliceWiz offers several novel filters via the \linkS4class{ASEFilter}
#' class. See [ASEFilter] for details.
#'
#' Once the \linkS4class{NxtSE} is annotated and filtered, differential
#' analysis is performed, using limma, DESeq2 or DoubleExpSeq wrappers.
#' These wrappers model isoform counts as log-normal, negative-binomial,
#' or beta-binomial distributions, respectively. See [ASE-methods] for details.
#'
#' Finally, SpliceWiz provides visualisation tools to illustrate alternative
#' splicing using coverage plots, including a novel method to normalise RNA-seq
#' coverage grouped by experimental condition. This approach accounts for
#' variations introduced by sequenced library size and gene expression. 
#' SpliceWiz efficiently computes and visualises means and variations in 
#' per-nucleotide coverage depth across alternate exons in genomic loci.
#'
#' The main functions are:
#'
#' * [Build-Reference-methods] - Prepares genome and gene annotation
#'   references from FASTA and GTF files and synthesizes the SpliceWiz reference
#'   for processing BAM files, collating the \linkS4class{NxtSE} object.
#' * [processBAM] - OpenMP/C++ based algorithm to analyse
#'   single or multiple BAM files.
#' * [collateData] - Collates an experiment based on multiple IRFinder outputs
#'   for individual samples, into one unified H5-based data structure.
#' * [makeSE] - Constructs a \linkS4class{NxtSE} (H5-based
#'   SummarizedExperiment) object, specialised to house measurements of retained
#'   introns and junction counts of alternative splice events.
#' * [applyFilters] - Use default or custom filters to remove alternative
#'   splicing or IR events pertaining to low-abundance genes and transcripts.
#' * [ASE-methods] - one-step method to perform differential alternate splice
#'   event (ASE) analysis on a NxtSE object using limma or DESeq2.
#' * [make_plot_data]: Functions that compile individual and group-mean percent
#'   spliced in (PSI) values of IR and alternative splice events; useful to
#'   produce scatter plots or heatmaps.
#' * [plotCoverage]: Generate RNA-seq coverage plots of
#'   individual samples or across samples grouped by user-specified conditions
#'
#' See the
#' [SpliceWiz Quick-Start](../doc/SW_QuickStart.html)
#' for worked examples on how to use SpliceWiz
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
