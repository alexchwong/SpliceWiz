#' View SpliceWiz Reference in read-able data frames
#'
#' These functions allow users to construct tables containing SpliceWiz's
#' reference of alternate splicing events, intron retention events, and
#' other relevant data
#'
#' @param reference_path The directory containing the SpliceWiz reference
#' @param directional (default `TRUE`) Whether to view IR events for stranded
#'   RNAseq `TRUE` or unstranded protocol `FALSE`
#' @return A data frame containing the relevant info. See details
#' @examples
#' ref_path <- file.path(tempdir(), "Reference")
#'
#' buildRef(
#'     reference_path = ref_path,
#'     fasta = chrZ_genome(),
#'     gtf = chrZ_gtf()
#' )
#'
#' df <- viewASE(ref_path)
#'
#' df <- viewIR(ref_path, directional = TRUE)
#'
#' df <- viewIntrons(ref_path)
#'
#' df <- viewIR_NMD(ref_path)
#'
#' df <- viewExons(ref_path)
#'
#' df <- viewGenes(ref_path)
#'
#' df <- viewProteins(ref_path)
#'
#' df <- viewTranscripts(ref_path)
#'
#' @seealso [Build-Reference-methods]
#' @name View-Reference-methods
#' @md
NULL

#' @describeIn View-Reference-methods Outputs summary of alternative splicing
#' events constructed by SpliceWiz
#' @export
viewASE <- function(reference_path) {
    .validate_reference(reference_path)
    
    targetFile <- file.path(reference_path, "fst", "Splice.fst")
    if(!file.exists(targetFile))
        .log(paste(targetFile, "does not exist!"))
        
    tmp <- read.fst(targetFile,
        columns = c(
            "EventType", "EventID", "EventName",
            "Event1a", "Event1b", "Event2a", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "transcript_id_a", "transcript_name_a",
            "transcript_id_b", "transcript_name_b"
        )
    )
    
    # Change names
    colnames(tmp) <- c(
        "EventType", "EventID", "EventName",
        "Inc_Event1", "Exc_Event1", "Inc_Event2", "Exc_Event2",
        "Inc_gene_id", "Exc_gene_id", "EventRegion",
        "Inc_transcript_id", "Inc_transcript_name",
        "Exc_transcript_id", "Exc_transcript_name"    
    )
    
    return(tmp)
}

#' @describeIn View-Reference-methods Outputs summary of assessed IRFinde-like
#' IR events, constructed by SpliceWiz
#' @export
viewIR <- function(reference_path, directional = TRUE) {
    .validate_reference(reference_path)
    
    if(directional) {
        targetFile <- file.path(reference_path, "fst", "Introns.Dir.fst")
    } else {
        targetFile <- file.path(reference_path, "fst", "Introns.ND.fst")    
    }
    if(!file.exists(targetFile))
        .log(paste(targetFile, "does not exist!"))
        
    if(directional) {
        tmp <- read.fst(targetFile, columns = c(
            "intron_id", "seqnames", "intron_start", "intron_end", "strand",
            "inclbases", "exclbases",
            "gene_name", "transcript_id",
            "known_exon_dir"
        ))
        colnames(tmp) <- c(
            "intron_id", "seqnames", "start", "end", "strand",
            "includedBases", "excludedBases",
            "gene_name", "transcript_id",
            "overlaps_known_exon"
        )
        return(tmp)
    } else {
        tmp <- read.fst(targetFile, columns = c(
            "intron_id", "seqnames", "intron_start", "intron_end", "strand",
            "inclbases", "exclbases",
            "gene_name", "transcript_id",
            "known_exon_nd", "antiover", "antinear"
        ))
        colnames(tmp) <- c(
            "intron_id", "seqnames", "start", "end", "strand",
            "includedBases", "excludedBases",
            "gene_name", "transcript_id",
            "overlaps_known_exon", "antisense_overlap", "antisense_nearby"
        )
        return(tmp)
    }
}

#' @describeIn View-Reference-methods Outputs summary of all introns from
#' the annotation, constructed by SpliceWiz
#' @export
viewIntrons <- function(reference_path) {
    .validate_reference(reference_path)
    
    targetFile <- file.path(reference_path, "fst", "junctions.fst")
    if(!file.exists(targetFile))
        .log(paste(targetFile, "does not exist!"))

    columns <- c(
        "intron_id", "Event", 
        "seqnames", "start", "end", "width", "strand",
        "intron_number", "gene_name", "gene_id", 
        "transcript_name", "transcript_id", "transcript_biotype",
        "splice_motif", "transcript_support_level", 
        "protein_id", "ccds_id"
    )

    tmp <- read.fst(targetFile)

    columns <- columns[columns %in% colnames(tmp)]
    tmp <- tmp[, columns]
    tmp
}

#' @describeIn View-Reference-methods Outputs information for every intron -
#' whether retention of the intron will convert the transcript to an NMD
#' substrate
#' @export
viewIR_NMD <- function(reference_path) {
    .validate_reference(reference_path)
    
    targetFile <- file.path(reference_path, "fst", "IR.NMD.fst")
    
    if(!file.exists(targetFile))
        .log(paste(targetFile, "does not exist!"))
        
    tmp <- read.fst(targetFile)
    
    tmp <- tmp[, c(
        "intron_id", "transcript_id", "intron_type",
        "splice_is_NMD", "IRT_is_NMD"
    )]
    colnames(tmp) <- c(
        "intron_id", "transcript_id", "intron_type",
        "spliced_transcript_is_NMD", "IR_transcript_is_NMD"
    )
    return(tmp)
}

#' @describeIn View-Reference-methods Outputs information for every exon 
#' from the annotation.
#' @export
viewExons <- function(reference_path) {
    .validate_reference(reference_path)
    
    targetFile <- file.path(reference_path, "fst", "Exons.fst")
    
    if(!file.exists(targetFile))
        .log(paste(targetFile, "does not exist!"))
        
    tmp <- read.fst(targetFile)
    excl_columns <- c(
        "score", "phase",
        "protein_id", "protein_version",
        "gene_group_stranded", "exon_group_stranded",
        "gene_group_unstranded", "exon_group_unstranded"
    )
    
    tmp <- tmp[, colnames(tmp)[!colnames(tmp) %in% excl_columns]]
    return(tmp)
}

#' @describeIn View-Reference-methods Outputs information for every gene 
#' from the annotation.
#' @export
viewGenes <- function(reference_path) {
    .validate_reference(reference_path)
    
    targetFile <- file.path(reference_path, "fst", "Genes.fst")
    
    if(!file.exists(targetFile))
        .log(paste(targetFile, "does not exist!"))
        
    tmp <- read.fst(targetFile)
    excl_columns <- c(
        "score", "phase",
        "transcript_id", "transcript_version", "transcript_name",
        "transcript_source", "transcript_biotype", "tag",                      "transcript_support_level", "ccds_id", 
        "exon_number", "exon_id", "exon_version", "protein_id",               "protein_version", "gene_display_name", 
        "gene_group_stranded", "gene_group_unstranded"
    )
    
    tmp <- tmp[, colnames(tmp)[!colnames(tmp) %in% excl_columns]]
    return(tmp)
}

#' @describeIn View-Reference-methods Outputs information for every 
#' protein-coding exon from the annotation.
#' @export
viewProteins <- function(reference_path) {
    .validate_reference(reference_path)
    
    targetFile <- file.path(reference_path, "fst", "Proteins.fst")
    
    if(!file.exists(targetFile))
        .log(paste(targetFile, "does not exist!"))
        
    tmp <- read.fst(targetFile)
    excl_columns <- c(
        "score", "exon_id", "exon_version"
    )
    
    tmp <- tmp[, colnames(tmp)[!colnames(tmp) %in% excl_columns]]
    return(tmp)
}

#' @describeIn View-Reference-methods Outputs information for every 
#' transcript from the annotation.
#' @export
viewTranscripts <- function(reference_path) {
    .validate_reference(reference_path)
    
    targetFile <- file.path(reference_path, "fst", "Transcripts.fst")
    
    if(!file.exists(targetFile))
        .log(paste(targetFile, "does not exist!"))
        
    tmp <- read.fst(targetFile)
    excl_columns <- c(
        "score", "phase",
        "exon_number", "exon_id", "exon_version", "protein_id",               "protein_version"
    )
    
    tmp <- tmp[, colnames(tmp)[!colnames(tmp) %in% excl_columns]]
    return(tmp)
}