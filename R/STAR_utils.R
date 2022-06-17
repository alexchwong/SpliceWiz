# Miscellaneous internal functions

# Creates a temporary FASTA file from locally-stored TwoBit
.STAR_get_FASTA <- function(reference_path) {
    genome.fa <- file.path(reference_path, "resource", "genome.fa")
    if (!file.exists(genome.fa)) {
        genome.fa <- paste0(genome.fa, ".temp")
        genome.2bit <- file.path(reference_path, "resource", "genome.2bit")
        if (!file.exists(genome.2bit)) {
            .log(paste(genome.2bit, "not found"))
        }
        .log("Extracting temp genome FASTA from TwoBit file", "message")
        rtracklayer::export(
            rtracklayer::import(TwoBitFile(genome.2bit)),
            genome.fa, "fasta"
        )
    }
    return(genome.fa)
}

# Creates a temporary unzipped GTF for STAR
.STAR_get_GTF <- function(reference_path) {
    transcripts.gtf <- file.path(reference_path, "resource", "transcripts.gtf")
    if (!file.exists(transcripts.gtf)) {
        if (!file.exists(paste0(transcripts.gtf, ".gz"))) {
            .log(paste(paste0(transcripts.gtf, ".gz"), "not found"))
        }
        .log("Extracting temp Annotation GTF from GZ file", "message")
        R.utils::gunzip(paste0(transcripts.gtf, ".gz"), remove = FALSE,
            overwrite = TRUE)
        file.rename(transcripts.gtf, paste0(transcripts.gtf, ".temp"))
        transcripts.gtf <- paste0(transcripts.gtf, ".temp")
    }
    return(transcripts.gtf)
}

.STAR_clean_temp_FASTA_GTF <- function(reference_path) {
    .log("Cleaning temp genome / gene annotation files", "message")
    genome.fa <- file.path(reference_path, "resource", "genome.fa.temp")
    transcripts.gtf <- file.path(reference_path,
        "resource", "transcripts.gtf.temp")
    if (file.exists(genome.fa)) file.remove(genome.fa)
    if (file.exists(transcripts.gtf)) file.remove(transcripts.gtf)
}
