test_that("buildRef fixes genome seqnames if contains descriptions", {
    genome <- rtracklayer::import(chrZ_genome(), "fasta")
    tmpFastaFile <- tempfile()
    names(GenomeInfoDb::seqinfo(genome)) <- "chrZ Description of chromosome"
    rtracklayer::export(genome, tmpFastaFile, "fasta")

    refPath <- file.path(tempdir(), "Reference_bR_sn")
    buildRef(
        fasta = tmpFastaFile, 
        gtf = chrZ_gtf(),
        reference_path = refPath
    )
    
    file.remove(tmpFastaFile)
    refFN <- file.path(refPath, "SpliceWiz.ref.gz")
    expect_equal(file.exists(refFN), TRUE)
})
