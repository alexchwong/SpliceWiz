test_that("SpliceWiz buildRef is behaving normally", {
    fasta <- "ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    gtf <- "ftp://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"
    
    message("Running buildRef with parameters")
    message("fasta = ", fasta)
    message("gtf = ", gtf)
    buildRef(
        fasta = fasta, 
        gtf = gtf,
        reference_path = file.path(tempdir(), "Reference"),
        genome_type = "hg38"
    )
    expect_equal(file.exists(file.path(tempdir(), "Reference", "SpliceWiz.ref.gz")),TRUE)
})
