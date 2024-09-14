test_that("SpliceWiz buildRef is behaving normally", {
    fasta <- "ftp://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    gtf <- "ftp://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz"
    ref_path <- file.path(tempdir(), "realRef")
    
    message("Running buildRef with parameters")
    message("fasta = ", fasta)
    message("gtf = ", gtf)
    buildRef(
        fasta = fasta, 
        gtf = gtf,
        reference_path = ref_path,
        genome_type = "hg38"
    )
    expect_equal(file.exists(file.path(ref_path, "SpliceWiz.ref.gz")),TRUE)
    
    fasta <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Gallus_gallus/latest_assembly_versions/GCF_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w/GCF_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w_genomic.fna.gz"
    gtf <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Gallus_gallus/latest_assembly_versions/GCF_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w/GCF_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w_genomic.gtf.gz"
    ref_path <- file.path(tempdir(), "realRef2")
    
    message("Running buildRef with parameters")
    message("fasta = ", fasta)
    message("gtf = ", gtf)
    buildRef(
        fasta = fasta, 
        gtf = gtf,
        reference_path = ref_path,
        genome_type = "hg38"
    )
    expect_equal(file.exists(file.path(ref_path, "SpliceWiz.ref.gz")),TRUE)
})
