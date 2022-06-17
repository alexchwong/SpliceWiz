test_that("SpliceWiz pipeline reproduces NxtSE object", {
    require(Rsubread)
    
    ref_path <- file.path(tempdir(), "refWithMapExcl")

    getResources(
        reference_path = ref_path,
        fasta = chrZ_genome(),
        gtf = chrZ_gtf()
    )

    generateSyntheticReads(
        reference_path = ref_path
    )

    setwd(ref_path)
    buildindex(basename = "./reference_index",
        reference = chrZ_genome())

    subjunc(
        index = "./reference_index",
        readfile1 = file.path(ref_path, "Mappability", "Reads.fa"),
        output_file = file.path(ref_path, "Mappability", "AlignedReads.bam"),
        useAnnotation = TRUE,
        annot.ext = chrZ_gtf(),
        isGTF = TRUE
    )

    calculateMappability(
        reference_path = ref_path,
        aligned_bam = file.path(ref_path, "Mappability", "AlignedReads.bam")
    )

    # Check 
    mappa_gr <- rtracklayer::import(
        file.path(ref_path, "Mappability", "MappabilityExclusion.bed.gz"),
        "bed"
    )
    
    expect_gr <- GenomicRanges::makeGRangesFromDataFrame(
        data.frame(
            seqnames = "chrZ",
            start = c(1, 101471),
            end = c(40, 101511),
            strand = "*"
        )
    )
    
    expect_identical(mappa_gr, expect_gr)
})
