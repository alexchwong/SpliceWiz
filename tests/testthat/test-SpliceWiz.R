test_that("SpliceWiz pipeline reproduces NxtSE object", {
    bams <- SpliceWiz_example_bams()
    chr_alias <- data.frame(old = "chrZ", new = "chrZ")

    buildRef(
        fasta = chrZ_genome(), 
        gtf = chrZ_gtf(),
        reference_path = file.path(tempdir(), "Reference"),
        chromosome_aliases = chr_alias
    )
    
    processBAM(bams$path, bams$sample,
        reference_path = file.path(tempdir(), "Reference"),
        output_path = file.path(tempdir(), "SpliceWiz_Output"),
        overwrite = TRUE,
        n_threads = 1
    )
    expr <- findSpliceWizOutput(file.path(tempdir(), "SpliceWiz_Output"))
    
    collateData(expr, 
        reference_path = file.path(tempdir(), "Reference"),
        output_path = file.path(tempdir(), "Collated_output"),
        novelSplicing = FALSE, overwrite = TRUE,
    )

    se <- makeSE(collate_path = file.path(tempdir(), "Collated_output"))
    
    # Test identical assays
    se_realized <- realize_NxtSE(se)
    
    se_compare <- SpliceWiz_example_NxtSE()
    
    expect_equal(
        assay(se_realized, "Included"), 
        assay(se_compare, "Included")
    )

    expect_equal(
        assay(se_realized, "Excluded"), 
        assay(se_compare, "Excluded")
    )

    expect_equal(
        assay(se_realized, "minDepth"), 
        assay(se_compare, "minDepth")
    )

    expect_equal(
        assay(se_realized, "Depth"), 
        assay(se_compare, "Depth")
    )

    expect_equal(
        assay(se_realized, "Coverage"), 
        assay(se_compare, "Coverage")
    )

    expect_equal(
        sampleQC(se_realized)[,-1], 
        sampleQC(se_compare)[,-1]
    )

    for(i in seq_len(ncol(se))) {
        cov1 <- getCoverage(covfile(se_realized)[i])
        cov2 <- getCoverage(covfile(se_compare)[i])
        expect_equal(cov1, cov2)

        expect_equal(
            tools::md5sum(covfile(se_realized)[i]), 
            tools::md5sum(covfile(se_compare)[i])
        )
    }


    expect_equal(
        up_inc(se_realized), 
        up_inc(se_compare)
    )

    expect_equal(
        up_exc(se_realized), 
        up_exc(se_compare)
    )

    expect_equal(
        down_inc(se_realized), 
        down_inc(se_compare)
    )

    expect_equal(
        down_exc(se_realized), 
        down_exc(se_compare)
    )

})
