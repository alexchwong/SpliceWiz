test_that("NxtIRF OpenMP produces same output regardless of threads", {
    n_threads = SpliceWiz:::Has_OpenMP()
    if(n_threads == 0) {
        message("NxtIRF was compiled without OpenMP")
        return(TRUE)
    }
    n_threads = min(2, parallel::detectCores())
    if(n_threads == 1) {
        message("Only 1 OpenMP processor available")
        return(TRUE)
    }
    
    if(!file.exists(file.path(tempdir(), "02H003.bam"))) {
        bams = NxtIRF_example_bams()
    } else {
        bams = Find_Bams(tempdir())
    }
    if(!file.exists(file.path(tempdir(), "Reference", "IRFinder.ref.gz"))) {
        BuildReference(
            fasta = chrZ_genome(), gtf = chrZ_gtf(),
            reference_path = file.path(tempdir(), "Reference")
        )
    }

    for(i in seq_len(n_threads)) {
        IRFinder(bams$path[1], paste0("thread_", i),
            reference_path = file.path(tempdir(), "Reference"),
            output_path = file.path(tempdir(), "IRFinder_test_threads"),
            overwrite = TRUE, n_threads = i, verbose = TRUE
        )
    }

    for(i in seq(2, n_threads)) {
        expect_equal(
            openssl::md5(file(file.path(tempdir(), "IRFinder_test_threads", paste0("thread_", 1, ".txt.gz")))), 
            openssl::md5(file(file.path(tempdir(), "IRFinder_test_threads", paste0("thread_", i, ".txt.gz"))))
        )
        expect_equal(
            openssl::md5(file(file.path(tempdir(), "IRFinder_test_threads", paste0("thread_", 1, ".cov")))), 
            openssl::md5(file(file.path(tempdir(), "IRFinder_test_threads", paste0("thread_", i, ".cov"))))
        )
    }
    
})