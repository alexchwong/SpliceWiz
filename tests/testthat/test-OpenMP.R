test_that("SpliceWiz OpenMP produces same output regardless of threads", {
    n_threads = SpliceWiz:::Has_OpenMP()
    if(n_threads == 0) {
        useOpenMP <- FALSE
        message("SpliceWiz was compiled without OpenMP")
        # return(TRUE)
    } else {
        useOpenMP <- TRUE
    }
    
    n_threads = min(2, parallel::detectCores())
    if(n_threads == 1) {
        message("Only 1 OpenMP processor available")
        return(TRUE)
    }
    
    if(!file.exists(file.path(tempdir(), "02H003.bam"))) {
        bams <- SpliceWiz_example_bams()
    } else {
        bams <- findBAMS(tempdir())
    }
    if(!file.exists(file.path(tempdir(), "Reference", "SpliceWiz.ref.gz"))) {
        buildRef(
            fasta = chrZ_genome(), gtf = chrZ_gtf(),
            reference_path = file.path(tempdir(), "Reference")
        )
    }

    for(i in seq_len(n_threads)) {
        processBAM(bams$path[1], paste0("thread_", i),
            reference_path = file.path(tempdir(), "Reference"),
            output_path = file.path(tempdir(), "pb_test_threads"),
            overwrite = TRUE, n_threads = i, verbose = TRUE,
            useOpenMP = useOpenMP
        )
    }

    for(i in seq(2, n_threads)) {
        expect_equal(
            tools::md5sum(file.path(tempdir(), "pb_test_threads", paste0("thread_", 1, ".txt.gz"))), 
            tools::md5sum(file.path(tempdir(), "pb_test_threads", paste0("thread_", i, ".txt.gz")))
        )
        expect_equal(
            tools::md5sum(file.path(tempdir(), "pb_test_threads", paste0("thread_", 1, ".cov"))), 
            tools::md5sum(file.path(tempdir(), "pb_test_threads", paste0("thread_", i, ".cov")))
        )
    }
    
})