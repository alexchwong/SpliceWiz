test_that("SpliceWiz is able to download stuff from Ensembl", {
    httplink <- "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/README"
    
    fn <- .parse_valid_file(httplink, force_download = TRUE)
    
    message("MD5 is ", tools::md5sum(fn))
    expect_equal(1,1)
})
