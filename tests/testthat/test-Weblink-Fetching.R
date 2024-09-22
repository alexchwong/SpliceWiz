test_that("SpliceWiz is able to download stuff from Ensembl", {
    httplink <- "https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/README"
    trueMD5 <- "815bad9eaf7d384766db42338d286b60" # as of 2024-09-14
    
    fn <- .parse_valid_file(httplink, force_download = TRUE)
    
    message("MD5 is ", tools::md5sum(fn))
    expect_equal(as.character(tools::md5sum(fn)), trueMD5)
})