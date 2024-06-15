test_that("alignRead works", {
    expect_error(alignRead(system.file("extdata", "ex_bam.bam", 
                                        package = "readAlign"),
                            "hg38", 16344))
    expect_equal(alignRead(system.file("extdata", "ex_bam.bam", 
                                        package = "readAlign"),
                            "hg38", 3261)$CIGAR[1], "18M1D132M")
    expect_equal(dim(alignRead(system.file("extdata", "ex_bam.bam", 
                                            package = "readAlign"),
                                "hg38", 3261))[2], 3)
    expect_error(alignRead(system.file("extdata", "ex_bam.bam", 
                                        package = "readAlign"),
                            "hg45", 3261))
    expect_error(alignRead("test.bam","hg38", 3261))
})

test_that("alignAllReads works", {
    expect_equal(dim(alignAllReads(system.file("extdata", "ex_bam.bam", 
                                                package = "readAlign"),
                                    "hg38")), c(5445,3))
    expect_error(alignAllReads(system.file("extdata", "ex_bam.bam", 
                                            package = "readAlign"),
                                "hg45"))
    expect_error(alignAllReads("test.bam","hg38"))
})
