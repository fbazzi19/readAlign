test_that("fileCheck works", {
    expect_equal(readAlign:::fileCheck(system.file("extdata", "ex_bam.bam", 
                                                    package = "readAlign")),
                "ex_bam.bam")
    expect_error(readAlign:::fileCheck("test"))
})

test_that("genomeCheck works", {
    expect_equal(readAlign:::genomeCheck("hg19"),
                BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
    expect_equal(readAlign:::genomeCheck("hg38"),
                BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
    expect_equal(readAlign:::genomeCheck("GRCh38"),
                BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38)
    expect_error(readAlign:::genomeCheck("hg36"))
})

test_that("cigarReader works", {
    expect_equal(readAlign:::cigarReader("2M1D3M"), 
                data.frame(lengths=c("2","1","3"),type=c("M","D","M")))
    expect_equal(readAlign:::cigarReader("100M"), 
                data.frame(lengths=c("100"),type=c("M")))
})
