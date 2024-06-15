## ----style, echo = FALSE, results = 'asis'------------------------------------
 library(BiocStyle)

## ----echo = FALSE-------------------------------------------------------------
 library(knitr)

## ----echo=FALSE---------------------------------------------------------------
 library(readAlign)

## -----------------------------------------------------------------------------
alignRead(system.file("extdata", "ex_bam.bam", package = "readAlign"),"hg38", 1)

## ----eval=FALSE---------------------------------------------------------------
#  alignAllReads(system.file("extdata", "ex_bam.bam", package = "readAlign"),"hg38")

## -----------------------------------------------------------------------------
utils::sessionInfo()

