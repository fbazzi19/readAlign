#' Align Read
#' 
#' This function aligns one specified read from a BAM or SAM file to 
#' the reference genome and prints out the resulting alignment. Function
#' will return an error in the case that the CIGAR string is not available
#' for the read. In the case of a SAM file, the function will first 
#' create a corresponding BAM file, which will be saved in the same 
#' location as the SAM file.
#' 
#' 
#' @usage alignRead(bamf, gnm, readidx)
#' @param bamf path to BAM or SAM file
#' @param gnm which reference genome to use, accepted options are "hg19",
#' "hg38", or "GRCh38"
#' @param readidx index number of read to align
#' @returns Nothing, prints the results of the alignment instead.
#' @author Fateema Bazzi\cr Politecnico di Milano\cr Maintainer: Fateema 
#' Bazzi\cr E-Mail: <fateemahani.bazzi@@mail.polimi.it>
#' @seealso \code{\link{alignAllReads}}\cr
#' @examples
#' alignRead(system.file("extdata", "ex_bam.bam", 
#'                         package = "readAlign"),"hg38", 1)
#' @importFrom Rsamtools asBam
#' @importFrom Rsamtools scanBam
#' @importFrom BiocGenerics as.data.frame
#' @export

alignRead <- function(bamf, gnm="hg19", readidx){
    #send file to function to ensure it is acceptable
    bf <- fileCheck(bamf)
  
    #get the proper genome from helper function
    g <- genomeCheck(gnm)
  
  
    #--preprocessing--
    bam <- Rsamtools::scanBam(bf)
    #make dataframe
    bamdf <- BiocGenerics::as.data.frame(bam[[1]])
  
    #check if user specified read is NA
    if (is.na(bamdf[readidx,]$cigar)) {
        stop("CIGAR string not available for specified read.")
    }
  
    #grab desired read
    read <- bamdf[readidx,]
  
    #read and corresponding reference sequence through the pipeline
    pipeline(read, gnm = g)
}

#' Align All Reads
#' 
#' This function aligns all reads from a BAM or SAM file that don't
#' contain NA to the reference genome and returns the resulting 
#' alignment in a DataFrame. In the case of a SAM file, the function 
#' will first create a corresponding BAM file, which will be saved in 
#' the same location as the SAM file.
#' 
#' 
#' @usage alignAllReads(bamf, gnm)
#' @param bamf path to BAM or SAM file
#' @param gnm which reference genome to use, accepted options are "hg19",
#' "hg38", or "GRCh38"
#' @returns Data frame containing the aligned reads where each row is 
#' a read alignment and the three columns are CIGAR, Ref., and Read.
#' @author Fateema Bazzi\cr Politecnico di Milano\cr Maintainer: Fateema 
#' Bazzi\cr E-Mail: <fateemahani.bazzi@@mail.polimi.it>
#' @seealso \code{\link{alignRead}}\cr
#' @examples
#' alignAllReads(system.file("extdata", "ex_bam.bam", 
#'                         package = "readAlign"),"hg38")
#' @importFrom Rsamtools asBam
#' @importFrom Rsamtools scanBam
#' @importFrom BiocGenerics as.data.frame
#' @importFrom stats na.omit
#' @export

alignAllReads <- function(bamf, gnm="hg19"){
    #send file to function to ensure it is acceptable
    bf <- fileCheck(bamf)
  
    #get the proper genome from helper function
    g <- genomeCheck(gnm)
  
    #--preprocessing--
    bam <- scanBam(bf)
    #make dataframe
    bamdf <- BiocGenerics::as.data.frame(bam[[1]])
  
    #remove NA
    bamdf <- bamdf[,!names(bamdf) %in% 
                    c("mrnm", "mpos")]
    bamdf <- stats::na.omit(bamdf)
    rownames(bamdf) <- seq_len(nrow(bamdf))
  
    #send each row through the pipeline
    apply(bamdf, 1, pipeline, makedf=TRUE, gnm=g)
  
    return(pkg.env$aligndf)
}