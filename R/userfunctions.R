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
#' @author Fateema Bazzi\cr Politecnico di Milano\cr Maintainer: Fateema 
#' Bazzi\cr E-Mail: <fateemahani.bazzi@@mail.polimi.it>
#' @seealso \code{\link{alignAllReads}}\cr
#' @examples
#' alignRead(system.file("extdata", "ex_bam.bam", 
#'                         package = "readAlign"),"hg38", 1)
#' @importFrom Rsamtools asBam
#' @importFrom Rsamtools scanBam
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom BSgenome.Hsapiens.NCBI.GRCh38 BSgenome.Hsapiens.NCBI.GRCh38
#' @importFrom BiocGenerics unlist
#' @importFrom Biostrings strsplit
#' @importFrom BiocGenerics as.data.frame
#' @importFrom utils tail
#' @export

alignRead <- function(bamf, gnm="hg19", readidx){
  #verify file is sam or bam. if sam, convert to bam
  file.name <- BiocGenerics::unlist(Biostrings::strsplit(bamf, '[.]'))
  if (utils::tail(file.name, n=1)=='bam'){
    bf <- bamf
  } else if (utils::tail(file.name, n=1)=='sam'){
    bf <- Rsamtools::asBam(bamf) #converts to bam file, new file is stored in same directory as sam
  } else {
    stop("Invalid file format. Accepted file types are BAM and SAM.")
  }
  
  #get the proper genome
  if (gnm=="hg19"){
    g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else if (gnm=="hg38"){
    g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  } else if (gnm=="GRCh38"){
    g <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
  } else {
    stop("Unsupported genome. Available genome options are hg19, hg38, and GRCh38.")
  }
  
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
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom BSgenome.Hsapiens.NCBI.GRCh38 BSgenome.Hsapiens.NCBI.GRCh38
#' @importFrom BiocGenerics unlist
#' @importFrom Biostrings strsplit
#' @importFrom BiocGenerics as.data.frame
#' @importFrom stats na.omit
#' @importFrom utils tail
#' @export

alignAllReads <- function(bamf, gnm="hg19"){
  #verify file is sam or bam. if sam, convert to bam
  file.name <- BiocGenerics::unlist(Biostrings::strsplit(bamf, '[.]'))
  if (utils::tail(file.name, n=1)=='bam'){
    bf <- bamf
  } else if (utils::tail(file.name, n=1)=='sam'){
    bf <- asBam(bamf) #converts to bam file, new file is stored in same directory as sam
  } else {
    stop("Invalid file format. Accepted file types are BAM and SAM.")
  }
  
  #get the proper genome
  if (gnm=="hg19"){
    g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else if (gnm=="hg38"){
    g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  } else if (gnm=="GRCh38"){
    g <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
  } else {
    stop("Unsupported genome. Available genome options are hg19, hg38, and GRCh38.")
  }
  
  #--preprocessing--
  bam <- scanBam(bf)
  #make dataframe
  bamdf <- BiocGenerics::as.data.frame(bam[[1]])
  
  #remove NA
  bamdf <- bamdf[,!names(bamdf) %in% 
                   c("mrnm", "mpos")]
  bamdf <- stats::na.omit(bamdf)
  rownames(bamdf) <- 1:nrow(bamdf)
  
  #send each row through the pipeline
  apply(bamdf, 1, pipeline, makedf=TRUE, gnm=g)
  
  return(pkg.env$aligndf)
}