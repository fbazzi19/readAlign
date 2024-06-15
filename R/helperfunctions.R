#global environment variables
pkg.env <- new.env(parent = emptyenv())
pkg.env$alignedread <- ""
pkg.env$alignedref <- ""
pkg.env$aligndf <- data.frame(CIGAR=character(),
                              Ref.=character(), 
                              Read=character())


#' CIGAR DataFrame
#' 
#' This function takes the CIGAR string and converts it into a dataframe
#' that is able to be read by downstream functions.
#' 
#' @usage cigarReader(cigar)
#' @param cigar CIGAR string associated with a read
#' @return DataFrame where each row contains the number of characters
#' a symbol applies to in the first column and the relevant symbol in 
#' the second column
#' @author Fateema Bazzi\cr Politecnico di Milano\cr Maintainer: Fateema 
#' Bazzi\cr E-Mail: <fateemahani.bazzi@@mail.polimi.it>
#' @seealso \code{\link{pipeline}}\cr
#' @importFrom Biostrings strsplit
#' @importFrom BiocGenerics cbind
#' @importFrom BiocGenerics lapply
#' @keywords internal


cigarReader <- function(cigar){
  #split cigar into sections
  #first split to only retain numbers
  lengths <- Biostrings::strsplit(cigar, '[^0-9]')[[1]]
  #then grab the letters
  type <- BiocGenerics::lapply(Biostrings::strsplit(cigar, '[0-9]'), function(z){ z[!is.na(z) & z != ""]})[[1]]
  #bind into one dataframe
  cigardf <- BiocGenerics::cbind(lengths, type)
  return(cigardf)
}


#' Read Alter-er
#' 
#' This function takes the CIGAR dataframe and a row in the
#' dataframe that indicates changes needed for alignment, and alters
#' the read and reference sequences as necessary based on the CIGAR.
#' 
#' @usage alterReads(idx, cdf)
#' @param idx relevant row in CIGAR dataframe
#' @param cdf CIGAR dataframe, output from cigarreader
#' @author Fateema Bazzi\cr Politecnico di Milano\cr Maintainer: Fateema 
#' Bazzi\cr E-Mail: <fateemahani.bazzi@@mail.polimi.it>
#' @seealso \code{\link{cigarReader}}\cr
#' \code{\link{alignment}}\cr
#' @importFrom BiocGenerics paste
#' @keywords internal


alterReads <- function(idx, cdf){
  n <- sum(as.integer(cdf[1:(idx-1),1])) #where alignment is referring to
  #string of dashes of the correct length
  dashes <- BiocGenerics::paste(replicate(as.integer(cdf[idx,1]), "-"), 
                                collapse = "")
  
  if (cdf[idx,2] =='I'){
    #insert dashes into ref sequence
    pkg.env$alignedref<- 
      BiocGenerics::paste(substr(pkg.env$alignedref, 1, n), dashes, 
                          substr(pkg.env$alignedref, (n+1), 
                                 (nchar(pkg.env$alignedref)-nchar(dashes))),
                          sep = "")
  } 
  else if (cdf[idx,2] =='H'){
    if (idx==1){
      #insert dashes into start of read
      pkg.env$alignedread <- BiocGenerics::paste(dashes, pkg.env$alignedread, 
                                                 sep = "")
    } else{
      #insert dashes at end of read
      pkg.env$alignedread <- BiocGenerics::paste(pkg.env$alignedread, dashes, 
                                                 sep = "")
    }
    
  }else {
    #insert dashes at relevant point in read
    pkg.env$alignedread<- 
      BiocGenerics::paste(substr(pkg.env$alignedread, 1, n), dashes, 
                          substr(pkg.env$alignedread, n+1, 
                                 nchar(pkg.env$alignedread)), sep = "")
  }
}

#' Alignment
#' 
#' This function sends the read and reference sequence to be aligned 
#' based on the CIGAR dataframe and either prints the alignment of
#' a single read or stores the alignments of the reads 
#' 
#' @usage alignment(cigar, cigardf, read, ref, makedf)
#' @param cigar CIGAR string associated with a read
#' @param cigardf dataframe of CIGAR string contents from cigarreader
#' @param read read sequence from BAM file
#' @param ref reference sequence from genome corresponding to read
#' @param makedf bool indicating if a dataframe of the results should be made
#' @author Fateema Bazzi\cr Politecnico di Milano\cr Maintainer: Fateema 
#' Bazzi\cr E-Mail: <fateemahani.bazzi@@mail.polimi.it>
#' @seealso \code{\link{cigarReader}}\cr
#' \code{\link{alterReads}}\cr
#' @importFrom BiocGenerics lapply
#' @keywords internal

alignment <- function(cigar, cigardf, read, ref, makedf){
  #M, =, X, S: leave as is
  #I: add dashes to the reference seq
  #D, P, N: add dashes to the read seq
  #H: add dashes to beginning/end of the read seq
  
  #define as global variables so they can be properly altered with lapply
  pkg.env$alignedread <-read
  pkg.env$alignedref <-ref
  
  #find indexes where dashes are being inserted
  idxchng <- which(cigardf[,2]!='M'&cigardf[,2]!='X'&
                     cigardf[,2]!='='&cigardf[,2]!='S')
  
  #apply alterations using the indexes that call for it
  BiocGenerics::lapply(idxchng, alterReads, cdf=cigardf)
  
  if (makedf){
    row <- data.frame(cigar, pkg.env$alignedref, pkg.env$alignedread)
    names(row) <- c("CIGAR", "Ref.", "Read")
    pkg.env$aligndf <- rbind(pkg.env$aligndf, row)
  } else{ #print
    #print aligned read
    print(paste0("CIGAR: ", cigar))
    print(paste0("Ref.: ", pkg.env$alignedref))
    print(paste0("Read: ", pkg.env$alignedread))
  }
  
}


#' Pipeline
#' 
#' This function gets the relevant information from the read, 
#' generates the CIGAR dataframe, finds the corresponding reference 
#' sequence off the reference genome, accounts for any hard clipping 
#' present at the start of a read, and sends all of this information
#' the the alignment function.
#' 
#' @usage pipeline(read, makedf, gnm)
#' @param read row from BAM file corresponding to the read being aligned
#' @param makedf bool indicating if a dataframe of the results should be made
#' @param gnm reference genome reads are aligned to
#' @author Fateema Bazzi\cr Politecnico di Milano\cr Maintainer: Fateema 
#' Bazzi\cr E-Mail: <fateemahani.bazzi@@mail.polimi.it>
#' @seealso \code{\link{cigarReader}}\cr
#' \code{\link{alignment}}\cr
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings getSeq
#' @keywords internal

pipeline<- function(read, makedf=FALSE, gnm){ #read is the entire row from bam
  cig <- read['cigar'][[1]] #CIGAR string
  readseq <- read['seq'][[1]] #read seq
  cigdf <- cigarReader(cig) #CIGAR string in readable format
  
  #Account for any beginning hard clipping offset
  hoffset<- 0
  if (cigdf[1,2]=='H'){
    hoffset <- as.integer(cigdf[1,1])
  }
  #get the sequence of the reference genome
  ir <- IRanges::IRanges((as.integer(read['pos'][[1]])-hoffset),
                               width = sum(as.integer(cigdf[,1])))
  range <- GenomicRanges::GRanges(read['rname'][[1]],ir)
  
  refseq <- Biostrings::getSeq(gnm,range)
  
  
  alignment(cig, cigdf, readseq, as.character(refseq), makedf)
  
}