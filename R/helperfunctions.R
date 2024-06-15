#global environment variables
pkg.env <- new.env(parent = emptyenv())
pkg.env$alignedread <- ""
pkg.env$alignedref <- ""
pkg.env$aligndf <- data.frame(CIGAR=character(),
                                Ref.=character(), 
                                Read=character())


#' File Checking
#' 
#' This function ensures that the file specified by the user exists
#' and is of an accepted type, and gives the file path. If SAM, the 
#' file path returned is for the newly created BAM file.
#' 
#' @usage fileCheck(filename)
#' @param filename string specifying the path to the user file.
#' @returns The resulting path to BAM file
#' @author Fateema Bazzi\cr Politecnico di Milano\cr Maintainer: Fateema 
#' Bazzi\cr E-Mail: <fateemahani.bazzi@@mail.polimi.it>
#' @importFrom Biostrings strsplit
#' @importFrom BiocGenerics unlist
#' @importFrom utils tail
#' @keywords internal

fileCheck <- function(filename){
    if (!file.exists(filename)){
        stop("Could not find specified BAM/SAM file.")
    }
    
    #verify file is sam or bam. if sam, convert to bam
    file.name <- BiocGenerics::unlist(Biostrings::strsplit(filename, '[.]'))
    if (utils::tail(file.name, n=1)=='bam'){
        bf <- filename
    } else if (utils::tail(file.name, n=1)=='sam'){
        bf <- asBam(filename) #converts to bam file, new file is stored in 
        #same directory as sam
    } else {
        stop("Invalid file format. Accepted file types are BAM and SAM.")
    }
    
    return(bf)
}

#' Genome Checking
#' 
#' This function ensures that the reference genome specified by the 
#' user exists and is of an accepted type, and creates the genome 
#' object in R.
#' 
#' @usage genomeCheck(gnm)
#' @param gnm string specifying the reference genome.
#' @returns The resulting genome object
#' @author Fateema Bazzi\cr Politecnico di Milano\cr Maintainer: Fateema 
#' Bazzi\cr E-Mail: <fateemahani.bazzi@@mail.polimi.it>
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom BSgenome.Hsapiens.NCBI.GRCh38 BSgenome.Hsapiens.NCBI.GRCh38
#' @keywords internal

genomeCheck <- function(gnm){
    if (gnm=="hg19"){
        g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    } else if (gnm=="hg38"){
        g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    } else if (gnm=="GRCh38"){
        g <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    } else {
        stop("Unsupported genome. Available genome options are 
            hg19, hg38, and GRCh38.")
    }
    
    return(g)
}

#' CIGAR DataFrame
#' 
#' This function takes the CIGAR string and converts it into a dataframe
#' that is able to be read by downstream functions.
#' 
#' @usage cigarReader(cigar)
#' @param cigar CIGAR string associated with a read
#' @returns DataFrame where each row contains the number of characters
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
    type <- BiocGenerics::lapply(Biostrings::strsplit(cigar, '[0-9]'), 
                                function(z){ z[!is.na(z) & z != ""]})[[1]]
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
#' @returns Nothing
#' @author Fateema Bazzi\cr Politecnico di Milano\cr Maintainer: Fateema 
#' Bazzi\cr E-Mail: <fateemahani.bazzi@@mail.polimi.it>
#' @seealso \code{\link{cigarReader}}\cr
#' \code{\link{alignment}}\cr
#' @importFrom BiocGenerics paste
#' @keywords internal


alterReads <- function(idx, cdf){
    n <- sum(as.integer(cdf[seq_len(idx-1),1]))
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
    } else if (cdf[idx,2] =='H'){
        if (idx==1){
            #insert dashes into start of read
            pkg.env$alignedread <- BiocGenerics::paste(dashes, 
                                                        pkg.env$alignedread, 
                                                        sep = "")
        } else{
            #insert dashes at end of read
            pkg.env$alignedread <- BiocGenerics::paste(pkg.env$alignedread, 
                                                        dashes, 
                                                        sep = "")
        }
    
    } else{
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
#' based on the CIGAR dataframe and stores the alignment of the read
#' in a DataFrame. 
#' 
#' @usage alignment(cigar, cigardf, read, ref, makedf)
#' @param cigar CIGAR string associated with a read
#' @param cigardf dataframe of CIGAR string contents from cigarreader
#' @param read read sequence from BAM file
#' @param ref reference sequence from genome corresponding to read
#' @returns Nothing
#' @author Fateema Bazzi\cr Politecnico di Milano\cr Maintainer: Fateema 
#' Bazzi\cr E-Mail: <fateemahani.bazzi@@mail.polimi.it>
#' @seealso \code{\link{cigarReader}}\cr
#' \code{\link{alterReads}}\cr
#' @importFrom BiocGenerics lapply
#' @keywords internal

alignment <- function(cigar, cigardf, read, ref){
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
    
    #add alignment to a dataframe
    row <- data.frame(cigar, pkg.env$alignedref, pkg.env$alignedread)
    names(row) <- c("CIGAR", "Ref.", "Read")
    pkg.env$aligndf <- rbind(pkg.env$aligndf, row)
    
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
#' @param gnm reference genome reads are aligned to
#' @returns Nothing
#' @author Fateema Bazzi\cr Politecnico di Milano\cr Maintainer: Fateema 
#' Bazzi\cr E-Mail: <fateemahani.bazzi@@mail.polimi.it>
#' @seealso \code{\link{cigarReader}}\cr
#' \code{\link{alignment}}\cr
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings getSeq
#' @keywords internal

pipeline<- function(read, gnm){ #read is the entire row from bam
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
    
    alignment(cig, cigdf, readseq, as.character(refseq))
    
}